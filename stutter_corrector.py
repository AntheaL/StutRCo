import h5py
import queue
import numpy as np
import logging
import queue
import pysam
import traceback
import multiprocessing

from operator import iadd
from multiprocessing import Manager, Process

from logger import start_process_logging, stop_process_logging
from reader import FilestreamBED

from read_util import (
    covers_str,
    expanded_seq,
    parse_changes,
    parse_insertions,
    fetch_indels,
    overlap_noise,
)

def str_to_regions(lines):
    res = []
    for line in lines:
        if line.startswith("#"):
            continue
        try:
            line = line.strip().replace("chr", "").upper().split()
            region = {
                "chrom": line[0],#.split('_')[0],
                "start": int(line[1]),
                "stop": int(line[2]),
                "motif_len": int(line[3]),
            }
            res.append(region)
        except (ValueError, IndexError):
            continue
    return res

class StutRCorrector(Process):
    def __init__(self, sites, bam, logs_dir, name, log_every=25000, batch_size=10000, pad_left=15, pad_right=20):
        super().__init__()
        multiprocessing.current_process().name = name
        self.sites = str_to_regions(sites)
        self.bam = pysam.AlignmentFile(bam)
        self.logs_dir = logs_dir
        self.log_every = log_every
        self.pad_left = pad_left
        self.pad_right = pad_right
        self.batch_size = batch_size
        self.iter = 0
        self.n_corr = 0
        self.n_umi = 0
        self.rm_reads = []
        self.info = dict(cell_freq=[], umi_freq=[], umi_count=[], chrom=[], pos=[])
        self.info["global"] = dict(n_umi=[], n_reads=[], chrom=[], pos=[])

    def run(self):
        log_handler=start_process_logging(self.logs_dir)
        try:
            for site in self.sites:
                self.iter += 1
                cell_dict = self.get_cell_families(site)
                self.process_site(site, cell_dict)
                if self.iter%self.log_every==0:
                    logging.info(
                        f"Iteration {self.iter}, removed {len(self.rm_reads)} reads across {self.n_umi} umi families."
                    )
        except Exception as e:
            tb = traceback.format_exc()
            logging.error(tb, e)
        stop_process_logging(log_handler)


    def get_cell_families(self, site):
        """
        Generate cell barcode clusters in a BAM file or within a region string.
        Reads are added to read_dict until a pair is found.
        """
        cell_dict = {}
        start = site["start"] - self.pad_left
        stop = site["stop"] + self.pad_right
        try:
            for read in self.bam.fetch(site["chrom"], max(start, 1), stop):
                if read.is_secondary or read.is_supplementary:
                    continue
                if not covers_str(read, site):
                    continue
                try:
                    cell_barcode = read.get_tag("CB")
                except:
                    continue
                cell_dict[cell_barcode] = cell_dict.get(cell_barcode, []) + [to_dict(read)]
        except ValueError:
            raise ValueError(f"Invalid site {site}")
        return cell_dict

    def process_site(self, site, cell_dict):
        n_umi = 0
        n_reads = 0
        cell_dict = get_umi_families(cell_dict)
        for cell_barcode, umi_families in cell_dict.items():
            self.n_umi += len(umi_families)
            self.parse_families(site, umi_families)
            n_umi += len(umi_families)
            self.info["global"]["n_umi"].append(n_umi)
            self.info["global"]["n_reads"].append(n_reads)
            self.info["global"]["chrom"].append(
                int(site["chrom"]) if site["chrom"].isdigit() else ord(site["chrom"])
            )
            self.info["global"]["pos"].append(site["start"])


    def parse_families(self, site, umi_families):

        indels_by_read = dict()
        indels_by_family = dict()
        n_by_indel = dict()
        n_reads = 0

        for UMI, family in umi_families.items():
            n_reads += len(family)
            indels = overlap_noise(family, site)
            uni, counts = np.unique(indels, return_counts=True)
            for u, c in zip(uni, counts):
                n_by_indel[u] = n_by_indel.get(u, 0) + c
            if len(uni) > 1:
                indels_by_read[UMI] = indels
                indels_by_family[UMI] = dict(zip(uni, counts))

        # finding most probable indel
        indel_by_umi = dict()
        for umi, v in indels_by_family.items():
            c = list(v.values())
            if len(c) > 1:
                indel_by_umi[umi] = max(v, key=lambda x: 10 * v[x] + n_by_indel[x])

        # assigning a matching read to each UMI and extending site
        start, stop = np.inf, 0

        n_reads = sum(list(n_by_indel.values()))
        for umi, x in indel_by_umi.items():
            i = next(i for i, v in enumerate(indels_by_read[umi]) if v == x)

            # TDOO: Extend read
            """
            start = np.min([read.pos for read in family])
            stop = np.max([read.pos+read.reference_len for read in family])
            rdict = family[i].to_dict()
            pad_read(rdict, start, stop)
            read = AlignedSegment.from_dict(rdict)
            """

            self.rm_reads += [r["name"] for j, r in enumerate(family) if j != i]
            self.info["cell_freq"].append(n_by_indel[x] / n_reads)
            self.info["umi_freq"].append(indels_by_family[umi][x] / len(indels_by_read[umi]))
            self.info["umi_count"].append(indels_by_family[umi][x])
            self.info["chrom"].append(
                int(site["chrom"]) if site["chrom"].isdigit() else ord(site["chrom"])
            )
            self.info["pos"].append(site["start"])


def get_umi_families(cell_dict):

    del_length_one(cell_dict)
    clean_dict = {}
    for cell_barcode, cell_reads in cell_dict.items():

        UMI_families = {}
        for i, read in enumerate(cell_reads):
            #logging.warning(f"Skipping read {i} which does not contain 'UB' tag.")
            umi = read["umi"]
            if umi is None:
                continue
            UMI_families[umi] = UMI_families.get(umi, []) + [read]
        del_length_one(UMI_families)

        if len(UMI_families) >= 1:
            clean_dict[cell_barcode] = UMI_families

    return clean_dict

def to_dict(read):
    return dict(
        name=read.query_name,
        ref_start=read.reference_start,
        ref_end=read.reference_end,
        sequence=read.query_alignment_sequence,
        cigartuples=read.cigartuples,
        umi = read.get_tag('UB') if read.has_tag('UB') else None
    )

def del_length_one(dict_):
    for key in [k for k in dict_.keys()]:
        if len(dict_[key]) <= 1:
            del dict_[key]
