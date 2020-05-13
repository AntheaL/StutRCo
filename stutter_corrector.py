import h5py
import sys
import os
import queue
import numpy as np
import logging
import queue
import pysam
import traceback
import multiprocessing
import psutil

from itertools import compress
from operator import iadd
from multiprocessing import Manager, Process

from logger import start_process_logging, stop_process_logging
from reader import FilestreamBED

from h5_util import save_h5
from read_util import (
    covers_str,
    expanded_seq,
    parse_changes,
    parse_insertions,
    fetch_indels,
    overlap_noise,
)


def str_to_regions(lines, chr_tag=""):
    while lines:
        line = lines.pop()
        if line.startswith("#"):
            continue
        line = line.strip().upper().replace("CHR", chr_tag).split()
        if not (line[0].isdigit() or line[0] in ["X", "Y"]):
            continue
        region = {
            "chrom": line[0],  # .split('_')[0],
            "start": int(line[1]),
            "stop": int(line[2]),
            "motif_len": int(line[3]),
        }
        yield region


class StutRCorrector(Process):
    def __init__(
        self,
        sites,
        bam,
        logs_dir,
        # send_end,
        name,
        log_every=25000,
        pad_left=15,
        pad_right=20,
        chr_tag="",
        **kwargs,
    ):
        super().__init__(name=name, **kwargs)
        self.sites = sites
        self.chr_tag = chr_tag
        self.bam = pysam.AlignmentFile(bam, "r")
        self.logs_dir = logs_dir
        self.log_every = int(log_every)
        # self.send_end = send_end
        self.pad_left = int(pad_left)
        self.pad_right = int(pad_right)
        self.iter = 0
        self.n_corr = 0
        self.n_umi = 0
        self.rm_reads = []
        self.info = dict()
        self.info["umi"] = dict(
            cfreq=[],
            freq=[],
            count=[],
            nucl=[],
            qual=[],
            avg_qual=[],
            motif_len=[],
            alleles=[],
            chrom=[],
            pos=[],
        )
        self.info["global"] = dict(
            n_umi=[],
            n_umi_corr=[],
            n_reads=[],
            n_reads_corr=[],
            nucl=[],
            motif_len=[],
            alleles=[],
            chrom=[],
            pos=[],
        )

    def run(self):
        log_handler = start_process_logging(self.logs_dir)
        try:
            for site in str_to_regions(self.sites):
                self.iter += 1
                start = site["start"] - self.pad_left
                stop = site["stop"] + self.pad_right
                cell_dict = self.get_cell_families(site)
                cell_dict = get_umi_families(cell_dict)
                for cell_barcode, umi_families in cell_dict.items():
                    self.parse_families(site, umi_families)
                if self.iter % self.log_every == 0:
                    logging.info(
                        "Iteration {}, removed {} reads across {} umi families out of {}. Memory usage: {}%.".format(
                            self.iter,
                            len(self.rm_reads),
                            self.n_corr,
                            self.n_umi,
                            psutil.virtual_memory().percent,
                        )
                    )
            logging.info(f"Saving {len(self.rm_reads)} reads to remove.")
            # self.send_end.send(self.rm_reads)
            rm_path = os.path.join(self.logs_dir, f"rm_{self.name}.txt")
            with open(rm_path, "w") as dst:
                dst.write("\n".join(self.rm_reads))
            dst_path = os.path.join(self.logs_dir, f"info_{self.name}.h5").replace(
                " ", "-"
            )
            logging.info(f"Saving info into temporary file {dst_path}.")
            save_h5(dst_path, self.info)
        except Exception as e:
            tb = traceback.format_exc()
            print(tb)
            raise e
        stop_process_logging(log_handler)

    def get_cell_families(self, site):
        """
        Generate cell barcode clusters in a BAM file or within a region string.
        Reads are added to read_dict until a pair is found.
        """
        cell_dict = {}
        start = site["start"] - self.pad_left
        stop = site["stop"] + self.pad_right
        for read in self.bam.fetch(site["chrom"], max(start, 1), stop):
            if read.is_secondary or read.is_supplementary:
                continue
            if not covers_str(read, site):
                continue
            try:
                cell_barcode = read.get_tag("CB")
            except:
                continue
            cell_dict[cell_barcode] = cell_dict.get(cell_barcode, []) + [read]

        return cell_dict

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
        n_umi_corr = 0
        for umi, v in indels_by_family.items():
            c = list(v.values())
            if len(c) > 1:
                n_umi_corr += 1
                indel_by_umi[umi] = max(v, key=lambda x: 10 * v[x] + n_by_indel[x])

        self.n_corr += n_umi_corr
        self.n_umi += len(umi_families)
        n_reads = sum(list(n_by_indel.values()))
        n_reads_corr = 0

        if not indel_by_umi:
            return

        nucls = np.unique(
            list(umi_families[umi][0].query_alignment_sequence), return_counts=True
        )
        nmax = np.argmax(nucls[1])
        nucl = (ord(nucls[0][nmax]) + nucls[1][nmax] / np.sum(nucls[1]),)

        for umi, x in indel_by_umi.items():
            # TDOO: Extend read
            """
            start = np.min([read.pos for read in family])
            stop = np.max([read.pos+read.reference_len for read in family])
            rdict = family[i].to_dict()
            pad_read(rdict, start, stop)
            read = AlignedSegment.from_dict(rdict)
            """

            is_concordant = [n == x for n in indels_by_read[umi]]
            n_reads_corr += sum(map(lambda x: 1 - x, is_concordant))
            self.rm_reads += [
                r.query_name for r, b in zip(umi_families[umi], is_concordant) if not b
            ]
            n_reads_umi = len(indels_by_read[umi])
            bq = map(lambda r: np.mean(r.query_qualities), umi_families[umi])

            umi_info = dict(
                cfreq=n_by_indel[x] / n_reads,
                freq=indels_by_family[umi][x] / n_reads_umi,
                count=indels_by_family[umi][x],
                alleles=len(indels_by_family[umi]),
                avg_qual=sum(bq) / n_reads_umi,
                qual=sum(compress(bq, is_concordant)) / sum(is_concordant),
                nucl=nucl,
                motif_len=site["motif_len"],
                chrom=int(site["chrom"])
                if site["chrom"].isdigit()
                else ord(site["chrom"]),
                pos=site["start"],
            )
            for k, v in umi_info.items():
                self.info["umi"][k].append(v)

        global_info = dict(
            n_umi=len(umi_families),
            n_umi_corr=n_umi_corr,
            n_reads=n_reads,
            n_reads_corr=n_reads_corr,
            nucl=nucl,
            motif_len=site["motif_len"],
            chrom=int(site["chrom"]) if site["chrom"].isdigit() else ord(site["chrom"]),
            alleles=len(n_by_indel),
            pos=site["start"],
        )
        for k, v in global_info.items():
            self.info["global"][k].append(v)


def get_umi_families(cell_dict):

    del_length_one(cell_dict)
    clean_dict = {}
    for cell_barcode, cell_reads in cell_dict.items():

        UMI_families = {}
        for i, read in enumerate(cell_reads):
            # logging.warning(f"Skipping read {i} which does not contain 'UB' tag.")
            umi = (read.get_tag("UB") if read.has_tag("UB") else None,)
            if umi is None:
                continue
            UMI_families[umi] = UMI_families.get(umi, []) + [read]
        del_length_one(UMI_families)

        if len(UMI_families) >= 1:
            clean_dict[cell_barcode] = UMI_families

    return clean_dict


def del_length_one(dict_):
    for key in [k for k in dict_.keys()]:
        if len(dict_[key]) <= 1:
            del dict_[key]
