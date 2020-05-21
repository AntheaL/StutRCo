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
from collections import OrderedDict
from read_util import (
    covers_str,
    expanded_seq,
    parse_changes,
    parse_insertions,
    fetch_indels,
    overlap_noise,
)


def str_to_regions(lines, site_keys=None, chr_tag=""):
    while lines:
        line = lines.pop()
        if line.startswith("#"):
            continue
        line = line.strip().upper().replace("CHR", chr_tag).split()
        if not (line[0].isdigit() or line[0] in ["X", "Y"]):
            continue
        region = {
            key: int(line[int(i)]) if line[int(i)].isdigit() else line[int(i)]
            for key, i in site_keys.items()
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
        cells=None,
        log_every=25000,
        pad_left=15,
        pad_right=20,
        site_keys=None,
        add_seq=False,
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
        self.add_seq = add_seq
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
            start=[],
            stop=[],
        )
        self.info["global"] = dict(
            n_umi=[],
            n_umi_corr=[],
            n_reads=[],
            n_reads_corr=[],
            nucl=[],
            motif_len=[],
            alleles=[],
            alleles_kept=[],
            chrom=[],
            start=[],
            stop=[],
        )
        self.sites_kept = set()
        if cells is not None:
            self.cells = set(
                map(lambda x: x.split()[0], open(cells, "r").read().splitlines())
            )
            self.has_cells = True
        else:
            self.cells = set()
            self.has_cells = False
        if site_keys is None:
            self.site_keys = OrderedDict(
                {v: i for i, v in enumerate(["chrom", "start", "stop"])}
            )
        else:
            self.site_keys = OrderedDict(site_keys)
        self.variant_f = open(self.pr_path.format("variants", "txt"), "w")
        self.variant_f.write(
            "chrom start cb a1 a2 n-umi-1 n-umi-2 n-reads-1 n-reads-2\n"
        )

    @property
    def pr_path(self):
        return os.path.join(
            self.logs_dir, "{}_" + str(self.name).replace(" ", "-") + ".{}"
        )

    def run(self):
        log_handler = start_process_logging(self.logs_dir)
        try:
            if self.has_cells:
                logging.info(f"The population has {len(self.cells)} cells.")
            for site in str_to_regions(self.sites, self.site_keys):
                self.iter += 1
                start = site["start"] - self.pad_left
                stop = site["stop"] + self.pad_right
                cell_dict = self.get_cell_families(site)
                cell_dict = get_umi_families(cell_dict)
                vars = []
                for cell_barcode, umi_families in cell_dict.items():
                    if not self.has_cells or (
                        self.has_cells and cell_barcode in self.cells
                    ):
                        self.n_umi += len(umi_families)
                        vars.append(
                            self.parse_families(site, cell_barcode, umi_families)
                        )
                indels = sum((var[1:3] for var in vars), [])
                if any([x != indels[0] for x in indels]):
                    self.variant_f.writelines(
                        [
                            " ".join(map(str, [site["chrom"], site["start"]] + var))
                            + "\n"
                            for var in vars
                        ]
                    )
                    site_list = [site[key] for key in self.site_keys]
                    if self.add_seq:
                        site_list.append(
                            next(iter(umi_families.values()))[
                                0
                            ].query_alignment_sequence
                        )
                    self.sites_kept.add(" ".join(map(str, site_list)) + "\n")
                    if not self.has_cells:
                        self.cells.update(cell_dict.keys())
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
            self.variant_f.close()
            rm_path = self.pr_path.format("rm", "txt")
            with open(rm_path, "w") as dst:
                dst.write("\n".join(self.rm_reads))
            site_path = self.pr_path.format("sites", "txt")
            with open(site_path, "w") as dst:
                names = list(self.site_keys.keys())
                if self.add_seq:
                    names.append("seq")
                dst.write(" ".join(names) + "\n")
                for site in sorted(self.sites_kept):
                    dst.write(site)
            bc_path = self.pr_path.format("cells", "txt")
            with open(bc_path, "w") as dst:
                dst.write("\n".join(self.cells))
            info_path = self.pr_path.format("info", "h5")
            logging.info(f"Saving info into temporary file {info_path}.")
            save_h5(info_path, self.info)

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
        for read in self.bam.fetch(str(site["chrom"]), max(start, 1), stop):
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

    def parse_families(self, site, cell_barcode, umi_families):

        indels_by_read = dict()
        indels_by_family = dict()
        n_by_indel = dict()
        n_reads = 0
        alleles_kept = []

        for UMI, family in umi_families.items():
            n_reads += len(family)
            indels = overlap_noise(family, site)
            uni, counts = np.unique(indels, return_counts=True)
            for u, c in zip(uni, counts):
                n_by_indel[u] = n_by_indel.get(u, 0) + c
            if len(uni) > 1:
                indels_by_read[UMI] = indels
                indels_by_family[UMI] = dict(zip(uni, counts))
            else:
                alleles_kept.append(uni[0])

        if len(n_by_indel) < 2:
            return [
                cell_barcode,
                u,
                u,
                len(umi_families),
                len(umi_families),
                n_by_indel[u],
                n_by_indel[u],
            ]

        # finding most probable indel
        indel_by_umi = dict()
        n_umi_corr = 0

        for umi, v in indels_by_family.items():
            c = list(v.values())
            n_umi_corr += 1
            indel_by_umi[umi] = max(v, key=lambda x: 10 * v[x] + n_by_indel[x])

        alleles_kept += list(indel_by_umi.values())
        uni, counts = np.unique(alleles_kept, return_counts=True)
        if len(uni) == 1:
            i2, i1 = 0, 0
        else:
            i2, i1 = np.argsort(counts)[-2:]

        self.n_corr += n_umi_corr
        n_reads = sum(list(n_by_indel.values()))
        n_reads_corr = 0

        nucls = np.unique(
            list(umi_families[UMI][0].query_alignment_sequence), return_counts=True
        )
        nmax = np.argmax(nucls[1])
        nucl = ord(nucls[0][nmax]) + nucls[1][nmax] / np.sum(nucls[1])

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
                " ".join((r.query_name, r.cigarstring))
                for r, b in zip(umi_families[umi], is_concordant)
                if not b
            ]
            n_reads_umi = len(indels_by_read[umi])
            bq = list(map(lambda r: np.mean(r.query_qualities), umi_families[umi]))

            umi_info = dict(
                cfreq=n_by_indel[x] / n_reads,
                freq=indels_by_family[umi][x] / n_reads_umi,
                count=indels_by_family[umi][x],
                alleles=len(indels_by_family[umi]),
                avg_qual=sum(bq) / n_reads_umi,
                qual=sum([x for x, b in zip(bq, is_concordant) if b])
                / sum(is_concordant),
                nucl=nucl,
                motif_len=site["motif_len"],
                chrom=site["chrom"]
                if isinstance(site["chrom"], int)
                else ord(site["chrom"]),
                start=site["start"],
                stop=site["stop"],
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
            chrom=site["chrom"]
            if isinstance(site["chrom"], int)
            else ord(site["chrom"]),
            alleles=len(n_by_indel),
            alleles_kept=len(alleles_kept),
            start=site["start"],
            stop=site["stop"],
        )
        for k, v in global_info.items():
            self.info["global"][k].append(v)

        return [
            cell_barcode,
            uni[i1],
            uni[i2],
            counts[i1],
            counts[i2],
            n_by_indel[uni[i1]],
            n_by_indel[uni[i2]],
        ]


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
