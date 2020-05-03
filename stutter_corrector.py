import h5py
import numpy as np
from read_util import (
    covers_str,
    expanded_seq,
    parse_changes,
    parse_insertions,
    fetch_indels,
    overlap_noise,
)


class StutRCorrector:
    def __init__(self, bed, bam, log_every=50000, pad_left=15, pad_right=20):
        self.bed = bed
        self.bam = bam
        self.log_every = log_every
        self.pad_left = pad_left
        self.pad_right = pad_right
        self.rm_reads = set()
        self.n_umi = 0
        self.n_corr = 0
        self.info = dict(cell_freq=[], umi_freq=[], umi_count=[], chrom=[], pos=[])
        self.info["global"] = dict(n_umi=[], n_reads=[], chrom=[], pos=[])


    def __call__(self):
        self.iter = 0
        for site in self.bed:
            self.process_site(site)

    def process_site(self, site):

        self.info["global"]["n_umi"].append(0)
        self.info["global"]["n_reads"].append(0)
        self.info["global"]["chrom"].append(int(site["chrom"]) if site["chrom"].isdigit() else ord(site["chrom"]))
        self.info["global"]["pos"].append(site["start"])

        self.iter += 1
        if self.log_every and self.iter % self.log_every == 0:
            print(
                f"Site {self.iter}, removed {len(self.rm_reads)} reads across {self.n_corr} UMI families out of {self.n_umi}."
            )
        cell_dict = self.get_site_umi_families(site)
        for cell_barcode, umi_families in cell_dict.items():
            self.parse_families(site, umi_families)
            self.info["global"]["n_umi"][-1] += len(umi_families)

    def save_info(self, dst_path):
        with h5py.File(dst_path, "w") as dst:
            for key, value in self.info.items():
                if isinstance(value, dict):
                    grp = dst.create_group(name=key)
                    for k,v in value.items():
                        grp.create_dataset(name=k, data=v)
                else:
                    dst.create_dataset(name=key, data=np.array(value))

    def parse_families(self, site, umi_families):

        indels_by_read = dict()
        indels_by_family = dict()
        n_by_indel = dict()

        for UMI, family in umi_families.items():
            self.info["global"]["n_reads"][-1] += len(family)
            indels = overlap_noise(family, site)
            uni, counts = np.unique(indels, return_counts=True)
            for u, c in zip(uni, counts):
                n_by_indel[u] = n_by_indel.get(u, 0) + c
            if len(uni) > 1:
                indels_by_read[UMI] = indels
                indels_by_family[UMI] = dict(zip(uni, counts))

        self.n_umi += len(umi_families)

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

            self.rm_reads.update([r.query_name for j, r in enumerate(family) if j != i])
            self.n_corr += 1
            self.info["cell_freq"].append(n_by_indel[x] / n_reads)
            self.info["umi_freq"].append(
                indels_by_family[umi][x] / len(indels_by_read[umi])
            )
            self.info["umi_count"].append(indels_by_family[umi][x])
            self.info["chrom"].append(int(site["chrom"]) if site["chrom"].isdigit() else ord(site["chrom"]))
            self.info["pos"].append(site["start"])

    def get_site_umi_families(self, site):
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

                if cell_barcode in cell_dict:
                    cell_dict[cell_barcode].append(read)
                else:
                    cell_dict[cell_barcode] = [read]
        except ValueError:
            raise ValueError(f"Invalid site {site}")
        del_length_one(cell_dict)

        clean_dict = {}
        for cell_barcode, cell_reads in cell_dict.items():

            UMI_families = {}
            for i, read in enumerate(cell_reads):
                try:
                    UMI_barcode = read.get_tag("UB")
                except KeyError:
                    print(f"Skipping read {i} which does not contain 'UB' tag.")
                    continue
                if UMI_barcode in UMI_families:
                    UMI_families[UMI_barcode].append(read)
                else:
                    UMI_families[UMI_barcode] = [read]
            del_length_one(UMI_families)

            if len(UMI_families) >= 1:
                clean_dict[cell_barcode] = UMI_families

        return clean_dict


def del_length_one(dict_):
    for key in [k for k in dict_.keys()]:
        if len(dict_[key]) <= 1:
            del dict_[key]
