import numpy as np
from read_util import (
    covers_str,
    expanded_seq,
    parse_changes,
    parse_insertions,
    fetch_indels,
    overlap_noise
)

class StutRCorrector:

    def __init__(self, vcf, bam, log_every=1000, pad_left=15, pad_right=20):
        self.vcf = vcf
        self.bam = bam
        self.log_every = log_every
        self.pad_left = pad_left
        self.pad_right = pad_right

    def __call__(self):

        rm_reads = set()
        n_corr = 0
        n_tot = 0
        log_i = 0

        for site in self.vcf:

            log_i += 1
            if log_i and log_i % self.log_every == 0:
                print(
                    f"Site {log_i}, removed {len(rm_reads)} reads across {n_corr} UMI families out of {n_tot}."
                )

            cell_dict = self.get_site_umi_families(site)

            for cell_barcode, UMI_families in cell_dict.items():

                indels_by_read = dict()
                indels_by_family = dict()
                n_by_indel = dict()

                for UMI, family in UMI_families.items():
                    indels = overlap_noise(family, site)
                    uni, counts = np.unique(indels, return_counts=True)
                    for u, c in zip(uni, counts):
                        n_by_indel[u] = n_by_indel.get(u, 0) + c
                    if len(uni) > 1:
                        indels_by_read[UMI] = indels
                        indels_by_family[UMI] = dict(zip(uni, counts))
                    n_tot += 1

                # finding most probable indel
                indel_by_umi = dict()
                for umi, v in indels_by_family.items():
                    c = list(v.values())
                    if len(c) > 1:
                        indel_by_umi[umi] = max(v, key=lambda x: 10 * v[x] + n_by_indel[x])

                # assigning a matching read to each UMI and extending site
                start, stop = np.inf, 0

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

                    rm_reads.update([r.query_name for j, r in enumerate(family) if j != i])
                    n_corr += 1

        return rm_reads


    def get_site_umi_families(self, site):
        """
        Generate cell barcode clusters in a BAM file or within a region string.
        Reads are added to read_dict until a pair is found.
        """
        cell_dict = {}
        start = site["pos"] - self.pad_left
        stop = site["pos"] + len(site["ref"]) + self.pad_right
        for read in self.bam.fetch(site["chrom"], max(start, 1), stop):
            if read.is_secondary or read.is_supplementary:
                continue
            if not covers_str(
                read, dict(start=site["pos"], stop=site["pos"] + len(site["ref"]))
            ):
                continue
            try:
                cell_barcode = read.get_tag("CB")
            except:
                continue

            if cell_barcode in cell_dict:
                cell_dict[cell_barcode].append(read)
            else:
                cell_dict[cell_barcode] = [read]

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
