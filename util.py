import numpy as np


def error_correction(vcf, bam, log_every=1000):

    rm_reads = set()
    n_corr = 0
    n_tot = 0
    log_i = 0

    for site in vcf:

        log_i += 1
        if log_i and log_i % log_every == 0:
            print(
                f"Site {log_i}, removed {len(rm_reads)} reads across {n_corr} UMI families out of {n_tot}."
            )

        cell_dict = get_site_umi_families(bam, site, 15, 20)

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


def get_site_umi_families(bam, site, pad_left, pad_right):
    """
	Generate cell barcode clusters in a BAM file or within a region string.
	Reads are added to read_dict until a pair is found.
	"""
    cell_dict = {}
    start = site["pos"] - pad_left
    stop = site["pos"] + len(site["ref"]) + pad_right
    for read in bam.fetch(site["chrom"], max(start, 1), stop):
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


def fetch_reads(bam, site, str_info, filter_dups=False, match_alt=False):
    for read in bam.fetch(site["chrom"], start=site["pos"] - 20, stop=site["pos"] + 20):
        if (
            (read.cigartuples is None)
            or read.is_unmapped
            or read.is_secondary
            or read.is_supplementary
            or (read.is_duplicate and filter_dups)
        ):
            continue
        start, end = read.reference_start, read.reference_end
        if covers_str(read, str_info):
            yield read


def get_real_start(read):
    start = read.reference_start
    cigar = read.cigartuples
    if cigar[0][0] == 4:
        start -= cigar[0][1]
    return start


def get_real_end(read):
    end = read.reference_end
    cigar = read.cigartuples
    if cigar[-1][0] == 4:
        end += cigar[-1][1]
    return end


def covers_str(read, str_info):
    return (
        read.reference_end > str_info["stop"]
        and read.reference_start < str_info["start"]
    )


def del_length_one(dict_):
    for key in [k for k in dict_.keys()]:
        if len(dict_[key]) <= 1:
            del dict_[key]


def collect_reference_STRs(bed):
    reference_STRs = {}
    for region in tqdm(
        bed,
        desc="Collecting STR Regions of Interest",
        total=len(bed),
        unit=" STR Regions",
    ):
        reference_STRs[region["STR_id"]] = region
    return reference_STRs


def expanded_seq(read, region):
    cigar_tup_dict = {i: s for i, s in enumerate(["M", "I", "D", "N", "S"])}
    start, end = read.reference_start, read.reference_end

    sequence = read.query_alignment_sequence
    cigar = read.cigartuples

    if cigar[0][0] == 4:
        cigar = cigar[1:]
    if cigar[-1][0] == 4:
        cigar = cigar[:-1]

    j = 0
    i = start
    seq_str = ""
    cig_str = ""

    k = 0

    for operator, length in cigar:
        operator = cigar_tup_dict[operator]
        k += length
        if operator == "D":
            seq_str += "X" * length
            i += length

        elif operator == "I":
            seq_str += sequence[j : j + length].lower()
            j += length
        else:
            seq_str += sequence[j : j + length]
            j += length
            i += length
        cig_str += operator * length

    front_clip, back_clip = (region[0] - start), (region[1] - end)
    if front_clip > 0:
        seq_str = seq_str[front_clip:]
        cig_str = cig_str[front_clip:]
    elif front_clip < 0:
        seq_str = abs(front_clip) * "_" + seq_str
        cig_str = abs(front_clip) * "_" + cig_str
    if back_clip < 0:
        seq_str = seq_str[:back_clip]
        cig_str = cig_str[:back_clip]
    elif back_clip > 0:
        seq_str += back_clip * "_"
        cig_str += back_clip * "_"

    return seq_str, cig_str


def parse_changes(cigar, C):
    changes = 0
    for c in cigar:
        if c == C:
            changes += 1
    return changes


def parse_insertions(f_cig, m_cig, e_cig):
    for c in f_cig[::-1]:
        if c != "I":
            break
        else:
            m_cig = c + m_cig
    for c in e_cig:
        if c != "I":
            break
        m_cig += c
    return parse_changes(m_cig, "I")


def fetch_indels(read, region):
    ref, cig = expanded_seq(read, region)
    deletions = parse_changes(cig[9:-11], "D")
    insertions = parse_insertions(cig[:9], cig[9:-11], cig[-11:])
    return deletions, insertions


def overlap_noise(family, site, unique=False):

    start, stop = site["pos"], site["pos"] + len(site["ref"])
    region = start - 10, stop + 10

    deletions, insertions = [], []
    for read in family:
        d, i = fetch_indels(read, region)
        deletions.append(d)
        insertions.append(i)

    if unique:
        uni_d, counts_d = np.unique(deletions, return_counts=True)
        uni_i, counts_i = np.unique(insertions, return_counts=True)
        return dict(
            uni=np.concatenate([uni_i, -uni_d]),
            counts=np.concatenate([counts_i, counts_d]),
        )

    return [i - d for i, d in zip(insertions, deletions)]
