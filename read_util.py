import numpy as np


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


def fetch_snvs(read, region):
    ref, cig = expanded_seq(read, region)
    return [i for i,c in enumerate(cigar) if c=='X']


def fetch_indels(read, region):
    ref, cig = expanded_seq(read, region)
    deletions = parse_changes(cig[9:-11], "D")
    insertions = parse_insertions(cig[:9], cig[9:-11], cig[-11:])
    return deletions, insertions


def overlap_snvs(famly, site):
    start, stop = site["pos"], site["pos"] + len(site["ref"])
    region = start - 10, stop + 10
    snvs = dict()
    for read in family:
        positions = fetch_snvs(read, region)
        for i in positions:
            snvs[i] = snvs.get(i, 0) + 1


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
