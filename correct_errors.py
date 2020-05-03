import sys
import os
import argparse
import pysam
import multiprocessing as mp
import subprocess

from reader import FilestreamBED
from stutter_corrector import StutRCorrector


def line_count(filename):
    return int(subprocess.check_output('wc -l {}'.format(filename), shell=True).split()[0])

def multiprocess_args(bam_path, bed_path, n_cores):
    n_lines = line_count(bam_path)
    n_per_core = n_lines//n_cores
    print(f"Dividing file with {n_lines} lines in {n_cores} chuncks of size {n_per_core}.")
    return zip(
        [bam_path]*n_cores,
        [bed_path]*n_cores,
        range(1,n_lines,n_per_core),
        [n_per_core]*(n_cores-1)+[None],
        [50000] + (n_cores-1)*[None]
    )

def wrapper(args):
    bam = pysam.AlignmentFile(args[0])
    bed = FilestreamBED(file=args[1], start=args[2], n_iters=args[3])
    stutter_corrector = StutRCorrector(bed, bam, log_every = args[-1])
    stutter_corrector()
    return stutter_corrector.rm_reads, stutter_corrector.info

def parse_result(data):

    rm_reads = set()
    for (_rm_reads, _info) in result:
        intersect = rm_reads.intersection(_rm_reads)
        if intersect:
            print(f"Caution! {len(intersect)} reads out of {len(_rm_reads)} are already present in main set of size {len(rm_reads)}.")
        rm_reads.update(_rm_reads)
    info = concatenate_dicts([x[1] for x in data])
    return rm_reads, info


def main():
    parser = argparse.ArgumentParser(description="Script for PCR stutter correction.")
    parser.add_argument("--logs-dir", type=str, required=True)
    parser.add_argument("--tag", type=str, required=True)
    parser.add_argument(
        "--bam", type=str, help="Path to Input Single Cell BAM", required=True
    )
    parser.add_argument("--write-log-every", type=int, required=False, default=5000000)
    parser.add_argument("--bed", type=str, help="Path to Bed with STR regions", required=True)
    parser.add_argument("--multiprocess", action='store_true')

    args = parser.parse_args()

    logs_dir = os.path.join(args.logs_dir, args.tag)
    os.mkdir(logs_dir)


    info_path = os.path.join(logs_dir, "info.h5")
    print(f"Saving info in {info_path}")
    print("Correcting errors.\n")

    if args.multiprocess:

        n_cores = mp.cpu_count()
        print(f"Using {n_cores} cores.")
        mp_args = multiprocess_args(args.bam, args.bed, n_cores)

        print(f"Running multiprocessing.")
        pool = mp.Pool(n_cores)
        result  = pool.map(wrapper, mp_args, 1)

        print(f"Combining results from different processes.")
        rm_reads, info = parse_result(result)
        save_h5(info_path, info)

    else:
        bam = pysam.AlignmentFile(args.bam)
        bed = FilestreamBED(args.bed)
        stutter_corrector = StutRCorrector(bed, bam, pool=pool)
        stutter_corrector()
        stutter_corrector.save_info(info_path)
        rm_reads = stutter_corrector.rm_reads

    bam_path = os.path.join(logs_dir, "output.bam")
    bam_out = pysam.AlignmentFile(bam_path, "w", template=bam)
    print(f"\n Saving new bam in {bam_path}")
    n_parsed = 0
    n_wrote = 0
    log_every = args.write_log_every
    for read in bam.fetch():
        if read.query_name not in rm_reads:
            bam_out.write(read)
            n_wrote += 1
        n_parsed += 1
        if log_every and n_parsed % log_every == 0:
            print(f"Wrote {n_wrote} reads, skipped {n_parsed-n_wrote}.")

    bam.close()
    bam_out.close()

if __name__ == "__main__":



    main()
