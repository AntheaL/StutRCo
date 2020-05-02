import sys
import os
import argparse
import pysam

from reader import FilestreamBED
from stutter_corrector import StutRCorrector


def main():
    parser = argparse.ArgumentParser(description="Script for PCR stutter correction.")
    parser.add_argument("--logs-dir", type=str, required=True)
    parser.add_argument("--tag", type=str, required=True)
    parser.add_argument(
        "--bam", type=str, help="Path to Input Single Cell BAM", required=True
    )
    parser.add_argument("--write-log-every", type=int, required=False, default=5000000)
    parser.add_argument("--bed", type=str, help="Path to Bed with STR regions", required=True)

    args = parser.parse_args()

    logs_dir = os.path.join(args.logs_dir, args.tag)
    os.mkdir(logs_dir)

    bam = pysam.AlignmentFile(args.bam)
    bed = FilestreamBED(args.bed)

    print("Correcting errors.\n")
    stutter_corrector = StutRCorrector(bed, bam)
    stutter_corrector()

    info_path = os.path.join(logs_dir, "info.h5")
    print(f"Saving info in {info_path}")
    stutter_corrector.save_info(info_path)

    bam_path = os.path.join(logs_dir, "output.bam")
    bam_out = pysam.AlignmentFile(bam_path, "w", template=bam)
    print(f"\n Saving new bam in {bam_path}")
    n_parsed = 0
    n_wrote = 0
    log_every = args.write_log_every
    for read in bam.fetch():
        if read.query_name not in stutter_corrector.rm_reads:
            bam_out.write(read)
            n_wrote += 1
        n_parsed += 1
        if log_every and n_parsed % log_every == 0:
            print(f"Wrote {n_wrote} reads, skipped {n_parsed-n_wrote}.")

    bam.close()
    bam_out.close()

if __name__ == "__main__":
    main()
