import sys
import os
import argparse
import pysam

from reader import FilestreamVCF
from util import error_correction


def main():
    parser = argparse.ArgumentParser(description="Script for PCR stutter correction.")
    parser.add_argument("--dst-path", type=str, required=True)
    parser.add_argument(
        "--bam", type=str, help="Path to Input Single Cell BAM", required=True
    )
    parser.add_argument("--write-log-every", type=int, required=False, default=5000000)
    parser.add_argument("--vcf", type=str, help="Path to Variant VCF", required=True)

    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam)
    vcf = FilestreamVCF(args.vcf)

    print("Correcting errors.\n")
    rm_reads = error_correction(vcf, bam)

    bam_out = pysam.AlignmentFile(args.dst_path, "w", template=bam)
    print(f"\n Saving new bam in {args.dst_path}")
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
