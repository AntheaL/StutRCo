import sys
import os
import argparse
import pysam
import logging

from snv_corrector import SNVCorrector

from logger import start_thread_logging, stop_thread_logging, config_root_logger

def main():
    parser = argparse.ArgumentParser(description="Script for PCR stutter correction.")
    parser.add_argument("--logs-dir", type=str, required=True)
    parser.add_argument("--tag", type=str, required=True)
    parser.add_argument(
        "--bam", type=str, help="Path to Input Single Cell BAM", required=True
    )
    parser.add_argument("--chr", type=str, required=True)

    args = parser.parse_args()

    logs_dir = os.path.join(args.logs_dir, args.tag)
    os.mkdir(logs_dir)

    print("Correcting errors.\n")
    snv_corrector = SNVCorrector(args.bam, chrom=args.chr, logs_dir=logs_dir)
    snv_corrector()

if __name__ == "__main__":
    main()
