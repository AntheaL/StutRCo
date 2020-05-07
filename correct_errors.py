import sys
import os
import argparse
import pysam
import logging
import queue
import random
import threading

from functools import reduce
from operator import iadd
from reader import FilestreamBED
from stutter_corrector import StutRCorrector
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager, Process, cpu_count
from logger import start_process_logging, stop_process_logging, config_root_logger

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script for PCR stutter correction.")
    parser.add_argument("--logs-dir", type=str, required=True)
    parser.add_argument("--tag", type=str, required=True)
    parser.add_argument(
        "--bam", type=str, help="Path to Input Single Cell BAM", required=True
    )
    parser.add_argument("--write-log-every", type=int, required=False, default=5000000)
    parser.add_argument("--bed", type=str, help="Path to Bed with STR regions", required=True)

    args = parser.parse_args()

    multiprocessing.current_process().name = "main"

    logs_dir = os.path.join(args.logs_dir, args.tag)
    os.mkdir(logs_dir)
    config_root_logger(logs_dir)

    #n_workers = cpu_count()
    n_workers = 16
    logging.info(f"Using {n_workers} cores.")

    m = Manager()
    queue = m.Queue(maxsize=cpu_count()-1)

    logging.info(f"Loading sites...")
    sites = [site for i, site in enumerate(open(args.bed, 'r')) if i<800000]
    random.shuffle(sites)

    n_sites = len(sites)
    bs = n_sites//n_workers
    idx = [i*bs for i in range(n_workers)]+[None]

    logging.info(f"Done. Each worker received {bs} sites.")

    processes = [
        StutRCorrector(
            sites[idx[i] : idx[i+1]],
            bam=args.bam,
            logs_dir=logs_dir,
            name="Process "+ str(i)
        )
        for i in range(n_workers)
    ]

    logging.info(f"Running processes.")
    [p.start() for p in processes]
    [p.join() for p in processes]

    """

    info_path = os.path.join(logs_dir, "info.h5")
    print(f"Saving info in {info_path}")
    stutrco.save_info(info_path)

    bam_path = os.path.join(logs_dir, "output.bam")
    bam_out = pysam.AlignmentFile(bam_path, "w", template=bam)
    print(f"\n Saving new bam in {bam_path}")
    n_parsed = 0
    n_wrote = 0
    log_every = args.write_log_every
    for read in bam.fetch():
        if read.query_name not in stutrco.rm_reads:
            bam_out.write(read)
            n_wrote += 1
        n_parsed += 1
        if log_every and n_parsed % log_every == 0:
            print(f"Wrote {n_wrote} reads, skipped {n_parsed-n_wrote}.")

    bam.close()
    bam_out.close()
    """
