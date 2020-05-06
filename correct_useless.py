import sys
import os
import argparse
import pysam
import logging
import threading

from time import time
from concurrent.futures import ThreadPoolExecutor
from snv_useless import SNVCorrector, correct_base
from multiprocessing import Manager

from logger import start_thread_logging, stop_thread_logging, config_root_logger


if __name__ == "__main__":

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
    bam = pysam.AlignmentFile(args.bam)
    executor = ThreadPoolExecutor(max_workers=8)
    manager = Manager()
    queue = manager.dict()
    queue.update(dict(
        add_new=set(), add_old=set(), rm_new=set(), rm_old=set()
    ))
    info = manager.dict()
    info.update(dict(step=0, n_families=0, n_umis=0, n_discordant=0, umi_freq=[]))
    sem = threading.Semaphore()

    snv_corrector = SNVCorrector(bam, logs_dir=logs_dir, queue=queue, sem=sem)

    threading.current_thread().name = "main"
    log_handler = start_thread_logging(logs_dir)
    logging.info("Identifying discordant UMIs.")
    t0 = time()

    logging.info("Fetching pileups....")
    # LOAD : Uniquement 100 pileups pour pas niquer ta mémoire
    i = 0
    pileups = []
    for x in bam.pileup(args.chr):
        pileups.append(x)
        i+=1
        if i==100:
            break
    i#pileups = [x for x in bam.pileup(args.chr)]
    logging.info(f"Done. Found {len(pileups)} pileups.")

    logging.info(f"Parsing chromosome {args.chr} in bam...")

    snv_corrector()
    for pileup in pileups:
        executor.submit(
            correct_base,
            (pileup, queue, info, sem, snv_corrector.save_every, snv_corrector.log_every)
        ).result()

    logging.info("Done. {} discordant families out of {}".format(
        info["n_discordant"], info["n_families"]
    ))
    logging.info(f"Elapsed time: {time()-t0}s")
    stop_thread_logging(log_handler)
    snv_corrector.bam.close()
