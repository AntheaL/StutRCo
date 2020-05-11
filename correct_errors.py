import sys
import glob
import os
import argparse
import pysam
import logging
import yaml
import queue
import random
import fileinput
import multiprocessing

from functools import reduce
from operator import iadd
from reader import FilestreamBED
from h5_util import concatenate_files, save_h5
from stutter_corrector import StutRCorrector
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Manager, Process, cpu_count
from logger import start_process_logging, stop_process_logging, config_root_logger

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script for PCR stutter correction.")
    parser.add_argument("--config", type=str, required=True)
    parser.add_argument("--tag", type=str, required=True)
    parser.add_argument(
        "--n-workers", type=int, help="Number of workers", required=False, default=2
    )

    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = yaml.load(f, Loader=yaml.BaseLoader)
        multiprocessing.current_process().name = "main"

    # n_workers = cpu_count()
    n_workers = args.n_workers
    logging.info(f"Using {n_workers} workers.")

    logs_dir = os.path.join(config["logs_dir"], args.tag)
    os.mkdir(logs_dir)
    config_root_logger(logs_dir)

    with open(os.path.join(logs_dir, "config.yml"), "w") as f:
        yaml.dump(config, f)

    m = Manager()
    queue = m.Queue(maxsize=cpu_count() - 1)

    logging.info(f"Loading sites...")
    sites = [site for i, site in enumerate(open(config["bed"], "r"))]
    random.shuffle(sites)

    n_sites = len(sites)
    bs = n_sites // n_workers
    idx = [i * bs for i in range(n_workers)] + [None]

    logging.info(f"Done. Each worker received {bs} sites.")
    logging.info(f"Running processes.")

    jobs = []
    pipe_list = []
    for i in range(n_workers):
        recv_end, send_end = multiprocessing.Pipe(False)
        p = StutRCorrector(
            sites[idx[i] : idx[i + 1]],
            bam=config["bam"],
            logs_dir=logs_dir,
            send_end=send_end,
            name="process-" + str(i),
            **config["corrector"],
        )
        jobs.append(p)
        pipe_list.append(recv_end)

    sites.clear()

    for p in jobs:
        p.start()

    rm_reads = [x.recv() for x in pipe_list]
    for proc in jobs:
        proc.join()

    info_path = os.path.join(logs_dir, "info.h5")
    info_files = glob.glob(os.path.join(logs_dir, "*.h5"))
    logging.info(f"Saving info in {info_path}...")
    info = concatenate_files(info_path, info_files,)
    logging.info("Removing temporary files...")
    for f in info_files:
        os.remove(f)

    log_path = os.path.join(logs_dir, "process.log")
    log_files = glob.glob(os.path.join(logs_dir, "*-*.log"))
    with open(log_path, "w") as file:
        input_lines = fileinput.input(log_files)
        file.writelines(input_lines)
    for f in log_files:
        os.remove(f)

    bam_path = os.path.join(logs_dir, "output.bam")
    bam_in = pysam.AlignmentFile(config["bam"])
    bam_out = pysam.AlignmentFile(bam_path, "w", template=bam_in)
    logging.info(f"Saving new bam in {bam_path}")
    n_parsed = 0
    n_wrote = 0
    log_every = int(config["write_log_every"])
    rm_reads = set(reduce(iadd, rm_reads))
    for read in bam_in.fetch():
        if read.query_name not in rm_reads:
            bam_out.write(read)
            n_wrote += 1
        n_parsed += 1
        if log_every and n_parsed % log_every == 0:
            logging.info(f"Wrote {n_wrote} reads, skipped {n_parsed-n_wrote}.")

    bam_in.close()
    bam_out.close()
