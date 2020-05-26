import sys
import glob
import os
import argparse
import pysam
import psutil
import logging
import yaml
import queue
import random
import fileinput
import multiprocessing

import numpy as np
from functools import reduce
from operator import iadd
from reader import FilestreamBED
from h5_util import concatenate_files, save_h5
from stutter_corrector import StutRCorrector
from multiprocessing import Manager, Process, cpu_count
from logger import start_process_logging, stop_process_logging, config_root_logger


def get_paths(logs_dir, itype, ext, sort=True):
    files = glob.glob(os.path.join(logs_dir, f"{itype}*.{ext}"))
    if sort:
        return sorted(files)
    return files


def merge_files(logs_dir, itype, ext="txt", sort=True, to_set=False, has_header=False):
    out_path = os.path.join(logs_dir, f"{itype}.{ext}")
    in_paths = get_paths(logs_dir, itype, ext, sort=sort)
    logging.info(f"Saving {itype} in {out_path}...")
    if to_set:
        lines = set()
        for path in in_paths:
            lines.update(open(path, "r").read().splitlines())
    else:
        lines = map(lambda x: x.rstrip("\n"), fileinput.input(in_paths))
    with open(out_path, "w") as dst:
        if has_header:
            header = next(iter(lines))
            dst.write(header + "\n")
            dst.write("\n".join(filter(lambda x: x != header, lines)))
        else:
            dst.write("\n".join(lines))
    logging.info("Removing temporary files...")
    for f in in_paths:
        os.remove(f)


def calls_to_array(logs_dir):

    src_path = os.path.join(logs_dir, "{}.txt")
    barcodes = {
        v: i
        for i, v in enumerate(open(src_path.format("cells"), "r").read().splitlines())
    }
    # sites = [s.split()[:2] for s in open(src_path.format('sites'), "r")][1:]
    # sites = {tuple(v): i for i, v in enumerate(sorted(sites))}
    n_bc = len(barcodes)
    logging.info(f"Number variant cells: {n_bc}.")

    dst_path = os.path.join(logs_dir, "variants_{}.csv")
    logging.info(f"Saving into {dst_path.format(1)}...")
    dst_1 = open(dst_path.format(1), "w")
    logging.info(f"saving into {dst_path.format(2)}...")
    dst_2 = open(dst_path.format(2), "w")

    logging.info("Scanning variants...")
    with open(src_path.format("variants"), "r") as src:
        next(src)
        s = None
        var_1, var_2 = [66] * n_bc, [66] * n_bc
        for l in src:
            l = l.strip("\n").split()
            i = barcodes[l[2]]
            if not s:
                s = (l[0], l[1])
            elif (l[0], l[1]) != s:
                dst_1.write(" ".join(map(str, var_1)) + "\n")
                dst_2.write(" ".join(map(str, var_2)) + "\n")
                var_1, var_2 = [66] * n_bc, [66] * n_bc
                s = (l[0], l[1])
            var_1[i] = int(l[3])
            var_2[i] = int(l[4])
    if any(x != 66 for x in var_1):
        dst_1.write(" ".join(map(str, var_1)) + "\n")
        dst_2.write(" ".join(map(str, var_2)) + "\n")
    dst_1.close()
    dst_2.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script for PCR stutter correction.")
    parser.add_argument("--config", type=str, required=True)
    parser.add_argument("--tag", type=str, required=True)
    parser.add_argument("--bam", type=str, required=False)
    parser.add_argument("--bed", type=str, required=False)
    parser.add_argument("--cells", type=str, required=False)
    parser.add_argument("--n-workers", type=int, required=False)

    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = yaml.load(f, Loader=yaml.BaseLoader)
        multiprocessing.current_process().name = "main"

    logs_dir = os.path.join(config["logs_dir"], args.tag)
    os.mkdir(logs_dir)
    config_root_logger(logs_dir)

    if args.n_workers:
        n_workers = args.n_workers
    else:
        cpu_affinity = psutil.Process().cpu_affinity()
        n_workers = len(cpu_affinity)
    logging.info(f"Using {n_workers} workers.")

    if args.bam:
        config["bam"] = args.bam
    assert "bam" in config, "Input bam file was not provided!"
    logging.info(f"Using bam from {config['bam']}")

    if args.bed:
        config["bed"] = args.bed
    assert "bed" in config, "Input bed file was not provided!"
    logging.info(f"Using bed from {config['bed']}")

    if args.cells:
        config["cells"] = args.cells

    with open(os.path.join(logs_dir, "config.yml"), "w") as f:
        yaml.dump(config, f)

    m = Manager()
    queue = m.Queue(maxsize=cpu_count() - 1)

    logging.info(f"Loading sites...")
    sites = [site for i, site in enumerate(open(config["bed"], "r"))]
    # random.shuffle(sites)

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
            cells=config.get("cells", None),
            logs_dir=logs_dir,
            # send_end=send_end,
            name="process-" + str(i).zfill(2),
            **config["corrector"],
        )
        jobs.append(p)
        pipe_list.append(recv_end)

    sites.clear()

    for p in jobs:
        p.start()

    # rm_reads = [x.recv() for x in pipe_list]

    are_alive = [True] * n_workers
    while any(are_alive):
        for i, proc in enumerate(jobs):
            if are_alive[i]:
                proc.join(timeout=0.1)
                if not proc.is_alive():
                    are_alive[i] = False
                    logging.info(f"Process {i} is not alive anymore!")

    info_path = os.path.join(logs_dir, "info.h5")
    info_files = get_paths(logs_dir, "info", "h5")
    logging.info(f"Saving info in {info_path}...")
    info = concatenate_files(info_path, info_files,)
    logging.info("Removing temporary files...")
    for f in info_files:
        os.remove(f)

    merge_files(logs_dir, "sites", has_header=True)
    merge_files(logs_dir, "cells", sort=False, to_set=True)
    merge_files(logs_dir, "variants", sort=False, has_header=True)
    merge_files(logs_dir, "process", "log")
    merge_files(logs_dir, "rm", sort=False)

    logging.info("Converting calls into arrays...")
    calls_to_array(logs_dir)

    rm_path = os.path.join(logs_dir, "rm.txt")
    rm_reads = set(open(rm_path).read().splitlines())
    head, ext = os.path.splitext(config["bam"])
    bam_name = os.path.basename(head) + "_out" + ext
    bam_path = os.path.join(logs_dir, bam_name)
    bam_in = pysam.AlignmentFile(config["bam"])
    bam_out = pysam.AlignmentFile(bam_path, "w", template=bam_in)
    logging.info(f"Saving new bam in {bam_path}")
    n_parsed = 0
    n_wrote = 0
    log_every = int(config["write_log_every"])
    for read in bam_in.fetch():
        if " ".join((read.query_name, read.cigarstring)) not in rm_reads:
            bam_out.write(read)
            n_wrote += 1
        n_parsed += 1
        if log_every and n_parsed % log_every == 0:
            logging.info(f"Wrote {n_wrote} reads, skipped {n_parsed-n_wrote}.")

    bam_in.close()
    bam_out.close()
