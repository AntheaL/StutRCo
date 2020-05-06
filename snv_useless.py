import os
import pysam
import threading
import queue
import logging

from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Manager

import numpy as np
from time import time

from logger import start_thread_logging, stop_thread_logging, config_root_logger

def correct_base(pileup_column, queue, info, sem, save_every, log_every=None):


    while True:
        continue
    print("FUCK FUCK FUCK")

    info["step"] += 1

    log_msg = "Step {}: position {}."
    bases = pileup_column.get_query_sequences()
    if all([b==bases[0] for b in bases]):
        queue["add_new"].update([r.alignment for r in pileup_column.pileups])
    families = get_families(bases, pileup_column.pileups)
    info["n_families"] += len(families)

    for v in families.values():
        uni, counts = np.unique(v['bases'], return_counts=True)

        if len(uni)>1:
            info["n_discordant"] += 1
            i = np.argmax(counts)
            info['umi_freq'].append(counts[i]/np.sum(counts))

            queue["add_new"].update([
                r.alignment for r,b in zip(v['reads'], v['bases']) if b==uni[i]
            ])
            queue["rm_new"].update([
                r.alignment for r,b in zip(v['reads'], v['bases']) if b!=uni[i]
            ])

        else:
            queue["add_new"].update([
                r.alignment for r,b in zip(v['reads'], v['bases']) if b==uni[i]
            ])

    if len(queue["add_new"]) >= save_every:
        sem.release()

    if info["step"]%log_every==0:
        logging.info(log_msg.format(
            self.step,
            pileup_column.reference_pos,
        ))


def get_families(bases, reads):
    families = dict()
    for base, read in zip(bases, reads):
        try:
            c, b = read.alignment.get_tag('CB'), read.alignment.get_tag('UB')
        except:
            continue
        if not (c,b) in families:
            families[(c,b)] = dict(bases=[], reads=[])
        families[(c,b)]['bases'].append(base)
        families[(c,b)]['reads'].append(read)
    return families


class SNVCorrector:

    def __init__(self, bam, logs_dir, queue, sem, log_every=50000, save_every=1000):
        self.bam = bam
        self.logs_dir = logs_dir
        self.log_every = log_every
        self.save_every = save_every
        self.queue = queue
        self.sem = sem
        self.wthread = threading.Thread(
            target=self.save_bam, name='writer', daemon=True, args=(self.queue,)
        )
        config_root_logger(logs_dir)

    def __call__(self):
        self.wthread.start()

    def save_bam(self, queue):
        thread_log_handler = start_thread_logging(self.logs_dir)
        save_path = os.path.join(self.logs_dir, "corrected.bam")
        logging.info(f"Saving new bam into {save_path}.")
        bam_out = pysam.AlignmentFile(save_path, "w", template=self.bam)
        log_every = self.log_every//self.save_every if self.log_every else None
        i = 0
        while self.sem.acquire():
            reads = queue["add_old"].difference(
                queue["rm_old"].union(queue["rm_new"])
            )
            queue["add_old"] = queue["add_new"]
            queue["add_new"] = set()
            queue["rm_old"] = queue["rm_new"]
            queue["rm_new"] =  set()

            t0 = time()
            for read in reads:
                bam_out.write(read)
            i += 1
            if log_every and i%log_every==0:
                logging.info(
                    f"Step {i}, wrote {len(reads)} new samples."
                )

        bam_out.close()
        stop_thread_logging(thread_log_handler)
