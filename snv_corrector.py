import os
import pysam
import threading
import logging
import queue

from concurrent.futures import ThreadPoolExecutor
from multiprocessing import Manager

import numpy as np
from time import time

from logger import start_thread_logging, stop_thread_logging, config_root_logger


def correct_base(reads, bases, save_every, log_every=None):  # queue, info, sem

    log_msg = "Step {}, position {}: {}"

    add_reads = set()
    rm_reads = set()
    umi_freq = []
    n_discordant = 0

    if all([b == bases[0] for b in bases]):
        add_reads.update(reads)
        families = []
    else:
        families = get_families(bases, reads)
        for v in families.values():
            uni, counts = np.unique(v["bases"], return_counts=True)
            if len(uni) > 1:
                n_discordant += 1
                i = np.argmax(counts)
                umi_freq.append(counts[i] / np.sum(counts))
                add_reads.update(
                    [r for r, b in zip(v["reads"], v["bases"]) if b == uni[i]]
                )
                rm_reads.update(
                    [r for r, b in zip(v["reads"], v["bases"]) if b != uni[i]]
                )
            else:
                add_reads.update(reads)

    """
    info["step"] += 1
    info["n_families"] += len(families)
    info["n_discordant"] += n_discordant
    queue["add_new"].update(add_reads)
    queue["rm_new"].update(rm_reads)
    if len(queue["add_new"]) >= save_every:
        sem.release()
    if info["step"]%log_every==0:
        logging.info(log_msg.format(
            self.step,
            pileup_column.reference_pos,
            info
        ))
    """


def get_families(bases, reads):
    families = dict()
    for base, read in zip(bases, reads):
        try:
            c, b = read.alignment.get_tag("CB"), read.alignment.get_tag("UB")
        except:
            continue
        if not (c, b) in families:
            families[(c, b)] = dict(bases=[], reads=[])
        families[(c, b)]["bases"].append(base)
        families[(c, b)]["reads"].append(read.alignment)
    return families


class SNVCorrector:
    def __init__(self, bam, logs_dir, chrom, log_every=50000, save_every=1000):
        self.bam = pysam.AlignmentFile(bam)
        self.chr = chrom

        self.logs_dir = logs_dir
        self.log_every = log_every
        self.save_every = save_every
        config_root_logger(logs_dir)

        manager = Manager()
        self.queue = manager.dict()
        self.queue.update(
            dict(add_new=set(), add_old=set(), rm_new=set(), rm_old=set())
        )
        self.info = manager.dict()
        self.info.update(
            dict(step=0, n_families=0, n_umis=0, n_discordant=0, umi_freq=[])
        )

        self.executor = ThreadPoolExecutor(max_workers=8)
        self.wthread = threading.Thread(
            target=self.save_bam, name="writer", args=(self.queue,), daemon=True
        )
        self.rqueue = queue.Queue(maxsize=2)
        self.rthread = threading.Thread(
            target=self.yield_reads, name="reader", args=(self.rqueue, 100000)
        )
        self.sem = threading.Semaphore()

    def __call__(self):

        # self.wthread.start()
        # self.rthread.start()
        self.bam_pileup = self.bam.pileup(self.chr)

        threading.current_thread().name = "main"
        log_handler = start_thread_logging(self.logs_dir)
        logging.info("Identifying discordant UMIs.")
        t0 = time()

        logging.info("Fetching pileups....")
        # LOAD : Uniquement 100 pileups pour pas niquer ta mémoire
        i = 0
        pileups = []
        bases = []
        """
        for x in self.bam.pileup(self.chr):
            pileups.append(x.pileups)
            bases.append(x.get_query_sequences())
            i+=1
            if i==5000000:
                break
        #pileups = [x for x in self.bam.pileup(self.chr)]
        logging.info(f"Done. Found {len(pileups)} pileups.")
        """

        logging.info(f"Parsing chromosome {self.chr} in bam...")
        i = 0
        # for pileup_column in self.bam.pileup(self.chr):
        while True:
            logging.info(i)
            # reads, bases = self.rqueue.get()
            reads, bases = self.yield_reads(100000)
            logging.info(i)
            self.executor.map(
                lambda x, y: correct_base(
                    x, y, self.save_every, self.log_every
                ),  # self.queue, self.info, self.sem, self.save_every, self.log_every),
                reads,
                bases,
                chunksize=2000,
            )
            i += 1
        # if list(ex)[0] is not None:
        #    print(list(ex))

        logging.info(
            "Done. {} discordant families out of {}".format(
                self.info["n_discordant"], self.info["n_families"]
            )
        )
        logging.info(f"Elapsed time: {time()-t0}s")
        stop_thread_logging(log_handler)
        self.bam.close()

    def yield_reads(self, chunksize=1):
        reads = []
        seqs = []
        i = 0
        while i < chunksize:
            pileup_column = next(self.bam_pileup)
            reads.append([r.alignment for r in pileup_column.pileups]),
            seqs.append(pileup_column.get_query_sequences())
            i += 1
        return reads, seqs

    def save_bam(self, queue):
        thread_log_handler = start_thread_logging(self.logs_dir)
        save_path = os.path.join(self.logs_dir, "corrected.bam")
        logging.info(f"Saving new bam into {save_path}.")
        bam_out = pysam.AlignmentFile(save_path, "w", template=self.bam)
        log_every = self.log_every // self.save_every if self.log_every else None
        i = 0
        while self.sem.acquire():
            reads = queue["add_old"].difference(queue["rm_old"].union(queue["rm_new"]))
            queue["add_old"] = queue["add_new"]
            queue["add_new"] = set()
            queue["rm_old"] = queue["rm_new"]
            queue["rm_new"] = set()

            t0 = time()
            for read in reads:
                bam_out.write(read)
            i += 1
            if log_every and i % log_every == 0:
                logging.info(f"Step {i}, wrote {len(reads)} new samples.")

        bam_out.close()
        stop_thread_logging(thread_log_handler)
