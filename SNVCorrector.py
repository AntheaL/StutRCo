import os
import pysam
import threading
import queue
import logging
import numpy as np
from time import time
from logger import start_thread_logging, stop_thread_logging, config_root_logger


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

    def __init__(self, bam, logs_dir, chrom, log_every=50000, save_every=1000):
        self.bam = pysam.AlignmentFile(bam)
        self.logs_dir = logs_dir
        self.log_every = log_every
        self.save_every = save_every
        self.chr = chrom
        self.rm_reads = set()
        config_root_logger(logs_dir)
        self.thread = threading.Thread(
            target=self.save_bam,
            name='Thread'
        )
        self.queue = queue.Queue()
        self.sem = threading.Semaphore()

    def __call__(self):

        self.info = dict(umi_freq=[])
        self.n_families = 0
        self.n_discordant = 0
        self.step = 0
        self.thread.start()
        self.parse_bam()
        self.bam.close()

    def save_bam(self):
        thread_log_handler = start_thread_logging(self.logs_dir)
        save_path = os.path.join(logs_dir, "bam_new.bam")
        logging.info(f"Saving new bam into {save_path}.")
        bam_out = pysam.AlignmentFile(save_path, "w", template=self.bam)
        log_every = self.log_every//self.save_every if self.log_every else None
        i = 0
        while self.sem.acquire():
            reads = self.queue.get()
            t0 = time()
            for read in reads:
                bam_out.write(read)
            i += 1
            if log_every and i%self.log_every==0:
                logging.info(
                    f"Step {i}, wrote {len(reads)} new samples."
                )
        bam_out.close()

        stop_thread_logging(thread_log_handler)

    def parse_bam(self):
        threading.current_thread().name = "main"
        log_handler = start_thread_logging(self.logs_dir)
        logging.info("Identifying discordant UMIs.")

        self.add_new = set()
        self.add_old = set()
        self.rm_new = set()
        self.rm_old = set()
        t0 = time()

        for pileup_column in self.bam.pileup(self.chr):
            self.correct_base(pileup_column)
            if self.step==200000:
                break

        logging.info(f"Done. {self.n_discordant} discordant families out of {self.n_families}")
        logging.info(f"Elapsed time: {time()-t0}s")
        stop_thread_logging(log_handler)

    def correct_base(self, pileup_column):
        log_msg = "Step {}: position {}, {} families out of {} reads."
        bases = pileup_column.get_query_sequences()
        if all([b==bases[0] for b in bases]):
            return
        families = get_families(bases, pileup_column.pileups)
        self.n_families += len(families)

        for v in families.values():

            uni, counts = np.unique(v['bases'], return_counts=True)
            if len(uni)>1:
                self.n_discordant += 1
                i = np.argmax(counts)
                self.info['umi_freq'].append(counts[i]/np.sum(counts))

                self.add_new.update([
                    r.alignment for r,b in zip(v['reads'], v['bases']) if b==uni[i]
                ])
                self.rm_new.update([
                    r.alignment for r,b in zip(v['reads'], v['bases']) if b!=uni[i]
                ])

        self.step += 1
        if self.step%self.save_every == 0:
            to_put = self.add_old.difference(self.rm_old.union(self.rm_new))
            self.queue.put(to_put)
            self.sem.release()
            self.add_old = self.add_new
            self.add_new = set()
            self.rm_old = self.rm_new
            self.rm_new = set()

        if self.step and self.step%self.log_every==0:
            logging.info(log_msg.format(
                self.step,
                pileup_column.reference_pos,
                len(families),
                len(bases)
            ))
