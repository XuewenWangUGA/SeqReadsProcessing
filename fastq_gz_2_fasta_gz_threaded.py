#!/usr/bin/env python3
''''
python fastq_to_fasta_multithreaded.py -i reads.fastq.gz -o filtered.fasta -l 5000 -t 16
later to do : pigz filtered.fasta will generate filtered.fasta.gz
__author__ ="Xuewen Wang"
'''
import argparse
import gzip
from Bio import SeqIO
import threading
import queue
from collections import Counter

# Thread-safe counters and storage
stats = Counter()
lengths = []
lengths_lock = threading.Lock()

# Queues
record_queue = queue.Queue(maxsize=1000)
result_queue = queue.Queue()
SENTINEL = object()

def reader_thread(input_path):
    """Read FASTQ records and push to queue."""
    with gzip.open(input_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            record_queue.put(record)
        for _ in range(num_workers):
            record_queue.put(SENTINEL)

def worker_thread(minlen):
    """Filter reads and push good ones to result queue."""
    while True:
        record = record_queue.get()
        if record is SENTINEL:
            result_queue.put(SENTINEL)
            break

        stats["total_reads"] += 1
        seqlen = len(record.seq)
        stats["total_bases"] += seqlen

        if seqlen >= minlen:
            stats["kept_reads"] += 1
            stats["kept_bases"] += seqlen
            with lengths_lock:
                lengths.append(seqlen)
            result_queue.put(record)

        record_queue.task_done()

def writer_thread(output_path):
    """Write good reads to fasta.gz."""
    done_count = 0
    #with gzip.open(output_path, "wt") as out_handle:
    with open(output_path, "wt") as out_handle:
        while True:
            item = result_queue.get()
            if item is SENTINEL:
                done_count += 1
                if done_count == num_workers:
                    break
                continue
            SeqIO.write(item, out_handle, "fasta")
            result_queue.task_done()

def compute_n50(lengths_list):
    if not lengths_list:
        return "N/A"
    sorted_lens = sorted(lengths_list, reverse=True)
    half = sum(sorted_lens) / 2
    acc = 0
    for l in sorted_lens:
        acc += l
        if acc >= half:
            return l
    return "N/A"

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ.GZ file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    parser.add_argument("-l", "--min-length", type=int, default=1, help="Min read length")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Worker threads")
    args = parser.parse_args()

    minlen = args.min_length
    num_workers = args.threads

    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Min Length: {minlen}")
    print(f"Threads: {num_workers}")
    print("-" * 40)

    # Start threads
    reader = threading.Thread(target=reader_thread, args=(args.input,))
    workers = [threading.Thread(target=worker_thread, args=(minlen,)) for _ in range(num_workers)]
    writer = threading.Thread(target=writer_thread, args=(args.output,))

    reader.start()
    for w in workers:
        w.start()
    writer.start()

    # Wait for all to complete
    reader.join()
    for w in workers:
        w.join()
    writer.join()

    # Print stats
    print("\n--- Summary Statistics ---")
    print(f"Total input reads: {stats['total_reads']}")
    print(f"Total input bases: {stats['total_bases']}")
    print(f"Filtered reads: {stats['kept_reads']}")
    print(f"Filtered bases: {stats['kept_bases']}")
    print(f"N50 of filtered reads: {compute_n50(lengths)}")
