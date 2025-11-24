import argparse
import gzip
import sys
import os
import math
from multiprocessing import Pool, cpu_count
from functools import partial

__AUTHOR__="Xuewen Wang"

def parse_size(size_str):
    """
    Parses genome size string to integer.
    Supports:
      - Suffixes: 1.5M, 3k, 2G (case-insensitive)
      - Raw Integers: 150345
    """
    if not size_str:
        return 0

    s = str(size_str).strip().lower()
    if not s:
        return 0

    units = {'k': 1_000, 'm': 1_000_000, 'g': 1_000_000_000}
    multiplier = 1

    # Check if the last character is a suffix (k, m, g)
    if s[-1] in units:
        multiplier = units[s[-1]]
        s = s[:-1] # Remove the suffix to get the number part

    try:
        # Convert to float first to handle decimals like "1.5", then integer
        return int(float(s) * multiplier)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid genome size: '{size_str}'. Use format like 500M, 3.2G, or a raw integer (e.g. 150345).")

def parse_args():
    parser = argparse.ArgumentParser(description="Summarize sequencing reads (FASTA/FASTQ) with parallel processing.")
    parser.add_argument("-r", "--reads", nargs='+', required=True, help="List of input files (supports .fasta, .fastq, .gz)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file path")
    parser.add_argument("-t", "--threads", type=int, default=min(4, cpu_count()), help="Number of threads for parallel processing")
    parser.add_argument("-g", "--genome_size", type=parse_size, default=0, help="Estimated genome size (bp) for coverage. Supports suffixes k/M/G (e.g. 500M) or raw integers (e.g. 150345).")
    return parser.parse_args()

def get_file_opener(filename):
    """Returns the appropriate open function based on file extension."""
    if filename.endswith(".gz"):
        return gzip.open
    return open

def analyze_file(file_path, genome_size=0):
    """
    Parses a FASTA or FASTQ file and calculates statistics.
    Returns a dictionary of stats.
    """
    stats = {
        "file": file_path,
        "format": "Unknown",
        "type": "DNA",  # Default assumption, checking first few bases
        "num_seqs": 0,
        "sum_len": 0,
        "min_len": float('inf'),
        "max_len": 0,
        "avg_len": 0,
        "Coverage (x)": "N/A"
    }

    try:
        opener = get_file_opener(file_path)
        # Open in text mode ('rt')
        with opener(file_path, 'rt') as f:
            first_char = f.read(1)
            f.seek(0)

            if not first_char:
                return stats # Empty file

            # Determine format
            if first_char == '>':
                stats["format"] = "FASTA"
                parse_func = parse_fasta
            elif first_char == '@':
                stats["format"] = "FASTQ"
                parse_func = parse_fastq
            else:
                stats["format"] = "Unknown"
                return stats

            # Process the file
            lengths, is_rna = parse_func(f)

            # Aggregate results
            if lengths:
                stats["num_seqs"] = len(lengths)
                stats["sum_len"] = sum(lengths)
                stats["min_len"] = min(lengths)
                stats["max_len"] = max(lengths)
                stats["avg_len"] = round(stats["sum_len"] / stats["num_seqs"], 1)

                if is_rna:
                    stats["type"] = "RNA"

                if genome_size > 0:
                    stats["Coverage (x)"] = round(stats["sum_len"] / genome_size, 2)
            else:
                stats["min_len"] = 0

    except Exception as e:
        sys.stderr.write(f"Error processing {file_path}: {str(e)}\n")
        return None

    return stats

def parse_fasta(handle):
    """Parses FASTA format, returns list of lengths and RNA detection bool."""
    lengths = []
    is_rna = False
    seq_len = 0
    checked_type = False

    # Simple state machine for FASTA
    for line in handle:
        line = line.strip()
        if not line: continue

        if line.startswith(">"):
            if seq_len > 0:
                lengths.append(seq_len)
            seq_len = 0
        else:
            seq_len += len(line)
            # Check for Uracil to detect RNA (heuristic on first few sequences)
            if not checked_type and len(lengths) < 100:
                if 'U' in line.upper():
                    is_rna = True
                    checked_type = True

    # Add the last sequence
    if seq_len > 0:
        lengths.append(seq_len)

    return lengths, is_rna

def parse_fastq(handle):
    """Parses FASTQ format, returns list of lengths and RNA detection bool."""
    lengths = []
    is_rna = False
    line_idx = 0
    checked_type = False

    # FASTQ is strictly 4 lines per record usually, but can handle wrapped lines if careful.
    # We will assume standard 4-line FASTQ for speed/simplicity or use a robust iterator.
    # This loop assumes standard 4-line FASTQ structure.

    try:
        while True:
            header = handle.readline()
            if not header: break # EOF

            seq = handle.readline().strip()
            handle.readline() # + separator
            handle.readline() # Quality scores

            if not seq: break # Premature end

            lengths.append(len(seq))

            if not checked_type and line_idx < 100:
                if 'U' in seq.upper():
                    is_rna = True
                    checked_type = True
            line_idx += 1

    except ValueError:
        pass

    return lengths, is_rna

def main():
    args = parse_args()

    # Prepare arguments for the worker function
    func = partial(analyze_file, genome_size=args.genome_size)

    # Run in parallel
    results = []
    with Pool(processes=args.threads) as pool:
        # Map returns results in order of input list
        for res in pool.imap(func, args.reads):
            if res:
                results.append(res)

    # Write Output
    if results:
        headers = ["file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len", "Coverage (x)"]

        try:
            with open(args.output, 'w') as out:
                # Write header
                out.write("\t".join(headers) + "\n")

                # Write rows
                for stat in results:
                    row = [str(stat.get(h, "")) for h in headers]
                    out.write("\t".join(row) + "\n")

            print(f"Successfully processed {len(results)} files. Output written to {args.output}")

        except IOError as e:
            print(f"Error writing output file: {e}")
    else:
        print("No results generated. Check input files.")

if __name__ == "__main__":
    main()
#usage: python readsSum.py -r genome.fq.gz -o hifi_reads.sum.txt -t 40
