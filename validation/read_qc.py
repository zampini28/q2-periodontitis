"""Read-level QC used to choose trimming and truncation before denoising.

    python -I read_qc.py <fastq.gz> [<fastq.gz> ...] [--reads N]

For a sample of reads from the given files, prints the read-length
distribution, the position of known 16S primers in the reads, and the
quartiles of the quality score every 10 positions. Nothing here looks at
taxonomy or indices, so the parameters it informs are fixed before any
index is computed.
"""
import argparse
import gzip
import re
from collections import Counter

import numpy as np

PRIMERS = {
    "27F": "AGAGTTTGATCMTGGCTCAG",
    "338R": "TGCTGCCTCCCGTAGGAGT",
    "341F": "CCTACGGGNGGCWGCAG",
    "338F": "ACTCCTACGGGAGGCAGCAG",
    "515F": "GTGYCAGCMGCCGCGGTAA",
    "534R": "ATTACCGCGGCTGCTGG",
    "785R": "GACTACHVGGGTATCTAATCC",
    "806R": "GGACTACNVGGGTWTCTAAT",
    "926R": "CCGYCAATTYMTTTRAGTTT",
}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
         "S": "[CG]", "W": "[AT]", "K": "[GT]", "M": "[AC]", "B": "[CGT]",
         "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}


def pattern(seq, prefix=14):
    """Regex for the 3' end of a primer (the part that must match)."""
    return re.compile("".join(IUPAC[b] for b in seq[-prefix:]))


def reads(paths, n):
    per_file = max(1, n // len(paths))
    for path in paths:
        with gzip.open(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i // 4 >= per_file:
                    break
                if i % 4 == 1:
                    seq = line.strip()
                elif i % 4 == 3:
                    yield seq, line.strip()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq", nargs="+")
    parser.add_argument("--reads", type=int, default=20000)
    args = parser.parse_args()
    patterns = {k: pattern(v) for k, v in PRIMERS.items()}
    lengths = Counter()
    hits = {k: Counter() for k in PRIMERS}
    quals = []
    for seq, qual in reads(args.fastq, args.reads):
        lengths[len(seq)] += 1
        window = seq[:80]
        for name, pat in patterns.items():
            m = pat.search(window)
            if m:
                hits[name][m.end()] += 1
        quals.append(np.frombuffer(qual.encode(), dtype=np.uint8) - 33)
    total = sum(lengths.values())
    print("reads sampled: %d" % total)
    q = np.percentile(list(lengths.elements()), [5, 25, 50, 75, 95])
    print("length percentiles 5/25/50/75/95: %s"
          % " ".join("%d" % v for v in q))
    print("primers found in the first 80 nt "
          "(fraction; most common end position):")
    for name, c in hits.items():
        n = sum(c.values())
        if n / total > 0.05:
            print("  %-5s %.2f  ends at %d" % (name, n / total,
                                               c.most_common(1)[0][0]))
    width = max(len(x) for x in quals)
    mat = np.full((len(quals), width), np.nan)
    for i, x in enumerate(quals):
        mat[i, :len(x)] = x
    print("quality quartiles (position: q25/median/q75, reads covering):")
    for pos in range(0, width, 10):
        col = mat[:, pos]
        col = col[~np.isnan(col)]
        if len(col) < 0.05 * total:
            break
        a, b, c = np.percentile(col, [25, 50, 75])
        print("  %3d: %2d/%2d/%2d  %.2f"
              % (pos + 1, a, b, c, len(col) / total))
    print("distinct quality values: %d"
          % len(np.unique(mat[~np.isnan(mat)])))


if __name__ == "__main__":
    main()
