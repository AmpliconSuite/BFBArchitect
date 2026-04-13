#!/usr/bin/env python3
"""Call CNVs from a BAM file using CNVkit.

Equivalent to the shell command:
cnvkit.py batch <bam> --seq-method wgs --drop-low-coverage --reference <cnn>
    --scatter --diagram -d <output_dir> -p <threads>
"""

from __future__ import annotations

import argparse
import subprocess


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Basic script employing cnvkit.py to call CNVs from alignments."
    )
    parser.add_argument("bam", help="Input BAM file")
    parser.add_argument("reference_cnn", help="CNVkit reference .cnn file")
    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("threads", type=int, help="Number of threads")
    return parser


def main() -> None:
    args = build_parser().parse_args()

    cmd = [
        "cnvkit.py",
        "batch",
        args.bam,
        "--seq-method",
        "wgs",
        "--drop-low-coverage",
        "--reference",
        args.reference_cnn,
        "--scatter",
        "--diagram",
        "-d",
        args.output_dir,
        "-p",
        str(args.threads),
    ]

    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()