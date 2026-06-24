"""Console entry point for CNVkit-based copy-number calling."""

from __future__ import annotations

import argparse
import subprocess
from contextlib import nullcontext
from importlib import resources
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Call CNVs from a BAM file using CNVkit."
    )
    parser.add_argument(
        "args",
        nargs="+",
        metavar="ARG",
        help=(
            "Either: <bam> <output_dir> <threads>, or legacy form: "
            "<bam> <reference_cnn> <output_dir> <threads>."
        ),
    )
    parser.add_argument(
        "--reference-cnn",
        default=None,
        help="CNVkit reference .cnn file. Defaults to the packaged hg38 5 kb reference.",
    )
    return parser


def _parse_args(argv):
    args = build_parser().parse_args(argv)
    if len(args.args) == 3:
        args.bam, args.output_dir, threads = args.args
    elif len(args.args) == 4:
        args.bam, legacy_reference, args.output_dir, threads = args.args
        if args.reference_cnn is None:
            args.reference_cnn = legacy_reference
    else:
        raise SystemExit(
            "Expected either <bam> <output_dir> <threads> or "
            "<bam> <reference_cnn> <output_dir> <threads>."
        )
    args.threads = int(threads)
    return args


def _reference_context(reference_cnn):
    if reference_cnn:
        return nullcontext(Path(reference_cnn))
    return resources.path("bfbarchitect.resources", "hg38full_ref_5k.cnn")


def main() -> None:
    args = _parse_args(None)

    with _reference_context(args.reference_cnn) as reference_cnn:
        cmd = [
            "cnvkit.py",
            "batch",
            args.bam,
            "--seq-method",
            "wgs",
            "--drop-low-coverage",
            "--reference",
            str(reference_cnn),
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
