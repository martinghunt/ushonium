#!/usr/bin/env python3

import argparse
import logging
import os

from ushonium import pipeline


def main(args=None):
    parser = argparse.ArgumentParser(
        description="make taxonium jsonl from fasta files",
        usage="ushonium [options] <samples_tsv> <outdir>",
    )
    parser.add_argument(
        "--ref_gb",
        help="Genbank file of reference [%(default)s]",
        default=pipeline.REF_GB,
        metavar="FILENAME",
    )
    parser.add_argument(
        "--matopt_opts",
        help="matOptimize options (don't use -T here because that's determined by the --cpus option to this script). [%(default)s]",
        metavar="STR",
        default="-m 0.000000001 -r 8",
    )
    parser.add_argument(
        "--metadata_tsv",
        help="TSV file of metadata from which to get annotations. If used, must also use --metacols",
        metavar="FILENAME",
    )
    parser.add_argument(
        "--metacols",
        help="Comma-separated list of columns to use from metadata_tsv file (see --metdata_tsv)",
        metavar="COL_x,COL_y,...",
    )
    parser.add_argument(
        "--metacol_name",
        help="Column name in metadata_tsv that has the name of each sample. Every name in --samples_tsv must be in this column. If not provided, then the default for usher_to_taxonium is used, which assumes 'strain'",
        metavar="STRING",
    )
    parser.add_argument(
        "--title",
        help="Title you end up seeing in taxonium browser [%(default)s]",
        default="title",
        metavar="STRING",
    )
    parser.add_argument(
        "--cpus",
        type=int,
        help="Number of CPUs [%(default)s]",
        default=1,
        metavar="INT",
    )
    parser.add_argument(
        "--indel_method",
        choices=["nothing", "as_ref", "N"],
        help="How to handle indels in the msa from mafft. nothing=do what mafft does which is keep '-'; as_ref: replace '-' with the ref call; N: replace '-' with N [%(default)s]",
        default="as_ref",
    )
    parser.add_argument(
        "--ref_start_end",
        help="Start and end coord of the ref to use, eg --ref_start_end 100,28000 would trim all MSAs to reference coords 100-28000 inclusive (1-based coords)",
        metavar="INT,INT",
    )
    parser.add_argument(
        "samples_tsv",
        help="TSV file of sample names (column 1) and FASTA filenames (column 2). No header line in file",
    )
    parser.add_argument("outdir", help="Output directory (will be created)")

    options = parser.parse_args()

    if options.metadata_tsv is not None:
        if options.metacols is None:
            raise Exception("If using --metadata_tsv, must also use --metacols")
        options.metadata_tsv = os.path.abspath(options.metadata_tsv)

    log = logging.getLogger()
    log.setLevel(logging.INFO)
    pipeline.run(options)


if __name__ == "__main__":
    main()
