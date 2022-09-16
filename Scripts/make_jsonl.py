#!/usr/bin/env python3

import argparse
from itertools import repeat
import logging
import multiprocessing
import os
import subprocess

import pyfastaq


def syscall(command, cwd=None, quiet=False):
    if not quiet:
        logging.info(f"Run: {command}")
    p = subprocess.run(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        cwd=cwd,
    )
    if p.returncode != 0:
        raise Exception(f"Command failed (exit code {p.returncode}): {command}")
    return p


def load_single_seq_fasta(filename):
    seqs = {}
    pyfastaq.tasks.file_to_dict(filename, seqs)
    assert len(seqs) == 1
    return list(seqs.values())[0]


def run_mafft_one_seq(name, to_align, ref_fa, ref_name, quiet):
    command = f"mafft --quiet --keeplength --add {to_align} {ref_fa}"
    process = syscall(command, quiet=quiet)
    # The aligned sequences are printed to stdout. We should have the
    # reference genome and the aligned genome in there (do not assume what order
    # they are output). We just want to extract the aligned seq.
    seqs = process.stdout.split(">")
    assert len(seqs) == 3
    seqs = [x for x in seqs if len(x) > 0 and not x.startswith(ref_name)]
    assert len(seqs) == 1
    # The aligned seq string we want will have linebreak characters in it,
    # which is fine because all we do is print it to a FASTA file later anyway
    return name, seqs[0].split("\n", maxsplit=1)[1].strip()


def load_samples_tsv(infile):
    samples = []
    sample_names = set()
    with open(options.samples_tsv) as f:
        for line in f:
            sample_name, filename = line.rstrip().split("\t")
            assert sample_name not in sample_names
            sample_names.add(sample_name)
            samples.append((sample_name, os.path.abspath(filename)))
    return samples


THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REF_GB = os.path.abspath(os.path.join(THIS_DIR, os.pardir, "data", "hu1.gb"))

parser = argparse.ArgumentParser(
    description="make taxonium jsonl from fasta files",
    usage="%(prog)s [options] <samples_tsv> <outdir>",
)
parser.add_argument(
    "--ref_gb",
    help="Genbank file of reference [%(default)s]",
    default=REF_GB,
    metavar="FILENAME",
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
    "--cpus", type=int, help="Number of CPUs [%(default)s]", default=1, metavar="INT"
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


logging.info(
    f"Load input TSV file of sample names and FASTA files {options.samples_tsv}"
)
samples = load_samples_tsv(options.samples_tsv)
logging.info(f"Got {len(samples)} samples from TSV file")


# At least one of the tools we're calling writes files in cwd.
# Change to output directory so we don't randomly write files somewhere else.
# Means using absolute paths for everything before changing cwd
options.ref_gb = os.path.abspath(options.ref_gb)
if options.metadata_tsv is not None:
    options.metadata_tsv = os.path.abspath(options.metadata_tsv)
os.mkdir(options.outdir)
os.chdir(options.outdir)

ref_fa = "00.ref.fa"
logging.info(
    f"Making reference fasta file from genbank file: {options.ref_gb} -> {ref_fa}"
)
pyfastaq.tasks.to_fasta(options.ref_gb, ref_fa)
ref_seq = load_single_seq_fasta(ref_fa)
msa_fa = "01.msa.fa"
f_out_msa = open(msa_fa, "w")
print(ref_seq, file=f_out_msa)


# Run the multiple alignments in batches. Each batch has size the same as
# number of cpus. At the end of each batch, append to the final FASTA file.
# Batching because otherwise would need to store all of the sequences in
# memory and write at the end. Could be a lot of sequences!
for i in range(0, len(samples), options.cpus):
    new_seqs = []
    if options.cpus == 1:
        for sample_name, sample_fa in samples[i : i + options.cpus]:
            logging.info(f"Running MSA on sample {sample_name}")
            new_seqs.append(
                run_mafft_one_seq(sample_name, sample_fa, ref_fa, ref_seq.id, False)
            )
    else:
        names = [x[0] for x in samples[i : i + options.cpus]]
        fastas = [x[1] for x in samples[i : i + options.cpus]]
        logging.info(
            f"Running batch of MSAs in parallel. Samples {','.join([x[0] for x in samples[i:i+options.cpus]])}"
        )
        with multiprocessing.Pool(processes=options.cpus) as p:
            new_seqs = p.starmap(
                run_mafft_one_seq,
                zip(names, fastas, repeat(ref_fa), repeat(ref_seq.id), repeat(True)),
            )

    for name, seq in new_seqs:
        print(f">{name}", seq, sep="\n", file=f_out_msa)

f_out_msa.close()
logging.info("Made MSA of all sequences")


vcf = "02.faToVcf.vcf"
logging.info(f"Making VCF file of all samples {vcf}")
syscall(f"faToVcf -includeNoAltN {msa_fa} {vcf}")


logging.info("Running usher")
empty_tree = "03.tmp.empty_tree"
with open(empty_tree, "w") as f:
    print("()", file=f)
tree_unoptimized = "03.usher.pb"
syscall(
    f"usher --sort-before-placement-3 --vcf {vcf} --tree {empty_tree} -T {options.cpus} --save-mutation-annotated-tree {tree_unoptimized} &> {tree_unoptimized}.log"
)


logging.info(f"Optimizing the tree from usher {tree_unoptimized}")
tree_optimized = "04.optimized.pb"
syscall(
    f"matOptimize -T {options.cpus} -r 8 -M 2 -i {tree_unoptimized} -o {tree_optimized} &> {tree_optimized}.log"
)
logging.info(f"Made optimized tree {tree_optimized}")


logging.info("Making taxonium jsonl file")
taxonium_out = "05.taxonium.jsonl.gz"
if options.metadata_tsv is not None:
    assert options.metacols is not None
    meta_opts = f"--metadata {options.metadata_tsv} --columns {options.metacols}"
    if options.metacol_name is not None:
        meta_opts += f" --key_column {options.metacol_name}"
else:
    meta_opts = ""
syscall(
    f'usher_to_taxonium --input {tree_optimized} --genbank {options.ref_gb} {meta_opts} --title "{options.title}" --output {taxonium_out}'
)
logging.info("Finished taxonium")
logging.info(f"All done. Final taxonium file: {taxonium_out}")
