from itertools import repeat
import logging
import multiprocessing
import os

import pyfastaq

from ushonium import mafft, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
REF_GB = os.path.abspath(os.path.join(this_dir, "data", "hu1.gb"))


def load_samples_tsv(infile):
    samples = []
    sample_names = set()
    with open(infile) as f:
        for line in f:
            sample_name, filename = line.rstrip().split("\t")
            assert sample_name not in sample_names
            sample_names.add(sample_name)
            samples.append((sample_name, os.path.abspath(filename)))
    return samples


def run(options):
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
    original_dir = os.getcwd()
    os.chdir(options.outdir)

    ref_fa = "00.ref.fa"
    logging.info(
        f"Making reference fasta file from genbank file: {options.ref_gb} -> {ref_fa}"
    )
    pyfastaq.tasks.to_fasta(options.ref_gb, ref_fa, line_length=0)
    ref_seq = utils.load_single_seq_fasta(ref_fa)
    msa_fa = "01.msa.fa"
    f_out_msa = open(msa_fa, "w")
    print(f">{ref_seq.id}", ref_seq.seq, sep="\n", file=f_out_msa)

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
                    mafft.run_mafft_one_seq(
                        sample_name,
                        sample_fa,
                        ref_fa,
                        ref_seq.id,
                        False,
                        options.indel_method,
                    )
                )
        else:
            names = [x[0] for x in samples[i : i + options.cpus]]
            fastas = [x[1] for x in samples[i : i + options.cpus]]
            logging.info(
                f"Running batch of MSAs in parallel. Samples {','.join([x[0] for x in samples[i:i+options.cpus]])}"
            )
            with multiprocessing.Pool(processes=options.cpus) as p:
                new_seqs = p.starmap(
                    mafft.run_mafft_one_seq,
                    zip(
                        names,
                        fastas,
                        repeat(ref_fa),
                        repeat(ref_seq.id),
                        repeat(True),
                        repeat(options.indel_method),
                    ),
                )

        for name, seq in new_seqs:
            print(f">{name}", seq, sep="\n", file=f_out_msa)

    f_out_msa.close()
    logging.info("Made MSA of all sequences")

    vcf = "02.faToVcf.vcf"
    logging.info(f"Making VCF file of all samples {vcf}")
    utils.syscall(f"faToVcf -includeNoAltN {msa_fa} {vcf}")

    logging.info("Running usher")
    empty_tree = "03.tmp.empty_tree"
    with open(empty_tree, "w") as f:
        print("()", file=f)
    tree_unoptimized = "03.usher.pb"
    utils.syscall(
        f"usher-sampled --sort-before-placement-3 --vcf {vcf} --tree {empty_tree} -T {options.cpus} --save-mutation-annotated-tree {tree_unoptimized}",
        log=f"{tree_unoptimized}.log",
    )

    logging.info(f"Optimizing the tree from usher {tree_unoptimized}")
    tree_optimized = "04.optimized.pb"
    utils.syscall(
        f"matOptimize {options.matopt_opts} --vcf {vcf} -T {options.cpus} -i {tree_unoptimized} -o {tree_optimized}",
        log=f"{tree_optimized}.log",
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
    utils.syscall(
        f'usher_to_taxonium --input {tree_optimized} --genbank {options.ref_gb} {meta_opts} --title "{options.title}" --output {taxonium_out}'
    )
    logging.info("Finished taxonium")
    logging.info(f"All done. Final taxonium file: {taxonium_out}")
    assert os.path.exists(taxonium_out)
    os.chdir(original_dir)
