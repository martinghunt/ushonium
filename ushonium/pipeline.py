from itertools import repeat
import logging
import multiprocessing
import os

import pyfastaq

from ushonium import mafft, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
REF_GB = os.path.abspath(os.path.join(this_dir, "data", "hu1.gb"))


class Pipeline:
    def __init__(self, options):
        # At least one of the tools we're calling writes files in cwd.
        # Change to output directory so we don't randomly write files somewhere else.
        # Means using absolute paths for everything before changing cwd
        self.ref_gb = utils.to_abs_path(options.ref_gb)
        self.outdir = utils.to_abs_path(options.outdir)
        self.samples_tsv = utils.to_abs_path(options.samples_tsv)
        self.fastas_fofn = utils.to_abs_path(options.fastas_fofn)
        self.matopt_opts = options.matopt_opts
        self.metadata_tsv = utils.to_abs_path(options.metadata_tsv)
        self.metacols = options.metacols
        self.metacol_name = options.metacol_name
        self.title = options.title
        self.cpus = options.cpus
        self.indel_method = options.indel_method
        self._set_ref_start_end(options.ref_start_end)
        self.multi_fa_mafft_chunk_size = 50

        self.ref_fa = "00.ref.fa"
        self.ref_name = None
        self.ref_seq = None
        self.tree_optimized = "04.optimized.pb"
        self.taxonium_out = "05.taxonium.jsonl.gz"

        if (self.fastas_fofn is None and self.samples_tsv is None) or (
            self.fastas_fofn is not None and self.samples_tsv is not None
        ):
            raise Exception(
                "Must provide either --fastas_fofn or --samples_tsv (but not both)"
            )

    def _set_ref_start_end(self, ref_start_end):
        if ref_start_end is None:
            self.ref_start = None
            self.ref_end = None
        else:
            try:
                self.ref_start, self.ref_end = ref_start_end.split(",")
                self.ref_start = int(self.ref_start) - 1
                self.ref_end = int(self.ref_end) - 1
            except:
                raise Exception(
                    f"Error parsing ref_start_end option: '{ref_start_end}'. Must be of the form START,END"
                )

    def load_samples_tsv(self):
        assert self.samples_tsv is not None
        samples = []
        sample_names = set()
        with open(self.samples_tsv) as f:
            for line in f:
                sample_name, filename = line.rstrip().split("\t")
                assert sample_name not in sample_names
                sample_names.add(sample_name)
                samples.append((sample_name, os.path.abspath(filename)))
        return samples

    def write_ref_fasta(self):
        logging.info(
            f"Making reference fasta file from genbank file: {self.ref_gb} -> {self.ref_fa}"
        )
        pyfastaq.tasks.to_fasta(self.ref_gb, self.ref_fa, line_length=0)
        ref_seq = utils.load_single_seq_fasta(self.ref_fa)
        self.ref_name = ref_seq.id.split()[0]
        self.ref_seq = ref_seq.seq

    def run_matOptimize(self, input_pb):
        logging.info(f"Optimizing the tree from usher {input_pb}")
        utils.syscall(
            f"matOptimize {self.matopt_opts} -T {self.cpus} -i {input_pb} -o {self.tree_optimized}",
            log=f"{self.tree_optimized}.log",
        )
        logging.info(f"Made optimized tree {self.tree_optimized}")

    def run_usher_to_taxonium(self):
        logging.info("Making taxonium jsonl file")
        if self.metadata_tsv is not None:
            assert self.metacols is not None
            meta_opts = f"--metadata {self.metadata_tsv} --columns {self.metacols}"
            if self.metacol_name is not None:
                meta_opts += f" --key_column {self.metacol_name}"
        else:
            meta_opts = ""
        utils.syscall(
            f'usher_to_taxonium --input {self.tree_optimized} --genbank {self.ref_gb} {meta_opts} --title "{self.title}" --output {self.taxonium_out}'
        )
        logging.info("Finished taxonium")
        logging.info(f"All done. Final taxonium file: {self.taxonium_out}")
        assert os.path.exists(self.taxonium_out)

    def run_using_samples_tsv(self):
        os.chdir(self.outdir)
        self.write_ref_fasta()

        logging.info(
            f"Load input TSV file of sample names and FASTA files {self.samples_tsv}"
        )
        samples = self.load_samples_tsv()
        logging.info(f"Got {len(samples)} samples from TSV file")

        msa_fa = "01.msa.fa"
        f_out_msa = open(msa_fa, "w")
        print(f">{self.ref_name}", self.ref_seq, sep="\n", file=f_out_msa)

        # Run the multiple alignments in batches. Each batch has size the same as
        # number of cpus. At the end of each batch, append to the final FASTA file.
        # Batching because otherwise would need to store all of the sequences in
        # memory and write at the end. Could be a lot of sequences!
        for i in range(0, len(samples), self.cpus):
            aln_seqs = {}
            names = []
            if self.cpus == 1:
                for sample_name, sample_fa in samples[i : i + self.cpus]:
                    logging.info(f"Running MSA on sample {sample_name}")
                    name, aln_seq = mafft.run_mafft_one_seq(
                        sample_name,
                        sample_fa,
                        self.ref_fa,
                        self.ref_name,
                        False,
                        self.indel_method,
                        self.ref_start,
                        self.ref_end,
                    )
                    names.append(sample_name)
                    aln_seqs[sample_name] = aln_seq
            else:
                names = [x[0] for x in samples[i : i + self.cpus]]
                fastas = [x[1] for x in samples[i : i + self.cpus]]
                logging.info(
                    f"Running batch of MSAs in parallel. Samples {','.join([x[0] for x in samples[i:i+self.cpus]])}"
                )
                with multiprocessing.Pool(processes=self.cpus) as p:
                    results = p.starmap(
                        mafft.run_mafft_one_seq,
                        zip(
                            names,
                            fastas,
                            repeat(self.ref_fa),
                            repeat(self.ref_name),
                            repeat(True),
                            repeat(self.indel_method),
                            repeat(self.ref_start),
                            repeat(self.ref_end),
                        ),
                    )
                aln_seqs = {x[0]: x[1] for x in results}

            for name in names:
                print(f">{name}", aln_seqs[name], sep="\n", file=f_out_msa)

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
            f"usher-sampled --sort-before-placement-3 --vcf {vcf} --tree {empty_tree} -T {self.cpus} --save-mutation-annotated-tree {tree_unoptimized}",
            log=f"{tree_unoptimized}.log",
        )

        self.run_matOptimize(tree_unoptimized)
        self.run_usher_to_taxonium()

    def load_and_check_fastas_fofn(self):
        filenames = []
        with open(self.fastas_fofn) as f:
            for filename in map(str.strip, f):
                if not os.path.exists(filename):
                    raise FileNotFoundError(f"FASTA file not found: {filename}")
                filenames.append(os.path.abspath(filename))
        return filenames

    def run_using_fastas_fofn(self):
        quiet = True
        fasta_input_files = self.load_and_check_fastas_fofn()
        os.chdir(self.outdir)
        self.write_ref_fasta()
        msa_fasta_dir = "01.msas"
        msa_fastas = [
            os.path.join(msa_fasta_dir, f"{i}.fa.gz")
            for i in range(len(fasta_input_files))
        ]

        # ------------- Make an MSA for each input fasta file -----------------
        msas_done = "01.msas.done"
        if not os.path.exists(msas_done):
            utils.syscall(f"rm -rf {msa_fasta_dir}")
            os.mkdir(msa_fasta_dir)
            with multiprocessing.Pool(processes=self.cpus) as p:
                p.starmap(
                    mafft.run_mafft_multi_fasta_chunked,
                    zip(
                        fasta_input_files,
                        repeat(self.ref_fa),
                        repeat(self.ref_name),
                        msa_fastas,
                        repeat(quiet),
                        repeat(self.indel_method),
                        repeat(self.ref_start),
                        repeat(self.ref_end),
                    ),
                )
            with open(msas_done, "w") as f:
                pass

        # ----------------------- Make VCF from each MSA ----------------------
        vcf_dir = "02.vcfs"
        vcfs_done = "02.vcfs.done"
        vcfs = [
            os.path.abspath(os.path.join(vcf_dir, f"{i}.vcf"))
            for i in range(len(fasta_input_files))
        ]
        if not os.path.exists(vcfs_done):
            utils.syscall(f"rm -rf {vcf_dir}")
            os.mkdir(vcf_dir)
            commands = [
                f"faToVcf -includeNoAltN {msa_fastas[i]} {vcfs[i]}"
                for i in range(len(fasta_input_files))
            ]
            with multiprocessing.Pool(processes=self.cpus) as p:
                p.map(utils.syscall, commands)
            with open(vcfs_done, "w") as f:
                pass

        # ------------------ iteratively make tree (unoptimized) --------------
        tree_file = None
        tree_iter_dir = "03.tree_iters"
        if not os.path.exists(tree_iter_dir):
            os.mkdir(tree_iter_dir)

        for i, vcf_file in enumerate(vcfs):
            iter_outdir = os.path.join(tree_iter_dir, f"iter.{i}")
            done_file = os.path.join(iter_outdir, "done")

            if not os.path.exists(done_file):
                utils.syscall(f"rm -rf {iter_outdir}")
                os.mkdir(iter_outdir)

                if i == 0:
                    logging.info(f"Tree iteration {i} start")
                    assert tree_file is None
                    tree_file = os.path.abspath(os.path.join(iter_outdir, "empty_tree"))
                    with open(tree_file, "w") as f:
                        print("()", file=f)
                    utils.syscall(
                        f"usher-sampled --sort-before-placement-3 --vcf {vcf_file} --tree {tree_file} -T {self.cpus} --save-mutation-annotated-tree tree.pb &> usher.stdouterr",
                        cwd=iter_outdir,
                    )
                else:
                    assert tree_file is not None
                    utils.syscall(
                        f"usher-sampled --sort-before-placement-3 --vcf {vcf_file} --load-mutation-annotated-tree {tree_file} -T {self.cpus} --save-mutation-annotated-tree tree.pb &> usher.stdouterr",
                        cwd=iter_outdir,
                    )

                with open(done_file, "w") as f:
                    pass
                logging.info(f"Tree iteration {i} finished")
            else:
                logging.info(
                    f"Tree iteration {i} skipped because already done {iter_outdir}"
                )

            tree_file = os.path.abspath(os.path.join(iter_outdir, "tree.pb"))
            assert os.path.exists(tree_file)

        # -------------- optimize tree and make taxonium file -----------------
        self.run_matOptimize(tree_file)
        self.run_usher_to_taxonium()

    def run(self):
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        original_dir = os.getcwd()

        if self.samples_tsv is not None:
            self.run_using_samples_tsv()
        elif self.fastas_fofn is not None:
            self.run_using_fastas_fofn()
        else:
            raise NotImplementedError()

        os.chdir(original_dir)


def run(options):
    p = Pipeline(options)
    p.run()
