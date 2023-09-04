from itertools import repeat
import multiprocessing
import subprocess

import pyfastaq

from ushonium import utils


def mafft_stdout_to_seqs(mafft_stdout, ref_name):
    # The aligned sequences are printed to stdout. We should have the
    # reference genome and the aligned genome in there (do not assume what order
    # they are output). We just want to extract the aligned seqs.
    seqs = mafft_stdout.split(">")
    # seqs looks like: ["", "ref\nACGT...", "to_align\nACGT...", ...]
    assert len(seqs) >= 3
    assert seqs[0] == ""
    ref_seq = None
    aln_seqs = {}
    seqs_to_parse = [x.split("\n", maxsplit=1) for x in seqs[1:]]
    for (name, seq_str) in seqs_to_parse:
        if name == ref_name:
            ref_seq = seq_str.replace("\n", "")
        else:
            assert name not in aln_seqs
            aln_seqs[name] = seq_str.replace("\n", "")

    return ref_seq, aln_seqs


def fix_indels(ref_seq, aln_seq, method):
    if method == "nothing":
        return aln_seq
    if method not in {"as_ref", "N"}:
        raise NotImplementedError(f"Not implemented {method}")

    aln_seq = list(aln_seq)
    for i, c in enumerate(aln_seq):
        if c == "-":
            aln_seq[i] = "n" if method == "N" else ref_seq[i]
    return "".join(aln_seq)


def force_ref_at_ends(ref_seq, aln_seq, ref_start, ref_end):
    new_aln_seq = []
    new_ref_seq = []
    ref_pos = -1
    for aln_start, c in enumerate(ref_seq):
        if c != "-":
            ref_pos += 1
        if ref_pos == ref_start:
            break
        if c != "-":
            new_aln_seq.append(c)
            new_ref_seq.append(c)

    while aln_start < len(ref_seq) and ref_seq[aln_start] == "-":
        aln_start += 1

    ref_pos = len([x for x in ref_seq if x != "-"])
    aln_end = len(ref_seq)
    for c in reversed(ref_seq):
        aln_end -= 1
        if c != "-":
            ref_pos -= 1
        if ref_pos == ref_end:
            break

    while aln_end > 0 and ref_seq[aln_end] == "-":
        aln_end -= 1

    new_end = ref_seq[aln_end + 1 :].replace("-", "")
    new_ref_seq.extend([ref_seq[aln_start : aln_end + 1], new_end])
    new_aln_seq.extend([aln_seq[aln_start : aln_end + 1], new_end])
    return "".join(new_ref_seq), "".join(new_aln_seq)


def replace_start_end_indels_with_N(seq):
    new_seq = list(seq)
    chars = {"n", "N", "-"}
    i = 0
    while i < len(new_seq) and new_seq[i] in chars:
        new_seq[i] = "N"
        i += 1
    i = len(new_seq) - 1
    while i >= 0 and new_seq[i] in chars:
        new_seq[i] = "N"
        i -= 1
    return "".join(new_seq)


def run_mafft(qry_seq, ref_seq, quiet, indel_method, ref_start, ref_end):
    assert ref_start == None == ref_end or None not in [ref_start, ref_end]
    assert indel_method in {"nothing", "as_ref", "N"}

    script = "\n".join(
        [
            f'''ref=">ref\n{ref_seq}"''',
            f'''qry=">qry\n{qry_seq}"''',
            """mafft --keeplength --add <(echo "$qry") <(echo "$ref")""",
        ]
    )
    p = subprocess.run(
        ["bash"],
        input=script,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if p.returncode != 0:
        raise Exception(
            f"Error running mafft. Stdout:\n{p.stdout}\n\nStderr:{p.stderr}"
        )

    ref_seq, aln_seqs = mafft_stdout_to_seqs(p.stdout, "ref")
    assert len(aln_seqs) == 1
    aln_seq = replace_start_end_indels_with_N(aln_seqs["qry"])

    if ref_start is not None:
        assert 0 <= ref_start < ref_end
        ref_seq, aln_seq = force_ref_at_ends(ref_seq, aln_seq, ref_start, ref_end)

    if indel_method != "nothing":
        aln_seq = fix_indels(ref_seq, aln_seq, indel_method)

    return aln_seq


def run_mafft_one_qry_fasta(qry_fa, ref_seq, quiet, indel_method, ref_start, ref_end):
    qry_seq = utils.load_single_seq_fasta(qry_fa)
    return run_mafft(qry_seq.seq, ref_seq, quiet, indel_method, ref_start, ref_end)


def run_mafft_multi_qry_fasta(
    qry_fa,
    ref_name,
    ref_seq,
    outfile,
    quiet,
    indel_method,
    ref_start,
    ref_end,
    cpus=1,
):
    f_out = pyfastaq.utils.open_file_write(outfile)
    print(f">{ref_name}", ref_seq, sep="\n", file=f_out)
    qry_names = []
    qry_seqs = []
    file_reader = pyfastaq.sequences.file_reader(qry_fa)
    for seq in file_reader:
        if cpus == 1:
            aln_seq = run_mafft(
                seq.seq, ref_seq, quiet, indel_method, ref_start, ref_end
            )
            print(f">{seq.id}", aln_seq, sep="\n", file=f_out)
            continue

        qry_names.append(seq.id)
        qry_seqs.append(seq.seq)
        if len(qry_names) >= cpus:
            with multiprocessing.Pool(processes=cpus) as pool:
                results = pool.starmap(
                    run_mafft,
                    zip(
                        qry_seqs,
                        repeat(ref_seq),
                        repeat(quiet),
                        repeat(indel_method),
                        repeat(ref_start),
                        repeat(ref_end),
                    ),
                )
            for name, aln_seq in zip(qry_names, results):
                print(f">{name}", aln_seq, sep="\n", file=f_out)
            qry_names = []
            qry_seqs = []

    if len(qry_names) > 0:
        with multiprocessing.Pool(processes=cpus) as pool:
            results = pool.starmap(
                run_mafft,
                zip(
                    qry_seqs,
                    repeat(ref_seq),
                    repeat(quiet),
                    repeat(indel_method),
                    repeat(ref_start),
                    repeat(ref_end),
                ),
            )
        for name, aln_seq in zip(qry_names, results):
            print(f">{name}", aln_seq, sep="\n", file=f_out)

    pyfastaq.utils.close(f_out)
