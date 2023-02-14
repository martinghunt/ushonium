from ushonium import utils


def mafft_stdout_to_seqs(mafft_stdout, ref_name):
    # The aligned sequences are printed to stdout. We should have the
    # reference genome and the aligned genome in there (do not assume what order
    # they are output). We just want to extract the aligned seqs.
    seqs = mafft_stdout.split(">")
    # seqs looks like: ["", "ref\nACGT...", "to_align\nACGT..."]
    assert len(seqs) == 3
    assert seqs[0] == ""
    seqs = [x.split("\n", maxsplit=1) for x in seqs[1:]]
    ref_seq = [x[1] for x in seqs if x[0].startswith(ref_name)]
    aln_seq = [x[1] for x in seqs if not x[0].startswith(ref_name)]
    assert len(ref_seq) == 1
    assert len(aln_seq) == 1
    return ref_seq[0].replace("\n", ""), aln_seq[0].replace("\n", "")


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


def run_mafft_one_seq(name, to_align, ref_fa, ref_name, quiet, indel_method):
    assert indel_method in {"nothing", "as_ref", "N"}
    command = f"mafft --quiet --keeplength --add {to_align} {ref_fa}"
    process = utils.syscall(command, quiet=quiet)
    ref_seq, aln_seq = mafft_stdout_to_seqs(process.stdout, ref_name)
    return name, fix_indels(ref_seq, aln_seq, indel_method)
