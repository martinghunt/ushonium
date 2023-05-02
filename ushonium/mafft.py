import os

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
        i =- 1
    return "".join(new_seq)


def run_mafft_one_seq(
    name, to_align, ref_fa, ref_name, quiet, indel_method, ref_start, ref_end
):
    assert ref_start == None == ref_end or None not in [ref_start, ref_end]
    assert indel_method in {"nothing", "as_ref", "N"}
    if to_align.endswith(".gz"):
        to_align_original = to_align
        to_align = f"tmp.{name}.fa"
        assert not os.path.exists(to_align)
        utils.syscall(f"gunzip -c {to_align_original} > {to_align}", quiet=quiet)
    else:
        to_align_original = None

    command = f"mafft --quiet --keeplength --add {to_align} {ref_fa}"
    process = utils.syscall(command, quiet=quiet)
    if to_align_original is not None:
        os.unlink(to_align)
    ref_seq, aln_seq = mafft_stdout_to_seqs(process.stdout, ref_name)
    aln_seq = replace_start_end_indels_with_N(aln_seq)
    if ref_start is not None:
        assert 0 <= ref_start < ref_end
        ref_seq, aln_seq = force_ref_at_ends(ref_seq, aln_seq, ref_start, ref_end)

    if indel_method != "nothing":
        aln_seq = fix_indels(ref_seq, aln_seq, indel_method)
    return name, aln_seq
