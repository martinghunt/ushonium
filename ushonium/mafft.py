import os
import tempfile

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


def run_mafft_one_fasta(
    to_align, ref_fa, ref_name, quiet, indel_method, ref_start, ref_end
):
    assert ref_start == None == ref_end or None not in [ref_start, ref_end]
    assert indel_method in {"nothing", "as_ref", "N"}
    # mafft can't take gzip files (and the <(zcat foo.gz) trick doesn't work)
    if to_align.endswith(".gz"):
        with tempfile.TemporaryDirectory() as tempdir:
            to_align_original = to_align
            to_align = os.path.join(tempdir, "tmp.fa")
            utils.syscall(f"gunzip -c {to_align_original} > {to_align}", quiet=quiet)
            command = f"mafft --quiet --keeplength --add {to_align} {ref_fa}"
            process = utils.syscall(command, quiet=quiet)
    else:
        command = f"mafft --quiet --keeplength --add {to_align} {ref_fa}"
        process = utils.syscall(command, quiet=quiet)

    ref_seq, aln_seqs = mafft_stdout_to_seqs(process.stdout, ref_name)
    for name in aln_seqs:
        aln_seqs[name] = replace_start_end_indels_with_N(aln_seqs[name])
        if ref_start is not None:
            assert 0 <= ref_start < ref_end
            ref_seq2, aln_seqs[name] = force_ref_at_ends(
                ref_seq, aln_seqs[name], ref_start, ref_end
            )
        else:
            ref_seq2 = ref_seq

        if indel_method != "nothing":
            aln_seqs[name] = fix_indels(ref_seq2, aln_seqs[name], indel_method)
    return aln_seqs


def run_mafft_one_seq(
    name, to_align, ref_fa, ref_name, quiet, indel_method, ref_start, ref_end
):
    aligned_seqs = run_mafft_one_fasta(
        to_align, ref_fa, ref_name, quiet, indel_method, ref_start, ref_end
    )
    assert len(aligned_seqs) == 1
    return name, list(aligned_seqs.values())[0]


def _write_fasta_chunk(outfile, names, seqs):
    with open(outfile, "w") as f:
        for name, seq in zip(names, seqs):
            print(f">{name}", seq, sep="\n", file=f)


def run_mafft_multi_fasta_chunked(
    to_align,
    ref_fa,
    ref_name,
    outfile,
    quiet,
    indel_method,
    ref_start,
    ref_end,
    chunk_size=50,
):
    chunk_number = 1
    file_reader = pyfastaq.sequences.file_reader(to_align)
    names = []
    seqs = []
    f_out = pyfastaq.utils.open_file_write(outfile)
    ref_seq = utils.load_single_seq_fasta(ref_fa)
    print(f">{ref_name}", ref_seq.seq, sep="\n", file=f_out)
    for seq in file_reader:
        names.append(seq.id)
        seqs.append(seq.seq)
        if len(names) >= chunk_size:
            tmp_fa = f"{outfile}.tmp.{chunk_number}.fa"
            _write_fasta_chunk(tmp_fa, names, seqs)
            aln_seqs = run_mafft_one_fasta(
                tmp_fa, ref_fa, ref_name, quiet, indel_method, ref_start, ref_end
            )
            for name in names:
                print(f">{name}", aln_seqs[name], sep="\n", file=f_out)
            os.unlink(tmp_fa)
            chunk_number += 1
            names = []
            seqs = []

    if len(names) > 0:
        tmp_fa = f"{outfile}.tmp.{chunk_number}.fa"
        _write_fasta_chunk(tmp_fa, names, seqs)
        aln_seqs = run_mafft_one_fasta(
            tmp_fa, ref_fa, ref_name, quiet, indel_method, ref_start, ref_end
        )
        for name in names:
            print(f">{name}", aln_seqs[name], sep="\n", file=f_out)
        os.unlink(tmp_fa)

    pyfastaq.utils.close(f_out)
