import os
import pytest

from ushonium import mafft, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "mafft")


@pytest.fixture(scope="session")
def test_data():
    outdir = "tmp.mafft_test_data"
    utils.syscall(f"rm -rf {outdir}")
    os.mkdir(outdir)

    start = "ACGCGAGCAGCGGTCAG"
    middle1 = "TCAGTCGGGC"
    middle2 = "AGACGCAAAC"
    end = "TGTGAGACGTCACACCATACACTGAACACACACGTACCGAGCTACGACGTCGTA"
    ref = start + middle1 + "TAT" + middle2 + "A" + end
    aln = start + "GT" + middle1 + middle2 + "T" + end
    expect_indel = start + middle1 + "---" + middle2 + "T" + end
    expect_indel_ref = start + middle1 + "TAT" + middle2 + "T" + end
    expect_indel_N = start + middle1 + "NNN" + middle2 + "T" + end

    data = {
        "ref_fa": os.path.join(outdir, "ref.fa"),
        "to_align_fa": os.path.join(outdir, "to_align.fa"),
        "expect_ignore_indel": expect_indel.lower(),
        "expect_indel_ref": expect_indel_ref.lower(),
        "expect_indel_N": expect_indel_N.lower(),
    }

    with open(data["ref_fa"], "w") as f:
        print(">ref", ref, sep="\n", file=f)
    with open(data["to_align_fa"], "w") as f:
        print(">to_align", aln, sep="\n", file=f)

    yield data
    utils.syscall(f"rm -rf {outdir}")


def test_mafft_stdout_to_seqs():
    mafft_stdout = "\n".join(
        [
            ">ref",
            "cag\nctaac",
            ">aln",
            "cat\nct--c",
        ]
    )
    got_ref, got_aln = mafft.mafft_stdout_to_seqs(mafft_stdout, "ref")
    assert got_ref == "cagctaac"
    assert got_aln == "catct--c"


def test_run_mafft_one_seq(test_data):
    name = "name"
    got_name, got_seq = mafft.run_mafft_one_seq(
        name,
        test_data["to_align_fa"],
        test_data["ref_fa"],
        "ref",
        False,
        "nothing",
        None,
        None,
    )
    assert got_name == name
    assert got_seq == test_data["expect_ignore_indel"]

    got_name, got_seq = mafft.run_mafft_one_seq(
        name,
        test_data["to_align_fa"],
        test_data["ref_fa"],
        "ref",
        False,
        "as_ref",
        None,
        None,
    )
    assert got_name == name
    assert got_seq == test_data["expect_indel_ref"]

    got_name, got_seq = mafft.run_mafft_one_seq(
        name,
        test_data["to_align_fa"],
        test_data["ref_fa"],
        "ref",
        False,
        "N",
        None,
        None,
    )
    assert got_name == name
    assert got_seq == test_data["expect_indel_N"]


def test_force_ref_at_ends():
    # refcoord 01 2345678 90 123
    # gapcoord 01234567890123456
    ref_seq = "AC-GTACGTG-GT-ATA"
    aln_seq = "-bcd-efghijklmn--"

    got_ref, got_aln = mafft.force_ref_at_ends(ref_seq, aln_seq, 0, 13)
    assert got_ref == "AC-GTACGTG-GT-ATA"
    assert got_aln == "-bcd-efghijklmn--"

    got_ref, got_aln = mafft.force_ref_at_ends(ref_seq, aln_seq, 1, 13)
    assert got_ref == "AC-GTACGTG-GT-ATA"
    assert got_aln == "Abcd-efghijklmn--"

    got_ref, got_aln = mafft.force_ref_at_ends(ref_seq, aln_seq, 2, 12)
    assert got_ref == "ACGTACGTG-GT-ATA"
    assert got_aln == "ACd-efghijklmn-A"

    got_ref, got_aln = mafft.force_ref_at_ends(ref_seq, aln_seq, 3, 11)
    assert got_ref == "ACGTACGTG-GT-ATA"
    assert got_aln == "ACG-efghijklmnTA"

    got_ref, got_aln = mafft.force_ref_at_ends(ref_seq, aln_seq, 4, 10)
    assert got_ref == "ACGTACGTG-GTATA"
    assert got_aln == "ACGTefghijklATA"

    got_ref, got_aln = mafft.force_ref_at_ends(ref_seq, aln_seq, 5, 9)
    assert got_ref == "ACGTACGTG-GTATA"
    assert got_aln == "ACGTAfghijkTATA"

    got_ref, got_aln = mafft.force_ref_at_ends(ref_seq, aln_seq, 6, 8)
    assert got_ref == "ACGTACGTGGTATA"
    assert got_aln == "ACGTACghiGTATA"

    got_ref, got_aln = mafft.force_ref_at_ends(ref_seq, aln_seq, 7, 8)
    assert got_ref == "ACGTACGTGGTATA"
    assert got_aln == "ACGTACGhiGTATA"

    ref_seq = "--ACGT-A--"
    aln_seq = "abcdefghi"
    got_ref, got_aln = mafft.force_ref_at_ends(ref_seq, aln_seq, 0, 4)
    assert got_ref == "ACGT-A"
    assert got_aln == "cdefgh"

    got_ref, got_aln = mafft.force_ref_at_ends(ref_seq, aln_seq, 1, 3)
    assert got_ref == "ACGTA"
    assert got_aln == "AdefA"
