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
        name, test_data["to_align_fa"], test_data["ref_fa"], "ref", False, "nothing"
    )
    assert got_name == name
    assert got_seq == test_data["expect_ignore_indel"]

    got_name, got_seq = mafft.run_mafft_one_seq(
        name, test_data["to_align_fa"], test_data["ref_fa"], "ref", False, "as_ref"
    )
    assert got_name == name
    assert got_seq == test_data["expect_indel_ref"]

    got_name, got_seq = mafft.run_mafft_one_seq(
        name, test_data["to_align_fa"], test_data["ref_fa"], "ref", False, "N"
    )
    assert got_name == name
    assert got_seq == test_data["expect_indel_N"]
