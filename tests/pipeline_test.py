import gzip
import json
import os
import pytest
from unittest import mock

from ushonium import pipeline, utils

this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "data", "pipeline")


def check_jsonl_contents(infile):
    with gzip.open(infile, "rb") as f:
        lines = [json.loads(x.rstrip()) for x in f]
    assert len(lines) == 6
    names = {x["name"] for x in lines if "name" in x}
    assert len(names) == 5  # expected samples plus nodes to make tree
    assert "sample_1" in names
    assert "sample_2" in names
    assert "sample_3" in names
    return lines


def test_pipeline():
    options = mock.Mock()
    options.ref_gb = pipeline.REF_GB
    options.metadata_tsv = None
    options.cpus = 1
    options.outdir = "tmp.test_pipeline"
    options.matopt_opts = "-m 0.000000001 -r 8"
    options.indel_method = "nothing"
    utils.syscall(f"rm -rf {options.outdir}")
    options.samples_tsv = "tmp.test_pipeline.samples.tsv"
    options.ref_start_end = None
    with open(options.samples_tsv, "w") as f:
        for i in range(1, 4):
            print(
                f"sample_{i}", os.path.join(data_dir, f"sample{i}.fa"), sep="\t", file=f
            )

    pipeline.run(options)
    expect_jsonl = os.path.join(options.outdir, "05.taxonium.jsonl.gz")
    assert os.path.exists(expect_jsonl)
    check_jsonl_contents(expect_jsonl)

    # rerun with 2 cpus and using metadata file
    utils.syscall(f"rm -r {options.outdir}")
    options.metadata_tsv = os.path.join(data_dir, "metadata.tsv")
    options.metacols = "col1"
    options.title = "TITLE"
    options.metacol_name = None
    pipeline.run(options)
    assert os.path.exists(expect_jsonl)
    json_dicts = check_jsonl_contents(expect_jsonl)
    # check that we got column1 from metadata for sample 1
    for d in json_dicts:
        if d.get("name", None) == "sample1":
            assert d["meta_col1"] == "s1.1"

    utils.syscall(f"rm -r {options.outdir}")
    os.unlink(options.samples_tsv)
