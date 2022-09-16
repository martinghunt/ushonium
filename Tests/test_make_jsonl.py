#!/usr/bin/env python3

import gzip
import json
import os
import pytest
import subprocess


this_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_dir, "test_make_jsonl_data")
scripts_dir = os.path.join(this_dir, os.pardir, "Scripts")


def check_jsonl_contents(infile):
    with gzip.open(infile, "rb") as f:
        lines = [json.loads(x.rstrip()) for x in f]
    assert len(lines) == 7
    names = {x["name"] for x in lines if "name" in x}
    assert len(names) == 6  # expected samples plus nodes to make tree
    assert "sample_1" in names
    assert "sample_2" in names
    assert "sample_3" in names
    return lines


def test_make_jsonl():
    outdir = "tmp.test_make_jsonl"
    subprocess.check_output(f"rm -rf {outdir}", shell=True)
    samples_tsv = "tmp.test_make_jsonl.samples.tsv"
    with open(samples_tsv, "w") as f:
        for i in range(1, 4):
            print(
                f"sample_{i}", os.path.join(data_dir, f"sample{i}.fa"), sep="\t", file=f
            )

    make_jsonl_py = os.path.join(scripts_dir, "make_jsonl.py")
    command = f"{make_jsonl_py} {samples_tsv} {outdir}"
    process = subprocess.run(command, shell=True)
    assert process.returncode == 0
    expect_jsonl = os.path.join(outdir, "05.taxonium.jsonl.gz")
    assert os.path.exists(expect_jsonl)
    check_jsonl_contents(expect_jsonl)

    # rerun with 2 cpus and using metadata file
    subprocess.check_output(f"rm -r {outdir}", shell=True)
    meta_tsv = os.path.join(data_dir, "metadata.tsv")
    command = f"{make_jsonl_py} --cpus 2 --metadata_tsv {meta_tsv} --metacols col1 {samples_tsv} {outdir}"
    process = subprocess.run(command, shell=True)
    assert process.returncode == 0
    expect_jsonl = os.path.join(outdir, "05.taxonium.jsonl.gz")
    assert os.path.exists(expect_jsonl)
    json_dicts = check_jsonl_contents(expect_jsonl)
    # check that we got column1 from metadata for sample 1
    for d in json_dicts:
        if d.get("name", None) == "sample1":
            assert d["meta_col1"] == "s1.1"

    subprocess.check_output(f"rm -r {outdir}", shell=True)
    os.unlink(samples_tsv)
