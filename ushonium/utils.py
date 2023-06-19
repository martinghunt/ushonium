import logging
import os
import subprocess
import pyfastaq


def syscall(command, cwd=None, quiet=False, log=None):
    if not quiet:
        logging.info(f"Run: {command}")
    p = subprocess.run(
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        cwd=cwd,
    )
    if log is not None:
        with open(log, "w") as f:
            print(p.stdout, file=f)
    if p.returncode != 0:
        raise Exception(
            f"Command failed (exit code {p.returncode}): {command}\n\n{p.stdout}\n"
        )
    return p


def load_single_seq_fasta(filename):
    seqs = {}
    pyfastaq.tasks.file_to_dict(filename, seqs)
    assert len(seqs) == 1
    return list(seqs.values())[0]


def to_abs_path(filename):
    if filename is None:
        return None
    else:
        return os.path.abspath(filename)
