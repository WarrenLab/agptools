from os.path import dirname, join
from unittest.mock import patch

import pytest
import screed

from agp import open_agp
from agp.agptools import main


def test_sanitize_help(capsys):
    with patch("sys.argv", ["agptools", "sanitize", "-h"]):
        with pytest.raises(SystemExit):
            main()

    out, err = capsys.readouterr()
    assert "path where sanitized contigs" in out


def test_sanitize_main(tmpdir):
    with patch(
        "sys.argv",
        [
            "agptools",
            "sanitize",
            join(dirname(__file__), "data", "test_sanitize.fa"),
            join(tmpdir, "test_out.fa"),
            join(dirname(__file__), "data", "test_sanitize.agp"),
            "-o",
            join(tmpdir, "test_out.agp"),
        ],
    ):
        main()

    test_out_agp = open_agp(join(tmpdir, "test_out.agp"))
    correct_out_agp = open_agp(join(dirname(__file__), "data", "sanitize_out.agp"))

    for row1, row2 in zip(test_out_agp, correct_out_agp):
        assert row1 == row2

    with screed.open(tmpdir / "test_out.fa") as test_out_fasta:
        with screed.open(
            join(dirname(__file__), "data", "sanitize_out.fa")
        ) as correct_out_fasta:
            for test_seq, correct_seq in zip(test_out_fasta, correct_out_fasta):
                assert test_seq == correct_seq
