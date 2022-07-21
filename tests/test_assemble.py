from os import devnull
from os.path import dirname, join
from unittest.mock import patch

import pytest
import screed

from agp import AgpRow
from agp.agptools import main
from agp.assemble import EmptyAgpError, NoSuchContigError, reverse_complement, run


def test_assemble(tmpdir):
    with patch(
        "sys.argv",
        [
            "agptools",
            "assemble",
            join(dirname(__file__), "data", "broken_contigs.fa"),
            join(dirname(__file__), "data", "test.agp"),
            "-o",
            join(tmpdir, "assemble_out.fa"),
        ],
    ):
        main()

    with screed.open(tmpdir / "assemble_out.fa") as test_out_fasta:
        with screed.open(
            join(dirname(__file__), "data", "assemble_out.fa")
        ) as assembled_fasta:
            for test_seq, correct_seq in zip(test_out_fasta, assembled_fasta):
                assert test_seq == correct_seq


def test_assemble_help(capsys):
    with patch("sys.argv", ["agptools", "assemble", "--help"]):
        with pytest.raises(SystemExit):
            main()

    out, err = capsys.readouterr()
    assert "AGP file assembling contigs into scaffolds" in out


def test_assemble_empty_agp():
    with patch(
        "sys.argv",
        [
            "agptools",
            "assemble",
            join(dirname(__file__), "data", "broken_contigs.fa"),
            join(dirname(__file__), "data", "empty.agp"),
        ],
    ):
        with pytest.raises(EmptyAgpError):
            main()


def test_reverse_complement():
    assert reverse_complement("ACTGcgTtCAGaaTTx") == "xAAttCTGaAcgCAGT"


def test_no_such_contig():
    agp_rows = [
        AgpRow("scaffold_18\t1930402\t2636956\t9\tW\ttig123\t1\t706555\t+"),
    ]
    with pytest.raises(NoSuchContigError) as err:
        with screed.open(
            join(dirname(__file__), "data", "broken_contigs.fa")
        ) as contigs_fasta:
            with open(devnull, "w") as outfile:
                run(contigs_fasta, outfile, agp_rows)

    assert "tig123" in str(err.value)
