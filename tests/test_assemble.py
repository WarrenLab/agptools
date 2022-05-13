from os import devnull
from os.path import dirname, join

import pytest
import screed

from agp import AgpRow, open_agp
from agp.assemble import NoSuchContigError, reverse_complement, run


def test_assemble(tmpdir):
    agp_rows = open_agp(join(dirname(__file__), "data", "test.agp"))
    with screed.open(
        join(dirname(__file__), "data", "broken_contigs.fa")
    ) as contigs_fasta:
        with open(tmpdir / "assemble_out.fa", "w") as outfile:
            run(contigs_fasta, outfile, agp_rows)

    with screed.open(tmpdir / "assemble_out.fa") as test_out_fasta:
        with screed.open(
            join(dirname(__file__), "data", "assemble_out.fa")
        ) as assembled_fasta:
            for test_seq, correct_seq in zip(test_out_fasta, assembled_fasta):
                assert test_seq == correct_seq


def test_reverse_complement():
    assert reverse_complement("ACTGcgTtCAGaaTT") == "AAttCTGaAcgCAGT"


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
