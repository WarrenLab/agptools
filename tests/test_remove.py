from os.path import dirname, join
from unittest.mock import patch

import pytest

from agp import AgpRow, is_string, open_agp
from agp.agptools import main
from agp.remove import ScaffoldNotFoundError, run


def test_remove_help(capsys):
    with patch("sys.argv", ["agptools", "remove", "--help"]):
        with pytest.raises(SystemExit):
            main()

    out, err = capsys.readouterr()
    assert "list of scaffolds to remove, one per line" in out


def test_remove(tmpdir):
    with patch(
        "sys.argv",
        [
            "agptools",
            "remove",
            join(dirname(__file__), "data", "scaffolds_to_remove.txt"),
            join(dirname(__file__), "data", "test.agp"),
            "-o",
            join(tmpdir, "remove_out.agp"),
        ],
    ):
        main()

    for row in open_agp(join(tmpdir, "remove_out.agp")):
        if not is_string(row):
            assert row.object not in ["scaffold_18", "scaffold_20"]


def test_run_remove(tmp_path):
    agp_rows = [
        AgpRow("scaffold_18\t1930402\t2636956\t9\tW\ttig123\t1\t706555\t+"),
        AgpRow("scaffold_19\t2636957\t2637456\t10\tN\t500\tscaffold\tyes\tna"),
        AgpRow("scaffold_20\t2637457\t3193955\t11\tW\ttig456\t1\t556499\t+"),
    ]
    with open(tmp_path / "test.agp", "w") as outfile:
        run({"scaffold_19"}, outfile, agp_rows)
    out_agp_rows = list(open_agp(tmp_path / "test.agp"))
    assert len(out_agp_rows) == 2
    assert not any(row.object == "scaffold_19" for row in out_agp_rows)


def test_remove_scaffold_not_found_error(tmp_path):
    agp_rows = [
        AgpRow("scaffold_18\t1930402\t2636956\t9\tW\ttig123\t1\t706555\t+"),
        AgpRow("scaffold_19\t2636957\t2637456\t10\tN\t500\tscaffold\tyes\tna"),
        AgpRow("scaffold_20\t2637457\t3193955\t11\tW\ttig456\t1\t556499\t+"),
    ]
    with pytest.raises(ScaffoldNotFoundError) as error:
        with open(tmp_path / "test.agp", "w") as outfile:
            run({"scaffold_21"}, outfile, agp_rows)
    assert "scaffold_21" in str(error.value)
