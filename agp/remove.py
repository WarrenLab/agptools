from typing import IO, Sequence, Set

from agp import AgpRow


class ScaffoldNotFoundError(Exception):
    pass


def scaffolds_list_type(filename: str) -> Set[str]:
    """Parse a list of scaffolds

    Given a filename pointing to a list of scaffolds, one per line,
    parse it and return a list of scaffold strings.

    Args:
        filename: path to a file containing a list of scaffolds

    Returns:
        the set of scaffold names in the file
    """
    scaffolds: Set[str] = set()
    with open(filename) as scaffolds_list_file:
        for line in scaffolds_list_file:
            scaffolds.add(line.strip())
    return scaffolds


def run(scaffolds_to_remove: Set[str], outfile: IO, agp_in: Sequence[AgpRow]):
    scaffolds_removed = set()
    for row in agp_in:
        if isinstance(row, str) or row.object not in scaffolds_to_remove:
            print(row, file=outfile)
        else:
            scaffolds_removed.add(row.object)
    remaining_scaffolds = scaffolds_to_remove - scaffolds_removed
    if remaining_scaffolds:
        raise ScaffoldNotFoundError(str(remaining_scaffolds))
