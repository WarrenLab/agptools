from typing import Iterable, List, Mapping, Optional, Sequence, TextIO

from agp import AgpRow
from agp.flip import reverse_rows


class FormatError(Exception):
    pass


class ScaffoldNotFoundError(Exception):
    pass


def renaming_file_type(filename: str) -> dict[str, tuple[str, str]]:
    """Parse a renaming file

    A renaming file has 2-3 columns:
    1. current scaffold name
    2. desired new scaffold name
    3. orientation (+/-)

    Parse one of these files and return it in dictionary format

    Args:
        filename: path to renaming file

    Returns:
        the file in dict format: d[old_name] = (new_name, orientation)
    """
    renaming_dict = {}
    with open(filename, "r") as renaming_file:
        for i, line in enumerate(renaming_file):
            fields = line.strip().split("\t")
            if len(fields) < 2:
                raise FormatError(
                    f"Line {i+1} of renaming file doesn't have enough columns."
                )
            if len(fields) == 2:  # default to + orientation
                renaming_dict[fields[0]] = (fields[1], "+")
            elif len(fields) >= 3:
                if fields[2] not in ["+", "-"]:
                    raise FormatError(
                        f"Line {i+1}: orientation column can only contain + or -"
                    )
                renaming_dict[fields[0]] = (fields[1], fields[2])

    return renaming_dict


def rename_rows(
    rows_to_rename: Iterable[AgpRow], new_name, orientation
) -> List[AgpRow]:
    """Rename a bunch of agp rows

    Rename a bunch of agp rows, and possibly reverse them too.

    Args:
        rows_to_rename:
            all of the rows of a scaffold that is being renamed
        new_name:
            the new name for this scaffold
        orientation:
            either "+" or "-", the latter if reverse orienting the
            scaffold is desired

    Returns:
        the renamed and perhaps reverse oriented rows
    """
    new_rows = []
    for row in rows_to_rename:
        row.object = new_name
        new_rows.append(row)
    if orientation == "-":
        return reverse_rows(new_rows)
    return new_rows


def run(
    renaming_map: Mapping[str, tuple[str, str]],
    outfile: TextIO,
    agp_rows: Sequence[AgpRow],
):
    """Run the rename module

    Args:
        renaming_map:
            maps the name of a contig to be renamed to a tuple
            containing the new name and orientation
        outfile: file where agp output should be sent
        agp_rows: all input AGP rows
    """
    rows_to_rename: List[AgpRow] = []
    current_object = "none"  # can't use None because dict needs a string
    renamed_objects: set[str] = set()
    for row in agp_rows:
        if isinstance(row, str):
            print(row, file=outfile)
        elif row.object in renaming_map:
            if rows_to_rename and current_object != row.object:
                # we're on a new scaffold that does need to be renamed,
                # but we haven't renamed and printed the last scaffold
                # yet, so let's do that first
                for row_to_print in rename_rows(
                    rows_to_rename, *renaming_map[current_object]
                ):
                    print(row_to_print, file=outfile)
                rows_to_rename = []
                renamed_objects.add(current_object)
            current_object = row.object
            rows_to_rename.append(row)
        else:
            if rows_to_rename and current_object != row.object:
                # we're on a new scaffold that does not need to be
                # renamed, but the previous scaffold does need to
                # be renamed, so let's rename and print it out first
                for row_to_print in rename_rows(
                    rows_to_rename, *renaming_map[current_object]
                ):
                    print(row_to_print, file=outfile)
                rows_to_rename = []
                renamed_objects.add(current_object)
                current_object = "none"
            # now we can print the new row
            print(row, file=outfile)

    if rows_to_rename:
        for row_to_print in rename_rows(rows_to_rename, *renaming_map[row.object]):
            print(row_to_print, file=outfile)
        renamed_objects.add(current_object)

    remaining_objects = renaming_map.keys() - renamed_objects
    if remaining_objects:
        raise ScaffoldNotFoundError(remaining_objects)
