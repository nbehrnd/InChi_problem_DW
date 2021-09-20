#!/usr/bin/env python3

# name:   table_stat.py
# author: nbehrnd@yahoo.com
# date:   2021-09-09 (YYYY-MM-DD)
# edit:
#
"""Write a sorted frequency count of similar/dissimilar InChi.

It looks like the probability of different InChi assigned by DataWarrior
and OpenBabel depends on the number of stereogenic centres and E/Z double
bonds.  Maybe there is some systematic pattern which may be extracted
in DataWarrior's exported list after running the DataWarrior macro (which
let DataWarrior compute SMILES, InChi, and InChiKeys and compare these
strings by their analogues assigned by OpenBabel).

After running, e.g.

python3 table_stat.py Random_Molecules_1000_check.txt

this script writes file frequency_list.tst of four tabulator separated
columns:
+ "stereo, E/Z, label diff_inchi:" labels the following one
+ a tuple which defines a type of molecules by the number of the chiral
  centres and E/Z double bonds count, and if their InChi assigned by
  OpenBabel and DataWarrior either is identical (1) or dissimilar (0)
+ "frequency:" labels the following one
+ the absolute number of molecules of this type

The report only includes types for molecules if there is at least one
matching the composite criterion "X Y 0" (or, "X Y 1").

Written for and initially tested in Linux Debian 12/bookworm (branch
testing); the script only relies on standard modules of Python 3.9.2.
"""
import argparse


def get_arguments():
    """Identify the input file to work with."""
    parser = argparse.ArgumentParser(description="""
    Yield a sorted frequency count of similar/dissimilar InChi DW/OpenBabel.
""")

    parser.add_argument(
        "file",
        help="DataWarrior's list file exported after running the macro.")
    args = parser.parse_args()

    data = args.file
    return data


def scrutiny(raw_data=""):
    """Analysis of DW's listing file."""
    # Per molecule in DW's listing file, extract count of stereo centres,
    # cistrans double bonds, and assigned diff_inchi label.
    survey = []
    with open(raw_data, mode="r") as source:
        for line in source:
            line = str(line).strip()
            data = line.split("\t")

            stereo = data[0]
            cistrans = data[1]
            diff_inchi = data[9]

            retain = str(f"{stereo} {cistrans} {diff_inchi}")
            survey.append(retain)

    # remove the header line:
    del survey[0]

    # for a frequency count, build a dictionary:
    counting = {}
    for instance in survey:
        counting.setdefault(instance, 0)
        counting[instance] = counting[instance] + 1

    # convert the dictionary into a list which may be sorted:
    listing = []
    for key, value in counting.items():
        retain = str(
            f"stereo, E/Z, label diff_inchi:\t{key}\tfrequency:\t{value}")
        listing.append(retain)
    listing.sort()

    # the eventual report:
    with open("frequency_list.txt", mode="w") as newfile:
        for element in listing:
            print(element)
            newfile.write(f"{element}\n")

    print("\nSee file 'frequency_list.txt' for a permanent record.")


def main():
    """Join the functions."""
    input_file = get_arguments()
    scrutiny(raw_data=input_file)


if __name__ == "__main__":
    main()
