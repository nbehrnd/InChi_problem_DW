#!/usr/bin/env python3

# name:   obabel_inchi.py
# author: nbehrnd@yahoo.com
# date:   2021-09-09 (YYYY-MM-DD)
# edit:
"""Prepare check of OpenBabel InChi and InChiKey with DW.

This is to probe if DataWarrior and OpenBabel already assign different
InChi (and not only InChiKeys) for a set of random molecules generated
by DataWarrior.

By a call like

python 3 obabel_inchi.py Random_Molecules.sdf

DataWarrior's exported .sdf of multiple molecules generated as a library
of random molecules is read.  The script then relays to OpenBabel to
count chiral centres, and E/Z double bonds; to assign SMILES, InChi, and
InChiKey.  Including a header line, these are stored as tabulator
separted columns in file 'analysis.txt'.

Copy-paste the content of 'analysis.txt' into DW (special paste with
header line) and let DW compute SMILES, InChi, and InChiKey.  Either
use macro ob_dw_check.dwam for the analysis, or 1) manually instruct
DW for these assignments and 2) check the identity of these strings by
the Boolean comparison, i.e.

str(obSMILES) == str(SMILES)

which either yields "1" (one) as .true. (identity) .or. "0" (zero) as
.false. (strings do not match fully).

For a quick summary per class of molecules (defined by number of chiral
centres, number of double bonds, matching/dissimilar InChi by OpenBabel
and DataWarrior),
1) export DataWarrior's result as a text file (File -> Save Special ->
   Textfile)
2) submit this new file to Python script table_stat.py.  For the use
   of table_stat.py, see the header of table_stat.py.

Written for and initially tested with Python 3.9.2 (e.g., provided by
Debian 12/bookworm, branch testing) and OpenBabel 3.1.1."""

import argparse

from openbabel import openbabel as ob
from openbabel import pybel


def get_arguments():
    """Comparison of InChiKeys by DataWarrior with those by OpenBabel."""
    parser = argparse.ArgumentParser(description="""
    Comparison of InChiKeys by DataWarrior with those by OpenBabel.""")

    parser.add_argument("dwfile", help="DataWarrior's .sdf file.")
    args = parser.parse_args()

    data = args.dwfile
    return data


def testing(raw_data=""):
    """Inspect each molecule in the .sdf DW has written.

    Writes to the CLI the number of chiral centres, E/Z bonds,
    OpenBabel SMILES, OpenBabel InChi, and OpenBabel InChiKeys."""
    # derived from
    # https://open-babel.readthedocs.io/en/latest/Stereochemistry/stereo.html
    num_cistrans = 0
    num_tetra = 0

    mol = pybel.readstring("smi", str(raw_data))
    m = mol.OBMol

    for genericdata in m.GetAllData(ob.StereoData):
        stereodata = ob.toStereoBase(genericdata)
        stereotype = stereodata.GetType()

        if stereotype == ob.OBStereo.CisTrans:
            cistrans = ob.toCisTransStereo(stereodata)
            if cistrans.IsSpecified():
                num_cistrans += 1

        elif stereotype == ob.OBStereo.Tetrahedral:
            tetra = ob.toTetrahedralStereo(stereodata)
            if tetra.IsSpecified():
                num_tetra += 1

    # compute an InChi
    inchi = str(mol.write("inchi")).strip()

    # compute an InChiKey
    inchikey = str(mol.write("inchikey")).strip()

    keep = str("{}\t{}\t{}\t{}\t{}".format(num_tetra, num_cistrans,
                                           str(mol).split()[0], inchi,
                                           inchikey))

    print(keep)
    with open("analysis.txt", mode="a") as record:
        record.write("{}\n".format(keep))


def main():
    """Join the functions."""
    heading_line = str("obstereo\tobEZ\tobSMILES\tobInChi\tobInChiKey")
    data = get_arguments()
    mols = list(pybel.readfile("sdf", data))
    print(heading_line)

    with open("analysis.txt", mode="w") as newfile:
        newfile.write(f"{heading_line}\n")
    for mol in mols:
        testing(mol)

    print("\nFile 'analysis.txt' contains a permanent record.")

if __name__ == "__main__":
    main()
