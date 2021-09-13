#!/bin/usr/env python3

# name:   rdkit_inchi.py
# author: nbehrnd@yahoo.com
# date:   2021-09-13 (YYYY-MM-DD)
# edit:
#
"""Report SMILES, InChi, and InChiKey for DW's .sdf.

Complementary to OpenBabel, a check with RDKit seems useful.  Starting
with DW's multi-molecule .sdf from the generation of random library of
drug-like molecules, run the following in the pattern of

python3 rdkit_inchi.py Random_Molecules_1000.sdf

to generate a listing of SMILES, InChi, and InChiKey as assigned by
RDKit.

Written for and initially tested on Linux Debian 12/bookworm (branch
testing).  Note: Python interpreter 3.9.7 and rdkit 2021.03.5 used
are now provided within a virtual environment of Miniconda."""
import argparse

import rdkit
from rdkit import Chem


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

    Writes RDKit's SMILES, InChi and InChiKey to CLI and report file."""

    with open("RDKit_analysis.txt", mode="w") as newfile:
        newfile.write("tetra\tez\trdSMILES\trdInChI\trdInChIKey\n")
        # inspired by
        # https://www.rdkit.org/docs/GettingStartedInPython.html#reading-sets-of-molecules
        molecules = Chem.SDMolSupplier(raw_data)
        for molecule in molecules:
            rd_smiles = str(Chem.MolToSmiles(molecule))
            rd_inchi = str(Chem.MolToInchi(molecule))
            rd_inchikey = str(Chem.MolToInchiKey(molecule))

            # Count of stereogenic centres, in reference to
            # https://www.rdkit.org/docs/source/rdkit.Chem.html?highlight=chiral#rdkit.Chem.FindMolChiralCenters
            tetra = Chem.FindMolChiralCenters(Chem.MolFromSmiles(rd_smiles),
                                              includeUnassigned=True,
                                              useLegacyImplementation=False)
            num_tetra = len(tetra)

            # Count of E/Z doublebonds not included by cyclic rings, e.g.
            # aiming to assign to cyclohexene = 0, and stilbene = 1.
            # Reference to
            # https://www.rdkit.org/docs/Cookbook.html?highlight=chem%20findpotentialstereo
            cistrans = Chem.FindPotentialStereo(molecule)
            num_cistrans = len(cistrans)

            retain = str(f"{num_tetra}\t{num_cistrans}\t{rd_smiles}\t{rd_inchi}\t{rd_inchikey}\n")
            newfile.write(retain)


def main():
    """Join the functions."""
    data = get_arguments()
    testing(data)

    print("\nFile 'RDKit_analysis.txt' contains a permanent record.")


if __name__ == "__main__":
    main()
