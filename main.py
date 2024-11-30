import urllib.request
import csv

from biotite.structure.io.pdb import PDBFile
from biotite.application import dssp
from biotite.structure import get_residue_count
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio.PDB import PDBParser
from Bio import PDB
import os


def main(load_pdb):
    if(load_pdb == True):
        ids = []
        with open("/root/Biology_project/OMBB_data.csv", "r") as file:
            reader = csv.reader(file)
            for row in reader:
                ids.append(row[0])
        for id in ids:
            # URL of the file
            url = f"http://prodata.swmed.edu/ecod/af2_pdb/structure?id={id}"

            # File to save the downloaded content
            output_file = f"/root/Biology_project/pdb_files/{id}.pdb"

            # Download the file
            try:
                urllib.request.urlretrieve(url, output_file)
                print(f"File saved as {output_file}")
            except Exception as e:
                print(f"Failed to download file: {e}")

    # Ensure DSSP is available
    dssp_exe = "/usr/bin/mkdssp"  # Path to the DSSP executable (if it's not in your PATH, use full path here)
    if not os.path.isfile(dssp_exe):
        raise FileNotFoundError(f"DSSP executable not found: {dssp_exe}")

    # Load the PDB file
    pdb_filename = "/root/Biology_project/pdb_files/e1af6A1.pdb"  # Replace with your PDB file path
    parser = PDB.PPBuilder()
    structure = PDB.PDBParser(QUIET=True).get_structure('protein', pdb_filename)

    # Run DSSP to get secondary structure information
    dssp = PDB.DSSP(structure[0], pdb_filename)

    # Print secondary structure assignments for all residues
    print("Secondary structure assignments:")
    for i, (index, ss) in enumerate(dssp):
        residue = index[1]  # The residue
        secondary_structure = ss  # The DSSP secondary structure code
        print(f"Residue {residue}: {secondary_structure}")

if __name__ == "__main__":
    load_pdb = False
    main(load_pdb)
