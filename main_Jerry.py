import urllib.request
import csv
import os
import argparse
from biotite.structure.io.pdb import PDBFile
from biotite.application import dssp
from biotite.structure import get_residue_count
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio.PDB import PDBParser
from Bio import PDB

#
def load_pdb_files(data_file, output_dir):
    
    assert data_file.endswith('.csv'), 'DataSet file must be csv'

    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)  # Create the directory
        except Exception as e:
            raise ValueError(f"Error creating directory '{output_dir}': {e}")
        
    ids = []
    with open(data_file, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            ids.append(row[0])
    
    for id in ids[1:]:#skip header
        # URL of the file
        url = f"http://prodata.swmed.edu/ecod/af2_pdb/structure?id={id}"

        # File to save the downloaded content
        #output_file = f"/root/Biology_project/pdb_files/{id}.pdb"
        output_file = f"{output_dir}/{id}.pdb"
        # Download the file
        try:
            urllib.request.urlretrieve(url, output_file)
            print(f"File saved as {output_file}")
        except Exception as e:
            print(f"Failed to download file: {e}")


def run(args):

    data_file = args.data_file
    load_pdb = args.load_pdb
    output_dir = args.output_dir

    if(load_pdb):
        load_pdb_files(data_file, output_dir)



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
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Strands Generation")
    # Add arguments
    parser.add_argument("--data_file", type=str, help="path for the dataset")
    parser.add_argument("--load_pdb", type=bool, default=False, help='load pdp file using wget or not')
    parser.add_argument("--output_dir", type=str, default='/root/Biology_project/pdb_files', help='dirctory to save pdb file in')

    args = parser.parse_args()
    run(args)
