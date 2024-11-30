import urllib.request
import csv
import os
import argparse
import mdtraj as md
import numpy as np
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
    #skip header
    ids = ids[1:]

    for id in ids:
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

    ids = []
    n_strands = []
    with open(data_file, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            ids.append(row[0])
            n_strands.append(row[1])

    # skip header
    ids = ids[1:]
    n_strands = n_strands[1:]

    for i, id in enumerate(ids):
        # Load PDB file
        pdb_file = f"{output_dir}/{id}.pdb"  # Replace with your PDB file path

        # Load the PDB file
        traj = md.load(pdb_file)

        # Compute DSSP to get secondary structure
        dssp = md.compute_dssp(traj)

        # Initialize a list to hold the beta strands as groups of residues
        beta_strands = []
        current_strand = []

        # Iterate over the secondary structure assignments
        for j, ss in enumerate(dssp[0]):
            if ss == 'E' or ss == 'B':  # 'B' for beta strand, 'E' for extended beta strand
                current_strand.append(traj.topology.atom(j).residue.index)
            else:
                if current_strand:  # When we reach a non-beta strand, finalize the current strand
                    beta_strands.append(current_strand)
                    current_strand = []  # Reset the current strand

        # If the last strand ends at the end of the chain, add it to the list
        if current_strand:
            beta_strands.append(current_strand)

        print(f'{i+2} : {len(beta_strands)} | {n_strands[i]}')





if __name__ == "__main__":
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Strands Generation")
    # Add arguments
    parser.add_argument("--data_file", type=str, help="path for the dataset")
    parser.add_argument("--load_pdb", type=bool, default=False, help='load pdp file using wget or not')
    parser.add_argument("--output_dir", type=str, default='/root/Biology_project/pdb_files', help='dirctory to save pdb file in')

    args = parser.parse_args()
    run(args)
