import urllib.request
import csv
import os
import argparse
from Bio.PDB import PDBParser, DSSP

# https://files.rcsb.org/download/1AF6.pdb
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
        domain = id[:4]
        url = f"https://files.rcsb.org/download/{domain}.pdb"

        # File to save the downloaded content
        #output_file = f"/root/Biology_project/pdb_files/{id}.pdb"
        output_file = f"{output_dir}/{domain}.pdb"
        # Download the file
        try:
            urllib.request.urlretrieve(url, output_file)
            print(f"File saved as {output_file}")
        except Exception as e:
            print(f"Failed to download file: {e}")

def group_residuals(beta_indcies):
    current_strand = [beta_indcies[0]]
    strands = []

    for i in beta_indcies[1:]:
        if i == current_strand[-1] + 1:
            current_strand.append(i)
        else:
            if len(current_strand) > 0 :
                strands.append(current_strand)
            current_strand = [i]
    if len(current_strand) > 0:
        strands.append(current_strand)
    return strands
    
def run():
 
    # Load the PDB file
    pdb_filename = "/root/Biology_project/pdb_files/1af6.pdb"  # Replace with your PDB file path
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_filename)
   
    # Use DSSP to calculate secondary structure
    model = structure[0]  # Select the first model (if multiple models exist)
    dssp = DSSP(model, pdb_filename)
    
    # Filter beta-strands (E) to identify potential beta-barrel regions
    beta_strands = []
    for chain in dssp.keys():
        residue, sec_struct, _ = chain, dssp[chain][2], dssp[chain][3]
        if sec_struct == 'E' or sec_struct == 'B':  # 'E' corresponds to beta-strands
            beta_strands.append(residue)
  
    for strand in beta_strands:
        print(strand)
    return

    indcies = [residue[1][1] for residue in beta_strands if residue[0] == 'B']

    strands = group_residuals(indcies)

    print(len(strands))


if __name__ == "__main__":

    run()
