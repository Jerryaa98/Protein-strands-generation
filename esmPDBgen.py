
import argparse
from huggingface_hub import login
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import torch 
import csv
import os

def run(args):
    # Will instruct you how to get an API key from huggingface hub, make one with "Read" permission.
    login(token = args.hf_token)
    database_file = "/root/Biology_project/data/OMBB_data.csv"
    pdb_dir = '/root/Biology_project/ESM_pdb_files_v2'
    if not os.path.exists(pdb_dir):
        os.makedirs(pdb_dir)
        
    model: ESM3InferenceClient = ESM3.from_pretrained("esm3-open").to("cuda") 
    with open(database_file, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header
        for row in reader:
            seq_id = row[0]
            seq_len = int(row[-1])
            seq = row[2]
            assert len(seq) == seq_len, 'Erorr'
            pdb_file = f'{pdb_dir}/{seq_id}.pdb'
            protein = ESMProtein(sequence=seq)
            num_steps = 1000
            #protein = model.generate(protein, GenerationConfig(track="sequence", num_steps=num_steps, temperature=0.7))
            protein = model.generate(protein, GenerationConfig(track="structure", num_steps=num_steps))
            protein.to_pdb(pdb_file)
            
     
if __name__ == "__main__":

    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Strands Generation")
    # Add arguments
    parser.add_argument("--hf_token", type=str, default="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")

    args = parser.parse_args()
    run(args)
