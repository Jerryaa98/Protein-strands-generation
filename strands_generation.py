
import argparse
import json
from huggingface_hub import login
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig
import csv
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import os


def run(args):
    # Will instruct you how to get an API key from huggingface hub, make one with "Read" permission.
    login(token = args.hf_token)
    model: ESM3InferenceClient = ESM3.from_pretrained("esm3-open").to("cuda") 
    strands_file = args.strands_file
    work_dir = args.work_dir

    generated_sequences_dir = f"{work_dir}/generated_sequences"
    if not os.path.exists(generated_sequences_dir):
        try:
            os.makedirs(generated_sequences_dir)  # Create the directory
        except Exception as e:
            raise ValueError(f"Error creating directory '{generated_sequences_dir}': {e}")
        
    with open(strands_file, 'r') as json_file:
        data = json.load(json_file)

    generated_data = {}
    for id in data.keys():
        seq = data[id]['seq']
        strands_indcies = data[id]['strands']
        state = []
        for i,strand in enumerate(strands_indcies):
            # mask
            masked_seq = ""
            for j, r in enumerate(seq):
                if j+1 in strand:
                    masked_seq += '_'
                else:
                    masked_seq += r

            protein = ESMProtein(sequence=masked_seq)

            # Generate and save seq
            num_steps = max(1, len(strand)//2)
            protein = model.generate(protein, GenerationConfig(track="sequence", num_steps=num_steps, temperature=0.7))
            output_seq = protein.sequence
            state.append(output_seq)
            # Generate and save pdb files 
            protein = model.generate(protein, GenerationConfig(track="structure", num_steps=num_steps))
            generated_sequences = f"{generated_sequences_dir}/{id}"
            # generated seq after each masking foreach seq (id)
            if not os.path.exists(generated_sequences):
                try:
                    os.makedirs(generated_sequences)  # Create the directory
                except Exception as e:
                    raise ValueError(f"Error creating directory '{generated_sequences}': {e}")
                
            protein.to_pdb(f"{generated_sequences}/generation{i}.pdb")

        generated_data[id] = state
        break

    with open(f'{generated_sequences_dir}/seq.json','w') as json_file:
        json.dump(generated_data, json_file)


if __name__ == "__main__":

    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Strands Generation")
    # Add arguments
    parser.add_argument("--strands_file", type=str, help='strands file json file')
    parser.add_argument("--work_dir", type=str, default='/root/Biology_project', help='dirctory to save pdb file in')
    parser.add_argument("--hf_token", type=str, default="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")
    args = parser.parse_args()
    run(args)
