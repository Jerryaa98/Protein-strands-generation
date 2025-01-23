import argparse
import json
import Bio
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser, Superimposer
import numpy as np
import os
import matplotlib.pyplot as plt
from tqdm import tqdm


def calculate_rmsd(pdb_file1, pdb_file2):

    parser = PDBParser(QUIET=True)
    structure1 = parser.get_structure("struct1", pdb_file1)
    structure2 = parser.get_structure("struct2", pdb_file2)
    
    # Extract Cα atoms from the first (and only) chain
    def get_ca_atoms_single_chain(structure):
        
        model = next(structure.get_models())  # Get the first model
        chain = next(model.get_chains())  # Get the first chain
        
        return [residue['CA'] for residue in chain if 'CA' in residue]
    
    ca_atoms1 = get_ca_atoms_single_chain(structure1)
    
    ca_atoms2 = get_ca_atoms_single_chain(structure2)
    
    print(pdb_file2)
    print(len(ca_atoms1), len(ca_atoms2))

    if len(ca_atoms1) != len(ca_atoms2):
        min_ca_atoms = min(len(ca_atoms2), len(ca_atoms1))
        ca_atoms1 = ca_atoms1[:min_ca_atoms]
        ca_atoms2 = ca_atoms2[:min_ca_atoms]
        #print("Structures have different numbers of Cα atoms, skipping this one")
        #return
    
    # Superimpose the structures
    super_imposer = Superimposer()
    super_imposer.set_atoms(ca_atoms1, ca_atoms2)
    super_imposer.apply(structure2.get_atoms())  # Align structure2 to structure1
    
    # Calculate RMSD
    rmsd = super_imposer.rms
    return rmsd



def plot_ss_comp(id, work_dir, generated_sequences_id_dir):
    generated_pdb_file = os.listdir(generated_sequences_id_dir)
    og_pdb_fil_path = f'{work_dir}/pdb_files/{id}.pdb'
    pdb_file = [f'{generated_sequences_id_dir}/{f}' for f in generated_pdb_file]
    ss_similiarty = []
    for generated_pdb in pdb_file:
        ss_similiarty.append(calculate_rmsd(og_pdb_fil_path, generated_pdb))

    return ss_similiarty

def plot_seq_comp(generated_sequences):
    sequences_identity = []
    sequences_similarity = []
    og_seq = generated_sequences[0]
    print(len(og_seq))
    for generated_seq in generated_sequences[1:]:
        print(len(generated_seq))
    return 

def run(args):

    work_dir = args.work_dir
    strategy = args.strategy
    generated_sequences_dir = f"{work_dir}/generated_sequences/{strategy}"
    sequences_file = f'{generated_sequences_dir}/seq.json'
    results_dir = f'{work_dir}/comparison_results/{strategy}_comparison'

    
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    with open(sequences_file,'r') as json_file:
        seq_data = json.load(json_file)

    print('Starting generation comparsion ...')

    id = 'e2wjrA1'

    generated_sequences_id_dir = f'{generated_sequences_dir}/{id}' 
    ss_similarity = plot_ss_comp(id, work_dir, generated_sequences_id_dir)
    print(ss_similarity)
    plt.plot(range(len(ss_similarity)), ss_similarity)
    plt.savefig('/root/Biology_project/test.png')        
    print('Done Comparing')



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Strands Generation")
    # Add arguments
    parser.add_argument("--work_dir", type=str, default='/root/Biology_project', help='dirctory to save pdb file in')
    parser.add_argument("--strategy",type=str, default='sequential',
                        choices=['sequential', 'reverse', 'random'])    
    args = parser.parse_args()
    run(args)


"""

[

[6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], 


[36,37,38,39,40,41,42, 43, 44 ,45, 46], 


[51, 52, 53, 54, 55, 56, 57, 58,59, 60, 61, 62, 63, 64, 65, 66, 67], 



[71, 72, 73, 74, 75, 76, 77, 78 ,79, 80, 81, 82, 82,, 83, 84, 85, 86, 87, 88], 


[96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114], 


[130, 131, 132, 133, 134, 136, 137, 138, 139, 140 ,141, 142, 143], 


[148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167], 


[174, 175, 176, 177, 178, 179 ,180 ,181, 182, 183, 184, 185, 186, 187, 188, 189, 190]
]

"""