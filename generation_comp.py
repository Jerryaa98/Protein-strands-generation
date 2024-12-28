import argparse
import json
import Bio
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser, Superimposer
import numpy as np
import os
import matplotlib.pyplot as plt
from tqdm import tqdm

def calculate_seq_similarity(seq1, seq2, matrix = substitution_matrices.load('BLOSUM62')):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length for this method.")
    
    similarity_score = 0
    valid_pairs = 0

    # Iterate through the sequence pairs
    for res1, res2 in zip(seq1, seq2):
        pair = (res1, res2)
        reverse_pair = (res2, res1)  # For cases where the matrix is symmetric

        # Check if the pair exists in the substitution matrix
        if pair in matrix:
            similarity_score += matrix[pair]
            valid_pairs += 1
        elif reverse_pair in matrix:
            similarity_score += matrix[reverse_pair]
            valid_pairs += 1
    
    # Normalize the similarity score by the length of the sequence
    if valid_pairs == 0:
        return 0  # Avoid division by zero
    normalized_similarity = similarity_score / valid_pairs
    return normalized_similarity

def calculate_seq_identity(seq1, seq2):
    # Calculate identity
    matches = sum(res1 == res2 for res1, res2 in zip(seq1, seq2))
    identity = (matches / max(len(seq1), len(seq2))) * 100
    return identity

def plot_seq_comp(generated_sequences):
    sequences_identity = []
    sequences_similarity = []
    for before_seq, after_seq in zip(generated_sequences[1:],generated_sequences[:-1]):
        sequences_identity.append(calculate_seq_identity(before_seq, after_seq))
        sequences_similarity.append(calculate_seq_similarity(before_seq, after_seq))
    
    return sequences_identity, sequences_similarity

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
    
    if len(ca_atoms1) != len(ca_atoms2):
        print("Structures have different numbers of Cα atoms, skipping this one")
        return
    
    # Superimpose the structures
    super_imposer = Superimposer()
    super_imposer.set_atoms(ca_atoms1, ca_atoms2)
    super_imposer.apply(structure2.get_atoms())  # Align structure2 to structure1
    
    # Calculate RMSD
    rmsd = super_imposer.rms
    return rmsd

def plot_ss_comp(id, work_dir, generated_sequences_id_dir):
    generated_pdb_file = sorted(os.listdir(generated_sequences_id_dir))
    og_pdb_fil_path = f'{work_dir}/pdb_files/{id}.pdb'
    pdb_file = [og_pdb_fil_path] + [f'{generated_sequences_id_dir}/{f}' for f in generated_pdb_file]
    ss_similiarty = []
    for before_pdb, after_pdb in zip(pdb_file[1:],pdb_file[:-1]):
        ss_similiarty.append(calculate_rmsd(before_pdb, after_pdb))

    return ss_similiarty

def plot_results(results_dir, id, sequences_identity, sequences_similarity, ss_similarity):
    x = range(1,len(sequences_similarity)+1)

    fig, axes = plt.subplots(1,3, figsize=(15,9))
    axes[0].plot(x, sequences_identity)
    axes[0].set_xlabel('Iteration Num')
    axes[0].set_ylabel('Sequence Identity (%)')
    axes[0].set_title(id)

    axes[1].plot(x, sequences_similarity)
    axes[1].set_xlabel('Iteration Num')
    axes[1].set_ylabel('Sequence Similarity (norm)')
    axes[1].set_title(id)

    axes[2].plot(x, ss_similarity)
    axes[2].set_xlabel('Iteration Num')
    axes[2].set_ylabel('Structural similarity')
    axes[2].set_title(id)

    plt.savefig(f'{results_dir}/{id}.png')
    plt.close()
    return 

def run(args):
    work_dir = args.work_dir
    generated_sequences_dir = f"{work_dir}/generated_sequences"
    sequences_file = f'{generated_sequences_dir}/seq.json'
    results_dir = f'{work_dir}/results'

    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    with open(sequences_file,'r') as json_file:
        seq_data = json.load(json_file)

    for id in tqdm(seq_data.keys()):
        generated_sequences = seq_data[id]
        sequences_identity, sequences_similarity = plot_seq_comp(generated_sequences)
        generated_sequences_id_dir = f'{generated_sequences_dir}/{id}' 
        ss_similarity = plot_ss_comp(id, work_dir, generated_sequences_id_dir)

        plot_results(results_dir, id, sequences_identity, sequences_similarity, ss_similarity)

    
if __name__ == "__main__":
     # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Strands Generation")
    # Add arguments
    parser.add_argument("--work_dir", type=str, default='/root/Biology_project', help='dirctory to save pdb file in')
    args = parser.parse_args()
    run(args)