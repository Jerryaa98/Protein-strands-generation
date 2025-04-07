import os
import csv 
from Bio.PDB import PDBParser, Superimposer
import matplotlib.pyplot as plt
import numpy as np

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
        #print("Structures have different numbers of Cα atoms, skipping this one")
        return -1
    
    # Superimpose the structures
    super_imposer = Superimposer()
    super_imposer.set_atoms(ca_atoms1, ca_atoms2)
    super_imposer.apply(structure2.get_atoms())  # Align structure2 to structure1
    
    # Calculate RMSD
    rmsd = super_imposer.rms
    return rmsd


if __name__ == "__main__":
    ECOD_PDB_dir = '/root/Biology_project/pdb_files'
    ESM_PDB_dir = '/root/Biology_project/ESM_pdb_files'
    ALPHA_PDB_dir = '/root/Biology_project/Alpha_Fold_model_0'
    database_file = "/root/Biology_project/data/OMBB_data.csv"
    comp_dict_esm = {}
    comp_dict_alpha_fold = {}
    with open(database_file, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header
        for row in reader:
            seq_id = row[0]
            ECOD_pdb_file = f'{ECOD_PDB_dir}/{seq_id}.pdb'
            ESM_pdb_file = f'{ESM_PDB_dir}/{seq_id}.pdb'
            ALPHA_pdb_file = f'{ALPHA_PDB_dir}/{seq_id}.pdb'
            if not os.path.exists(ALPHA_pdb_file):
                rmsd_alpha_fold = -1
            else:   
                rmsd_alpha_fold = calculate_rmsd(ECOD_pdb_file, ALPHA_pdb_file)
            if not os.path.exists(ESM_pdb_file):
                rmsd_esm = -1
            else:
                rmsd_esm = calculate_rmsd(ECOD_pdb_file, ESM_pdb_file)
            comp_dict_alpha_fold[seq_id] = rmsd_alpha_fold
            comp_dict_esm[seq_id] = rmsd_esm


    x = np.arange(len(comp_dict_esm.keys()))  # positions for Y-axis (sequence IDs)
    width = 0.4  # width of the bars

    # Prepare values
    esm_values = list(comp_dict_esm.values())
    alpha_fold_values = list(comp_dict_alpha_fold.values())
    labels = list(comp_dict_esm.keys())

    fig, ax = plt.subplots(figsize=(12, 16))  # taller figure for vertical sequence IDs
    mean_esm = np.mean([x for x in comp_dict_esm.values() if x > -1])
    mean_alpha_fold = np.mean([x for x in comp_dict_alpha_fold.values() if x > -1])
    # Add mean ± std as text inside the plot
    mean_text = f"AVG RMSD (ESM): {mean_esm:.2f}\nAVG RMSD (AlphaFold): {mean_alpha_fold:.2f}\n"
    ax.text(0.95, 0.05, mean_text, transform=ax.transAxes, fontsize=14,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='black'))

    # Plot horizontal bars, shift by width to avoid overlap
    bars1 = ax.barh(x - width/2, esm_values, height=width, label='ESM RMSD', color='b', alpha=0.7)
    bars2 = ax.barh(x + width/2, alpha_fold_values, height=width, label='AlphaFold RMSD', color='r', alpha=0.7)

    # Labels, title, and legend
    ax.set_xlabel('RMSD')
    ax.set_ylabel('Seq ID')
    ax.set_title('Comparison of RMSD for Each Seq ID')

    ax.set_yticks(x)
    ax.set_yticklabels(labels)

    ax.legend()
    plt.tight_layout()
    plt.savefig('/root/Biology_project/Protein-strands-generation/images/COMP.png')
    

