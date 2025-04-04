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
    comp_dict = {}
    with open(database_file, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header
        for row in reader:
            seq_id = row[0]
            ECOD_pdb_file = f'{ECOD_PDB_dir}/{seq_id}.pdb'
            ESM_pdb_file = f'{ESM_PDB_dir}/{seq_id}.pdb'
            ALPHA_pdb_file = f'{ALPHA_PDB_dir}/{seq_id}.pdb'
            if not os.path.exists(ALPHA_pdb_file):
                print(seq_id)
                continue
            rmsd_alpha = calculate_rmsd(ESM_pdb_file, ALPHA_pdb_file)
            #rmsd_esm = calculate_rmsd(ECOD_pdb_file, ESM_pdb_file)
            #print(f'{seq_id} | ESM : {rmsd_esm} | ALPHA FOLD : {rmsd_alpha}')
            comp_dict[seq_id] = rmsd_alpha





    x = np.arange(len(comp_dict.keys()))  # X-axis positions
    width = 0.4  # Width of bars

    # Compute mean and standard deviation
    mean = np.mean([x for x in comp_dict.values() if x > -1])
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(16, 9))

    # Add mean ± std as text inside the plot
    mean_text = f"MEAN: {mean:.2f} "
    ax.text(0.95, 0.05, mean_text, transform=ax.transAxes, fontsize=12,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='black'))

    # Labels, title, and legend
    bars1 = ax.bar(x - width/2, list(comp_dict.values()), width, label='RMSD', color='b', alpha=0.7)
    ax.set_xlabel('Seq ID')
    ax.set_ylabel('RMSD')
    ax.set_title('Comparison of RMSD for Each Seq ID ')
    ax.set_xticks(x)
    ax.set_xticklabels(comp_dict.keys(), rotation=90, ha='right')
    ax.legend()

    # Show plot
    plt.tight_layout()
    plt.savefig('/root/Biology_project/Protein-strands-generation/images/AlphaFold_ESM_COMP.png')
    

