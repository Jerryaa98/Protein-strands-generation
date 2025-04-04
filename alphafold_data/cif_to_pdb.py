from Bio.PDB import MMCIFParser, PDBIO
import os

for idx in range(5):
    # Set your input and output directories
    input_dir = f'/root/Biology_project/alphafold_data/model_{idx}'
    output_dir = f'/root/Biology_project/alphafold_data/model_{idx}/model_{idx}_pdb_conversions'
    os.makedirs(output_dir, exist_ok=True)

    parser = MMCIFParser(QUIET=True)
    io = PDBIO()

    # Loop through all .cif files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith('.cif'):
            cif_path = os.path.join(input_dir, filename)
            pdb_id = os.path.splitext(filename)[0]  # Remove ".cif"
            pdb_path = os.path.join(output_dir, pdb_id + '.pdb')

            try:
                structure = parser.get_structure(pdb_id, cif_path)
                io.set_structure(structure)
                io.save(pdb_path)
                print(f"✅ Converted Model {idx}: {filename} → {pdb_id}.pdb")
            except Exception as e:
                print(f"❌ Failed to convert Model {idx} {filename}: {e}")
