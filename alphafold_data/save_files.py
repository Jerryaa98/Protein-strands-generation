import os
import shutil
import re

# Define the source directory to search and destination directories
root_dir = '/root/Biology_project/alphafold_data/data'
model_0_dir = '/root/Biology_project/alphafold_data/model_0'
model_1_dir = '/root/Biology_project/alphafold_data/model_1'
model_2_dir = '/root/Biology_project/alphafold_data/model_2'
model_3_dir = '/root/Biology_project/alphafold_data/model_3'
model_4_dir = '/root/Biology_project/alphafold_data/model_4'

# Regex pattern to match CIF filenames
pattern = re.compile(r"fold_seq_(\d{1,3})_model_([01234])\.cif$")

# Make sure destination folders exist
os.makedirs(model_0_dir, exist_ok=True)
os.makedirs(model_1_dir, exist_ok=True)
os.makedirs(model_2_dir, exist_ok=True)
os.makedirs(model_3_dir, exist_ok=True)
os.makedirs(model_4_dir, exist_ok=True)

# Traverse the root directory
for dirpath, _, filenames in os.walk(root_dir):
    for filename in filenames:
        match = pattern.match(filename)
        if match:
            k = int(match.group(1))
            if 1 <= k <= 102:
                model = match.group(2)
                new_filename = f"seq_{k}.cif"
                src_file = os.path.join(dirpath, filename)

                if model == '0':
                    dst_file = os.path.join(model_0_dir, new_filename)
                elif model == '1':
                    dst_file = os.path.join(model_1_dir, new_filename)
                elif model == '2':
                    dst_file = os.path.join(model_2_dir, new_filename)
                elif model == '3':
                    dst_file = os.path.join(model_3_dir, new_filename)
                elif model == '4':
                    dst_file = os.path.join(model_4_dir, new_filename)
                else:
                    print(f"Unknown model {model} for file {filename}. Skipping.")
                    continue

                shutil.copy2(src_file, dst_file)
                print(f"Copied: {src_file} -> {dst_file}")
