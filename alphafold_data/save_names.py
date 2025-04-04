import os
import pandas as pd

# Paths
csv_path = '/root/Biology_project/data/OMBB_data.csv'
model_0_dir = '/root/Biology_project/alphafold_data/model_0'
model_1_dir = '/root/Biology_project/alphafold_data/model_1'
model_2_dir = '/root/Biology_project/alphafold_data/model_2'
model_3_dir = '/root/Biology_project/alphafold_data/model_3'
model_4_dir = '/root/Biology_project/alphafold_data/model_4'

# Make sure destination folders exist
os.makedirs(model_0_dir, exist_ok=True)
os.makedirs(model_1_dir, exist_ok=True)
os.makedirs(model_2_dir, exist_ok=True)
os.makedirs(model_3_dir, exist_ok=True)
os.makedirs(model_4_dir, exist_ok=True)


# Load the CSV
df = pd.read_csv(csv_path)

# Loop over the rows (assuming 1-based indexing)
for k in range(1, 103):  # 1 to 102 inclusive
    row = df.iloc[k - 1]  # Adjust for 0-based indexing
    file_id = row['id']

    old_name = f"seq_{k}.cif"
    new_name = f"{file_id}.cif"

    old_path_0 = os.path.join(model_0_dir, old_name)
    new_path_0 = os.path.join(model_0_dir, new_name)

    old_path_1 = os.path.join(model_1_dir, old_name)
    new_path_1 = os.path.join(model_1_dir, new_name)

    old_path_2 = os.path.join(model_2_dir, old_name)
    new_path_2 = os.path.join(model_2_dir, new_name)

    old_path_3 = os.path.join(model_3_dir, old_name)
    new_path_3 = os.path.join(model_3_dir, new_name)

    old_path_4 = os.path.join(model_4_dir, old_name)
    new_path_4 = os.path.join(model_4_dir, new_name)

    if os.path.exists(old_path_0):
        os.rename(old_path_0, new_path_0)
        print(f"Renamed in model_0: {old_name} -> {new_name}")
    if os.path.exists(old_path_1):
        os.rename(old_path_1, new_path_1)
        print(f"Renamed in model_1: {old_name} -> {new_name}")
    if os.path.exists(old_path_2):
        os.rename(old_path_2, new_path_2)
        print(f"Renamed in model_2: {old_name} -> {new_name}")
    if os.path.exists(old_path_3):
        os.rename(old_path_3, new_path_3)
        print(f"Renamed in model_3: {old_name} -> {new_name}")
    if os.path.exists(old_path_4):
        os.rename(old_path_4, new_path_4)
        print(f"Renamed in model_4: {old_name} -> {new_name}")
