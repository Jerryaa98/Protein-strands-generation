import re

pdb_path = '/root/Biology_project/pdb_files/e2wjrA1.pdb'
# Read the PDB file
with open(pdb_path, 'r') as file:
    lines = file.readlines()

# Filter out lines where the 6th column (occupancy) is missing for ATOM/HETATM lines
filtered_lines = []
for line in lines:
    # Check if the line starts with ATOM 
    if line.startswith("ATOM"):
        # Use regex to split by one or more spaces
        columns = re.split(r'\s+', line.strip())
        # Include the line if there are 11 or more columns
        if len(columns) == 11:
            filtered_lines.append(line)
    else:
        # Include all other lines (non-ATOM/HETATM)
        filtered_lines.append(line)

# Write the filtered lines back to a new PDB file
with open('/root/Biology_project/pdb_files/filtered_e2wjrA1.pdb', 'w') as output_file:
    output_file.writelines(filtered_lines)

print("Filtered PDB file created: 'filtered_file.pdb'")
