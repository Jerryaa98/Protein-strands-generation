# converts alphafold data from all models cif files to pdb with the relative names

python save_files.py # cuts the relative info from alphafold's data
python save_names.py # saves the files with the correct names
python cif_to_pdb.py # converts all files to pdb