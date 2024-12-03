import urllib.request
import csv
import os
import argparse
import numpy as np

# https://files.rcsb.org/download/1AF6.pdb
def load_pdb_files(data_file, output_dir):
    
    assert data_file.endswith('.csv'), 'DataSet file must be csv'

    ids = []
    with open(data_file, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            ids.append(row[0])
    #skip header
    ids = ids[1:]

    for id in ids:
        # URL of the file
        domain = id[1:5]
        url = f"https://files.rcsb.org/download/{domain}.pdb"

        # File to save the downloaded content
        #output_file = f"/root/Biology_project/pdb_files/{id}.pdb"
        output_file = f"{output_dir}/{domain}.pdb"
        # Download the file
        try:
            urllib.request.urlretrieve(url, output_file)
            print(f"File saved as {output_file}")
        except Exception as e:
            print(f"Failed to download file: {e}")


def run(args):

    data_file = args.data_file
    load_pdb = args.load_pdb
    work_dir = args.work_dir
    pdb_dir = f"{work_dir}/pdb_files"
    if not os.path.exists(pdb_dir):
        try:
            os.makedirs(pdb_dir)  # Create the directory
        except Exception as e:
            raise ValueError(f"Error creating directory '{pdb_dir}': {e}")
    dssp_dir = f"{work_dir}/dssp_files"
    if not os.path.exists(dssp_dir):
        try:
            os.makedirs(dssp_dir)  # Create the directory
        except Exception as e:
            raise ValueError(f"Error creating directory '{dssp_dir}': {e}")

    if(load_pdb):
        load_pdb_files(data_file, pdb_dir)

    domains = os.listdir(pdb_dir)
    for domain_ in domains[0:1]:
        domain = domain_.split('.')[0]

        dssp_output_file = f"{dssp_dir}/{domain}.txt"
        try :
            os.system(f'mkdssp {pdb_dir}/{domain_} >> {dssp_output_file}')
        except Exception as e:
            print(f"Failed to excute dssp : {e} at file {dssp_output_file}")
    
    
    
        



        







if __name__ == "__main__":
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Strands Generation")
    # Add arguments
    parser.add_argument("--data_file", type=str, help="path for the dataset")
    parser.add_argument("--load_pdb", type=bool, default=False, help='load pdp file using wget or not')
    parser.add_argument("--work_dir", type=str, default='/root/Biology_project', help='dirctory to save pdb file in')

    args = parser.parse_args()
    run(args)
