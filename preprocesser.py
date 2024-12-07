import urllib.request
import csv
import os
import argparse
import numpy as np
import json 
import time
import subprocess
# ECOD ss download link : http://prodata.swmed.edu/ecod/af2_pdb/structure?id=e2iahA3
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

        url = f"http://prodata.swmed.edu/ecod/af2_pdb/structure?id={id}"

        # File to save the downloaded content
        #output_file = f"/root/Biology_project/pdb_files/{id}.pdb"
        output_file = f"{output_dir}/{id}.pdb"
        # Download the file
        try:
            urllib.request.urlretrieve(url, output_file)
            print(f"File saved as {output_file}")
        except Exception as e:
            print(f"Failed to download file: {e}")

def get_strands(ss_file):
    residue_indcies = []
    with open(ss_file, 'r') as file:
        for line in file :
            ss_type = line[24]
            index = int(line[17:21].strip())
            if ss_type.strip() == 'E' or ss_type.strip() == 'B':
                residue_indcies.append(index)

    # grouping
    strands = []
    current_strand = [residue_indcies[0]]
    for residue in residue_indcies[1:] :
        if residue == current_strand[-1] + 1:
            current_strand.append(residue)
        else:
            strands.append(current_strand)
            current_strand = [residue]
    return current_strand

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
        
    ss_dir = f"{work_dir}/ss_files"
    if not os.path.exists(ss_dir):
        try:
            os.makedirs(ss_dir)  # Create the directory
        except Exception as e:
            raise ValueError(f"Error creating directory '{ss_dir}': {e}")

    if(load_pdb):
        load_pdb_files(data_file, pdb_dir)

    ids = os.listdir(pdb_dir)
    faulty_ids = []

    for id in ids:
        id = id.split('.')[0]
        ss_output_file = f"{ss_dir}/{id}.txt"
        
        if os.path.exists(ss_output_file):
            continue

        try :
            # cliping the stands
            command = f"stride {pdb_dir}/{id}.pdb | grep '^ASG' >> {ss_output_file}"
            with open(ss_output_file, 'w') as output_file:
                subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            faulty_ids.append(id)
            print(f"Error executing stride: {e}")
    # print(faulty_ids)

    # Save faulty IDs to a JSON file
    fault_file = f"{work_dir}/faulty_ids.json"
    try:
        with open(fault_file, 'w') as json_file:
            json.dump(faulty_ids, json_file, indent=4)
        print(f"Faulty IDs saved to {fault_file}")
    except Exception as e:
        print(f"Error saving faulty IDs: {e}")


    # save meta/data in json file
    data = {}
    with open(data_file, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            id = row[0]
            seq = row[2]

            ss_output_file = f"{ss_dir}/{id}.txt"
            strands = get_strands(ss_output_file)

            state = {'id': id, 'seq': seq, 'strands':strands}
            data[id] = state

    with open(f'{work_dir}/strands.json','w') as json_file:
        json.dump(data, json_file)

if __name__ == "__main__":
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Strands Generation")
    # Add arguments
    parser.add_argument("--data_file", type=str, default='./OMBB_data.csv', help="path for the dataset")
    parser.add_argument("--load_pdb", type=bool, default=False, help='load pdp file using wget or not')
    parser.add_argument("--work_dir", type=str, default='/root/Biology_project', help='dirctory to save pdb file in')
    args = parser.parse_args()
    run(args)
