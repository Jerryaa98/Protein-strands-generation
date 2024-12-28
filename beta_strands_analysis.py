import argparse
import json
import numpy as np
import matplotlib.pyplot as plt 

def run(args):
    work_dir = args.work_dir
    strands_file = args.strands_file
    with open(strands_file, 'r') as file:
        data = json.load(file)

    strands_lens = []
    for id in data.keys():
        for strand in data[id]['strands']:
            strands_lens.append(len(strand))

    uniqe_lens, counts = np.unique(strands_lens, return_counts=True)
    plt.bar(uniqe_lens, counts, color='skyblue', edgecolor='black')
    plt.xlabel('Beta Strand length')
    plt.ylabel('Frequency')
    plt.savefig(f'{work_dir}/BetaStrandLength.png')



if __name__ == "__main__":
         # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Strands Generation")
    # Add arguments
    parser.add_argument("--work_dir", type=str, default='/root/Biology_project', help='dirctory to save pdb file in')
    parser.add_argument("--strands_file", type=str)
    args = parser.parse_args()
    run(args)