import json
import argparse
from tqdm import tqdm
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
def calculate_seq_identity(seq1, seq2):
    # Calculate identity
    matches = sum(res1 == res2 for res1, res2 in zip(seq1, seq2))
    identity = (matches / max(len(seq1), len(seq2))) * 100
    return identity


def plotHeapMap(results_dir, strategy, id, param, param_range, data, iter_num):
    seq_comp = [[0]*len(data) for _ in data]
    for i, seq1 in enumerate(data):
        for j, seq2 in enumerate(data):
            seq_comp[i][j] = calculate_seq_identity(seq1, seq2)

    # Create the heatmap
    seq_comp = np.array(seq_comp)
    plt.figure(figsize=(12, 10))
    sns.heatmap(seq_comp, annot=False, cmap='coolwarm', cbar=True, xticklabels=param_range, yticklabels=param_range)
    # Adjust the axis labels
    plt.xticks(rotation=45, ha="right", fontsize=8)  # Rotate x-axis labels
    plt.yticks(fontsize=8)  # Adjust font size for y-axis labels
    plt.xlabel('Seed')
    plt.ylabel('Seed')
    # Add a title
    plt.title(f"Generated sequences in the {iter_num}th iteration comparsion using {strategy} strategy and different {param}s")

    # Save the figure
    curr_dir = f"{results_dir}/{strategy}/{param}/{id}"
    if not os.path.exists(curr_dir):
        try:
            os.makedirs(curr_dir)  # Create the directory
        except Exception as e:
            raise ValueError(f"Error creating directory '{curr_dir}': {e}")
        
    plt.savefig(f"{curr_dir}/generation_{iter_num}.png", dpi=100)  


def run(args):
    strands_file = args.strands_file
    work_dir = args.work_dir
    results_dir = f'{work_dir}/hyperparam_comp_results'
    gen_seq_ana_dir = args.gen_seq_ana_dir
    strategy = args.strategy
    MIN_STRAND_LEN = args.MIN_STRAND_LEN
    MAX_STRAND_LEN = args.MAX_STRAND_LEN
    gen_seq_ana_strategy_dir =  f'{gen_seq_ana_dir}/{strategy}'
    param = args.param

    param_range = []
    # generated randomly 
    if args.param == 'seed' :
        param_range = [1055, 1491, 1779, 2152, 2280, 
                       2795, 3017, 3188, 3829, 4604, 
                       4762, 4772, 6507, 6600, 7631, 
                       7659, 8357, 8531, 8956, 9213]
        
    elif args.param == 'temperature' :
        param_range = [ 0.01322636856325421, 0.04699817656041905, 0.07464317162580314, 0.12045644716760062, 
                        0.13916569878527563, 0.15008777501514126, 0.21204356514687728, 0.21947129622082506, 
                        0.3432534576778067, 0.38764365536787626, 0.3883917640946134, 0.3895152114525644, 
                        0.42611710474265585, 0.4502129715137575, 0.4748397666642211, 0.5754839007807748, 
                        0.5759104997976984, 0.6851649929442852, 0.7112923885903855, 0.8936197308074946]
    else :
        raise ValueError('Param analysis isn\'t supported')
    
    with open(strands_file, 'r') as json_file:
        data = json.load(json_file)

    generated_seq = dict()
    for curr_param in param_range:
        seq_file = f'{gen_seq_ana_strategy_dir}/{param}_{curr_param}/seq.json'
        with open(seq_file, 'r') as json_file :
            gen_seq = json.load(json_file)
            generated_seq[curr_param] = gen_seq

    for id in tqdm(data.keys()):
        gen_iter_num = len([strand for strand in data[id]['strands'] if len(strand) >= MIN_STRAND_LEN and len(strand) <= MAX_STRAND_LEN])
        gen_iter_num += 1 # for the orginal seq
            
        for i in range(gen_iter_num):
            seq_ith_gen_with_param = []
            for curr_param in param_range:
                seq_ith_gen_with_param.append(generated_seq[curr_param][id][i])
            # save comparsion betweeen the generated sequences in the ith iteration between all params
            plotHeapMap(results_dir, strategy, id, param, param_range, seq_ith_gen_with_param, i)

    return 

if __name__ == "__main__" :
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="HyperParam generation  analysis")
    # Add arguments
    parser.add_argument("--strands_file", type=str, default='/root/Biology_project/data/strands.json', help='strands file json file')
    parser.add_argument("--work_dir", type=str, default='/root/Biology_project', help='dirctory to save pdb file in')
    parser.add_argument("--gen_seq_ana_dir", type=str, default='/root/Biology_project/generated_sequences_analysis')
    parser.add_argument("--MIN_STRAND_LEN",type=int,default=5)
    parser.add_argument("--MAX_STRAND_LEN",type=int,default=25)
    parser.add_argument("--strategy",type=str,
                        choices=['sequential', 'reverse', 'random'],
                        default='sequential')
    parser.add_argument("--param", type=str,
                        choices=['seed', 'temperature'],
                        default='seed')
    args = parser.parse_args()
    run(args)

