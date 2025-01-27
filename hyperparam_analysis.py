
import argparse
import json
from huggingface_hub import login
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig
import csv
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import os

from tqdm import tqdm
from strandsloader import SequentialStrandLoader, ReverseStrandLoader, RandomStrandLoader

import random
import torch

def init_loader(strategy, data, seed):
    if strategy == 'sequential':
        return SequentialStrandLoader(data)

    elif strategy == 'reverse':
        return ReverseStrandLoader(data)
    
    elif strategy == 'random':
        return RandomStrandLoader(data, seed)
    
    else :
        raise ValueError('Strategy not supported!')

def analyze_seeds(args, param_range):
    # Will instruct you how to get an API key from huggingface hub, make one with "Read" permission.
    login(token = args.hf_token)
    model: ESM3InferenceClient = ESM3.from_pretrained("esm3-open").to("cuda") 
    strands_file = args.strands_file
    work_dir = args.work_dir
    strategy = args.strategy
    seed = args.seed

    MIN_STRAND_LEN = args.MIN_STRAND_LEN
    MAX_STRAND_LEN = args.MAX_STRAND_LEN


    generated_sequences_dir = f"{work_dir}/generated_sequences_analysis/{strategy}"
    if not os.path.exists(generated_sequences_dir):
        try:
            os.makedirs(generated_sequences_dir)  # Create the directory
        except Exception as e:
            raise ValueError(f"Error creating directory '{generated_sequences_dir}': {e}")

    with open(strands_file, 'r') as json_file:
        data = json.load(json_file)


    print('Analyzing the afficet of different seed on the generation')
    print('Starting Strand Generation ...')

    generated_data = {}
    for curr_seed in param_range :
        # set current seed
        random.seed(curr_seed)
        torch.manual_seed(curr_seed)  # Set the seed for CPU
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(curr_seed)  # Set the seed for all GPUs (if using CUDA)
        # directory per seed
        curr_generated_sequences_dir = f'{generated_sequences_dir}/seed_{curr_seed}'
        if not os.path.exists(curr_generated_sequences_dir):
            try:
                os.makedirs(curr_generated_sequences_dir)  # Create the directory
            except Exception as e:
                raise ValueError(f"Error creating directory '{curr_generated_sequences_dir}': {e}")
        
        for id in tqdm(data.keys()):
            seq = data[id]['seq']
            strands_indcies = data[id]['strands']
            state = [seq]
            
            strand_loader = init_loader(strategy, strands_indcies, seed)
            for i, strand in enumerate(strand_loader):
                if len(strand) < MIN_STRAND_LEN or len(strand) > MAX_STRAND_LEN:
                    continue
                # mask
                masked_seq = ""
                for j, r in enumerate(state[-1]):
                    if j+1 in strand:
                        masked_seq += '_'
                    else:
                        masked_seq += r

                protein = ESMProtein(sequence=masked_seq)

                # Generate and save seq
                num_steps = max(1, len(strand)//2)
                protein = model.generate(protein, GenerationConfig(track="sequence", num_steps=num_steps, temperature=args.temperature))
                output_seq = protein.sequence
                state.append(output_seq)
                # Generate and save pdb files 
                protein = model.generate(protein, GenerationConfig(track="structure", num_steps=num_steps))
                generated_sequences = f"{curr_generated_sequences_dir}/{id}"
                # generated seq after each masking foreach seq (id)
                if not os.path.exists(generated_sequences):
                    try:
                        os.makedirs(generated_sequences)  # Create the directory
                    except Exception as e:
                        raise ValueError(f"Error creating directory '{generated_sequences}': {e}")
                    
                protein.to_pdb(f"{generated_sequences}/generation{i}.pdb")

            generated_data[id] = state

        with open(f'{curr_generated_sequences_dir}/seq.json','w') as json_file:
            json.dump(generated_data, json_file, indent=4)
    
        print(f'Finished Strand Generation for seed : {curr_seed}!')
    print('Finished Strand Generation !')

def analyze_temperature(args, param_range):
    # Will instruct you how to get an API key from huggingface hub, make one with "Read" permission.
    login(token = args.hf_token)
    model: ESM3InferenceClient = ESM3.from_pretrained("esm3-open").to("cuda") 
    strands_file = args.strands_file
    work_dir = args.work_dir
    strategy = args.strategy
    seed = args.seed

    MIN_STRAND_LEN = args.MIN_STRAND_LEN
    MAX_STRAND_LEN = args.MAX_STRAND_LEN


    generated_sequences_dir = f"{work_dir}/generated_sequences_analysis/{strategy}"
    if not os.path.exists(generated_sequences_dir):
        try:
            os.makedirs(generated_sequences_dir)  # Create the directory
        except Exception as e:
            raise ValueError(f"Error creating directory '{generated_sequences_dir}': {e}")

    with open(strands_file, 'r') as json_file:
        data = json.load(json_file)

    print('Analyzing the afficet of different temperature on the generation')
    print('Starting Strand Generation ...')

    generated_data = {}
    for curr_temperature in param_range :
        # set a fixed seed
        random.seed(args.seed)
        torch.manual_seed(args.seed)  # Set the seed for CPU
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(args.seed)  # Set the seed for all GPUs (if using CUDA)

        # directory per temperature
        curr_generated_sequences_dir = f'{generated_sequences_dir}/temperature_{curr_temperature}'
        if not os.path.exists(curr_generated_sequences_dir):
            try:
                os.makedirs(curr_generated_sequences_dir)  # Create the directory
            except Exception as e:
                raise ValueError(f"Error creating directory '{curr_generated_sequences_dir}': {e}")
        
        for id in tqdm(data.keys()):
            seq = data[id]['seq']
            strands_indcies = data[id]['strands']
            state = [seq]
            
            strand_loader = init_loader(strategy, strands_indcies, seed, args.NxLoop)
            for i, strand in enumerate(strand_loader):
                if len(strand) < MIN_STRAND_LEN or len(strand) > MAX_STRAND_LEN:
                    continue
                # mask
                masked_seq = ""
                for j, r in enumerate(state[-1]):
                    if j+1 in strand:
                        masked_seq += '_'
                    else:
                        masked_seq += r

                protein = ESMProtein(sequence=masked_seq)

                # Generate and save seq
                num_steps = max(1, len(strand)//2)
                protein = model.generate(protein, GenerationConfig(track="sequence", num_steps=num_steps, temperature=curr_temperature))
                output_seq = protein.sequence
                state.append(output_seq)
                # Generate and save pdb files 
                protein = model.generate(protein, GenerationConfig(track="structure", num_steps=num_steps))
                generated_sequences = f"{curr_generated_sequences_dir}/{id}"
                # generated seq after each masking foreach seq (id)
                if not os.path.exists(generated_sequences):
                    try:
                        os.makedirs(generated_sequences)  # Create the directory
                    except Exception as e:
                        raise ValueError(f"Error creating directory '{generated_sequences}': {e}")
                    
                protein.to_pdb(f"{generated_sequences}/generation{i}.pdb")

            generated_data[id] = state

        with open(f'{curr_generated_sequences_dir}/seq.json','w') as json_file:
            json.dump(generated_data, json_file, indent=4)
    
        print(f'Finished Strand Generation for temperature : {curr_temperature}!')
    print('Finished Strand Generation !')

def run(args):
    param_range = []
    # generated randomly 
    if args.param == 'seed' :
        param_range = [1055, 1491, 1779, 2152, 2280, 
                       2795, 3017, 3188, 3829, 4604, 
                       4762, 4772, 6507, 6600, 7631, 
                       7659, 8357, 8531, 8956, 9213]
        
        analyze_seeds(args, param_range)

    elif args.param == 'temperature' :
        param_range = [ 0.01322636856325421, 0.04699817656041905, 0.07464317162580314, 0.12045644716760062, 
                        0.13916569878527563, 0.15008777501514126, 0.21204356514687728, 0.21947129622082506, 
                        0.3432534576778067, 0.38764365536787626, 0.3883917640946134, 0.3895152114525644, 
                        0.42611710474265585, 0.4502129715137575, 0.4748397666642211, 0.5754839007807748, 
                        0.5759104997976984, 0.6851649929442852, 0.7112923885903855, 0.8936197308074946]

        analyze_temperature(args, param_range)
    else :
        raise ValueError('Param analysis isn\'t supported')


if __name__ == "__main__":

    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(description="Strands Generation")
    # Add arguments
    parser.add_argument("--strands_file", type=str, default='/root/Biology_project/data/strands.json', help='strands file json file')
    parser.add_argument("--work_dir", type=str, default='/root/Biology_project', help='dirctory to save pdb file in')
    parser.add_argument("--hf_token", type=str, default="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")
    parser.add_argument("--MIN_STRAND_LEN",type=int,default=5)
    parser.add_argument("--MAX_STRAND_LEN",type=int,default=25)
    parser.add_argument("--seed",type=int,default=42)
    parser.add_argument("--temperature",type=float,default=0.7)
    parser.add_argument("--strategy",type=str,
                        choices=['sequential', 'reverse', 'random'],
                        default='sequential')
    parser.add_argument("--param", type=str,
                        choices=['seed', 'temperature'],
                        default='seed')
    parser.add_argument("--NxLoop", type=int, default=1, help="Number of (Loop) generatens per protein") 

    args = parser.parse_args()
    run(args)
