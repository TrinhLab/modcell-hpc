#!/usr/bin/env python3

"""
Formats population using the table format compatible with other ModCell tools. It only keeps non-dominated individuals.Ignores metadata in the .pop file.

Usage examples:
    pop2csv.py problem_path population_path
    pop2csv.py problem_path population_path -o ouput_path
"""

import os, argparse
import pandas as pd
import numpy as np
from libdomi import fdom

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('problem_path', help='Path to the problem that the population is built from')
    parser.add_argument('pop_path', help='Path to the .pop file to be converted')
    default_output = "<pop_path>.csv"
    parser.add_argument('-o','--output_path', help='output file name', default=default_output, type=str)
    parser.add_argument('-w','--write_all', help='writes all solutions in the population instead of only dominated ones', action='store_true')
    default_alpha = 10000
    parser.add_argument('-a','--alpha', help='Will check for solutions that violate alpha and discard them', default=default_alpha, type=int)
    args = parser.parse_args()

    if args.output_path == default_output:
        output_path = args.pop_path.replace('.pop', '.csv')
    else:
        output_path = args.output_path

    if args.alpha == default_alpha:
        print("Warning: Alpha constraint is not strictly enforced in .pop due to the use of a penalty function. If you do not specify --alpha for your problem, solutions that have more deletions than intended will not be discarded in the output .csv file")

    # Mapping dictionaries for ids
    mdf = pd.read_csv(os.path.join(args.problem_path, "modelidmap.csv"))
    rdf = pd.read_csv(os.path.join(args.problem_path, "rxnidmap.csv"))
    m_map = mdf.set_index('new_model_ids')['og_model_ids'].to_dict()
    r_map = rdf.set_index('new_ids')['all_ids'].to_dict()

    # Parse file into table
    in_individual= False
    in_deletions = False
    in_modules = False
    in_objectives = False
    line_n = 0

    individuals = []
    indv = None

    def append_indv(indv):
        if indv['Deletion_id']: # Drop individuals without deletions
            indv['Deletion_id'] = ', '.join(indv['Deletion_id'])
            individuals.append(indv)

    with open(args.pop_path, 'r') as f:
        for rline in f:
            line = rline.strip()
            # Set state // Note: This assumes the lines are encountered in a given order so it is not the most robust solution
            if line == "#INDIVIDUAL":
                in_objectives = False
                if indv:
                    append_indv(indv)
                indv = {}
                indv['Deletion_id'] = []
                continue
            if line == "#DELETIONS":
                in_deletions = True
                in_modules = False
                continue
            if line == "#MODULES":
                in_deletions = False
                in_modules = True
                continue
            if line == "#OBJECTIVES":
                in_modules = False
                in_deletions = False
                in_objectives = True
                continue
            if line == "#ENDFILE":
                append_indv(indv)
                break

            if in_deletions:
                indv['Deletion_id'].append(r_map[line])

            elif in_modules:
                ml = line.split(',')
                model_id = m_map[ml[0]]
                module_rxns = ml[1:]
                indv["{}(module)".format(model_id)] = ', '.join([r_map[rid] for rid in module_rxns])

            elif in_objectives:
                ml = line.split(',')
                model_id = m_map[ml[0]]
                objective_value = ml[1]
                indv["{}(objective)".format(model_id)] = objective_value

    df = pd.DataFrame(individuals)

    # Check alpha
    sizes = df['Deletion_id'].map(lambda x: len(x.split(',')))
    if any(sizes > args.alpha):
        df = df.loc[sizes <= args.alpha,:]
        print(f"Designs that meet alpha={args.alpha}:\t {sum(sizes>0) - sum(sizes>args.alpha)}/{sum(sizes>0)}")

    # Reorganize columns
    model_ids = list(m_map.values())
    df = df[['Deletion_id'] + ['{}(module)'.format(x) for x in model_ids] + ['{}(objective)'.format(x) for x in model_ids]]

    # Remove duplicated solutions
    obj_idx = df.columns.str.contains(r'\(objective\)')
    prior_size = df.shape[0]
    df = df.drop_duplicates(subset=df.columns[np.logical_not(obj_idx)])  # Do  not consider objectives in case of FP error
    print(f"Unique solutions:\t\t {df.shape[0]}/{prior_size}")

    if not args.write_all:
        # Only keep non-dominated solutions
        a = df.iloc[:,obj_idx].values
        non_dominated_idx = fdom(a)
        print(f"Non-dominated solutions:\t {np.sum(non_dominated_idx)}/{df.shape[0]}")
        if (np.sum(non_dominated_idx) < df.shape[0]):
            df = df.iloc[non_dominated_idx.flatten(),:]
            df.reset_index(inplace=True, drop=True)

    # Add solution index
    df.insert(loc=0, column='Solution index', value=[_ for _ in range(df.shape[0])])

    #Save
    df.to_csv(output_path, index=False)


if __name__ == '__main__':
    main()
