#!/usr/bin/env python3

"""
Formats population using the table format compatible with other ModCell tools. It only keeps non-dominated individuals.

Notes:
    - Add option to keep dominated individuals?
"""

import os, argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('problem_path', help='Path to the problem that the population is built from')
    parser.add_argument('pop_path', help='Path to the .pop file to be converted')
    default_output = "<pop_path>.csv"
    parser.add_argument('-o','--output_path', help='output file name', default=default_output, type=str)
    args = parser.parse_args()

    if args.output_path == default_output:
        output_path = args.pop_path.replace('.pop', '.csv')

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

    with open(args.pop_path, 'r') as f:
        for rline in f:
            line = rline.strip()
            # Set state // Note: This assumes the lines are encountered in a given order so it is not the most robust solution
            if line == "#INDIVIDUAL":
                in_objectives = False
#                in_individual = True
                if indv:
                    indv['Deletion_id'] = ', '.join(indv['Deletion_id'])
                    individuals.append(indv)
                indv = {}
                indv['Deletion_id'] = []
                continue
            if line == "#DELETIONS":
                in_deletions = True
                continue
            if line == "#MODULES":
                in_deletions = False
                in_modules = True
                continue
            if line == "#OBJECTIVES":
                in_modules = False
                in_objectives = True
                continue
            if line == "#ENDFILE":
                continue

            if in_deletions:
                indv['Deletion_id'].append(r_map[line])

            if in_modules:
                ml = line.split(',')
                model_id = m_map[ml[0]]
                module_rxns = ml[1:]
                indv["{}(module)".format(model_id)] = ', '.join([r_map[rid] for rid in module_rxns])

            if in_objectives:
                ml = line.split(',')
                model_id = m_map[ml[0]]
                objective_value = ml[1]
                indv["{}(objective)".format(model_id)] = objective_value

    df = pd.DataFrame(individuals)

    # Reorganize columns
    model_ids = list(m_map.values())
    df = df[['Deletion_id'] + ['{}(module)'.format(x) for x in model_ids] + ['{}(objective)'.format(x) for x in model_ids]]

    # Only keep non-dominated solutions
    obj_idx = df.columns.str.contains(r'\(objective\)')
    idx_to_drop = []
    for idx, row in df.iterrows():
        for idx_b, row_b in df.iterrows():
            if idx_b != idx:
                comp = row[obj_idx] <= row_b[obj_idx]
                if comp.all():
                    idx_to_drop.append(idx)
                    break
    print("Non-dominated solutions: {}/{}".format(df.shape[0] - len(idx_to_drop), df.shape[0]))
    df = df.drop(idx_to_drop)

    # Add solution index
    df.insert(loc=0, column='Solution index', value=[_ for _ in range(df.shape[0])])

    #Save
    df.to_csv(output_path, index=False)


if __name__ == '__main__':
   main()
