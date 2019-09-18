#!/usr/bin/env python3

"""
Does the reverse of pop2csv.py. Useful to provide an interesting input population.

Usage examples:
    csv2pop.py problem_path pf_path
    csv2pop.py problem_path pf_path -o ouput_path
"""

import os, argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('problem_path', help='Path to the problem that the population is built from')
    parser.add_argument('pf_path', help='Path to the .csv file to be converted')
    default_output = "<pf_path>_pf.pop"
    parser.add_argument('-o','--output_path', help='output file name', default=default_output, type=str)
    args = parser.parse_args()

    if args.output_path == default_output:
        output_path = "{}_pf.pop".format(args.pf_path.replace('.csv', ''))
    else:
        output_path = args.output_path

    # Mapping dictionaries for ids
    mdf = pd.read_csv(os.path.join(args.problem_path, "modelidmap.csv"))
    rdf = pd.read_csv(os.path.join(args.problem_path, "rxnidmap.csv"))
    m_map = mdf.set_index('og_model_ids')['new_model_ids'].to_dict()
    r_map = rdf.set_index('all_ids')['new_ids'].to_dict()

    pf_df = pd.read_csv(args.pf_path)

    def write_design(row):
        """ Takes a row of the pf.csv dataframe an dwrites it to the currently open file
        """
        f.write("#INDIVIDUAL\n")
        f.write("#DELETIONS\n")
        for rxn in row['Deletion_id'].split(', '):
            f.write(r_map[rxn] + "\n")
        f.write("#MODULES\n")
        for model_id in m_map.keys():
            item = str(row['{}(module)'.format(model_id)])
            new_rxns=""
            if item != "nan":
                new_rxns = ",".join([r_map[rxn] for rxn in item.split(", ")])
            if new_rxns != "":
                f.write('{},{}\n'.format(m_map[model_id], new_rxns))
            else:
                f.write('{}\n'.format(m_map[model_id]))

        f.write("#OBJECTIVES\n")

    with open(output_path, 'w') as f:
        f.write("#METADATA\n")
        f.write("population_size={}\n".format(len(pf_df.index)))
        f.write("alpha=0\n")
        f.write("beta=0\n")

        for row in pf_df.iterrows():
            write_design(row[1])

    print("Output written to: ", output_path)


if __name__ == '__main__':
    main()
