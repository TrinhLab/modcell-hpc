"""Functions common to multiple methods"""

import seaborn as sns
import pandas as pd
import csv


def setup_plot_look(args):
    if args.font_size:
        fs = args.font_size
    else:
        fs = 12

    if args.width:
        sizet = tuple([args.width, args.height])
    else:
        sizet = tuple([4, 4])

    sns.set(style='whitegrid',
            font=args.font, font_scale=fs/12, # this is not the exat font size in pt but nothing else works
            rc={'figure.figsize':sizet})


def square_plot(ax):
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    ax.set_aspect(abs(x1-x0)/abs(y1-y0))


def find_n_prods(pf_path):
    df = pd.read_csv(pf_path)
    return df.columns.str.contains("objective").sum()


def get_prod_id2name(args):
    if args.id2name_path:
        with open(args.id2name_path) as f:
            return {d['id']: d['name'] for d in csv.DictReader(f, fieldnames=('id', 'name'))}
    else:
            return None
