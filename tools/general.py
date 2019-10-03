# Functions used by multiple methods

import seaborn as sns
import pandas as pd
import csv

def setup_plot_look(args):
    sns.axes_style({'font.family': ['sans-serif'],
        'font.sans-serif': [args.font],
        'xtick.labelsize' : 12,
        'ytick.labelsize' : 12,
        'axes.labelsize' : 12,
        });


def setup_plot_size(args):
    sns.set(rc={'figure.figsize':tuple([args.width, args.height])})


def find_n_prods(pf_path):
    df = pd.read_csv(pf_path)
    return df.columns.str.contains("objective").sum()


def square_plot(ax):
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    ax.set_aspect(abs(x1-x0)/abs(y1-y0))


def get_prod_id2name(args):
    if args.id2name_path:
        with open(args.id2name_path) as f:
            return {d['id']: d['name'] for d in csv.DictReader(f, fieldnames=('id', 'name'))}
    else:
            return None
