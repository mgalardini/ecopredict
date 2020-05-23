#!/usr/bin/env python

__author__ = "Marco Galardini"
__version__ = '0.1.0'


def get_options():
    import argparse

    # create the top-level parser
    description = "Get pairwise phylogenetic distance"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('tree', action='store',
                        help='tree file')
    return parser.parse_args()

if __name__ == "__main__":
    import sys
    import numpy as np
    import pandas as pd
    from Bio import Phylo
    import itertools

    options = get_options()

    t = Phylo.read(options.tree, 'newick')
    ref = None
    for x in t.get_terminals():
        x.name = x.name.split('_')[0]
        if 'genome' in x.name:
            x.name = 'NT12001'
            ref = x
    strains = {x.name for x in t.get_terminals()}

    res = []
    for s1, s2 in itertools.combinations(strains, 2):
        res.append((s1, s2, t.distance(s1, s2)))
    r = pd.DataFrame(res)
    r.columns = ['strain1', 'strain2',
                 'phylo']
    r.dropna(inplace=True)
    r.to_csv(sys.stdout, sep='\t', index=False)
