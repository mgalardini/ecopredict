
# coding: utf-8

# In[1]:

# Plotting imports

import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('white')

plt.rc('font', family='sans-serif')
plt.rc('font', serif='Arial')
plt.rc('text', usetex='false')

try:
    plt.style.use('../custom.mplstyle')
except IOError:
    plt.rc('font', size=10)
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)
    plt.rc('axes', labelsize=12, titlesize=12)
    plt.rc('legend', fontsize=8)


# In[2]:

import itertools
import numpy as np
import pandas as pd
from scipy import stats
from Bio import Phylo


# In[3]:

from nadist.pure.spatial import euclidean_pdist
from scipy.spatial.distance import pdist, squareform


# In[4]:

from sklearn import linear_model


# In[5]:

# In[20]:

import numpy as np
import pandas as pd
from scipy import stats
from ete3 import Tree

import itertools
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.colors as colors
from statsmodels.sandbox.stats.multicomp import multipletests
from ete3 import TreeStyle, NodeStyle, Tree, RectFace


# In[28]:

def get_tree(infile):
    tree = Tree(infile)

    for x in tree.traverse():
        if not x.is_leaf():
            continue
        x.name = x.name.replace("'", '').split('.')[0]
        if x.name == 'genome':
            x.name = 'NT12001_189'
    strains = {x.name.split('_')[0]
               for x in tree.traverse()
               if x.is_leaf()}
    for s in strains:
        nodes = sorted([x
                        for x in tree.traverse()
                        if x.name.startswith(s)],
                       key=lambda x: x.name)
        if len(nodes) == 1:
            continue
        for node in nodes[1:]:
            node.delete()
    for x in tree.traverse():
        if not x.is_leaf():
            continue
        x.name = x.name.split('_')[0]
    tree.set_outgroup(tree.get_midpoint_outgroup())

    return tree


# In[29]:

def draw_tree(pstrains, names=False, subset=None):
    evol = {x.rstrip()
            for x in open('input/evolution_experiment.txt')}
    
    tree = get_tree('input/tree.nwk')

    strains = {x.name for x in tree.traverse()
               if x.name != ''}
    if subset is None:
        subset = strains
    evol_tree = {x.name
                 for x in tree.traverse()
                 if x.name in evol}

    commensals = {}
    for x in strains:
        if pstrains[x] == 'Commensal strain':
            commensals[x] = colors.cnames['blue']
        elif pstrains[x] == 'Pathogenic strain':
            commensals[x] = colors.cnames['red']
        else:
            commensals[x] = colors.cnames['white']

    ref = NodeStyle()
    ref['fgcolor'] = sns.xkcd_rgb['light red']
    inner = NodeStyle()
    inner['size'] = 0
    for n in tree.traverse():
        if not n.is_leaf():
            n.set_style(inner)
            continue
        if n.name not in subset:
            continue
        if not names:
            r = RectFace(10, 3,
                         commensals[n.name],
                         commensals[n.name])
            n.add_face(r, 0,
                       position="aligned")
        ev = NodeStyle()
        ev["fgcolor"] = "black"
        ev['size'] = 3
        n.set_style(ev)

    circular_style = TreeStyle()
    circular_style.mode = "c"
    if not names:
        circular_style.scale = 1300
    else:
        circular_style.scale = 7300
    circular_style.show_leaf_name = names
    return tree, circular_style


# In[30]:

strains = pd.read_table('input/strains.tsv')
strains.set_index('Strain Identifier', inplace=True)


# In[31]:

sdict = strains['Phenotype'].to_dict()


# In[31]:

tree, circular_style = draw_tree(sdict,
                                 False)
tree.render('notebooks/sfig3c.svg',
            h=7.5,
            units="in",
            tree_style=circular_style,
            dpi=90,
            );
