"""
orca_fig.py

garrett parzych and sean bergen

this file creates a dataset from multiple runs of the ORCA software, and
then uses that dataset to create figres
"""
import pandas as pd
import matplotlib.pyplot as plt
import os
import cmcrameri.cm as cmc

# we can get the path to where scripts, orca, and figs are as follows:
scripts_path = os.path.dirname(os.path.abspath(__file__))
figs_path = os.join(os.path.dirname(scripts_path), "figures")
orca_path = os.join(os.path.dirname(scripts_path), "ORCA")

######## dataset creation ########



######## figure  creation ########

# data we are reading in is in the format
# [mol_name, n_edges, n_heavyatoms, ai, search_space, time_to_cmp, sgi_time, 
#  search_time, states_searched]
#
# mol_name    : name of the molecule
# n_edges     : number of edges in the molecule's graph
# n_heavyatoms: number of heavy atoms in the molecule
# ai          : computed assembly index of the molecule
# search_space: number of duplicate structures in the molecule
# time_to_cmp : total time to complete a full run, averaged over x runs
# sgi_time    : time to compute subgraph-isomorphisms, averaged over x runs
# search_time : time to compute assembly index post sgi, averaged over x runs
# states_searched: number of states searched when computing assembly index

# at this point, we should have three pandas dataframes:
#   - naive
#   - logbound
#   - seetbound
#
# the plot we are making is a scatter plot that shows for each of these
# three dataframes, color of a point indicating what df it comes from
#   x:search space
#   y:time to compute
fig, ax = plt.subplots(tight_layout=True)

ax.scatter('search_space', 'time_to_cmp', data=naive, color=cmc.batlowS(0),
        label="Naive")
ax.scatter('search_space', 'time_to_cmp', data=logbound, color=cmc.batlowS(1),
        label="Log-bound")
ax.scatter('search_space', 'time_to_cmp', data=seetbound, color=cmc.batlowS(2),
        label="Seet-bound")

ax.set_xlabel("Search Space")
ax.set_ylabel("Time to Compute")

fig.legend(loc='outside right center', fontsize='small')

fig.savefig(figs_path + "/method_comparison.svg")
