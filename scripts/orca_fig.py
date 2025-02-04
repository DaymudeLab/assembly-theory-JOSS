"""
orca_fig.py

garrett parzych and sean bergen

this file creates a dataset from multiple runs of the ORCA software, and
then uses that dataset to create figres
"""
import pandas as pd
import matplotlib.pyplot as plt

######## dataset creation ########



######## figure  creation ########

# data we are reading in is a csv in the format
# [mol_name, n_edges, n_heavyatoms, ai, search_space, sgi_time, search_time, 
#  states_searched]
#
# mol_name    : name of the molecule
# n_edges     : number of edges in the molecule's graph
# n_heavyatoms: number of heavy atoms in the molecule
# ai          : computed assembly index of the molecule
# search_space: number of duplicate structures in the molecule
# sgi_time    : time to compute subgraph-isomorphisms, averaged over x runs
# search_time : time to compute assembly index post sgi, averaged over x runs
# states_searched: number of states searched when computing assembly index


