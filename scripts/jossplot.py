"""
jossplot: Plot molecules' number of disjoint isomorphic subgraph pairs vs.
          their average assembly-theory assembly index calculation time
"""

import argparse
import cmcrameri.cm as cmc
import csv
import json
import matplotlib.pyplot as plt
import os
import os.path as osp
import pandas as pd

if __name__ == "__main__":
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description=__doc__)
    def_crit_path = osp.join("assembly-theory", "target", "criterion", "jossplot")
    parser.add_argument('-C', '--crit_path', type=str, default=def_crit_path,
                        help='Path to jossplot criterion output directory')
    parser.add_argument('-F', '--figs_path', type=str, default='figures',
                        help='Path to figures directory')
    parser.add_argument('-D', '--datasets', nargs='+', default=['gdb17_800'],
                        help='Space-separated list of datasets to include')
    args = parser.parse_args()

    # Load molecule isomorphic subgraph pair counts.
    iso_sub_pairs = {}
    with open(osp.join(args.crit_path, "jossplot.csv"), 'r') as f:
        csv_reader = csv.reader(f)
        for row in csv_reader:
            iso_sub_pairs[row[0]] = int(row[1])

    # Load assembly index calculation times per molecule per bound option.
    mol_dfs = {}
    for bound in ["naive", "logbound", "intbound", "allbounds"]:
        # Set up results dict to turn into a DataFrame.
        results = {"mol_name": [], "iso_sub_pairs": [], "time": []}

        # Gather per molecule data.
        for mol_name, mol_path in [(d.name, d.path)
                                   for d in os.scandir(args.crit_path)
                                   if d.name.split('-')[0] in args.datasets]:
            # Set molecule name and isomorphic subgraph pair count.
            results["mol_name"].append(mol_name)
            results["iso_sub_pairs"].append(iso_sub_pairs[mol_name])

            # Retrieve and set assembly index calculation time in seconds.
            with open(osp.join(mol_path, bound, "new", "estimates.json")) as f:
                mean_time = float(json.load(f)["mean"]["point_estimate"])
                results["time"].append(mean_time / 1e9)

        # Convert results dict to DataFrame.
        mol_dfs[bound] = pd.DataFrame(data=results)
        mol_dfs[bound].sort_values(by="iso_sub_pairs", inplace=True)

    # Plot molecules' number of duplicate isomorphic subgraphs vs. their mean
    # assembly index calculation time as a scatter plot.
    fig, ax = plt.subplots(dpi=300, facecolor='w', tight_layout=True)
    ax.scatter("iso_sub_pairs", "time", data=mol_dfs["naive"],
               s=3, color=cmc.batlow(0.2), label="bb-naive")
    ax.scatter("iso_sub_pairs", "time", data=mol_dfs["logbound"],
               s=3, color=cmc.batlow(0.4), label="bb-logbound")
    ax.scatter("iso_sub_pairs", "time", data=mol_dfs["intbound"],
               s=3, color=cmc.batlow(0.6), label="bb-intbound")
    ax.scatter("iso_sub_pairs", "time", data=mol_dfs["allbounds"],
               s=3, color=cmc.batlow(0.8), label="bb-allbounds")
    ax.set(xlabel="# Disjoint Isomorphic Subgraph Pairs", yscale='log',
           ylabel="Assembly Index Calculation Time (seconds, log scale)")
    ax.legend(loc='best', fontsize='small')
    fig.savefig(osp.join(args.figs_path, "jossplot.pdf"))
