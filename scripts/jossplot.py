"""
jossplot: Plot molecules' number of duplicate isomorphic subgraphs vs. their
          average ORCA assembly index calculation time
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
    def_crit_path = osp.join("ORCA", "target", "criterion", "jossplot")
    parser.add_argument('-C', '--crit_path', type=str, default=def_crit_path,
                        help='Path to jossplot criterion output directory')
    parser.add_argument('-F', '--figs_path', type=str, default='figures',
                        help='Path to figures directory')
    args = parser.parse_args()

    # Load molecule duplicate isomorphic subgraph counts.
    dup_iso_subs = {}
    with open(osp.join(args.crit_path, "jossplot.csv"), 'r') as f:
        csv_reader = csv.reader(f)
        for row in csv_reader:
            dup_iso_subs[row[0]] = int(row[1])

    # Load assembly index calculation times per molecule per bound option.
    mol_dfs = {}
    for bound in ["naive", "logbound", "addbound"]:
        # Set up results dict to turn into a DataFrame.
        results = {"mol_name": [], "dup_iso_subs": [], "time": []}

        # Gather per molecule data.
        for mol_name, mol_path in [(d.name, d.path) for d in
                                   os.scandir(osp.join(args.crit_path, bound))
                                   if d.name[:5] == 'gdb17']:
            # Set molecule name and duplicate isomorphic subgraph count.
            results["mol_name"].append(mol_name)
            results["dup_iso_subs"].append(dup_iso_subs[mol_name])

            # Retrieve and set assembly index calculation time in seconds.
            with open(osp.join(mol_path, "new", "estimates.json")) as f:
                mean_time = float(json.load(f)["mean"]["point_estimate"])
                results["time"].append(mean_time)
                results["time"][-1] /= 1e9

        # Convert results dict to DataFrame.
        mol_dfs[bound] = pd.DataFrame(data=results)
        mol_dfs[bound].sort_values(by="dup_iso_subs", inplace=True)

    # Plot molecules' number of duplicate isomorphic subgraphs vs. their mean
    # assembly index calculation time as a scatter plot.
    fig, ax = plt.subplots(dpi=300, facecolor='w', tight_layout=True)
    ax.scatter("dup_iso_subs", "time", data=mol_dfs["naive"],
               s=3, color=cmc.batlow(0.2), label="Naive")
    ax.scatter("dup_iso_subs", "time", data=mol_dfs["logbound"],
               s=3, color=cmc.batlow(0.5), label="Log. Bound")
    ax.scatter("dup_iso_subs", "time", data=mol_dfs["addbound"],
               s=3, color=cmc.batlow(0.8), label="Int. Add. Bound")
    ax.set(xlabel="# Duplicate Isomorphic Subgraphs", yscale='log',
           ylabel="ORCA Assembly Index Calculation Time (seconds, log scale)")
    ax.legend(loc='best', fontsize='small')
    fig.savefig(osp.join(args.figs_path, "jossplot.png"))
