"""
bench_stats: Compute benchmark statistics from Rust (criterion) and Go outputs
"""

import argparse
from collections import defaultdict
from itertools import product
import json
import numpy as np
import os.path as osp
import pandas as pd
import scipy.stats as sps

if __name__ == "__main__":
    # Parse command line arguments.
    parser = argparse.ArgumentParser(description=__doc__)
    def_ago_path = osp.join('assembly_go', 'cmd', 'app', 'datasets_bench.tsv')
    parser.add_argument('-A', '--ago_path', type=str, default=def_ago_path,
                        help='Path to assembly_go benchmark output')
    def_crit_path = osp.join('assembly-theory', 'target', 'criterion', 'datasets')
    parser.add_argument('-C', '--crit_path', type=str, default=def_crit_path,
                        help='Path to criterion datasets benchmark output')
    args = parser.parse_args()

    # Define fixed constants.
    datasets = ['gdb13_1201', 'gdb17_800', 'checks', 'coconut_220']
    at_algs = ['naive', 'logbound', 'intbound', 'allbounds']

    # Create results dict of the form {algo: {dataset: (mean, 95% conf.)}}.
    results = defaultdict(dict)

    # Load assembly_go benchmark outputs.
    ago_df = pd.read_csv(args.ago_path, sep='\t', header=None,
                         names=['bench', 'reps', 'time'], engine='python',
                         skiprows=4, skipfooter=2).dropna()
    ago_df['dataset'] = ago_df.bench.apply(
            lambda x: x.strip().split('/')[1].split('-')[0])
    ago_df.time = ago_df.time.apply(lambda x: float(x.split()[0]))

    # Compute assembly_go means and 95% confidence intervals.
    for dataset in datasets:
        # Skip very slow benchmarks that we don't do multiple samples for.
        if dataset == 'coconut_220':
            continue

        times = np.array(ago_df.loc[ago_df['dataset'] == dataset].time)
        mean_time = times.mean()
        conf = sps.t.interval(0.95, len(times)-1, loc=mean_time,
                              scale=sps.sem(times))
        conf_perc = (conf[1] - conf[0]) / 2 / mean_time * 100
        results['assembly_go'][dataset] = (mean_time, conf_perc)

    # Do the same for assembly-theory's three algorithm variants.
    for dataset, at_alg in product(datasets, at_algs):
        # Skip very slow benchmarks that we don't do multiple samples for.
        if dataset == 'coconut_220' and at_alg in ['naive', 'logbound']:
            continue

        with open(osp.join(args.crit_path, dataset, at_alg, 'new',
                           'estimates.json'), 'r') as f:
            stats = json.load(f)
            mean_time = stats['mean']['point_estimate']
            conf_low = stats['mean']['confidence_interval']['lower_bound']
            conf_high = stats['mean']['confidence_interval']['upper_bound']
            conf_perc = (conf_high - conf_low) / 2 / mean_time * 100
            results[f'at-{at_alg}'][dataset] = (mean_time, conf_perc)

    # Print results.
    results_df = pd.DataFrame(results)
    results_df = results_df.map(lambda x: f'{x[0] / 1e9:.3f} s ± {x[1]:.2f}%'
                                if not pd.isna(x) else 'NaN')
    print(results_df)
