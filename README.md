# ORCA-JOSS

This repository contains the source files for our [Journal of Open Source Science](https://joss.theoj.org/) manuscript on `ORCA`: Open, Reproducible Calculation of Assembly Indices.
This is a collaboration in the Biodesign Center for Biocomputing, Security and Society at Arizona State University involving&mdash;in alphabetical order by last name, with PIs at the end&mdash;Sean Bergen, Devendra Parkar, Garrett Parzych, Olivia Smith, Devansh Vimal, Joshua J. Daymude, and Cole Mathis.

The associated codebase is contained in this repository as a submodule and its current version can be found [on GitHub](https://github.com/DaymudeLab/ORCA).

### Benchmarking `ORCA`

Set up the `ORCA` benchmark by copying the rust benchmark file into the submodule and then going to the appropriate directory:

```
cp scripts/benchmark.rs ORCA/benches/
cd ORCA
```

Then run the benchmark with

```
cargo bench gdb13
cargo bench gdb17
```

### Genrating Plots for `ORCA`

Set up the `ORCA` benchmark by copying the rust benchmark file into the submodule and then going to the appropriate directory:

```
cp scripts/benchmark.rs ORCA/benches/
cd ORCA
```

Then run the benchmark with

```
cargo bench plot
```

Finally run the python script using following command to generate a comparision scatter plot for different molecules and methods used as `method_comparison.svg` in  `figures` folder:

```
python3 scripts/orca_fig.py
```
