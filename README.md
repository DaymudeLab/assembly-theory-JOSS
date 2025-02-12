# ORCA-JOSS

This repository contains the source files for our [Journal of Open Source Science](https://joss.theoj.org/) manuscript on `ORCA`: Open, Reproducible Calculation of Assembly Indices (see the [GitHub repository](https://github.com/DaymudeLab/ORCA)).
This is a collaboration in the Biodesign Center for Biocomputing, Security and Society at Arizona State University involving&mdash;in alphabetical order by last name, with PIs at the end&mdash;Sean Bergen, Devendra Parkar, Garrett Parzych, Olivia Smith, Devansh Vimal, Joshua J. Daymude, and Cole Mathis.


## Using This Repository

The primary function of this repository is to house the source files for our JOSS manuscript: `paper.md`, `paper.bib`, and `figures/`.
Instructions for compiling a preview of the paper can be found in the [JOSS docs](https://joss.readthedocs.io/en/latest/paper.html#locally).

If you additionally wish to reproduce our benchmarks and associated figures (instructions below), you will need Rust, Go, Python, and Git LFS.
After cloning this repository, run the following to also pull the submodules:

```shell
git submodule init
git submodule update
```


## Benchmarking Method and Instructions for Reproduction

This paper includes a small benchmark of our `ORCA` assembly index calculations against those of [`assembly_go`](https://github.com/croningp/assembly_go), an earlier implementation written in Go ([Jirasek et al., 2024](https://doi.org/10.1021/acscentsci.4c00120)).
For the purposes of reproducibilty, this repository includes the versions of `ORCA` and `assembly_go` that we benchmarked as submodules.
The benchmarks measure only the runtime required to compute assembly pathways and indices, but exclude all setup and teardown (e.g., loading `.mol` files into internal molecule/graph representations).
The molecule datasets used for benchmarking are described in `paper.md`.


### Benchmarking `ORCA`

Set up the `ORCA` benchmark by copying the Rust benchmark file into the appropriate submodule and then going to the corresponding directory:

```shell
cp scripts/benchmark.rs ORCA/benches/
cd ORCA
```

Then run the benchmark with

```shell
cargo bench datasets
```

The results can be viewed in the report created at `ORCA/target/criterion/datasets/report/index.html`.


### Benchmarking `assembly_go`

Set up the `assembly_go` benchmark by copying the Go benchmark file into the appropriate submodule and then going to the corresponding directory:

```shell
cp scripts/main_test.go assembly_go/cmd/app/
cd assembly_go/cmd/app
```

Then run the benchmark with

```shell
go test -bench=. -cpu=1 -count=<iters> -timeout=0 > assembly_go_bench.txt
```

where `<iters>` is replaced by the number of iterations you want to run the benchmark and average the times over.
For our paper, we used `-count=100`.

To compute statistics, install `benchstat` with:

```shell
go install golang.org/x/perf/cmd/benchstat@latest
```

Then run:

```shell
benchstat assembly_go_bench.txt
```

> [!TIP]
> If this fails, you probably forgot to [set your `GOPATH`](https://go.dev/wiki/SettingGOPATH).

> [!NOTE]
> `benchstat` computes the median benchmark time and the 95% confidence interval for the median.
> This is not directly comparable to what we report in the paper, which is the *mean* benchmark time and 95% confidence interval for the mean.


## Generating Plots for `ORCA`

Our manuscript includes a scatterplot of molecules' numbers of duplicate isomorphic subgraphs (a rough estimate of their complexity) vs. their mean `ORCA` assembly index calculation time.
Instructions for reproducing that figure on your own hardware are below.

We use the [`uv`](https://docs.astral.sh/uv/) tool to manage Python environments.
[Install it](https://docs.astral.sh/uv/getting-started/installation/) and then run the following to get all dependencies:

```shell
uv sync
```

Copy the Rust benchmark file into the appropriate submodule (if you haven't already) and then go to the corresponding directory:

```shell
cp scripts/benchmark.rs ORCA/benches/
cd ORCA
```

Then run the benchmark generating the data for the plot with

```shell
cargo bench jossplot
```

Finally, come back to this directory and run the Python plotting script to generate `figures/jossplot.png`:

```shell
cd ..
uv run scripts/jossplot.py
```
