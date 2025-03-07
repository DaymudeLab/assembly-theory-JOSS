---
title: 'ORCA: Open, Reproducible Calculation of Assembly Indices'
tags:
  - assembly theory
  - biochemistry
  - astrobiology
  - Rust
authors:
  - name: Sean Bergen
    orcid: 0009-0004-3570-5120
    affiliation: "1, 2"
  - name: Devendra Parkar
    orcid: 0009-0009-0133-8875
    affiliation: "1, 2"
  - name: Garrett Parzych
    orcid: 0009-0008-4789-9603
    affiliation: "1, 2"
  - name: Olivia M. Smith
    orcid: 0009-0004-2299-3522
    affiliation: "1, 3"
  - name: Devansh Vimal
    orcid: 0009-0006-2794-8995
    affiliation: 1
  - name: Joshua J. Daymude
    orcid: 0000-0001-7294-5626
    affiliation: "1, 2"
  - name: Cole Mathis
    orcid: 0000-0001-8424-9169
    corresponding: true
    affiliation: "1, 3"
affiliations:
  - name: Biodesign Center for Biocomputing, Security and Society, Arizona State University, United States
    index: 1
  - name: School of Computing and Augmented Intelligence, Arizona State University, United States
    index: 2
  - name: School of Complex Adaptive Systems, Arizona State University, United States
    index: 3
date: 6 March 2025
bibliography: paper.bib
---

# Summary

We present `ORCA` (**O**pen, **R**eproducible **C**omputation of **A**ssembly Indices), a Rust package for computing *assembly indices* of covalently bonded molecular structures.
This is a key complexity measure of *assembly theory*, a recent theoretical framework qunatifying selection across diverse systems, most importantly chemistry [@Walker2024-experimentallymeasured; @Sharma2023-assemblytheory].
`ORCA` is designed for researchers and practitioners alike, providing (i) extensible, high-performance implementations of assembly index calculation algorithms, (ii) comprehensive benchmarks against which current and future algorithmic improvements can be tested, and (iii) Python bindings and `RDKit`-compatible data loaders to support integration with existing computational pipelines.



# Background

*Assembly theory* (AT) is a recently developed body of theoretical and empirical work focused on characterizing selection in diverse physical systems, most importantly in chemistry [@Sharma2023-assemblytheory; @Walker2024-experimentallymeasured].
Objects are defined in AT as entities that are finite, distinguishable, decomposable, and persistent in time.
AT characterizes objects based on their *assembly index*, the minimum number of recursive subcontructions required to construct the object starting from a given set of building blocks [@Jirasek2024-investigatingquantifying; @Seet2024-rapidcomputation].
The most commonly studied application domain of AT to date is molecular chemistry, where bonds act as the basic building blocks and the quantity of interest is the *molecular assembly index* (MA); see \autoref{fig:assemblyindex} for an example.
It has previously been shown that MA can be measured for covalently-bonded molecules using standard analytical techniques such as tandem mass spectrometry as well as infrared and nuclear magnetic resonance spectroscopy [@Jirasek2024-investigatingquantifying], enabling a novel approach to life detection based on AT [@Marshall2021-identifyingmolecules].
Beyond life detection, AT and MA have been proposed in methods to generate novel therapeutic drugs, identify environmental pollutants, and gain new insights into evolutionary history by inferring relationships directly from metabolomic data [@Liu2021-exploringmapping; @Kahana2024-constructingmolecular].

![*Assembly Pathways for Anthracene*. Starting with bonds as building blocks (yellow), a joining operation yields progressively larger structures by combining any two compatible structures that have already been constructed (arrows). These intermediate structures must obey valence rules but otherwise do not have to be physically accessible or chemically synthesizable. There may be many assembly pathways from building blocks to a target structure&mdash;in this case, Anthracene (green)&mdash;but the length of any shortest such pathway (blue) is that structure's assembly index.\label{fig:assemblyindex}](figures/anthracene.pdf){ width=100% }



# Statement of Need

Despite AT's promising applications, computing MA efficiently remains a challenge.
In general, exact MA calculation is an NP-hard problem [@Kempes2024-assemblytheory]; i.e., the necessary computing resources are likely to grow exponentially with a molecule's number of bonds.
Previous software to compute MA have been approximate, closed-source, platform-dependent, or written in languages rarely used by the broader scientific community.
For example, the original software to compute a split-branch approximation of MA (an upper bound on the exact value) was written in C++ and depended on the MSVC compiler, making it difficult to deploy to non-Windows machines [@Marshall2021-identifyingmolecules].
Machine learning methods only provide approximate MA values [@Gebhard2022-inferringmolecular].
The more recent `AssemblyGo` implementation computes MA exactly but is written in Go, yielding worse performance and posing an accessibility barrier for most scientists who are unfamiliar with the language [@Jirasek2024-investigatingquantifying].
Finally, the latest `AssemblyCPP` implementation has achieved significant performance milestones through a branch-and-bound approach. It is again written in C++ but is not available for public use, prohibiting its use and verification by the community [@Seet2024-rapidcomputation].

With `ORCA`, we provide a high-performance, cross-platform Rust package for fast MA calculation while also providing Python bindings for key functionality, offering the best efficiency without sacrificing accessibility.
We chose Rust for its advantages of cross-platform support, memory-safety, performant runtime, convenient parallelism, and integrated testing and documentation [@Perkel2020-whyscientists].
By including test and benchmark suites, we also lay a foundation for fair, reproducible comparisons of future algorithmic improvements and new techniques.



# Design and Current Algorithms

`ORCA` is not a single algorithmic implementation of assembly index calculations; rather, it is a framework and source of ground truth within which a diversity of algorithmic approaches can be validated and compared.
We purposely designed `ORCA` with a modular algorithm interface and data structures that can be easily extended to handle new algorithmic developments introduced as AT matures.

Currently, `ORCA` implements several top-down, branch-and-bound algorithm variants in which a molecule is recursively fragmented and the MA of smaller fragments are used to determine the MA of their parents [@Marshall2021-identifyingmolecules; @Jirasek2024-investigatingquantifying; @Seet2024-rapidcomputation].
We briefly summarize these algorithms below, but emphasize that `ORCA` is not limited to this top-down, recursive approach.

- `ORCA`-naive fully enumerates and searches all non-duplicate assembly pathways in an efficient order.
- `ORCA`-logbound improves over the naive method by eliminating any assembly pathways longer than $\log_2b$, where $b$ is the molecule's number of bonds [@Jirasek2024-investigatingquantifying].
- `ORCA`-intbound improves over the logarithmic bound by eliminating any assembly pathways longer than a bound provided by an integer addition chain [@Seet2024-rapidcomputation].
- `ORCA`-allbounds simultaneously applies the previous integer addition chain bound and a novel bound provided by a vector addition chain.



# Functionality and Usage Examples

`ORCA` can be used to compute assembly indices as a standalone executable, as a library imported by other Rust code, or via a Python interface.
Here, we provide usage examples of each; in the next section, we demonstrate testing and benchmarking functionality.


## Building and Running the Executable

Rust provides the `cargo` build system and package manager for dependency management, compilation, packaging, and versioning.
To build the standalone executable, run:

```shell
cargo build --release
```

This creates an optimized, portable, standalone executable named `target/release/orca`.
It takes as input a path to a `.mol` file and returns that molecule's integer assembly index:

```shell
> ./target/release/orca data/checks/anthracene.mol
6
```

Running with the `--verbose` flag provides additional information, including the input molecule's *number of disjoint, isomorphic subgraph pairs* (i.e., the number of times any molecular substructure is repeated inside the molecule) and the size of the top-down algorithm's *search space* (i.e., its total number of recursive calls).

```shell
> ./target/release/orca data/checks/anthracene.mol --verbose
Assembly Index: 6
Duplicate subgraph pairs: 406
Search Space: 3143
```

By default, `ORCA` uses its fastest algorithm for assembly index calculation (currently `ORCA`-allbounds, see the previous section).
To use a specific bound or disable bounds altogether, set the `--bounds` or `--no-bounds` flags:

```shell
# ORCA-naive, no bounds
./target/release/orca <molpath> --no-bounds

# ORCA-logbound, only logarithmic bound (Jirasek et al., 2024)
./target/release/orca <molpath> --bounds log

# ORCA-intbound, only integer addition chain bound (Seet et al., 2024)
./target/release/orca <molpath> --bounds int-chain

# ORCA-allbounds, both integer and vector addition chain bounds
./target/release/orca <molpath> --bounds int-chain vec-chain
```

Finally, the `--molecule-info` flag prints the molecule's graph representation as a vertex and edge list, the `--help` flag prints a guide to this command line interface, and the `--version` flag prints the current `ORCA` version.


## Installing and using PyORCA

The python library uses `maturin` as a build tool. This needs to be run in a virtual environment. Use the following commands to build and install the library:
```shell
pip install maturin
maturin develop
```

PyORCA computes the assembly index of molecules using RDKit's `Mol` class. Here's a basic example:

```python
import pyorca
from rdkit import Chem

anthracene = Chem.MolFromSmiles("c1ccc2cc3ccccc3cc2c1")
pyorca.molecular_assembly(anthracene)  # 6
```

## Core Functions  

`pyorca` provides three main functions:

- **`molecular_assembly(mol: Chem.Mol, bounds: set[str] = None, no_bounds: bool = False, timeout: int = None, serial: bool = False) -> int`**  
  Computes the assembly index of a given molecule.
  - `timeout` (in seconds) sets a limit on computation time, raising a `TimeoutError` if exceeded.  
  - `serial=True` forces a serial execution mode, mainly useful for debugging.


- **`molecular_assembly_verbose(mol: Chem.Mol, bounds: set[str] = None, no_bounds: bool = False, timeout: int = None, serial: bool = False) -> dict`**  
  Returns additional details, including the number of duplicated isomorphic subgraphs (`duplicates`) and the size of the search space (`space`).  
  - `timeout` (in seconds) sets a limit on computation time, raising a `TimeoutError` if exceeded.  
  - `serial=True` forces a serial execution mode, mainly useful for debugging.

- **`molecule_info(mol: Chem.Mol) -> str`**  
  Returns a string representation of the molecule’s atom and bond structure for debugging.


# Tests and Benchmarks

`ORCA` includes test and benchmark suites for software validation and performance evaluation, respectively.
Both suites are backed by curated reference datasets representing different classes of molecules, arranged roughly in order of increasing molecular size and complexity:

- `gdb13_1201`: 1,201 small, organic molecular structures sampled from GDB-13, a database of enumerated chemical structures containing Carbon, Hydrogen, Nitrogen, Oxygen, Sulfur, and Chlorine that are constrained only by valence rules and quantum mechanics but may not be chemically stable or synthesizable [@Reymond2015-chemicalspace].
Our sample includes all 201 molecules in GDB-13 with 4&ndash;5 heavy atoms and 200 randomly sampled molecules for each number of heavy atoms from 6&ndash;10.
These molecules' MA range from 2&ndash;9.
- `gdb17_800`: 800 organic molecular structures sampled from the larger GDB-17 database, which includes additional nuclei beyond GDB-13 such as the halogens Flourine and Iodine [@Reymond2015-chemicalspace].
Compared to GDB-13, these molecules are typically larger and represent more structural diversity.
Our sample includes 200 randomly sampled molecules for each number of heavy atoms from 14&ndash;17.
These molecules' MA range from 5&ndash;15.
- `checks`: 15 named molecules (e.g., anthracene, aspirin, caffeine, morphine) primarily used for rapid testing.
These molecules' number of heavy atoms range from 5&ndash;28 and have MA from 3&ndash;18.
- `coconut_220`: 220 natural products sampled from the COCONUT database [@Sorokina2021-coconutonline].
Natural products (or secondary metabolites) are a rich source of evolved chemical complexity, often exhibiting drug-like properties.
Subsets of this database were used to benchmark recent algorithmic progress in [@Seet2024-rapidcomputation]. 
Our sample includes 20 randomly sampled molecules for each number of heavy atoms from 15&ndash;25.
These molecules' MA range from 5&ndash;20.

Ground truth MA values were calculated using the publicly avaiable `AssemblyGo` implementation [@Jirasek2024-investigatingquantifying].
We curated these reference datasets for their structural diversity and approachable runtime on commodity hardware.
Larger, more demanding datasets will be easily added as needed.

The `ORCA` test suite contains unit tests validating internal functionality and database tests checking the calculation of correct assembly indices for all molecules in any of our reference datasets.
Each reference dataset contains an `ma-index.csv` file with ground truth assembly indices.
Incorrect calculations are flagged for developer review.

Our benchmark suite evaluates `ORCA` performance by running repeated assembly index calculations over individual molecules or entire reference datasets.
We leverage the `criterion` package for Rust to automatically collect detailed timing statistics, charts, and estimates of performance improvements and regressions.
As an example, \autoref{tab:benchtimes} shows `ORCA` performance across our three reference datasets against that of `AssemblyGo` [@Jirasek2024-investigatingquantifying], another recent implementation written in Go.
Depending on the dataset and choice of `ORCA` algorithm, `ORCA` outperforms `AssemblyGo` by one to three orders of magnitude.
The 8.9&ndash;**TODO**x speedup of `ORCA`-logbound over `AssemblyGo` most clearly represents the efficiency of Rust over Go, since both use the same branch-and-bound approach with a logarithmic bound.
Further algorithmic improvements such as the integer addition chain bound [@Seet2024-rapidcomputation] and our novel vector addition chain bound yield more dramatic speedups for larger molecules, like those up to 965x for `coconut_220`.
**TODO**: If `ORCA`-naive is slower than `AssemblyGo`, explain that here.
This internal comparison showcases `ORCA` as a framework capable of comparing multiple algorithmic approaches on equal footing, free of differences in underlying datasets or language-specific efficiency issues.

: \label{tab:benchtimes} Mean benchmark execution times for `AssemblyGo` [@Jirasek2024-investigatingquantifying] vs. `ORCA` across reference datasets.
The benchmark times the MA calculation of all molecules in a given dataset in sequence, excluding the time required to parse and load `.mol` files into internal molecular graph representations.
`AssemblyGo` uses its default parameters while `ORCA` tests each of its algorithm variants independently.
Each benchmark was run on a Linux machine with a 5.7 GHz Ryzen 9 7950X CPU (16 cores) and 64 GB of memory.
Means are reported over 20 samples per software&ndash;dataset pair, except those marked with an $\ast$ which were prohibitively expensive to run multiple times.

|               | `AssemblyGo`  | `ORCA`-naive  | `ORCA`-logbound | `ORCA`-intbound | `ORCA`-allbounds     |
| --------- | --------: | --------: | -----------: | -----------: | -----------: | 
| `gdb13_1201`  |       1.048 s |       0.118 s |         0.117 s |         0.110 s |          **0.109 s** |
| `gdb17_800`   |     128.396 s |       7.041 s |         5.644 s |         4.476 s |          **4.331 s** |
| `checks`      |     296.584 s |      13.964 s |         2.787 s |         2.191 s |          **2.001 s** |
| `coconut_220` | 6.09 h$^\ast$ |   TODO$^\ast$ |     TODO$^\ast$ |        30.123 s |         **22.911 s** |

If finer-grained timing insights are needed, `ORCA` can also benchmark assembly index calculations for each individual molecule in a reference dataset.
For example, \autoref{fig:timescatter} shows the calculation time of each molecule in `gdb17_800` for three different algorithm settings.
This is useful for teasing out which molecules are "hard" and characterizing where algorithmic improvements make the largest impact.

![*Per-Molecule Benchmark Times*. The mean assembly index calculation time across 100 samples for each molecule (dot) in `gdb17_800` as a function of the molecule's number of duplicate isomorphic subgraphs, a measure roughly correlated with the molecule's size and complexity. The same three `ORCA ` algorithm settings from \autoref{tab:benchtimes} are shown here.\label{fig:timescatter}](figures/jossplot.pdf){ width=75% }



# Availability and Governance

`ORCA` source code and documentation are openly available on [GitHub](https://github.com/DaymudeLab/ORCA).
Following the standard practice for Rust packages, `ORCA` is dual-licensed under the MIT and Apache-2.0 licenses.
External feedback and code contributions are handled through the usual Issues and Pull Request interfaces; guidelines for contributions are listed in `HACKING.md`.
The project's *maintainers* (initially Vimal, Daymude, and Mathis) will govern the project using the committee model: high-level decisions about the project's direction require maintainer consensus, major code changes require majority approval, hotfixes and patches require at least one approval, new maintainers may be added by unanimous decision of the existing maintainers, and existing maintainers may step down with advance notice.



# Author Contributions

GP, DV, and CM formalized the branch-and-bound algorithm design.
GP and SB formalized the integer and vector addition chain bounds.
DV was the primary software developer (architecture, command line interface, molecule representations, unit tests, parallelism, performance engineering).
GP implemented all bound calculations.
DP and DV implemented the `.mol` file parser and dataset-based benchmarks.
CM implemented the Python interface.
OMS curated all reference datasets and assembly index ground truths with input from CM.
SB and JJD wrote the `AssemblyGo` benchmarks.
JJD conducted and analyzed the benchmarks shown in \autoref{tab:benchtimes}.
DP, SB, GP, and JJD produced \autoref{fig:timescatter}.
JJD and CM wrote the paper.



# Acknowledgements

JJD and GP are supported in part by NSF award CCF-2312537.
DV is supported by the ASU Biodesign Institute.
**TODO**: Other acks?



# References
