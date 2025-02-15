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
date: 3 February 2025
bibliography: paper.bib
---

# Summary

We present `ORCA` (**O**pen, **R**eproducible **C**omputation of **A**ssembly Indices), a Rust package for computing *assembly indices* of biochemical structures.
This is a key complexity measure of *assembly theory*, a recent theoretical framework qunatifying selection across diverse systems, most importantly chemistry [@Walker2024-experimentallymeasured; @Sharma2023-assemblytheory].
`ORCA` is designed for researchers and practitioners alike, providing (i) extensible, high-performance implementations of assembly index calculation algorithms, (ii) comprehensive benchmarks against which current and future algorithmic improvements can be tested, and (iii) Python bindings and `RDKit`-compatible data loaders to support integration with existing computational pipelines.



# Background

*Assembly theory* (AT) is a recently developed body of theoretical and empirical work focused on characterizing selection in chemical systems [@Sharma2023-assemblytheory; @Walker2024-experimentallymeasured].
Objects are defined in AT as entities that are finite, distinguishable, decomposable, and persistent in time.
AT characterizes objects based on their *assembly index*, the minimum number of recursive subcontructions required to construct the object starting from a given set of building blocks [@Jirasek2024-investigatingquantifying; @Seet2024-rapidcomputation].
The most commonly studied application domain of AT to date is molecular chemistry, where bonds act as the basic building blocks and the quantity of interest is the *molecular assembly index* (MA); see \autoref{fig:assemblyindex} for an example.
It has previously been shown that MA can be measured for covalently-bonded molecules using standard analytical techniques such as tandem mass spectrometry as well as infrared and nuclear magnetic resonance spectroscopy [@Jirasek2024-investigatingquantifying], enabling a novel approach to life detection based on AT [@Marshall2021-identifyingmolecules].
Beyond life detection, AT and MA have been proposed in methods to generate novel therapeutic drugs, identify environmental pollutants, and gain new insights into evolutionary history by inferring relationships directly from metabolomic data [@Liu2021-exploringmapping; @Kahana2024-constructingmolecular].

![*Assembly Pathways for Anthracene*. Starting with bonds as building blocks (yellow), a joining operation yields progressively larger structures by combining any two compatible structures that have already been constructed (arrows). These intermediate structures must obey valence rules but otherwise do not have to be physically accessible or chemically synthesizable. There may be many assembly pathways from building blocks to a target structure&mdash;in this case, Anthracene (green)&mdash;but the length of any shortest such pathway (blue) is that structure's assembly index.\label{fig:assemblyindex}](figures/anthracene.pdf){ width=100% }



# Statement of Need

Despite AT's promising applications, computing MA efficiently remains a challenge.
In general, exact MA calculation is an NP-hard problem [@Kempes2024-assemblytheory]; i.e., the necessary computing resources are likely to grow exponentially with a molecule's number of bonds.
Previous software to compute MA have been closed-source, platform-dependent, or written in languages rarely used by the broader scientific community.
For example, the original software to compute a split-branch approximation of MA (an upper bound on the exact value) was written in C++ and depended on the MSVC compiler, making it difficult to deploy to non-Windows machines [@Marshall2021-identifyingmolecules].
The more recent `AssemblyGo` implementation computes MA exactly, but is written in Go, yielding worse performance than alternatives and posing an accessibility barrier for most scientific practitioners who are unfamiliar with the language [@Jirasek2024-investigatingquantifying].
Finally, the latest `AssemblyCPP` implementation is again written in C++ but is closed-source, prohibiting its use and verification by the community [@Seet2024-rapidcomputation].

With `ORCA`, we provide a high-performance, cross-platform Rust package for fast MA calculation while also providing Python bindings for key functionality, offering the best efficiency without sacrificing accessibility.
We chose Rust for its advantages of cross-platform support, memory-safety, performant runtime, convenient parallelism, and integrated testing and documentation [@Perkel2020-whyscientists].
By including test and benchmark suites, we also lay a foundation for fair, reproducible comparisons of future algorithmic improvements and new techniques.



# Design

**TODO**: CM & JD are leaning toward folding this section into the previous paragraph. One issue about the ML sentence, which is the main new idea here, is that it suggests the presence of a lot of data to train on, which we're not exactly providing here.

`ORCA` is not a single algorithmic implementation of assembly index calculations; rather, it is a library that can be used to implement a diversity of algorithmic approaches.
As AT matures, we expect new algorithmic implementations will develop. 
The design philosophy behind `ORCA` is to provide a source of ground truth and robust comparison for future implementations.
These could include novel methods for exact calculation of assembly indices, or they may be approximation methods that leverage advances in machine learning [@Gebhard2022-inferringmolecular; @Marshall2021-identifyingmolecules].



# Functionality and Examples

`ORCA` provides a stand-alone executable that can be compiled via `cargo`, as well as libraries for Rust and Python. 
The primary function of ORCA is to enable users to compute assembly indices for molecules of interest, and benchmark algorithmic changes against a standard suite of molecules.
Here we provide examples of how to use the exectuable, how to use the Python library, and how to run tests and benchmarks.

## Building and running the executable
Building the executable is handled by `cargo`. Inside the main repository you simply need to run:
```shell
cargo run tests/mol/aspirin.mol
```
This will build the exectuable and compute the assembly index for aspirin (the given input). 
For repeated use you do not need to rebuild, instead the exectuable can built using 
```shell
cargo build -r
``` 
This will generate a binary `target/release/orca` which can be used. For example 
```shell
./target/release/orca testt/input/aspirin.mol 
```

## Installing and using the Python library
[Install instructions]

Once the library is installed you can use it to compute assembly indices directly on RDKIT `Mol` objects in the following way:
```python
from rdkit import Chem as Chem

aspirin_mol = Chem.MolFromSmiles("O=C(C)Oc1ccccc1C(=O)O")
orca.compute_assembly_index(aspirin_mol) # 8
```


# Tests and Benchmarks

`ORCA` includes test and benchmark suites for software validation and performance evaluation, respectively.
Both suites are backed by curated reference datasets representing different classes of molecules, arranged roughly in order of increasing molecular size and complexity:

- `checks`: 15 named molecules (e.g., anthracene, aspirin, caffeine, morphine) primarily used for rapid testing and sanity checking.
These molecules' number of heavy atoms range from 5&ndash;28 and have MA from 3&ndash;18.
- `gdb13_1201`: 1,201 small, organic molecular structures sampled from GDB-13, a database of enumerated chemical structures containing Carbon, Hydrogen, Nitrogen, Oxygen, Sulfur, and Chlorine that are constrained only by valence rules and quantum mechanics but may not be chemically stable or synthesizable [@Reymond2015-chemicalspace].
Our sample includes all 201 molecules in GDB-13 with 4&ndash;5 heavy atoms and 200 randomly sampled molecules for each number of heavy atoms from 6&ndash;10.
These molecules' MA range from 2&ndash;9.
- `gdb17_800`: 800 organic molecular structures sampled from the larger GDB-17 database, which includes additional nuclei beyond GDB-13 such as the halogens Flourine and Iodine [@Reymond2015-chemicalspace].
Compared to GDB-13, these molecules are typically larger and represent more structural diversity.
Our sample includes 200 randomly sampled molecules for each number of heavy atoms from 14&ndash;17.
These molecules' MA range from 5&ndash;15.
- `coconut_220`: 220 natural products sampled from the COCONUT database [@Sorokina2021-coconutonline].
Natural products (or secondary metabolites) are a rich source of evolved chemical complexity, often exhibiting drug-like properties.
Subsets of this database were used to benchmark recent algorithmic progress in [@Seet2024-rapidcomputation]. 
Our sample includes 20 randomly sampled molecules for each number of heavy atoms from 15&ndash;25.
These molecules' MA range from 5&ndash;20.

We curated these reference datasets for their structural diversity and approachable runtime on commodity hardware.
Larger, more demanding datasets can be easily added as needed.

The `ORCA` test suite contains unit tests validating internal functionality and database tests checking the calculation of correct assembly indices for all molecules in any of our reference datasets.
Each reference dataset contains an `ma-index.csv` file with ground truth assembly indices.
Incorrect calculations are flagged for developer review.

Our benchmark suite evaluates `ORCA` performance by running repeated assembly index calculations over individual molecules or entire reference datasets.
We leverage the `criterion` package for Rust to automatically collect detailed timing statistics, charts, and estimates of performance improvements and regressions.
As an example, \autoref{tab:benchtimes} shows `ORCA` performance across our three reference datasets against that of `AssemblyGo` [@Jirasek2024-investigatingquantifying], another recent implementation written in Go.
Depending on the dataset and choice of `ORCA` algorithm, `ORCA` outperforms `AssemblyGo` by 16.9&ndash;213.8x.
The 16.9&ndash;18.1x speedup on the `gdb13_1201` dataset most clearly represents the efficiency of Rust over Go, since those molecules are so small that they barely benefit from algorithmic improvements.
Algorithmic improvements such as branch-and-bound with an integer addition chain bound [@Seet2024-rapidcomputation] over the trivial logarithmic bound [@Jirasek2024-investigatingquantifying] or no bound at all ("naive") yield more dramatic speedups for larger molecules, like those in `gdb17_800`.
This internal comparison showcases `ORCA` as a framework capable of comparing multiple algorithmic approaches on equal footing, free of differences in underlying datasets or language-specific efficiency issues.

: \label{tab:benchtimes} Benchmark execution times for `AssemblyGo` [@Jirasek2024-investigatingquantifying] vs. `ORCA`.
`AssemblyGo` uses its default parameters.
`ORCA` has three algorithm settings: "naive" which fully enumerates all non-duplicate assembly pathways; "logbound" which improves over "naive" by eliminating any assembly pathways longer than $\log_2b$, where $b$ is the molecule's number of bonds [@Jirasek2024-investigatingquantifying]; and "addbound" which improves over "logbound" by eliminating any assembly pathways longer than a bound provided by an integer addition chain [@Seet2024-rapidcomputation].
The benchmark times the sequential MA calculation of all molecules in a given dataset, excluding the time required to parse and load `.mol` files into internal molecular graph representations.
We repeated the benchmark 100 times on a single CPU for each software&ndash;dataset pair.
All results are reported as mean runtime $\pm$ 95% confidence interval.

| Dataset       | `AssemblyGo`        | `ORCA`-naive        | `ORCA`-logbound     | `ORCA`-addbound     |
| ------- | ----------: | ----------: | ----------: | ----------: |
| `checks`      | TODO                | TODO                | TODO                | TODO                |
| `gdb13_1201`  | 1.943 s $\pm$ 3.28% | 0.115 s $\pm$ 0.07% | 0.114 s $\pm$ 0.07% | 0.107 s $\pm$ 0.02% |
| `gdb17_800`   |  1239 s $\pm$ 0.20% | 38.06 s $\pm$ 0.35% | 19.40 s $\pm$ 0.57% | 5.796 s $\pm$ 0.37% |
| `coconut_220` | TODO                | TODO                | TODO                | TODO                |

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
DV was the primary software developer (architecture, command line interface, molecule representations, unit tests, performance engineering).
GP implemented all bound calculations.
DP and DV implemented the `.mol` file parser and dataset-based benchmarks.
CM implemented the Python interface.
OMS curated all reference datasets and assembly index ground truths with input from CM.
SB and JJD wrote the `AssemblyGo` benchmarks.
JJD conducted and analyzed the benchmarks shown in \autoref{tab:benchtimes}.
DP, SB, and GP produced the molecule-based visualization in \autoref{fig:timescatter}.
JJD and CM wrote the paper.



# Acknowledgements

JJD and GP are supported in part by NSF award CCF-2312537.
DV is supported by the ASU Biodesign Institute.
**TODO**: Other acks?



# References
