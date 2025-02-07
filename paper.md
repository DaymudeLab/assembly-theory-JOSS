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
    equal-contrib: true
    affiliation: "1, 2"
  - name: Devendra Parkar
    orcid: 0009-0009-0133-8875
    equal-contrib: true
    affiliation: "1, 2"
  - name: Garrett Parzych
    orcid: 0009-0008-4789-9603
    equal-contrib: true
    affiliation: "1, 2"
  - name: Olivia Smith
    orcid: 0009-0004-2299-3522
    equal-contrib: true
    affiliation: "1, 3"
  - name: Devansh Vimal
    orcid: 0009-0006-2794-8995
    equal-contrib: true
    affiliation: 1
  - name: Joshua J. Daymude
    orcid: 0000-0001-7294-5626
    equal-contrib: true
    affiliation: "1, 2"
  - name: Cole Mathis
    orcid: 0000-0001-8424-9169
    corresponding: true
    equal-contrib: true
    affiliation: "1, 3"
affiliations:
  - name: Biodesign Center for Biocomputing, Security and Society, Arizona State University, United States
    index: 1
  - name: School of Computing and Augmented Intelligence, Arizona State University, United States
    index: 2
  - name: School of Complex Adaptive Systems, Arizona State University, United States
    index: 3
date: 3 February 2025
bibliography: ref.bib
---

# Summary

We present `ORCA` ([O]{.ul}pen, [R]{.ul}eproducible [C]{.ul}omputation of [A]{.ul}ssembly Indices), a Rust package for computing *assembly indices* and *minimum assembly pathways* of biochemical structures.
These are key complexity measures of *assembly theory*, a recent theoretical framework characterizing how selection occurs across diverse systems, most importantly chemistry `[@Walker2024-experimentallymeasured; @Sharma2023-assemblytheory]`.
`ORCA` is designed for researchers and practitioners alike, providing (*i*) extensible, high-performance implementations of assembly index calculation algorithms, (*ii*) comprehensive benchmarks against which current and future algorithmic improvements can be tested, and (*iii*) Python bindings and `RDKit`-compatible data loaders to support integration with existing computational pipelines.



# Background

*Assembly theory* (AT) is a recently developed body of theoretical and empirical work focused on characterizing selection in chemical systems `[@Sharma2023-assemblytheory; @Walker2024-experimentallymeasured]`.
Objects are defined in AT as entites that are finite, distinguishable, persist in time and decomposable. 
AT characterizes objects based on their *assembly index* (AI). The most interesting application area of AT is molecules and the assembly index of molecules, the *molecular assembly index* is abbreviated (MA). 
The AI of an object is defined as the minimum number of recursive subcontructions required to construct a target structure starting from a given set of building blocks (e.g., bonds for molecules).`[@Jirasek2024-investigatingquantifying; @Seet2024-rapidcomputation]`; see Figure \autoref{fig:assemblyindex} for an example.
It has previously been shown that MA can be measured for covalently-bonded molecules using standard analytical techniques such as tandem mass spectrometry as well as infrared and nuclear magnetic resonance spectroscopy `[@Jirasek2024-investigatingquantifying]`, enabling a novel approach to life detection based on AT `[@Marshall2021-identifyingmolecules]`.
Beyond life detection, AT and MA have been proposed in methods to generate novel therapeutic drugs, identify environmental pollutants, and gain new insights into evolutionary history by inferring relationships directly from metabolomic data `[@Liu2021-exploringmapping; @Kahana2024-constructingmolecular]`.

![*Assembly Pathways for Anthracene*. Starting with bonds as building blocks (yellow), a joining operation yields progressively larger structures by combining any two compatible structures that have already been constructed (arrows). These intermediate structures must obey valence rules but otherwise do not have to be physically accessible or chemically synthesizable. There may be many assembly pathways from building blocks to a target structure&mdash;in this case, Anthracene (green)&mdash;but the length of any shortest such pathway (blue) is that structure's assembly index.\label{fig:assemblyindex}](figures/anthracene.pdf){ width=80% }



# Statement of Need

Despite AT's promising applications, computing MA efficiently remains a challenge.
In general, exact MA calculation is an NP-hard problem `[@Kempes2024-assemblytheory]`; i.e., the necessary computing resources are likely to grow exponentially with the target structure's size.
Previous software to compute assembly indices have been closed-source, platform-dependent, or written in languages rarely used by the broader scientific community.
For example, the original software to compute a split-branch approximation of MA (an upper bound on the exact value) was written in C++ and depended on the MSVC compiler, making it difficult to deploy to non-Windows machines `[@Marshall2021-identifyingmolecules]`.
The more recent `AssemblyGo` implementation computes MA exactly, but is written in Go, yielding worse performance than alternatives and posing an accessibility barrier for most scientific practitioners who are unfamiliar with the language `[@Jirasek2024-investigatingquantifying]`.
Finally, the latest `AssemblyCPP` implementation is again written in C++ but is closed-source, prohibiting its use and verification by the community `[@Seet2024-rapidcomputation]`.

With `ORCA`, we provide a high-performance, cross-platform Rust package for fast MA calculation while also providing Python bindings for key functionality, offering the best efficiency without sacrificing accessibility.
By including test and benchmark suites, we additionally lay the foundation for fair, reproducible comparisons of future algorithmic improvements and new techniques.

# Design
`ORCA` is not a single algorithmic implementation of assembly index calculations. 
Rather it is a library that can be used to implement a diversity of algorithmic approaches. 
As AT matures, we expect new algorithmic implementations will develop. 
The design philosophy behind `ORCA` is to provide a source of ground truth, and robust comparison for future implementations.
These could include novel methods for exact calculation of assembly indices, or they may be approximation methods that leverage advances in machine learning `[@Gebhard2022-inferringmolecular, @Marshall2021-identifyingmolecules]`.


# Contributing and Governance



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

## Running tests and benchmarks

**TODO**


# Tests and Benchmarks

`ORCA` includes test and benchmark suites for software validation and performance evaluation, respectively.
Both suites are backed by curated reference datasets representing different classes of molecules, arranged roughly in order of increasing molecular size and complexity:

- `gdb13_1200`: 1,200 small, organic molecular structures sampled from GDB-13, a database of enumerated chemical structures containing Carbon, Hydrogen, Nitrogen, Oxygen, Sulfur, and Chlorine that are constrained only by valence rules and quantum mechanics but may not be chemically stable or synthesizable `[@Reymond2015-chemicalspace]`.
Our sample includes all XX molecules in GDB-13 with XX&ndash;XX heavy atoms and 200 randomly sampled molecules for each number of heavy atoms between XX and XX.
These data contain molecules with assembly indices between XX&ndash;XX.
- `gdb17_1000`: 1,000 organic molecular structures sampled from the larger GDB-17 database, which includes additional nuclei beyond GDB-13, such as the halogens Flourine and Iodine `[@Reymond2015-chemicalspace]`.
Compared to GDB-13, these molecules are typically larger and represent more structural diversity.
Our sample includes **TODO** explain heavy atoms/distribution.
These data contain molecules with assembly indices between XX&ndash;XX.
- `coconut_50`: 50 natural products sampled from the COCONUT database `[@Sorokina2021-coconutonline]`.
Natural products (or secondary metabolites) are a rich source of evolved chemical complexity, often exhibiting drug-like properties.
We selected XX ... [Get info from Olivia]
Subsets of this database were used to benchmark recent algorithmic progress in `[@Seet2024-rapidcomputation]`. 
These data contain molecules with assembly indices between XX&ndash;XX.

We curated these reference datasets for their structural diversity and approachable runtime on commodity hardware.
Larger, more complicated datasets can be easily added as needed.

The `ORCA` test suite contains unit tests validating internal functionality and database tests checking the calculation of correct assembly indices for all molecules in any of our reference datasets.
Each reference dataset contains an `ma-index.csv` file with ground truth assembly indices.
Incorrect calculations are flagged for developer review.

Our benchmark suite evaluates `ORCA` performance by running repeated assembly index calculations over individual molecules or entire reference datasets.
We leverage the `criterion` package for Rust to automatically collect detailed timing statistics, charts, and estimates of performance improvements and regressions.
As an example, \autoref{tab:benchtimes} shows `ORCA` performance across our three reference datasets against that of `AssemblyGo` `[@Jirasek2024-investigatingquantifying]`, another recent implementation written in Go.
As an aside, this showcases `ORCA` as not just one algorithm's implementation, but as a framework capable of comparing multiple algorithmic approaches on equal footing, free of differences in underlying datasets or language-specific efficiency issues.

: \label{tab:benchtimes}**TODO**: Caption for table.

| Dataset         | `AssemblyGo` [-@Jirasek2024-investigatingquantifying] | `ORCA`-nobounds | `ORCA`-logbound | `ORCA`-seetbound |
| :-------------- | ----------------------------------------------------: | --------------: | --------------: | ---------------: |
| `gdb13_1200`    | 2.409 ± 2% s                                          |                 |                 |                  |
| `gdb17_1000`    |                                                       |                 |                 |                  |
| `coconut_50`    |                                                       |                 |                 |                  |



# Formatting Considerations

$\LaTeX$ formatting works as usual, with single dollar signs for inline math, double dollar signs for self-standing equations, and begin-end equation environments for labeled equations.

Citations look like this:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }



# Acknowledgements

J.J.D. and G.P. are supported in part by NSF award CCF-2312537.
**TODO**: Other acks?


# References
