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
Among AT's key measures of structural complexity, the *(molecular) assembly index* (MA) quantifies the length of the shortest *assembly pathway* for a given structure, i.e., the minimum number of recursive subcontructions required to construct a target structure starting from a given set of building blocks (e.g., bonds for molecules) `[@Jirasek2024-investigatingquantifying; @Seet2024-rapidcomputation]`; see Figure \autoref{fig:assemblyindex} for an example.
It has previously been shown that MA can be measured for covalently-bonded molecules using standard analytical techniques such as tandem mass spectrometry as well as infrared and nuclear magnetic resonance spectroscopy `[@Jirasek2024-investigatingquantifying]`, enabling a novel approach to life detection based on AT `[@Marshall2021-identifyingmolecules]`.
Beyond life detection, AT and MA have been proposed in methods to generate novel therapeutic drugs, identify environmental pollutants, and gain new insights into phylogenomics by inferring evolutionary relationships directly from metabolomic data `[@Liu2021-exploringmapping; @Kahana2024-constructingmolecular]`.

![*Assembly Pathways for Anthracene*. Starting with bonds as building blocks (yellow), a joining operation yields progressively larger structures by combining any two compatible structures that have already been constructed (arrows). These intermediate structures must obey valence rules but otherwise do not have to be physically accessible or chemically synthesizable. There may be many assembly pathways from building blocks to a target structure&mdash;in this case, Anthracene (green)&mdash;but the length of any shortest such pathway (blue) is that structure's assembly index.\label{fig:assemblyindex}](figures/anthracene.pdf){ width=80% }



# Statement of Need

Despite AT's promising applications, computing MA efficiently remains a challenge.
In general, exact MA calculation is an NP-hard problem `[@Kempes2024-assemblytheory]`; i.e., the necessary computing resources are likely to grow exponentially with the target structure's size.
Previous software to compute assembly indices have been closed-source, platform-dependent, or written in languages rarely used by the broader scientific community.
For example, the original software to compute a split-branch approximation of MA (an upper bound on the exact value) was written in C++ and depended on the MSVC compiler, making it difficult to deploy to non-Windows machines `[@Marshall2021-identifyingmolecules]`.
The more recent `AssemblyGo` implementation computes MA exactly but is written in Go, yielding worse performance than alternatives and posing an accessibility barrier for most scientific practitioners who are unfamiliar with the language `[@Jirasek2024-investigatingquantifying]`.
Finally, the latest `AssemblyCPP` implementation is again written in C++ but is closed-source, prohibiting its use and verification by the community `[@Seet2024-rapidcomputation]`.

With `ORCA`, we provide a high-performance, cross-platform Rust package for fast MA calculation while also providing Python bindings for key functionality, offering the best efficiency without sacrificing accessibility.
By including test and benchmark suites, we additionally lay the foundation for fair, reproducible comparisons of future algorithmic improvements and new techniques.



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

TODO: Fill these in, if we have any.



# References
