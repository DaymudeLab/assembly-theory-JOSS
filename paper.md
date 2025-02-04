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



# Statement of Need

TODO: Why is the software necessary?



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
