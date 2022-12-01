# Lattice QCD for Novices

## Introduction

This repository contains the source codes and worked experiments of the exercises proposed in "Lattice QCD for Novices" by G.P. Lepage.
The original paper can be downloaded at the following link on [arXiv](https://arxiv.org/abs/hep-lat/0506036v1). The codes are mainly written in C++,
with a small component in Python (to produce plots) and a few Bash scripts and Makefiles to build and run the applications. The documentation is provided in html format.

## Structure

The repository is organized in the following main directories:
1. Vegas_Integration: source code for the first exercise of the paper on the Vegas integration of a 1D quantum harmonic oscillator
2. 1D_Path_Integration: source code for the second set of exercises on the Metropolis algorithm of a 1D quantum system
3. QCD: source code for the third set of exercises on the Metropolis algorithm for gluonic path integrals
4. Experiments: this directory contains the worked exercises using the codes described in the points from 1 to 3

## Documentation

More detailed information on the code may be found in the code directories (Vegas_Integration, 1D_Path_Integration, QCD) using the documentation.html link as a html, Doxygen-generated documentation.
A PowerPoint presentation is included in this directory, which summarises the more general theoretical aspects as presented in Lepage's article. In the presentation, the exercises and their solutions
are described, too.

---

## Bibliography

Lepage, G. Peter, Lattice QCD for Novices, 2005, arXiv, [arXiv:hep-lat/0506036](https://arxiv.org/abs/hep-lat/0506036v1)
