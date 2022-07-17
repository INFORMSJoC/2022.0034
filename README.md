[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# CacheTest

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[This is a Template](https://doi.org/10.1287/ijoc.2019.0934) by Aharon Ben-Tal and Ernst Roos

## Cite

To cite this software, please cite the [paper](https://doi.org/10.1287/ijoc.2019.0934) using its DOI and the software itself, using the following DOI.

[![DOI](https://zenodo.org/badge/285853815.svg)](https://zenodo.org/badge/latestdoi/285853815)

Below is the BibTex for citing this version of the code.

```
@article{comax2022,
  author =        {Ben-Tal, Aharon and Roos, Ernst},
  publisher =     {INFORMS Journal on Computing},
  title =         {An Algorithm for maximizing a convex function based on its minimum},
  year =          {2022},
  doi =           {10.5281/zenodo.3977566},
  url =           {https://github.com/INFORMSJoC/JoCTemplate},
}  
```

## Description

The goal of this software is to maximize convex functions over a convex feasible set.

The software is written in MATLAB and uses the publically available YALMIP package 
(https://yalmip.github.io/), and the MOSEK solver (https://www.mosek.com/). For the latter,
free academic licenses are available.

## Usage

The software can be ran via the "RunInstances" function that takes a directory name as
input. It will solve all instances in the specified folder. Instances should be provided
as ".mat" files with two variables with the following specification:
* objective - Structure that describes the objective through these fields:
  * f         - Function handle that evaluates the objective function at a given input x.
  * grad      - Function handle that evaluates the gradient of the objective function at a given input x.
* feasible_set - Structure that describes the feasible set through these fields:
  * max_linear        - Function that maximizes a linear function with the first argument as coefficients over the feasible set and returns the optimal solution. 
  * random_boundary   - Function that returns a random point on the boundary of the feasible set.
  * n                 - Dimension of the feasible set.
  
## Replicating

To replicate the results in any of the tables in the paper, simply run the "RunInstances"  
function on the relevant sub directory in the instances directory.

## Support

For support in using this software, submit an
[issue](https://github.com/tkralphs/JoCTemplate/issues/new).
