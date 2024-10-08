[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# An Algorithm for maximizing a convex function based on its minimum

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[*An Algorithm for maximizing a convex function based on its minimum*](https://doi.org/10.1287/ijoc.2022.1238) by Aharon Ben-Tal and Ernst Roos

## Cite

To cite this material, please cite this repository, using the following DOI: [10.1287/ijoc.2022.1238.cd](https://doi.org/10.1287/ijoc.2022.1238.cd)

Below is the BibTex for citing this version of the code.

```
@article{comax2022,
  author =        {Ben-Tal, Aharon and Roos, Ernst},
  publisher =     {INFORMS Journal on Computing},
  title =         {An Algorithm for maximizing a convex function based on its minimum},
  year =          {2022},
  doi =           {10.1287/ijoc.2022.1238.cd},
  url =           {https://github.com/INFORMSJoC/2022.0034},
}  
```

## Description

The goal of this software is to maximize convex functions over a convex feasible set.

The software is written in MATLAB and uses the publically available YALMIP package 
(https://yalmip.github.io/), and the MOSEK solver (https://www.mosek.com/). For the latter,
free academic licenses are available.

## Usage

The software can be ran via the `RunInstances` function that takes a directory name as
input. It will solve all instances in the specified folder. Instances should be provided
as ".mat" files with two variables with the following specification:
* objective - Structure that describes the objective through these fields:
  * f         - Function handle that evaluates the objective function at a given input x.
  * grad      - Function handle that evaluates the gradient of the objective function at a given input x.
* feasible_set - Structure that describes the feasible set through these fields:
  * max_linear        - Function that maximizes a linear function with the first argument as coefficients over the feasible set and returns the optimal solution. 
  * random_boundary   - Function that returns a random point on the boundary of the feasible set.
  * n                 - Dimension of the feasible set.
  
For example, running 
`RunInstances("..\Data\Table 1\")`
will run all instances used to generate the results in Table 4, which will output the following:
```
Processing Instance 1
------------------------------------------------
Best Value found: 394.751 
Total time spent: 2.241 
------------------------------------------------
Processing Instance 2
------------------------------------------------
Best Value found: 884.751 
Total time spent: 2.162 
------------------------------------------------
Processing Instance 3
------------------------------------------------
Best Value found: 4674.677 
Total time spent: 1.782 
------------------------------------------------
Processing Instance 4
------------------------------------------------
Best Value found: 175705.589 
Total time spent: 2.026 
------------------------------------------------
Processing Instance 5
------------------------------------------------
Best Value found: 692613.047 
Total time spent: 3.802 
------------------------------------------------
Processing Instance 6
------------------------------------------------
Best Value found: 6020787.416 
Total time spent: 5.186 
------------------------------------------------
Processing Instance 7
------------------------------------------------
Best Value found: 1855739.984 
Total time spent: 6.196 
------------------------------------------------
```
as well as the content of Table 4.
  
## Replicating

To replicate the results in any of the tables in the paper, simply run the "RunInstances"  
function on the relevant sub directory in the instances directory.
