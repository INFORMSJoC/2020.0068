[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Projective Cutting-Planes for robust linear programming and Cutting-Stock problems

## Cite

The final version of this repository, with updated bibliographical information, is available at [GitHub](https://github.com/INFORMSJoC/2020.0068).

## Description

The goal of this software is to demonstrate the efficiency of the proposed method (Projective Cutting Planes) on two problems

1. Robust Linear Programming
2. Cutting Stock with multiple lengths
    
The code for the first problem is in the `robust-lp` folder.
The code for the second one is in the `cut-stock` folder.

A part of the code is actually shared. The shared source code files can be found
in the ``src_shared'' folder.

## Building

For both problems, you can build the executable by typing:
```
make 
```

Be sure to make clean before building a different version of the code.
```
make clean
```

## Results and replicating

For both programs, it is enough to type `./main` to see the command line
options. Two examples are provided below, first for robust optimization and then for
cutting stock:

1. The command `./main prj instances/maros.txt -m -tabularLatexOnly` will execute Projective Cutting Planes (because of argument `prj`) and output
only the tabular data that can be easily put in a latex table.
`` Ratio (robustobj-nominalObj)/nominalObj:  12.11 LOWGAP_ITERS       52 LOWGAP_TIME    0.2167 ITERS   53 TIME   0.8324 MULTICUTS 10015``
This is how the results from Table 2 (with gamma=10) have been generated.

2. The command `./main instances/m1M100n100.1bp` will simply execute Projective Cutting Planese on the very
first instance from the `wascher.txt' benchmark set. The very last printed line 

## Ongoing Development

This code is being developed on an on-going basis at a private 
github. Requests for copies of the latest code source may be addressed to
daniel.porumbel@cnam.fr.

