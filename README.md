[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)]([https://pubsonline.informs.org/journal/ijoc])

# [Projective Cutting-Planes for robust linear programming and Cutting-Stock problems](https://doi.org/10.1287/ijoc.2022.1160)

## Cite

To cite this material, please cite this repository, using the following DOI.

[![DOI](https://zenodo.org/badge/424937237.svg)](https://zenodo.org/badge/latestdoi/424937237)

Below is the BibTex for citing this version of the code.

```
@article{PCP4RLPaCSP2021,
  author =        {D. Porumbel},
  publisher =     {INFORMS Journal on Computing},
  title =         {Projective Cutting-Planes for robust linear programming and Cutting-Stock Problems},
  year =          {2021},
  doi =           {10.5281/zenodo.5745335},
  url =           {[https://github.com/INFORMSJoC/2020.0068](https://pubsonline.informs.org/doi/10.1287/ijoc.2022.1160)}
}  
```

## Description and folder structure

The goal of this software is to demonstrate the efficiency of the proposed method (Projective Cutting Planes) on two problems

1. Robust Linear Programming
2. Cutting Stock with multiple lengths
    
The source code for the first problem is provided in the `robust-lp` folder.
The source code for the second one is provided in the `cut-stock` folder.

In fact, a part of the code is actually shared by the two pieces of software. The shared source code files can be found in the `src_shared` folder.

Inside both the `robust-lp` and the `cut-stock` folder, there is `src` folder with the source code 
file and a folder `instances` folder containing the benchmark data set.

The file `CODE_GUIDELINES` from the `src_shared` folder describe the practices and standards used to write the whole software.
Such information may help one more easily understand the source code.

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
options. Two command line examples are provided below, first for robust 
optimization and then for cutting stock:

1. The command `./main prj instances/maros.txt -m -tabularLatexOnly` will execute Projective Cutting Planes (because of argument `prj`) and output
only the tabular data that can be easily modifed and compiled using `pdflatex` to obtain a pdf file.  This is how the results from Table 2 have been generated.

``` Ratio (robustobj-nominalObj)/nominalObj:  12.11 LOWGAP_ITERS       52 LOWGAP_TIME    0.2167 ITERS   53 TIME   0.8324 MULTICUTS 10015```


2. The command `./main instances/m1M100n100.1bp` will simply execute Projective Cutting Planese on the very
first instance from the `wascher.txt` benchmark set. The very last printed line provides the tabular data
that can be integrated into a latex table to generate a pdf document.

## Ongoing Development

This code is being developed on an on-going basis using a private github. Requests for copies of the latest code source may be addressed to daniel.porumbel@cnam.fr.

