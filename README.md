[![DOI](https://zenodo.org/badge/73704566.svg)](https://zenodo.org/badge/latestdoi/73704566)

#### Experimentation framework for the manuscript:

> Miklós Homolya, Lawrence Mitchell, Fabio Luporini, and David
> A. Ham. "TSFC: a structure-preserving form compiler." (2017).

#### Usage

The installation required to run the measurements with the modified
version of FFC that supports non-affine geometries is incompatible
with the installation that is required to run the rest of the
measurements.  I suggest using separate `virtualenv`s.

For the former one needs the correct versions of
COFFEE [[10.5281/zenodo.573235](https://doi.org/10.5281/zenodo.573235)],
FFC [[10.5281/zenodo.573237](https://doi.org/10.5281/zenodo.573237)],
FIAT [[10.5281/zenodo.573238](https://doi.org/10.5281/zenodo.573238)],
Instant [[10.5281/zenodo.573255](https://doi.org/10.5281/zenodo.573255)]
and UFL [[10.5281/zenodo.573236](https://doi.org/10.5281/zenodo.573236)]
installed, then run measurements with:

    python run.py ffc-bendy >>data.txt

The latter needs recent versions of
COFFEE [[10.5281/zenodo.573267](https://doi.org/10.5281/zenodo.573267)],
dijitso [[10.5281/zenodo.573287](https://doi.org/10.5281/zenodo.573287)],
FFC [[10.5281/zenodo.573270](https://doi.org/10.5281/zenodo.573270)],
FIAT [[10.5281/zenodo.573269](https://doi.org/10.5281/zenodo.573269)],
FInAT [[10.5281/zenodo.573266](https://doi.org/10.5281/zenodo.573266)],
TSFC [[10.5281/zenodo.573271](https://doi.org/10.5281/zenodo.573271)],
and UFL [[10.5281/zenodo.573268](https://doi.org/10.5281/zenodo.573268)],
and run measurements with:

    python run.py current >>data.txt

Once the data is collected, call

    python plot.py data.txt

to produce the violin plots.  This step requires `matplotlib`, and
creates a file called `example.pgf`.  This file shall be copied next
to the LaTeX source code of the paper before building the paper.

The `data.txt` file used for the submitted paper is also uploaded into
this repository.
