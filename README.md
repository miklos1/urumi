#### Experimentation framework for the manuscript:

> Miklós Homolya, Lawrence Mitchell, Fabio Luporini, and David
> A. Ham. "TSFC: a structure-preserving form compiler." (2017).

#### Usage

The installation required to run the measurements with the modified
version of FFC that supports non-affine geometries is incompatible
with the installation that is required to run the rest of the
measurements.  I suggest using separate `virtualenv`s.

For the former one needs the correct versions of COFFEE, FFC, FIAT,
Instant and UFL installed, then run measurements with:

    python run.py ffc-bendy >>data.txt

The latter needs COFFEE, dijitso, FFC, FIAT, FInAT, TSFC, and UFL, and
run measurements with:

    python run.py current >>data.txt

Once the data is collected, call

    python plot.py data.txt

to produce the violin plots.  This step requires `matplotlib`, and
creates a file called `example.pgf`.  This file shall be copied next
to the LaTeX source code of the paper before building the paper.

The `data.txt` file used for the submitted paper is also uploaded into
this repository.
