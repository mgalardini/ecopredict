ecopredict
=========

Scripts and pipeline to predict conditional growth phenotypes in
the *Escherichia coli* strain reference panel (EcoRef).
Also, the repository can generate all the figures of the related manuscript.

Note
----

The pipeline and scripts come with limited documentation.
Please do get in touch with the author (Marco Galardini, marco@ebi.ac.uk) if you need any guidance.

Usage
-----

The input data for this repository are strongly dependent on the
[pangenome_variation](https://github.com/mgalardini/pangenome_variation) and
[screenings](https://github.com/mgalardini/screenings) pipelines, as well as the
chemical genomics data.
Since the necessary input data is too big in size to be hosted in a git
repository, [it can be downloaded from here](http://www.ebi.ac.uk/~marco/inputs.7zip) and uncompressed in the base directory.
Input data contains both the *E. coli* genotypes and phenotypes, as well
as the precomputed predicted impact of non-synonymous substitutions and related
information.

The makefile contains the various bits of the pipeline:
* `make constraints` will analyse evolutionary constraints in functional and structural regions of the preoteome
* `make sickness` will compute the gene level probabilities of functional loss
* `make score` will compute the conditional growth predictions based on the chemical genomics data
* `make auc` will assess the predictions computed in the previous step
* `make bootstraps` will run three shuffling strategies to verify the robustness of the predictions (might need several Gb of disk, so beware)
* `make collect` will collect the results of the bootstraps
* `make associations` will run gene presence/absence associations with the phenotypes
* `make followup` will simulate the effects of gene complementations and compare them with actual experimental data
* `make plots` will generate individual plots
* `make figures` will combine all the plots and generate the figures as they appear in the manuscript (minus a few manual edits) 

(minimum) prerequisites
-------------

* Scoary (the version used by this repository is present as a submodule)
* bedtools
* python (2.7+ AND 3.3+), plus the following libraries:
    * jupyter
    * notebook
    * nbconvert
    * biopython
    * ete3
    * fastcluster
    * netwrokx
    * requests
    * numpy
    * pandas
    * matplotlib
    * seaborn
    * svgutils
* a working internet connection might be needed to download some data from uniprot and other sources

Copyright
---------

Copyright (C) <2015> EMBL-European Bioinformatics Institute

This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

Neither the institution name nor the name ecopredict
can be used to endorse or promote products derived from
this software without prior written permission.
For written permission, please contact <marco@ebi.ac.uk>.

Products derived from this software may not be called ecopredict
nor may ecopredict appear in their names without prior written
permission of the developers. You should have received a copy
of the GNU General Public License along with this program.
If not, see <http://www.gnu.org/licenses/>.
