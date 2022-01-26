# Contents of this directory
This directory contains files from a run of the SSI Influenza phylogeny pipeline.
Read more in on the Phylogeny documentation page in InfluenzaDocumentation
Contact Jakob Nybo Nissen for details.

### genotypes.txt
This lists consensus sequences in three categories. Currently only consensus sequences with at least one segment passing quality control will be listed.

The three groups are:
    * Known: The genotype is unambiguously resolvable
    * Indeterminate: The sample could belong to two or more genotypes due to missing segments
    * New genotypes: The sample does not fit in any known genotype

### phylogeny_versions.txt
This file displays the version of the pipeline and the Julia language used to generate the files in this directory.

### trees
This directory contains phylogenetic trees (as PDF a Newick files) for a subset of the processed segments.
In the PDF file, samples are highlighted in red, the others are reference sequences. 

### tmp
This directory contains files that are considered internal to the pipeline.
That means the content is subject to change and mostly not documented.
Currently, the `tmp/blast` dir contains BLAST results to references, this can be useful to check why some segments were or were not assigned to a particular clade.
