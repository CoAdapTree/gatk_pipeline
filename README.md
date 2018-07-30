# coadaptree
=====
## Assumed requirements
1. install an anaconda (not miniconda) environment with python 2.7 (eg: conda create -n py27 python=2.7)
    1. export anaconda path in $HOME/.bashrc (this should automatically be done when installing anaconda)
    1. source anaconda env within $HOME/.bashrc (source activate py27)
1. install bwa and export path in $HOME/.bashrc
1. install samtools with anaconda (after activating environment: conda install -c bioconda samtools)
1. clone the pipeline repo to the server and create a symlink in $HOME so that it can be accessed via $HOME/pipeline

-----

## Using the pipeline
- To kick off the pipeline, just run 01a_trim-fastq.py from the home node, and it will run the rest of the preprocessing pipeline automatically by serially sbatching jobs (through 01b). I can implement a script within 01a or 01b to check to see when all of the final bam files are created at the end of 01b.py which could then kick off the 02.py script automatically, but I haven't done this yet (as of now we'd have to kick off the 02.py script once all bam files are ready to be merged).
