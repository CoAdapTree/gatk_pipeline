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
