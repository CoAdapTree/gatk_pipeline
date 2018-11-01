# GATK individual and poolseq pipeline
-----
## Assumed requirements
1. install an anaconda (not miniconda) environment with python 3.7 (eg: `conda create -n py3 python=3.7)`)
    1. export anaconda path in `$HOME/.bashrc` (this should automatically be done when installing anaconda)
    1. source anaconda env within `$HOME/.bashrc` on the last line of the file (`source activate py3`)
1. install bwa and export path in `$HOME/.bashrc`
1. install samtools with anaconda (after activating environment: `conda install -c bioconda samtools`)
1. install picardtools with anaconda (after activating env: `conda install -c bioconda picard`)
1. install gatk or make sure there is a module (eg `module load gatk/4.0.8.1`)
1. copy the following into your `$HOME/.bashrc` file so that the `def-someuser` reflects your compute canada account
    ```
    export SLURM_ACCOUNT=def-someuser  
    export SBATCH_ACCOUNT=$SLURM_ACCOUNT  
    export SALLOC_ACCOUNT=$SLURM_ACCOUNT
    ```
1. clone the pipeline repo to the server and create a symlink in `$HOME` so that it can be accessed via `$HOME/pipeline`

-----

## Using the pipeline
- To kick off the pipeline, first source your bashrc (`source ~/.bashrc`), `cd ~/pipeline`, and run `00_start-pipeline.py` from the home node, and it will run the rest of the preprocessing pipeline automatically by serially sbatching jobs (through `02_scatter-gvcf.py`). See example datatable.txt file needed for `00_start-pipeline.py`. Within `scheduler.py`, the user can specify the total number of jobs they wish to have in their queue and the pipeline with automatically schedule new jobs as others finish (as of this commit, the scheduler only controls the queue for `HaplotypeCaller` in `02_scatter-gvcf.py`).

`(py3) [user@host pipeline]$ python 00_start-pipeline.py /path/to/folder/with/fastq.gzfiles`
