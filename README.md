






[![](https://img.shields.io/badge/CoAdapTree-UBC-green.svg?style=popout&logo=data:image/png;base64,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)](http://coadaptree.forestry.ubc.ca)


# GATK individual and poolseq pipeline
-----
## Assumed requirements
1. Access to an HPC with a scheduler (e.g., slurm, SGE, PBS, TORQUE) - this pipeline assumes slurm
1. install an anaconda (not miniconda) or virtual environment with python 3.7 (eg: `conda create -n py3 python=3.7` or `virtualenv --no-download ~/py3`)
    1. source env within `$HOME/.bashrc` on the last line of the file (`source /path/to/conda/bin/activate py3` or `source ~/py3/bin/activate`)
1. install bwa and export path in `$HOME/.bashrc` or load  module (eg `module load bwa/0.7.17`)
1. install samtools with anaconda (after activating environment: `conda install -c bioconda samtools`) or load module (eg `module load samtools/1.9`)
1. install picardtools with anaconda (after activating env: `conda install -c bioconda picard`) or load module (eg `module load picard/2.18.9`)
1. install gatk or load module (eg `module load gatk/4.0.8.1`)
1. install bcftools load module (eg `module load bcftools/1.9`)
1. copy the following into your `$HOME/.bashrc` file so that the `def-someuser` reflects your non-RAC compute canada account
    ```
    export SLURM_ACCOUNT=def-someuser  
    export SBATCH_ACCOUNT=$SLURM_ACCOUNT  
    export SALLOC_ACCOUNT=$SLURM_ACCOUNT
    ```
1. clone the pipeline repo to the server and create a symlink in `$HOME` so that it can be accessed via `$HOME/gatk_pipeline`

-----

## Using the pipeline
- First create interval files (see `02_scatter-gvcf.py`) and install python modules from pipeline imports with `pip install <module>`. 
- See example datatable.txt file needed for `00_start-pipeline.py`. 
- Within `scheduler.py`, the user can specify the total number of `HaplotypeCaller` jobs they wish to have in their queue and the pipeline with automatically schedule new jobs as others finish (scheduler.py only controls the queue for `HaplotypeCaller` in `02_scatter-gvcf.py` - this queue threshold is independent of the threshold for genotyping).
- Within `genotyping_scheduler.py`, the user can specify the total number of `GenotypeGVCF` jobs they wish to have in their queue and the pipeline will automatically schedule new jobs as others finish (`genotyping_scheduler.py` only controls queue for `GenotypeGVCFs` in `03_combine...py` - this queue threshold is independent from the threshold for `HaplotypeCaller`) 

- To kick off the pipeline, source your bashrc (`source ~/.bashrc`) to activate the python env, `cd ~/gatk_pipeline`, and run `00_start-pipeline.py` from the home node, and it will run the rest of the preprocessing pipeline automatically by serially sbatching jobs (through `02_scatter-gvcf.py`).

`(py3) [user@host gatk_pipeline]$ python 00_start-pipeline.py /path/to/folder/with/fastq.gzfiles`
