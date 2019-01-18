







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
1. clone the pipeline repo to the server and create a symlink in `$HOME` so that it can be accessed via `$HOME/pipeline`

-----

## Using the pipeline
- First create interval files (see `02_scatter-gvcf.py`) and install python modules from pipeline imports with `pip install <module>`. 
- See example datatable.txt file needed for `00_start-pipeline.py`. 
- Within `scheduler.py`, the user can specify the total number of `HaplotypeCaller` jobs they wish to have in their queue and the pipeline with automatically schedule new jobs as others finish (scheduler.py only controls the queue for `HaplotypeCaller` in `02_scatter-gvcf.py` - this queue threshold is independent of the threshold for genotyping).
- Within `genotyping_scheduler.py`, the user can specify the total number of `GenotypeGVCF` jobs they wish to have in their queue and the pipeline will automatically schedule new jobs as others finish (`genotyping_scheduler.py` only controls queue for `GenotypeGVCFs` in `03_combine...py` - this queue threshold is independent from the threshold for `HaplotypeCaller`) 

- To kick off the pipeline, source your bashrc (`source ~/.bashrc`) to activate the python env, `cd ~/pipeline`, and run `00_start-pipeline.py` from the home node, and it will run the rest of the preprocessing pipeline automatically by serially sbatching jobs (through `02_scatter-gvcf.py`).

`(py3) [user@host pipeline]$ python 00_start-pipeline.py /path/to/folder/with/fastq.gzfiles`
