### usage
# python 03_combine-gvcfs.py /path/to/pool/fastqfolder/ /path/to/ref.fa sampname
### 


### imports
import sys
import os
from os import path as op
from os import listdir
import pickle
def ls(DIR):
    return sorted([f for f in listdir(DIR)])
def fs (DIR):
    return sorted([op.join(DIR,f) for f in ls(DIR)])
### 

### args
thisfile, fqdir, ref, samp = sys.argv
###



gvcfdir = op.join(fqdir,'vcfs')
comdir  = op.join(fqdir,'combined_gvcfs')

if not op.exists(comdir):
    os.makedirs(comdir)
    
files = [f for f in fs(gvcfdir) if samp in f]
print 'len(files) =',len(files)

text = '''#!/bin/bash
#SBATCH --time=02:59:00
#SBATCH --nodes=1
#SBATCH --mem=8000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=%s%s%s
#SBATCH --export=all
#SBATCH --output=gvcf%s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu

source $HOME/.bashrc
module load gatk/4.0.0.0

gatk CombineGVCFs 
'''



