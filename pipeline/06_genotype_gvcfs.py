### usage
# python 06_genotype_gvcfs.py /path/to/parentdir/used/in/00_start-pipeline/command 
### 


### imports
import sys
import os
from os import path as op
from os import listdir
import pickle
import numpy as np
def uni(mylist):
    return (np.unique(mylist).tolist())
def ls(DIR):
    return sorted([f for f in listdir(DIR)])
def fs (DIR):
    return sorted([op.join(DIR,f) for f in ls(DIR)])
### 

### args
thisfile, parentdir = sys.argv
if parentdir.endswith("/"):
    parentdir = parentdir[:-1]
###

poolref = pickle.load(open('/home/lindb/scratch/alldata/poolref.pkl','rb'))
combdir = op.join(parentdir,'combined_gvcfs/combined_final')
shdir   = op.join(parentdir,'shfiles/genotyped_gvcfs')
if not op.exists(shdir):
    os.makedirs(shdir)
vcfs    = [v for v in fs(combdir) if v.endswith('gz') and 'genotypedgvcf.vcf.gz' not in v]

shfiles = []
for v in vcfs:
    poolkey = [x for x in poolref.keys() if x in v][0] # get the first key to determine ref
    label = op.basename(v).split("_combined")[0]
    ref = poolref[poolkey]
    out = v.replace(".vcf.gz","_genotypedgvcf.vcf.gz")
    text = '''#!/bin/bash
#SBATCH --time=7-0:0:0
#SBATCH --nodes=1
#SBATCH --mem=150000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=genotype-%(label)s
#SBATCH --export=all
#SBATCH --output=genotype_gvcf-%(label)s-%%j.out

source $HOME/.bashrc
module load gatk/4.0.8.1

gatk --java-options "-Xmx4g" GenotypeGVCFs -R %(ref)s -V %(v)s -O %(out)s

''' % locals()
    filE = op.join(shdir,'%s.sh' % label)
    with open(filE,'w') as o:
        o.write("%s" % text)
    shfiles.append(filE)
    print (filE)
    
for s in shfiles:
    os.chdir(op.dirname(s))
    os.system('sbatch %s' % s)
    
    
