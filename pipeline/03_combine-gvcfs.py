### usage
# python 03_combine-gvcfs.py /path/to/pool/fastqfolder/ /path/to/ref.fa sampname
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

### reqs
samp2pool = pickle.load(open(op.join(parentdir,'samp2pool.pkl'),'rb'))
rginfo    = pickle.load(open(op.join(parentdir,'rginfo.pkl'),'rb'))
poolsamps = pickle.load(open(op.join(parentdir,'poolsamps.pkl'),'rb')) #key=pool val=list(sampnames)
poolref   = pickle.load(open(op.join(parentdir,'poolref.pkl'),'rb'))   #key=pool val=/path/to/ref.fa
pools = uni(list(samp2pool.values()))
samps = list(rginfo.keys())
###

pooldirs = []
for p in pools:
    pooldir = op.join(parentdir,p)
    assert op.exists(pooldir)
    pooldirs.append(pooldir)

shfiles = []
for pdir in pooldirs:
    pool    = op.basename(pdir)
    ref     = poolref[pool]
    gvcfdir = op.join(pdir,'vcfs')
    assert op.exists(gvcfdir)
    comdir  = op.join(pdir,'combined_gvcfs')
    shdir   = op.join(pdir,'shfiles/combine_gvcfs')
    for x in [comdir,shdir]:
        if not op.exists(x):
            os.makedirs(x)
        
    # get list of gvcf files for each samp
    for samp in poolsamps[pool]:
        sampfiles  = [f for f in fs(gvcfdir) if samp in f and f.endswith('.gz')]
        outfile    = op.join(comdir,'%s_combined.g.vcf.gz' % samp)
        variantcmd = '--variant ' + ' --variant '.join([x for x in sorted(sampfiles)])
        filE = op.join(shdir,'%s.sh' % samp)
        text       = '''#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --nodes=1
#SBATCH --mem=30000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=combine_%(samp)s
#SBATCH --export=all
#SBATCH --output=combine-gvcf_%(samp)s_%%j.out 

# for debugging
echo $0
echo %(filE)s

source $HOME/.bashrc
module load gatk/4.0.8.1

gatk CombineGVCFs -R %(ref)s %(variantcmd)s -O %(outfile)s

''' % locals()
        with open(filE,'w') as o:
            o.write("%s" % text)
        shfiles.append(filE)
        print(filE)
        

for s in shfiles:
    os.chdir(op.dirname(s))
    os.system('sbatch %s' % s)

        



