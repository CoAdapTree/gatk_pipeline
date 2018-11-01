### usage
# python 03b_combine_and_genotype_within_pool.py /path/to/parentdir/used/in/00_start-pipeline/command 
### 

### purpose
# instead of calling snps when combining some of our pools, call snps for indiviual pools
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
def createdirs(dirs):
    for d in dirs:
        if not op.exists(d):
            os.makedirs(d)
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
ploidy    = pickle.load(open(op.join(parentdir,'ploidy.pkl'),'rb'))    #not used yet, but maybe if pooled needs more time
pools = uni(list(samp2pool.values()))
samps = list(rginfo.keys())
###

# get a list of subdirectory pool dirs created earlier in pipeline
pooldirs = []
for p in pools:
    pooldir = op.join(parentdir,p)
    assert op.exists(pooldir)
    pooldirs.append(pooldir)

# create some dirs
outdir = op.join(parentdir,'snps')
shdir  = op.join(parentdir,'shfiles/select_variants_single_pools')
createdirs([outdir,shdir])
        
# get a list of files and partition by those that belong to true pools
poolfiles = {}
for d in pooldirs:
    check = op.basename(d).split("_")[1]
    if check.startswith("p"):
        poolfiles[op.basename(d)] = []
        vcfdir = op.join(d,'vcfs')
        [poolfiles[op.basename(d)].append(f) for f in fs(vcfdir) if f.endswith('.gz')]

# make sh files
shfiles = []
for pool in poolfiles:
    for f in sorted(poolfiles[pool]):
        #no need to combine files by scaff, just create symlinks
        scaff    = op.basename(f).split("_scaff")[-1].replace(".g.vcf.gz","")
        combfile = op.join(outdir,'%s--%s_combined.vcf.gz' % (pool,scaff))
        tbifile  = f.replace(".gz",".gz.tbi")
        tbilink  = combfile.replace(".gz",".gz.tbi")
        try:
            os.symlink(f,combfile)
            os.symlink(tbifile,tbilink)
        except:
            print("symlink already exists")
            continue # this will skip re-doing any poolseq samps that were done solo from 03a_.py
        gfile    = combfile.replace("_combined.vcf.gz","_genotyped.vcf.gz")
        snpfile  = gfile.replace("_genotyped.vcf.gz","_snps.vcf.gz")
        cmd = '''gatk --java-options "-Xmx4g" GenotypeGVCFs -R %(ref)s -V %(combfile)s -O %(gfile)s

gatk SelectVariants -R %(ref)s -V %(gfile)s --select-type-to-include SNP -O %(snpfile)s

''' % locals()

        filE = op.join(shdir,'%(pool)s--%(scaff)s.sh' % locals())
        print('\t','\t',filE)
        scaff = op.basename(f).split("scaff")[1].split(".g.")[0]
        text = '''#!/bin/bash
#SBATCH --time=23:59:59
#SBATCH --nodes=1
#SBATCH --mem=30000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=%(pool)s--%(scaff)s
#SBATCH --export=all
#SBATCH --output=%(pool)s--%(scaff)s---%%j.out 

source $HOME/.bashrc
module load gatk/4.0.8.1

%(cmd)s

''' % locals()
        with open(filE,'w') as o:
            o.write("%s" % text)
        shfiles.append(filE)
print(shdir)
print(shfiles)