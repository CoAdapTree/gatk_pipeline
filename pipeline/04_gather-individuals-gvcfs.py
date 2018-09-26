### purpose
# to combine individual sample gvcfs that were previously combined across intervals
###

### usage
# python 04_gather-individual-gvcfs.py /path/to/parentdir/used/in/00_start-pipeline/command 
###

### imports
import sys
import os
import shutil
import pickle
from os import path as op
from os import listdir
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
print(pooldirs)

outdir = op.join(parentdir,'combined_gvcfs')
if not op.exists(outdir):
    os.makedirs(outdir)
    
shfiles = []
for p in pooldirs:
    pooldir = op.dirname(p)
    cohort = op.basename(p)
    ref = poolref[cohort]  # /path/to/ref.fa
    shdir = op.join(pooldir,'shfiles/combine_gvcfs_again')
    if not op.exists(shdir):
        os.makedirs(shdir)
    combdir = op.join(p,'combined_gvcfs')
    gzfiles = [f for f in fs(combdir) if f.endswith('vcf.gz')]
    outfile = op.join(outdir,'%s_all_combined.vcf.gz' % cohort)
    variantcmd = '--variant ' + ' --variant '.join([x for x in sorted(gzfiles)])
    filE = op.join(shdir,'combine_cohort_%s.sh' % cohort)
    if len(gzfiles) == 0:
        print('found no files in %s' % combdir)
        continue
    if len(gzfiles) == 1: # if it's a poolseq file, just move it since no need to recombine again
        shutil.move(gzfiles[0],outfile)
        tbifile = gzfiles[0].replace(".gz",".gz.tbi")
        shutil.move(tbifile,op.join(outdir,op.basename(tbifile)))
        print('moved %s to %s' % (gzfiles[0],outfile))
    else:    
        text       = '''#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --nodes=1
#SBATCH --mem=30000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=%(cohort)s-cohort_combine
#SBATCH --export=all
#SBATCH --output=cohort_combine_gvcf-%(cohort)s-%%j.out 

# for debugging
cat $0
echo %(filE)s

source $HOME/.bashrc
module load gatk/4.0.8.1

gatk CombineGVCFs -R %(ref)s %(variantcmd)s -O %(outfile)s

''' % locals()
        with open(filE,'w') as o:
            o.write("%s" % text)
        shfiles.append(filE)
        print(filE)

# # # sbatch shfiles
for s in shfiles:
    os.chdir(op.dirname(s))
    os.system('sbatch %s' % s)



