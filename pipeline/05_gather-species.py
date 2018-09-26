### purpose
# to combine species-level sample gvcfs that were previously combined across intervals and library
###

### usage
# python 05_gather-species-gvcfs.py /path/to/parentdir/used/in/00_start-pipeline/command 
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

poolref = pickle.load(open(op.join(parentdir,'poolref.pkl'),'rb'))

# make some dirs
combdir = op.join(parentdir,'combined_gvcfs') # created in 04.py
finaldir = op.join(combdir,'combined_final')
shdir    = op.join(parentdir,'shfiles/combine_gvcfs_final')
for x in [shdir,finaldir]:
    if not op.exists(x):
        os.makedirs(x)

# of the gvcfs combined by library, get a list by species' pool
spp = {}
files = [f for f in fs(combdir) if 'tbi' not in f and f.endswith('.gz')]
for f in files:
    splits = op.basename(f).split("_")
    sp = splits[0]
    if not sp in spp:
        spp[sp] = {}
    if splits[1].startswith('p'):
        kind = 'pooled'
    else:
        kind = 'unpooled'
    if not kind in spp[sp]:
        spp[sp][kind] = []
    spp[sp][kind].append(f)
        
# make sh files
shfiles = []
for sp in spp:
    print(sp)
    for kind in spp[sp]:
        writesh = False
        print('\t',kind)
        outbasename = "_".join([op.basename(f).split("_all_combined")[0] for f in spp[sp][kind]]) + '_combined.vcf.gz'
        outfile = op.join(finaldir,outbasename)
        if len(spp[sp][kind]) > 1:
            print('\t','\t','combining variant commands')
            key = op.basename(spp[sp][kind][0]).split("_all")[0]
            ref = poolref[key]
            variantcmd = '--variant ' + ' --variant '.join([x for x in sorted(spp[sp][kind])])
            writesh = True
        elif len(spp[sp][kind]) == 1:
            # no need to combine, just move with tbi file
            f = spp[sp][kind][0]
            shutil.move(f,outfile)
            print('\t','\t','moving',f,'to',outfile)
            ftbi = f.replace(".gz",".gz.tbi")
            otbi = outfile.replace(".gz",".gz.tbi")
            print('\t','moving',ftbi,'to',otbi)
            shutil.move(ftbi,otbi)
        else:
            print (sp,kind,'has no files in',combdir)
        
        if writesh == True:
            text = '''#!/bin/bash
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --mem=30000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=%(sp)s_%(kind)s-spp_combine
#SBATCH --export=all
#SBATCH --output=%(sp)s_%(kind)s-spp_combine-%%j.out 

source $HOME/.bashrc
module load gatk/4.0.8.1

gatk CombineGVCFs -R %(ref)s %(variantcmd)s -O %(outfile)s


''' % locals()
            filE = op.join(shdir,'%(sp)s_%(kind)s-spp_combine.sh' % locals())
            print('\t','\t',filE)
            with open(filE,'w') as o:
                o.write("%s" % text)
            shfiles.append(filE)

# for s in shfiles:
#     os.chdir(op.dirname(s))
#     os.system('sbatch %s' % s)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        