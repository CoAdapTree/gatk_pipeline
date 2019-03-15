### usage
# python 00_start-pipeline.py /path/to/folder/with/all/fastq/files/ 
### 

### imports
from __future__ import division
import sys
import os
import pandas as pd
import numpy as np
from os import path as op
from os import listdir as ls
import pickle
import math
import time

def fs(DIR):
    return sorted([op.join(DIR,f) for f in ls(DIR)])
def sbatch(DIR):
    os.system("cd %s" % DIR)
    files = [f for f in fs(DIR) if '.sh' in f]
    for f in files:
        os.system('sbatch %s' % f)
        time.sleep(2)
def pkldump(obj,f):
    with open(f,'wb') as o:
        pickle.dump(obj,o,protocol=pickle.HIGHEST_PROTOCOL)
def uni(mylist):
    return (np.unique(mylist).tolist())
def luni(mylist):
    return (len(uni(mylist)))
os.system('source $HOME/.bashrc')
###

### args
thisfile, parentdir = sys.argv
###


# read in the datatable, save rginfo for later
if parentdir.endswith("/"):
    parentdir = parentdir[:-1]
print ('reading datatable, getting rginfo')
datatable = op.join(parentdir,'datatable.txt')
try:
    assert op.exists(datatable)
except AssertionError:
    print ('the datafile is not in the necessary path: %s' % datatable)
    sys.exit(3)
data      = pd.read_csv(datatable,sep='\t')
rginfo    = {} #key=sampname vals=rginfo
f2pool    = {} #key=filename val=pool
samp2pool = {} #key=samp val=pool
poolref   = {} #key=pool val=ref.fa
ploidy    = {} #key=pool val=ploidy
poolsamps = {} #key=pool val=sampnames
for row in data.index:
    if type(data.loc[row,'pool_name']) == float: # if pool_name is blank/nan, it is a poolseq
        data.loc[row,'pool_name'] = data.loc[row,'sample_name']
    samp = data.loc[row,'sample_name']
    pool = data.loc[row,'pool_name']
    if not pool in poolsamps:
        poolsamps[pool] = []
    if not samp in poolsamps[pool]:
        poolsamps[pool].append(samp)
    samp2pool[samp] = pool
    df = data[data['pool_name'] == pool].copy()
    try:
        assert luni(df['ploidy']) == 1
    except AssertionError as e:
        print ("the ploidy values for some elements with pool name '%s' are not the same" % pool)
        sys.exit(1)
    if not pool in ploidy:
        ploidy[pool] = data.loc[row,'ploidy']
    if pool in poolref:
        try:
            assert poolref[pool] == data.loc[row,'ref']
        except AssertionError as e:
            print ("ref genome for samples in %s pool seems to have different paths in datatable.txt" % pool)
            print ("samples in a pool should have the same ref genome")
            sys.exit(1)
    else:
        poolref[pool] = data.loc[row,'ref']
    rginfo[samp] = {}
    for col in ['rgid','rglb','rgpl','rgpu','rgsm']: # rg info columns
        rginfo[samp][col] = data.loc[row,col]    
    for f in [data.loc[row,'file_name_r1'],data.loc[row,'file_name_r2']]:
        f2pool[f] = pool
pkldump(rginfo,op.join(parentdir,'rginfo.pkl'))
pkldump(ploidy,op.join(parentdir,'ploidy.pkl'))
pkldump(f2pool,op.join(parentdir,'samp2pool.pkl'))
pkldump(poolsamps,op.join(parentdir,'poolsamps.pkl'))
pkldump(poolref,op.join(parentdir,'poolref.pkl'))

# make pool dirs
print ("making pool dirs")
pools = uni(data['pool_name'].tolist())
pooldirs = []
for p in pools:
    DIR = op.join(parentdir,p)
    if not op.exists(DIR):
        os.makedirs(DIR)
    if op.exists(DIR):
        pooldirs.append(DIR)

# get list of files from datatable, make sure they exist in parentdir, create symlinks in /parentdir/<pool_name>/
print ('creating symlinks')
files = [f for f in fs(parentdir) if 'fastq' in f and 'md5' not in f]
datafiles = data['file_name_r1'].tolist()
[datafiles.append(x) for x in data['file_name_r2'].tolist()]
for f in datafiles:
    try:
        src = op.join(parentdir,f)
        assert op.exists(src)
        pooldir = op.join(parentdir,f2pool[f])
        dst = op.join(pooldir,f)
        print ("creating fastq symlinks")
        if not op.exists(dst):
            os.symlink(src,dst)
        pklsrc = op.join(parentdir,'ploidy.pkl') # no need to assert 
        pkldst = op.join(pooldir,'ploidy.pkl')
        print ("creating ploidy symlink")
        if not op.exists(pkldst):
            os.symlink(pklsrc,pkldst)
        rgsrc = op.join(parentdir,'rginfo.pkl')  # no need to assert 
        rgdst = op.join(pooldir,'rginfo.pkl')
        print ("creating rginfo symlink")
        if not op.exists(rgdst):
            os.symlink(rgsrc,rgdst)
    except AssertionError as e:
        print ('could not find %s in %s' % (f,parentdir))
        sys.exit(1)

# create sh files
print ('writing sh files')
shdir = op.join(parentdir,'shfiles')
if not op.exists(shdir):
    os.makedirs(shdir)
for fqdir in pooldirs:
    pool = op.basename(fqdir)
    ref = poolref[pool]
    text = '''#!/bin/bash
#SBATCH --job-name=%(pool)s_start
#SBATCH --time=02:59:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=%(pool)s_start_%%j.out
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

source $HOME/.bashrc
cd $HOME/gatk_pipeline

python 01a_trim-fastq.py %(fqdir)s %(ref)s
''' % locals()
    filE = op.join(shdir,"%s.sh" % pool)
    with open(filE,'w') as o:
        o.write("%s" % text)

# sbatch jobs
print ('sbatching sh files')
# os.chdir(shdir)
# sbatch(shdir)





