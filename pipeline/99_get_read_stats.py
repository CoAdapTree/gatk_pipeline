### usage
# put this in an SBATCH file and submit to queue with 32 engines
# python 99_get_read_stats.py /path/to/parentdir/used/in/00_start-pipeline/command engines<int>
### 
##!/bin/bash
##SBATCH --time=12:0:0
##SBATCH --nodes=1
##SBATCH --ntasks=32
##SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=500M
##SBATCH --job-name=99_get_read_stats
##SBATCH --output=99_get_read_stats---%j.out
##source $HOME/.bashrc
##cd $HOME/pipeline
##module load samtools/1.9
##cd $HOME/pipeline
##python 99_get_read_stats.py $1 32

### imports
import sys
import os
import pandas as pd
from os import path as op
from os import listdir
import pickle
import numpy as np
import json
from collections import OrderedDict
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
thisfile, parentdir, engines = sys.argv
if parentdir.endswith("/"):
    parentdir = parentdir[:-1]
###

### reqs
os.system('echo getting reqs')
samp2pool = pickle.load(open(op.join(parentdir,'samp2pool.pkl'),'rb'))
pools     = uni(list(samp2pool.values()))
###

# get a list of subdirectory pool dirs created earlier in pipeline
os.system('echo getting pooldirs')
pooldirs = []
for p in pools:
    pooldir = op.join(parentdir,p)
    assert op.exists(pooldir)
    pooldirs.append(pooldir)

######## TRIMMING DATA ########
# get the json data from trimming
os.system('echo getting trim data')
data = {}
count = 0
for p in pooldirs:
    trimdir = op.join(p,'trimmed')
    jsons = [f for f in fs(trimdir) if f.endswith('.json')]
    count += len(jsons)
    for j in jsons:
        with open(j,'r') as f:
            data[op.basename(j)] = json.load(f)

# put data into a dataframe, and sort columns
os.system('echo reading data')
readinfo = OrderedDict()
samps    = []
for j in sorted(data):
    samp = j.replace(".json","").split(".")[-1].split("__")[0]
    for when in ['before_filtering','after_filtering']:
        for which in ['total_reads','total_bases','q20_bases','q30_bases']:
            what = "%s-%s" % (which,when.replace("filtering","trimming"))
            if what not in readinfo:
                readinfo[what] = OrderedDict()
            readinfo[what][samp] = data[j]['summary'][when][which]
    if not 'trim_command' in readinfo:
        readinfo['trim_command'] = OrderedDict()
    readinfo['trim_command'][samp] = data[j]['command']
    samps.append(samp)

# get counts from downstream
os.system('echo getting bam counts')
key = ['mapped_bamfile','filtered_bamfile','dedup_bamfile']
for k in key:
    readinfo[k] = OrderedDict()
for p in pooldirs:
    os.system('echo %s' % p)
    for i,d in enumerate(['sorted_bamfiles','filtered_indexed_sorted_bamfiles','dedup_rg_filtered_indexed_sorted_bamfiles']):
        os.system('echo %s' % d)
        DIR = op.join(p,d)
        assert op.exists(DIR)
        bams = [f for f in fs(DIR) if f.endswith('.bam')]
        for b in bams:
            if not d == 'dedup_rg_filtered_indexed_sorted_bamfiles':
                samp = op.basename(b).split("_R1R2")[0].split(".")[-1]
            else:
                samp = op.basename(b).replace("_rd.bam","")
            assert samp in samps
#             readinfo[key[i]][samp] = 'somenum' # for testing
            num = os.popen('samtools view -@ %(engines)s -c %(b)s' % locals()).read().replace("\n","")
            readinfo[key[i]][samp] = num
            
# make the dataframe
os.system('echo creating dataframe')
df = pd.DataFrame(readinfo)
df['samp'] = list(df.index)
df.index = range(len(df.index))
order = ['samp',
         'total_reads-before_trimming', 
         'total_reads-after_trimming', 
         'total_bases-before_trimming', 
         'total_bases-after_trimming', 
         'q30_bases-before_trimming', 
         'q30_bases-after_trimming', 
         'q20_bases-before_trimming', 
         'q20_bases-after_trimming', 
         'mapped_bamfile', 
         'filtered_bamfile', 
         'dedup_bamfile',
         'trim_command']
df = df[[x for x in order]].copy()
file = op.join(parentdir,'readinfo.txt')
df.to_csv(file,sep='\t',index=False)
print('created file:',file)








