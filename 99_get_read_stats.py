"""Convert json file from fastp into data table.
Count reads in bamfiles from throughout pipeline.

### usage
# module load samtools/1.9
# python 99_get_read_stats.py parentdir 32
###

### purpose
# get counts for trimmed, bams
###
"""

# imports
import os, sys, json, pandas as pd
from tqdm import tqdm
from os import path as op
from collections import OrderedDict
from coadaptree import fs, uni, pklload


# args
thisfile, parentdir, engines = sys.argv
if parentdir.endswith("/"):
    parentdir = parentdir[:-1]


# reqs
print('getting reqs')
samp2pool = pklload(op.join(parentdir, 'samp2pool.pkl'))
pools = uni(list(samp2pool.values()))


# get a list of subdirectory pool dirs created earlier in pipeline
print('getting pooldirs')
pooldirs = []
for p in pools:
    pooldir = op.join(parentdir, p)
    pooldirs.append(pooldir)


# TRIMMING DATA
# get the json data from trimming
print('getting trim data')
data = {}
count = 0
for p in pooldirs:
    trimdir = op.join(p, '01_trimmed')
    jsons = [f for f in fs(trimdir) if f.endswith('.json')]
    count += len(jsons)
    for j in jsons:
        with open(j,'r') as f:
            data[op.basename(j)] = json.load(f)


# put data into a dataframe, and sort columns
print('reading data')
readinfo = OrderedDict()
samps = []
for j in sorted(data):
    jsplit = j.split("__trimmed")[0]
    splits = jsplit.split(".")
    samp = '.'.join([splits[-1]] + splits[:-1])
    for when in ['before_filtering', 'after_filtering']:
        for which in ['total_reads', 'total_bases', 'q20_bases', 'q30_bases']:
            what = "%s-%s" % (which, when.replace("filtering", "trimming"))
            if what not in readinfo:
                readinfo[what] = OrderedDict()
            readinfo[what][samp] = data[j]['summary'][when][which]
    if not 'trim_command' in readinfo:
        readinfo['trim_command'] = OrderedDict()
    readinfo['trim_command'][samp] = data[j]['command']
    samps.append(samp)


# BAM DATA
# get counts from downstream
print('getting bam counts')
key = ['mapped_bamfile', 'dedup_bamfile', 'realigned_bamfile']
for k in key:
    readinfo[k] = OrderedDict()
for p in pooldirs:
    print(p)
    for i,d in enumerate(['02c_sorted_bamfiles',
                          '03_dedup_rg_filtered_indexed_sorted_bamfiles']):
        print(d)
        DIR = op.join(p, d)
        bams = [f for f in fs(DIR) if f.endswith('.bam')]
        for b in tqdm(bams):
            if d == '02c_sorted_bamfiles':
                splits = op.basename(b).split("_R1R2")[0].split(".")
                samp = '.'.join([splits[-1]] + splits[:-1])
            else:
                samp = op.basename(b).replace("_rd.bam", "")
            num = os.popen('samtools view -@ %(engines)s -c %(b)s' % locals()).read().replace("\n", "")
            readinfo[key[i]][samp] = num


# make the dataframe
print('creating dataframe')
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
         'dedup_bamfile',
         'trim_command']
df = df[[x for x in order]].copy()
file = op.join(parentdir, 'readinfo.txt')
df.to_csv(file, sep='\t', index=False)
print('created file:', file)
