"""
### execution
# python 01a_trim-fastq.py /path/to/pooldir /path/to/ref.fa
###
"""


import os
import sys
import time
import shutil
import subprocess
from os import path as op
from coadaptree import fs, pklload, pkldump, get_email_info

# args
thisfile, pooldir, ref = sys.argv
parentdir = op.dirname(pooldir)
pool = op.basename(pooldir)
f2samp = pklload(op.join(parentdir, 'f2samp.pkl'))
adaptors = pklload(op.join(parentdir, 'adaptors.pkl'))
for arg, path in [('pooldir', pooldir), ('ref', ref)]:
    if not op.exists(path):
        print("The argument does not exist in the specified path:\narg = %s\npath =%s" % (arg, path))
        sys.exit(1)


# make some dirs
shdir = op.join(pooldir, 'shfiles')
shtrimDIR = op.join(shdir, '01_trimmed_shfiles')  # cmd.sh files
trimDIR = op.join(pooldir, '01_trimmed')          # outfiles
for d in [shtrimDIR, trimDIR]:
    if not op.exists(d):
        os.makedirs(d)
mfile = op.join(parentdir, 'msgs.txt')
###


def writetomfile(text):
    with open(mfile, 'a') as m:
        m.write("%s\n" % text)


# get the fastq.gz files
os.chdir(pooldir)
gzfiles = [f for f in fs(pooldir) if 'R1' in f]
lgz = len(gzfiles)
text = 'found %(lgz)s R1 fastq.gz files in %(pooldir)s' % locals()
print('\t%s' % text)
writetomfile(text)

# match seq pairs to samp, alert if pair not found
seq_pairs = {}
for f in gzfiles:
    samp = f2samp[f]
    if samp not in seq_pairs:
        seq_pairs[samp] = []
    read2 = f.replace("_R1", "_R2")
    if op.exists(read2):
        seq_pairs[samp].append((op.abspath(f), op.abspath(read2)))
    else:
        text = '\nWARNING: no pair for %s\n' % f
        writetomfile(text)
text = "found %s R1/R2 seq pairs\n" % str(len([f for samp, files in seq_pairs.items() for f in files]))
print('\t%s' % text)
writetomfile(text)

# write sh files
shfiles = []
email_text = get_email_info(parentdir, '01')
samp2_r1r2out = {}
for samp, pairs in seq_pairs.items():
    samp2_r1r2out[samp] = []
    header = '''#!/bin/bash
#SBATCH --job-name=%(pool)s-%(samp)s-trim
#SBATCH --time=02:59:00
#SBATCH --mem=5000M
#SBATCH --cpus-per-task=16
#SBATCH --output=%(pool)s-%(samp)s-trim_%%j.out
%(email_text)s

source $HOME/.bashrc
export PYTHONPATH="${PYTHONPATH}:$HOME/pipeline"
export SQUEUE_FORMAT="%%.8i %%.8u %%.12a %%.68j %%.3t %%16S %%.10L %%.5D %%.4C %%.6b %%.7m %%N (%%r)"

module load fastp/0.19.5

''' % locals()
    
    newtext = ''''''
    for r1, r2 in pairs:
        r1adaptor, r2adaptor = list(adaptors[samp].values())
        r1out = op.join(trimDIR, op.basename(r1).split(".fastq")[0] + "_trimmed.fastq.gz")
        r2out = op.join(trimDIR, op.basename(r2).split(".fastq")[0] + "_trimmed.fastq.gz")
        html = r1out.replace("R1", "").replace(".fastq.gz", "_R1_R2_stats")
        json = r1out.replace("R1", "").replace(".fastq.gz", "_R1_R2")
        logfile = r1out.replace("R1", "").replace(".fastq.gz", "_R1_R2_stats.log")
        samp2_r1r2out[samp].append((r1out, r2out))

        text = '''fastp -i %(r1)s -o %(r1out)s -I %(r2)s -O %(r2out)s --disable_quality_filtering \
-g --cut_window_size 5 --cut_mean_quality 30 --n_base_limit 20 --length_required 75 \
-h %(html)s.html --cut_by_quality3 --thread 16 --json %(json)s.json \
--adapter_sequence %(r1adaptor)s --adapter_sequence_r2 %(r2adaptor)s > %(logfile)s

''' % locals()
        newtext = newtext + text
        
    suffix = '''# once finished, map using bwa mem 
python $HOME/pipeline/02_bwa-map_view_sort_index_flagstat.py %(parentdir)s %(samp)s

''' % locals()
    
    text = header + newtext + suffix
    
    filE = op.join(shtrimDIR, '%(pool)s-%(samp)s-trim.sh' % locals())
    shfiles.append(filE)
    with open(filE, 'w') as o:
        o.write("%s" % text)
pkldump(samp2_r1r2out, op.join(pooldir, 'samp2_r1r2out.pkl'))


print('\tshcount =', len(shfiles))
print('\tshdir = ', shtrimDIR)
# qsub the files
for sh in shfiles:
    os.chdir(op.dirname(sh))     # want sbatch outfiles in same folder as sh file
    print('\tshfile=', sh)
    subprocess.call([shutil.which('sbatch'), sh])
    # os.system('sbatch %s' % sh)
    time.sleep(2)
