####
# assumes fastq.gz files
###

###
# execution:
# python 01a_trim-fastq.py /path/to/fastq.gz-files/ /path/to/ref.fa <int-number-of-shfiles>
###

###
# imports
from __future__ import division
import sys
import os
from os import path as op
from os import listdir as ls
import math
def fs(DIR):
    return sorted([op.join(DIR,f) for f in ls(DIR)])
os.system('source $HOME/.bashrc')
###


# args
fqdir   = sys.argv[1] # path to raw fastq files
ref     = sys.argv[2] # path to reference genome used for mapping (/home/lindb/scratch/ptaeda.v1.01.reduced.pseudo.fasta)
lensh   = int(sys.argv[3])
print "lensh =",lensh
for i,arg in enumerate([ref,fqdir]):
    # make sure the args exist
    try:
        assert op.exists(arg)
    except AssertionError as e:
        print "The %s'th argument does not exist in the specified path" % str(i)
        sys.exit(1)

# more imports and aliases (don't think I need this any more)
pipedir = op.dirname(op.abspath(sys.argv[0]))         # this is where the git repo (pipeline) is cloned


# make some dirs
msgdir    = op.join(fqdir,'messages')
shdir     = op.join(fqdir,'shfiles')
shtrimDIR = op.join(shdir,'trimmed_shfiles') # will make shdir below
trimDIR   = op.join(fqdir,'trimmed')         # outfiles
for d in [shtrimDIR,trimDIR,msgdir]:
    if not op.exists(d):
        os.makedirs(d)
mfile = op.join(fqdir,'messages/msgs.txt')
        
###



###
        
# get the fastq.gz files
os.chdir(fqdir)
gzfiles = [f for f in fs(fqdir) if 'R1' in f]
lgz     = len(gzfiles)
# !echo 'found '$lgz' gz files in '$fqdir >> $fqdir'/messages/msgs.txt' # only works in jupyter notebooks :'(
# instead of ^, do (lame/boring):
text = 'found %s R1 fastq.gz files in %s' % (lgz, fqdir)
print text
with open(mfile,'w') as o:
    o.write("%s\n" % text)

# match seq pairs, alert if pair not found
seq_pairs = []
for f in gzfiles:
    assert 'R1' in f
    read2 = f.replace("_R1","_R2")
    if op.exists(read2):
        seq_pairs.append((f,read2))
#     else:
#         !echo 'no pair for '$f >> $fqdir'/messages/msgs.txt' # #iheartjupyter
    else:
        text = 'no pair for %s' % f
        with open(mfile,'a') as o:
            o.write("%s\n" % text)
print "found %s R1/R2 seq pairs" % str(len(seq_pairs))
print "type(len(seq_pairs)) =",type(len(seq_pairs))
print "type(lensh) =",type(lensh)
print "len(seq_pairs) <= lensh?", len(seq_pairs) <= lensh
# determine how many commands per sh file
if len(seq_pairs) <= lensh:
    # one command per sh file
    ceil = 1
else:
    # multiple commands per sh file
    ceil = math.ceil(len(seq_pairs)/lensh)
print "ceil =",ceil
    
shcount = 0
fcount  = 0
tcount  = 0
text    = ''''''
for s in seq_pairs:
#     print s
    r1    = op.abspath(s[0])
    if r1.endswith("fastq"):
        r1out = op.join(trimDIR,op.basename(r1).replace('.fastq','_trimmed.fastq'))
    else:
        r1out = op.join(trimDIR,op.basename(r1).replace('.fastq.gz','_trimmed.fastq'))
    r2    = op.abspath(s[1])
    if r2.endswith("fastq"):
        r2out = op.join(trimDIR,op.basename(r2).replace('.fastq','_trimmed.fastq'))
    else:
        r2out = op.join(trimDIR,op.basename(r2).replace('.fastq.gz','_trimmed.fastq'))
    html = r1out.replace("R1","").replace(".fastq","_R1_R2_stats")
    json = r1out.replace("R1","").replace(".fastq","_R1_R2")
    logfile = r1out.replace("R1","").replace(".fastq","_R1_R2_stats.log")
    cmd  = '''fastp -i %s -o %s -I %s -O %s -g --cut_window_size 5 --cut_mean_quality 30 --qualified_quality_phred 30 --unqualified_percent_limit 20 --n_base_limit 5 --length_required 75 -h %s.html --cut_by_quality3 --thread 16 --json %s.json --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT > %s


cd $HOME/pipeline
# once finished, map using bwa mem 
python 01b_bwa-map.py %s %s %s %s %s           
''' % (r1  , r1out,
       r2  , r2out,
       html, json ,  logfile,
       ref , r1out,  r2out  , shdir, str(tcount).zfill(3)
      )
# (r1,  r1out,
#        r2,  r2out,
#        r1out.replace("R1_trimmed.fastq","R1_R2_trimmed_stats"),
#        r1out.replace("R1_trimmed.fastq","R1_R2_trimmed"),
#        ref,  r1out,  r2out,  shdir,  str(tcount).zfill(3)
#       )
    text = text + cmd
    
    fcount += 1
    tcount += 1
    if fcount == ceil or tcount == len(seq_pairs):
        shz = str(shcount).zfill(3)
        header = '''#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --job-name=trim%s
#SBATCH --export=all
#SBATCH --time=11:59:00
#SBATCH --mem=1000mb
#SBATCH --cpus-per-task=16
''' % (shz)
        text = header + text
        filE = op.join(shtrimDIR,'trim_%s.sh' % shz)
        with open(filE,'w') as o:
            o.write("%s" % text)
        shcount += 1
        fcount = 0
        text = ''''''
print 'shcount =',shcount
        
# qsub the files
shs = fs(shtrimDIR)
os.chdir(shtrimDIR) # want sbatch outfiles in same folder as sh file
for f in shs:
#     !sbatch $f # jupyter es el mejor
    os.system('sbatch %s' % f)
#####
