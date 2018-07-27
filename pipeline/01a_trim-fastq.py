from __future__ import division
####
# assumes fastq.gz files
###

#####
# execution:
# python 01a_trim-fastq.py /path/to/fastq.gz-files/ /path/to/ref.fa 
#####

#####
# changable
lensh = 950 # number of sh files that can be sbatched
#####

#####
# imports
import sys
from os import path as op

# args
fqdir   = sys.argv[1] # path to raw fastq files
ref     = sys.argv[2] # path to reference genome used for mapping (/home/lindb/scratch/ptaeda.v1.01.reduced.pseudo.fasta)
for i,arg in enumerate([ref,fqdir]):
    # make sure the args exist
    try:
        assert op.exists(arg)
    except AssertionError as e:
        print "The %s'th argument does not exist in the specified path" % str(i)
        sys.exit(1)

# more imports and aliases
pipedir = op.dirname(op.abspath(sys.argv[0]))         # this is where the git repo (pipeline) is pulled
try:
    imp = open(op.join(pipedir,'pythonimports.py')).read()
except IOError as e:
    print "Error: couldn't find pythonimports.py file in %s. Check implementation of directory tree." % pipedir
    sys.exit(1)
exec(imp)

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
cd(fqdir)
gzfiles = [op.abspath(f) for f in fs(fqdir) if 'R1' in f]
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

# determine how many commands per sh file
if len(seq_pairs) <= lensh:
    # one command per sh file
    ceil = 1
else:
    # multiple commands per sh file
    ceil = math.ceil(len(seq_pairs)/lensh)

shcount = 0
fcount  = 0
tcount  = 0
text    = ''''''
for s in seq_pairs:
#     print s
    r1    = op.abspath(s[0])
    r1out = op.join(trimDIR,op.basename(r1).replace('.fastq.gz','_trimmed.fastq'))
    r2    = op.abspath(s[1])
    r2out = op.join(trimDIR,op.basename(r2).replace('.fastq.gz','_trimmed.fastq'))
    cmd   = '''fastp -i %s -o %s -I %s -O %s -g --cut_window_size 5 --cut_mean_quality 30 --qualified_quality_phred 30 --unqualified_percent_limit 20 --n_base_limit 5 --length_required 75 -h %s.html --cut_by_quality3 --thread 16 --json %s.json


cd $HOME/pipeline
# once finished, map using bwa mem 
python 01b_bwa-map.py %s %s %s %s %s           
''' % (r1,  r1out,
       r2,  r2out,
       r1out.replace("R1_trimmed.fastq","R1_R2_trimmed_stats"),
       r1out.replace("R1_trimmed.fastq","R1_R2_trimmed"),
       ref,  r1out,  r2out,  shdir,  str(tcount).zfill(3)
      )
    text = text + cmd
    
    fcount += 1
    tcount += 1
    if fcount == ceil or tcount == len(seq_pairs):
        shz = str(shcount).zfill(3)
        header = '''#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --job-name=trim%s
#SBATCH --export=all
#SBATCH --time=01:00:00
#SBATCH --mem=1000mb
#SBATCH --cpus-per-task=16
''' % (shz)
        text = header + text
        filE = op.join(shtrimDIR,'trim_%s.sh' % str(shcount).zfill(3))
        with open(filE,'w') as o:
            o.write("%s" % text)
        shcount += 1
        fcount = 0
        text = ''''''

# qsub the files
shs = fs(shtrimDIR)
os.system('cd %s' % shtrimDIR)
cd(shtrimDIR) # just in case os.system(cd) doesn't work, we want outfiles in same directory as sh files
for f in shs:
#     !qsub $f # jupyter es el mejor
    os.system('sbatch %s' % f)
#####