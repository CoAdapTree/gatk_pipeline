###
# assumes:
# bwa index -a bwtsw ref.fa
# samtools faidx ref.fa
# java -jar picard.jar CreateSequenceDictionary \
#    REFERENCE=ref.fa \ 
#    OUTPUT=ref.dict
###

import sys 
import os
import random
from os import path as op

# get argument inputs
thisfile,ref,r1out,r2out,shdir,tcount = sys.argv

# create a dir for stdout and stderr
bwashdir = op.join(shdir,'bwa_shfiles')
if not op.exists(bwashdir):
    os.makedirs(bwashdir)
stdout = op.join(bwashdir,'bwa_%s.o' % tcount)
stderr = op.join(bwashdir,'bwa_%s.e' % tcount)

print 'basename(r1out) =',op.basename(r1out)
sam = op.basename(r1out).replace("R1_trimmed.fastq","R1_R2_trimmed.sam")
print 'sam =',sam
samdir = op.join(op.dirname(shdir),'samfiles')
print 'samdir =',samdir
samfile = op.join(samdir,sam)
print 'samfile =',samfile
text = '''#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --cpus-per-task=16
#SBATCH --job-name=bwa%s
#SBATCH --export=all
#SBATCH --time=01:00:00
#SBATCH --mem=50000mb
#SBATCH --output=%%x-%%j.out # needs two %% because of text replace in python

source $HOME/.bashrc

bwa mem -t 16 -M %s %s %s > %s
''' % (str(tcount).zfill(3),
       ref,
       r1out,
       r2out,
       samfile
      )

# !cat $text | qsub    # so much easier with jupyter
qsubfile = op.join(bwashdir,'bwa_%s.sh' % str(random.randint(1,1000000000)) ) # to avoid 2+ files written at same time
with open(qsubfile,'w') as o:
    o.write("%s" % text)
os.system('cd %s' % op.dirname(qsubfile))
os.chdir(op.dirname(qsubfile)) # redundant to ^
os.system("sbatch %s" % qsubfile)
os.remove(qsubfile) # there will be way too many, n = n_reads

    
    
    
    
    
