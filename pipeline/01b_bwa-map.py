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


ref    = sys.argv[1]
r1out  = sys.argv[2]
r2out  = sys.argv[3]
shdir  = sys.argv[4]
tcount = sys.argv[5]

# create a dir for stdout and stderr
bwashdir = op.join(shdir,'bwa_shfiles')
if not op.exists(bwashdir):
    os.makedirs(bwashdir)
stdout = op.join(bwashdir,'bwa_%s.o' % tcount)
stderr = op.join(bwashdir,'bwa_%s.e' % tcount)

samfile = r1out.replace("R1_trimmedandfiltered.fastq.gz","R1_R2_trimmedandfiltered.sam")
text = '''#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --job-name=bwa%s
#SBATCH --export=all
#SBATCH --time=0:10:00
#SBATCH --mem=600mb

bwa mem -t 2 -M %s %s %s > %s
''' % (str(tcount).zfill(3),
       stdout,
       stderr,
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
cd(op.dirname(qsubfile)) # redundant to ^
os.system("sbatch %s" % qsubfile)
os.remove(qsubfile) # there will be way too many, n = n_reads

    
    
    
    
    
