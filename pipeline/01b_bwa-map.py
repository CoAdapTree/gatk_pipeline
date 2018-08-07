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


### create dirs and filenames
bwashdir  = op.join(shdir,'bwa_shfiles')
fqdir     = op.dirname(shdir)
#bwa: fastq -> sam
sam       = op.basename(r1out).replace("R1_trimmed.fastq","R1_R2_trimmed.sam")
samdir    = op.join(fqdir,'samfiles')
samfile   = op.join(samdir,sam)
#samtools view: sam -> bam
bam       = op.basename(samfile).replace('.sam','.bam')
bamdir    = op.join(fqdir,'bamfiles')
bamfile   = op.join(bamdir,bam)
#samtools sort: bamfile -> sortfile
sort      = op.basename(bamfile).replace('.bam','_sorted.bam')
sortdir   = op.join(fqdir,'sorted_bamfiles')
sortfile  = op.join(sortdir,sort)
#samtools index: sortfile -> indexfile
# index     = op.basename(sortfile).replace('.bam','.bai')
# indexfile = op.join(sortdir,index)
#samtools view: indexfile -> filtfile
filt      = op.basename(sortfile).replace('.bam','_filtered.bam')
filtdir   = op.join(fqdir,'filtered_indexed_sorted_bamfiles')
filtfile  = op.join(filtdir,filt)

for d in [bamdir,samdir,sortdir,filtdir,bwashdir]:
    if not op.exists(d):
        os.makedirs(d)

# send it off
text = '''#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --cpus-per-task=32
#SBATCH --job-name=bwa%s
#SBATCH --export=all
#SBATCH --time=03:00:00
#SBATCH --mem=50000mb
#SBATCH --output=%%x-%%j.out # needs two %% because of text replace in python

source $HOME/.bashrc

bwa mem -t 32 -M %s %s %s > %s
samtools view -@ 32 -Sb %s > %s
samtools sort -@ 32 %s > %s
samtools index %s

#filter
samtools view -@ 32 -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b %s > %s
samtools index %s
''' % (str(tcount).zfill(3),
       ref,  r1out,  r2out,  samfile,
       samfile,  bamfile,
       bamfile,  sortfile,
       sortfile,
       sortfile, filtfile,
       filtfile
      )

# !cat $text | qsub    # so much easier with jupyter
qsubfile = op.join(bwashdir,'bwa_%s.sh' % str(random.randint(1,1000000000)) ) # to avoid 2+ files written at same time
with open(qsubfile,'w') as o:
    o.write("%s" % text)
# os.system('cd %s' % op.dirname(qsubfile))
os.chdir(op.dirname(qsubfile)) # redundant to ^
os.system("sbatch %s" % qsubfile)
# os.remove(qsubfile) 

    
    
    
    
    

