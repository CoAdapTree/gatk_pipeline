###
# assumes:
# bwa index -a bwtsw ref.fa
# samtools faidx ref.fa
# java -jar picard.jar CreateSequenceDictionary \
#    REFERENCE=ref.fa \ 
#    OUTPUT=ref.dict
###

###
# usage: 01b_bwa-map_rginfo_mark_build.py /path/to/ref.fa /path/to/trimmedR1.fastq /path/to/trimmedR2.fastq /sbatch/dir/ <int>
###

import sys 
import os
import random
from os import path as op
import pickle

# get argument inputs
thisfile,ref,r1out,r2out,shdir,tcount = sys.argv


### create dirs and filenames
bwashdir  = op.join(shdir,'bwa_shfiles')
fqdir     = op.dirname(shdir)
rginfo    = pickle.load(open(op.join(fqdir,'rginfo.pkl')))
count = 0
for samp in rginfo:
    if samp in r1out:
        print samp
        RG = rginfo[samp]
        count += 1
try: # make sure rginfo makes sense
    assert count == 1
except AssertionError as e:
    print "count = %i" % count
    print "if count = 0, could not find samp's rginfo. If count > 1, samp names are not consistent"
    sys.exit(1)
print RG

#bwa: fastq -> sam
sam       = op.basename(r1out).replace("R1","").replace("_trimmed.fastq","R1R2_trimmed.sam")
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
# nothing
#samtools view: indexfile -> filtfile
filt      = op.basename(sortfile).replace('.bam','_filtered.bam')
filtdir   = op.join(fqdir,'filtered_indexed_sorted_bamfiles')
filtfile  = op.join(filtdir,filt)
#AddOrReplaceReadGroups
rg        = op.basename(filtfile).replace('.bam','_rg.bam')
rgdir     = op.join(fqdir,'rg_filtered_indexed_sorted_bamfiles')
rgfile    = op.join(rgdir,rg)
#MarkDuplicates
dupdir  = op.join(fqdir,'dedup_rg_filtered_indexed_sorted_bamfiles')
pool    = op.basename(op.dirname(op.dirname(rgout)))
samp    = op.basename(rfile).split("---")[1].split("_R1R2")[0].split(".")[1]
mout    = op.join(rgdir,"%s.bam" % samp)
dupfile = op.join(dupdir,"%s_rd.bam" % samp)
dupstat = op.join(dupdir,"%s_rd_dupstat.txt" % samp)
#BuildBamIndex
# nothing

for d in [bamdir,samdir,sortdir,filtdir,bwashdir,rgdir,dupdir]:
    if not op.exists(d):
        os.makedirs(d)
        
# send it off
text = '''#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --cpus-per-task=32
#SBATCH --job-name=bwa%s
#SBATCH --export=all
#SBATCH --time=02:59:00
#SBATCH --mem=30000mb
#SBATCH --output=bwa%s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

source $HOME/.bashrc

# map, sam to bam, sort by coordinate, index
bwa mem -t 32 -M %s %s %s > %s
samtools view -@ 32 -Sb %s > %s
samtools sort -@ 32 %s > %s
samtools index %s

# filter
samtools view -@ 32 -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b %s > %s
samtools index %s

# add rginfo
picard AddOrReplaceReadGroups RGID=%s RGLB=%s RGPL=%s RGPU=%s RGSM=%s I=%s O=%s
samtools index %s

# remove dups
picard MarkDuplicates I=%s O=%s M=%s VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true

# Build bam index for GATK
picard BuildBamIndex I=%s

# call GVCF
cd $HOME/pipeline
python 02_scatter-gvcf.py %s %s %s
''' % (str(tcount).zfill(3), str(tcount).zfill(3),
       ref,  r1out,  r2out,  samfile,
       samfile,  bamfile,
       bamfile,  sortfile,
       sortfile,
       sortfile, filtfile,
       filtfile,
       RG['rgid'], RG['rglb'], RG['rgpl'], RG['rgpu'], RG['rgsm'], filtfile, rgfile,
       rgfile,
       rgfile, dupfile, dupstat,
       dupfile,
       rgfile, fqdir, ref, str(tcount).zfill(3)
      )

# !cat $text | qsub    # so much easier with jupyter
qsubfile = op.join(bwashdir,'bwa_%s.sh' % tcount) # to avoid 2+ files written at same time
with open(qsubfile,'w') as o:
    o.write("%s" % text)
# os.system('cd %s' % op.dirname(qsubfile))
os.chdir(op.dirname(qsubfile)) # redundant to ^
os.system("sbatch %s" % qsubfile)

