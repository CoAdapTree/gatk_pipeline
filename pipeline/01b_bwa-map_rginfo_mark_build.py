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
print('fqdir=',fqdir)
print(op.join(fqdir,'rginfo.pkl'))
rginfo    = pickle.load(open(op.join(fqdir,'rginfo.pkl'),'rb'))
count = 0
for samp in rginfo:
    if samp in r1out:
        print (samp)
        RG = rginfo[samp]
        count += 1
try: # make sure rginfo makes sense
    assert count == 1
except AssertionError as e:
    print ("count = %i" % count)
    print ("if count = 0, could not find samp's rginfo. If count > 1, samp names are not consistent")
    sys.exit(1)
print ("RG=",RG)
rgid = RG['rgid']
rglb = RG['rglb']
rgpl = RG['rgpl']
rgpu = RG['rgpu']
rgsm = RG['rgsm']

#bwa: fastq -> sam
sam      = op.basename(r1out).replace("R1","").replace("_trimmed.fastq","R1R2_trimmed.sam")
samdir   = op.join(fqdir,'samfiles')
samfile  = op.join(samdir,sam)
#samtools view: sam -> bam
bam      = op.basename(samfile).replace('.sam','.bam')
bamdir   = op.join(fqdir,'bamfiles')
bamfile  = op.join(bamdir,bam)
#samtools sort: bamfile -> sortfile
sort     = op.basename(bamfile).replace('.bam','_sorted.bam')
sortdir  = op.join(fqdir,'sorted_bamfiles')
sortfile = op.join(sortdir,sort)
#samtools index: sortfile -> indexfile
# nothing
#samtools view: indexfile -> filtfile
filt     = op.basename(sortfile).replace('.bam','_filtered.bam')
filtdir  = op.join(fqdir,'filtered_indexed_sorted_bamfiles')
filtfile = op.join(filtdir,filt)
#AddOrReplaceReadGroups
rg      = op.basename(filtfile).replace('.bam','_rg.bam')
rgdir   = op.join(fqdir,'rg_filtered_indexed_sorted_bamfiles')
rgfile  = op.join(rgdir,rg)
#MarkDuplicates
dupdir  = op.join(fqdir,'dedup_rg_filtered_indexed_sorted_bamfiles')
pool    = op.basename(op.dirname(op.dirname(rgfile)))
samp    = op.basename(rgfile).split("---")[1].split("_R1R2")[0].split(".")[1]
print ('samp=',samp)
mout    = op.join(rgdir,"%s.bam" % samp)
dupfile = op.join(dupdir,"%s_rd.bam" % samp)
dupstat = op.join(dupdir,"%s_rd_dupstat.txt" % samp)
#BuildBamIndex
# nothing

for d in [bamdir,samdir,sortdir,filtdir,bwashdir,rgdir,dupdir]:
    if not op.exists(d):
        os.makedirs(d)

strt = str(tcount).zfill(3)
# send it off
text = '''#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --job-name=bwa_%(samp)s
#SBATCH --time=02:59:00
#SBATCH --mem=30000M
#SBATCH --output=bwa_%(samp)s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

source $HOME/.bashrc
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.18.9

# map, sam to bam, sort by coordinate, index
bwa mem -t 32 -M %(ref)s %(r1out)s %(r2out)s > %(samfile)s
samtools view -@ 32 -Sb %(samfile)s > %(bamfile)s
samtools sort -@ 32 %(bamfile)s > %(sortfile)s
samtools index %(sortfile)s

# filter
samtools view -@ 32 -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b %(sortfile)s > %(filtfile)s
samtools index %(filtfile)s

# add rginfo
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups RGID=%(rgid)s RGLB=%(rglb)s RGPL=%(rgpl)s RGPU=%(rgpu)s RGSM=%(rgsm)s I=%(filtfile)s O=%(rgfile)s
samtools index %(rgfile)s

# remove dups
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=%(rgfile)s O=%(dupfile)s M=%(dupstat)s VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true

# Build bam index for GATK
java -jar $EBROOTPICARD/picard.jar BuildBamIndex I=%(dupfile)s

# call GVCF
cd $HOME/pipeline
python 02_scatter-gvcf.py %(rgfile)s %(fqdir)s %(ref)s %(strt)s
''' % locals()

qsubfile = op.join(bwashdir,'bwa_%s.sh' % tcount) # to avoid 2+ files written at same time
with open(qsubfile,'w') as o:
    o.write("%s" % text)
os.chdir(op.dirname(qsubfile)) # redundant to ^
os.system("sbatch %s" % qsubfile)
print ('finished!')
