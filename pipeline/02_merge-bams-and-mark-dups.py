import sys
import os
from os import path as op
from os import listdir as ls
def fs (DIR):
    return (sorted([op.join(DIR,f) for f in ls(DIR)]))

###
# execution: 02_merge-bams-and-mark-dups.py /path/to/fastq.gz-files/
###

thisfile, fqdir = sys.argv

# create dirs
mergedir = op.join(fqdir,'filtered_indexed_sorted_bamfiles')
shdir    = op.join(fqdir,'shfiles')
for d in [mergedir,shdir]:
    assert op.exists(d)
mergeoutdir = op.join(fqdir,'merged_filtered_indexed_sorted_bamfiles')
dupdir = op.join(fqdir,'dedup_merged_filtered_indexed_sorted_bamfiles')
mshdir = op.join(shdir,'merge_shfiles')
for d in [mergeoutdir,dupdir,mshdir]:
    if not op.exists(d):
        os.makedirs(d)
    
# create filenames
mfiles  = fs(mergedir) # bam files to be merged
out     = ".".join([x for x in op.basename(mfiles[0]).split(".")[:3]]) # name by lane info: eg 'paired_HI.0748.006'
mout    = op.join(mergeoutdir,"%s.bam" % out)
dupfile = op.join(dupdir,"%s_rd.bam" % out)
dupstat = op.join(dupdir,"%s_rd_dupstat.txt" % out)

# run it
text = '''#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --cpus-per-task=32
#SBATCH --job-name=mergebams
#SBATCH --export=all
#SBATCH --time=03:00:00
#SBATCH --mem=50000mb
#SBATCH --output=%%x-%%j.out 

# merge and index
samtools merge -@ 32 -f %s %s
samtools index %s

# remove dups
picard MarkDuplicates.jar I=%s O=%s M=%s VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true

# create read groups?
''' % (mout, " ".join([x for x in mfiles]),
       mout,
       mout, dupfile, dupstat
      )
filE = op.join(mshdir,"merge-dedup.sh")
with open(filE,'w') as o:
    o.write("%s" % text)
os.chdir(mshdir)
os.system("sbatch %s" % filE)