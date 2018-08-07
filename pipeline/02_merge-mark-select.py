import sys
import os
from os import path as op
from os import listdir as ls
def fs (DIR):
    return (sorted([op.join(DIR,f) for f in ls(DIR)]))

###
# execution: 02_merge-bams-and-mark-dups.py /path/to/fastq.gz-folder/ ploidy<int>
###

thisfile, fqdir, ref, ploidy = sys.argv

# create dirs
mergedir = op.join(fqdir,'filtered_indexed_sorted_bamfiles')
shdir    = op.join(fqdir,'shfiles')
for d in [mergedir,shdir]:
    assert op.exists(d)
mergeoutdir = op.join(fqdir,'merged_filtered_indexed_sorted_bamfiles')
dupdir      = op.join(fqdir,'dedup_merged_filtered_indexed_sorted_bamfiles')
mshdir      = op.join(shdir,'merge_shfiles')
gatkdir     = op.join(fqdir,'gatkdir')
vcfdir    = op.join(fqdir,'vcfs')
for d in [mergeoutdir,dupdir,mshdir,gatkdir,vcfdir]:
    if not op.exists(d):
        os.makedirs(d)
    
# create filenames
mfiles    = [f for f in fs(mergedir) if f.endswith('bam')] # bam files to be merged
out       = ".".join([x for x in op.basename(mfiles[0]).split(".")[:3]]) # name by lane info: eg 'paired_HI.0748.006'
mout      = op.join(mergeoutdir,"%s.bam" % out)
dupfile   = op.join(dupdir,"%s_rd.bam" % out)
dupstat   = op.join(dupdir,"%s_rd_dupstat.txt" % out)
gatklist  = op.join(gatkdir,'realignment_targets.list')
realigned = op.join(gatkdir,'realigned_reads.bam')
rawvcf    = op.join(vcfdir,'raw_%s.vcf' % out)
snpvcf    = op.join(vcfdir,'snps_%s.vcf' % out)
indelvcf  = op.join(vcfdir,'indel_%s.vcf' % out)

# run it
text = '''#!/bin/bash
#SBATCH --account=def-saitken
#SBATCH --cpus-per-task=32
#SBATCH --job-name=mergebams
#SBATCH --export=all
#SBATCH --time=03:00:00
#SBATCH --mem=500000mb
#SBATCH --output=%%x-%%j.out 

source $HOME/.bashrc

# merge and index
samtools merge -@ 32 -f %s %s
samtools index %s

# remove dups
picard MarkDuplicates I=%s O=%s M=%s VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true

# add read groups?

# Build bam index for GATK
picard BuildBamIndex I=%s

# Realign around INDELs
module load gatk/3.8
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 32 -R %s -I %s -o %s
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s

# call variants
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T HaplotypeCaller -ploidy %s -nt 32 -R %s -I %s -o %s
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SelectVariants -R %s -V %s -selectType SNP -o %s
java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T SelectVariants -R %s -V %s -selectType INDEL -o %s
''' % (mout, " ".join([x for x in mfiles]),
       mout,
       mout, dupfile,  dupstat,
       dupfile,
       ref,  dupfile,  gatklist,
       ref,  dupfile,  gatklist,  realigned,
       ploidy,  ref,  realigned,  rawvcf,
       ref,  rawvcf,  snpvcf,
       ref,  rawvcf,  indelvcf
      )
filE = op.join(mshdir,"merge-dedup.sh")
with open(filE,'w') as o:
    o.write("%s" % text)
os.chdir(mshdir)
os.system("sbatch %s" % filE)