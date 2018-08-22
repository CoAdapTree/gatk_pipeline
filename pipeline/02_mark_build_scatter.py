###
# execution: 02_mark_build_scatter.py /path/to/rgout /path/to/fastq.gz-folder/ /path/to/ref.fa <int passed from 01b.py>
###

###
# when using samtools from anaconda, change default_jvm_mem_opts to "-Xms512m -Xmx4g" in <which picard> file
###

### imports 
import sys
import os
from os import path as op
from os import listdir
import pickle
def ls(DIR):
    return sorted([f for f in listdir(DIR)])
def fs (DIR):
    return sorted([op.join(DIR,f) for f in ls(DIR)])
###

### args
thisfile, rgout, fqdir, ref, tcount = sys.argv
### 

# create dirs
rgdir = op.join(fqdir,'rg_filtered_indexed_sorted_bamfiles')
shdir    = op.join(fqdir,'shfiles')
for d in [rgdir,shdir]:
    assert op.exists(d)
# mergeoutdir = op.join(fqdir,'merged_rg_filtered_indexed_sorted_bamfiles')
dupdir  = op.join(fqdir,'dedup_rg_filtered_indexed_sorted_bamfiles')
mshdir  = op.join(shdir,'mark_shfiles')
gatkdir = op.join(fqdir,'gatkdir')
vcfdir  = op.join(fqdir,'vcfs')
for d in [dupdir,mshdir,gatkdir,vcfdir]:
    if not op.exists(d):
        os.makedirs(d)

    
# create filenames
pool      = op.basename(op.dirname(op.dirname(rgout))) 
samp      = op.basename(rgout).split("---")[1].split("_R1R2")[0].split(".")[1]
mout      = op.join(rgdir,"%s.bam" % samp)
dupfile   = op.join(dupdir,"%s_rd.bam" % samp)
dupstat   = op.join(dupdir,"%s_rd_dupstat.txt" % samp)
gatklist  = op.join(gatkdir,'realignment_targets.list')
realigned = op.join(gatkdir,'realigned_reads.bam')
rawvcf    = op.join(vcfdir,'raw_%s.g.vcf.gz' % samp)
snpvcf    = op.join(vcfdir,'snps_%s.g.vcf.gz' % samp)
indelvcf  = op.join(vcfdir,'indel_%s.g.vcf.gz' % samp)

#get ploidy and rginfo
PLOIDY = pickle.load(open(op.join(fqdir,'ploidy.pkl')))
ploidy = int(PLOIDY[pool])
print 'ploidy =',ploidy
rginfo = pickle.load(open(op.join(fqdir,'rginfo.pkl')))

# create sh files
shfiles = []
if ploidy > 2: #poolseq
    print "this is a poolseq file"
    # run it
    text = '''#!/bin/bash
#SBATCH --time=14-00:00 # 14 days
#SBATCH --nodes=1
#SBATCH --mem=140000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=%smark%s
#SBATCH --export=all
#SBATCH --output=mark%s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

source $HOME/.bashrc
module load gatk/4.0.0.0

# remove dups
picard MarkDuplicates I=%s O=%s M=%s VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true

# Build bam index for GATK
picard BuildBamIndex I=%s

# call variants
module load gatk/4.0.0.0
gatk HaplotypeCaller --sample-ploidy %s -R %s --genotyping-mode DISCOVERY -ERC GVCF -I %s -O %s  


''' % (pool,
       str(tcount).zfill(3),
       str(tcount).zfill(3),
       rgout, dupfile,  dupstat,
       dupfile,
       ploidy,  ref,  dupfile,  rawvcf,
      )
# ''' % ("%s_%s" % (scaff,str(tcount).zfill(3)),
#        "%s_%s" % (scaff,str(tcount).zfill(3)),
#        rgout, dupfile,  dupstat,
#        dupfile,
#        ploidy,  ref,  dupfile,  rawvcf,
#       )
    filE = op.join(mshdir,'%s.sh' % samp)
    with open(filE,'w') as o:
        o.write("%s" % text)
    shfiles.append(filE)
else: # non poolseq
    print "this is an individual's file"
    text = '''#!/bin/bash
#SBATCH --time=7-00:00 # 7 days
#SBATCH --nodes=1
#SBATCH --mem=140000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=%smark%s
#SBATCH --export=all
#SBATCH --output=mark%s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

source $HOME/.bashrc
module load gatk/4.0.0.0

# remove dups
picard MarkDuplicates I=%s O=%s M=%s VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true

# Build bam index for GATK
picard BuildBamIndex I=%s

# call variants
module load gatk/4.0.0.0
gatk HaplotypeCaller --sample-ploidy %s -R %s --genotyping-mode DISCOVERY -ERC GVCF -I %s -O %s  


''' % (pool,
       str(tcount).zfill(3),
       str(tcount).zfill(3),
       rgout, dupfile,  dupstat,
       dupfile,
       ploidy,  ref,  dupfile,  rawvcf,
      )



    filE = op.join(mshdir,'%s.sh' % samp)
    with open(filE,'w') as o:
        o.write("%s" % text)
    shfiles.append(filE)

print shfiles

for s in shfiles:
    os.chdir(op.dirname(s))
    os.system('sbatch %s' % s)







