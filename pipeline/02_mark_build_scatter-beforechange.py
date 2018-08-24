###
# execution: 02_mark_build_scatter.py /path/to/rgout /path/to/fastq.gz-folder/ /path/to/ref.fa <int passed from 01b.py>
###

###
# when using samtools from anaconda, change default_jvm_mem_opts to "-Xms512m -Xmx4g" in <which picard> file
###

###
# note that you will need to change the location of the -L interval files (idir) or (pdir)
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
[os.remove(f) for f in fs(mshdir) if f.endswith('.sh')] # don't want to submit the wrong files
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
shcount = 0
if ploidy > 2: #poolseq
    print "this is a poolseq file"
    # run it
    pdir = '/scratch/lindb/testdata/intervals/pooled'
    pfiles = [f for f in fs(pdir) if f.endswith('.list')]
    for scaff in pfiles:
        s = "scaff%s" % scaff.split(".list")[0].split("scaff_")[1]
        text = '''#!/bin/bash
#SBATCH --time=02:59:00
#SBATCH --nodes=1
#SBATCH --mem=140000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=%s%s%s
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
gatk HaplotypeCaller --sample-ploidy %s -R %s --genotyping-mode DISCOVERY -ERC GVCF -I %s -O %s -L %s

''' % (s, pool, str(tcount).zfill(3),
       str(tcount).zfill(3),
       rgout, dupfile,  dupstat,
       dupfile,
       ploidy,  ref,  dupfile,  rawvcf.replace(".g.vcf.gz","_%s.g.vcf.gz" % s),
       scaff
      )
        filE = op.join(mshdir,'%s_%s.sh' % (samp,s))
        with open(filE,'w') as o:
            o.write("%s" % text)
        shfiles.append(filE)
else: # non poolseq
    print "this is an individual's file"
    idir = '/scratch/lindb/testdata/intervals/individual'
    ifiles = [f for f in fs(idir) if f.endswith('.list')]
    print 'len(ifiles)=',len(ifiles)
    for scaff in ifiles:
        s = "scaff%s" % scaff.split(".list")[0].split("scaff_")[1]
        text = '''#!/bin/bash
#SBATCH --time=02:59:00
#SBATCH --nodes=1
#SBATCH --mem=140000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=%s%s%s
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
gatk HaplotypeCaller --sample-ploidy %s -R %s --genotyping-mode DISCOVERY -ERC GVCF -I %s -O %s -L %s

''' % (s,  pool,  str(tcount).zfill(3),
       str(tcount).zfill(3),
       rgout,  dupfile,  dupstat,
       dupfile,
       ploidy,  ref,  dupfile,  rawvcf.replace(".g.vcf.gz","_%s.g.vcf.gz" % s),
       scaff
      )
        filE = op.join(mshdir,'%s_%s.sh' % (samp,s))
        print filE
        with open(filE,'w') as o:
            o.write("%s" % text)
        shfiles.append(filE)

print shfiles

for s in shfiles:
    os.chdir(op.dirname(s))
    os.system('sbatch %s' % s)







