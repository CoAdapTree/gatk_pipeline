###
# execution: 02_mark_build_scatter.py dupfile /path/to/rgout.bam /path/to/fastq.gz-folder/ /path/to/ref.fa <int from 01b.py>
###

###
# when using samtools from anaconda, change default_jvm_mem_opts to "-Xms512m -Xmx4g" in `which picard` file
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

print ('fqdir =',fqdir)

# create dirs
pardir = op.dirname(fqdir)
rgdir  = op.join(fqdir,'rg_filtered_indexed_sorted_bamfiles')
shdir  = op.join(fqdir,'shfiles')
for d in [rgdir,shdir]:
    assert op.exists(d)
gvcfdir  = op.join(shdir,'gvcf_shfiles')
gatkdir  = op.join(fqdir,'gatkdir')
vcfdir   = op.join(fqdir,'vcfs')
scheddir = op.join(op.dirname(fqdir),'shfiles/gvcf_shfiles')
for d in [gvcfdir,gatkdir,vcfdir,scheddir]:
    if not op.exists(d):
        os.makedirs(d)
    
# create filenames
pool    = op.basename(op.dirname(op.dirname(rgout))) 
samp    = op.basename(rgout).split("---")[1].split("_R1R2")[0].split(".")[1]
dupdir  = op.join(fqdir,'dedup_rg_filtered_indexed_sorted_bamfiles')
dupfile = op.join(dupdir,"%s_rd.bam" % samp)
assert op.exists(dupfile)
rawvcf    = op.join(vcfdir,'raw_%s.g.vcf.gz' % samp)


#get ploidy 
PLOIDY = pickle.load(open(op.join(fqdir,'ploidy.pkl'),'rb'))
ploidy = int(PLOIDY[pool])
assert type(ploidy) == int
print ('ploidy =',ploidy)

# create sh files
shfiles = []
shcount = 0
if ploidy > 2: #poolseq
    print ("this is a poolseq file")
    scafdir = op.join(pardir,'intervals/pooled')
else:
    print ("this is an individual's file")
    scafdir = op.join(pardir,'intervals/individual')
    
scaffiles = [f for f in fs(scafdir) if f.endswith('.list')]
assert len(scaffiles) > 0
for scaff in scaffiles:
    s    = "scaff%s" % scaff.split(".list")[0].split("scaff_")[1]
    filE = op.join(gvcfdir,'%s_%s.sh' % (samp,s))
    shz  = str(tcount).zfill(3)
    vcf  = rawvcf.replace(".g.vcf.gz","_%s.g.vcf.gz" % s)
    text = '''#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30000M
#SBATCH --job-name=%(s)s%(pool)s%(shz)s
#SBATCH --export=all
#SBATCH --output=gvcf%(shz)s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu

# for debugging 
cat $0 
echo %(filE)s

source $HOME/.bashrc
module load gatk/4.0.8.1

# resubmit jobs with errors
cd $HOME/pipeline
python rescheduler.py %(fqdir)s

# fill up the queue
cd $HOME/pipeline
python scheduler.py %(fqdir)s

# call variants
gatk HaplotypeCaller --sample-ploidy %(ploidy)s -R %(ref)s --genotyping-mode DISCOVERY -ERC GVCF -I %(dupfile)s -O %(vcf)s -L %(scaff)s --minimum-mapping-quality 20

# keep running jobs until time runs out
echo 'getting help from gvcf_helper'
cd $HOME/pipeline
python gvcf_helper.py %(fqdir)s


''' % locals()
    with open(filE,'w') as o:
        o.write("%s" % text)
    # now create a symlink in scheddir
    dst = op.join(scheddir,op.basename(filE))
    if not op.exists(dst):
        os.symlink(filE,dst)


# submit to scheduler
pipedir = os.popen('echo $HOME/pipeline').read().replace("\n","")
os.system('python %s %s' % (op.join(pipedir,'scheduler.py'),fqdir))




