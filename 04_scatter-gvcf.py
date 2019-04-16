"""
### usage
# python 04_scatter-gvcf.py dupfile pooldir ref samp
###

###
# note that you will need to change the location of the -L interval files (idir) or (pdir)
###
"""

import sys
from coadaptree import *

### args
thisfile, dupfile, pooldir, samp = sys.argv
### 

print ('pooldir =', pooldir)

# info
parentdir = op.dirname(pooldir)
pool = op.basename(pooldir) 
ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]

# create dirs
shdir = op.join(pooldir, 'shfiles')
gvcfdir = op.join(shdir, '04_gvcf_shfiles')
vcfdir = op.join(pooldir, 'vcfs')
scheddir = op.join(parentdir, 'shfiles/gvcf_shfiles')
for d in [shdir, gvcfdir, vcfdir, scheddir]:
    makedir(d)
    
# create filenames
dupdir = op.dirname(dupfile)
rawvcf = op.join(vcfdir, 'raw_%s.g.vcf.gz' % samp)

#get ploidy 
ploidy = pklload(op.join(parentdir, 'ploidy.pkl'))[pool]
print ('ploidy =', ploidy)

# create sh files
shfiles = []
shcount = 0
if ploidy > 2: #poolseq
    print ("this is a poolseq file")
else:
    print ("this is an individual's file")
scafdir = op.join(op.dirname(ref), 'intervals')
    
scaffiles = [f for f in fs(scafdir) if f.endswith('.list')]
os.system('echo found %s intervalfiles' % str(len(scaffiles)))
email_text = get_email_info(parentdir, '04')
for scaff in scaffiles:
    s = "scatter-%s" % scaff.split(".list")[0].split("batch_")[1]
    filE = op.join(gvcfdir, f'{pool}-{samp}-{s}.sh')
    vcf = rawvcf.replace(".g.vcf.gz", "_%s.g.vcf.gz" % s)
    tbi = vcf.replace(".gz", ".gz.tbi")
    text = f'''#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --job-name={pool}-{samp}-{s}
#SBATCH --output={pool}-{samp}-{s}_%j.out 
{email_text}

# for debugging 
cat $0 
echo {filE}

source $HOME/.bashrc

# resubmit jobs with errors
python $HOME/gatk_pipeline/rescheduler.py {pooldir}

# fill up the queue
python scheduler.py {pooldir}

# call variants
module load gatk/4.1.0.0
gatk HaplotypeCaller --sample-ploidy {ploidy} -R {ref} --genotyping-mode DISCOVERY -ERC GVCF -I {dupfile} -O {vcf} -L {scaff} --minimum-mapping-quality 20
module unload gatk

# keep running jobs until time runs out
echo 'getting help from gvcf_helper'
python $HOME/gatk_pipeline/gvcf_helper.py {pooldir} {tbi}

'''
    with open(filE, 'w') as o:
        o.write("%s" % text)
    # now create a symlink in scheddir
    dst = op.join(scheddir, op.basename(filE))
    if not op.exists(dst):
        os.symlink(filE, dst)

# submit to scheduler
os.system('python $HOME/gatk_pipeline/scheduler.py %s' % pooldir)
