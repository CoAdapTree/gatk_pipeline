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

print 'fqdir =',fqdir

# create dirs
rgdir = op.join(fqdir,'rg_filtered_indexed_sorted_bamfiles')
shdir    = op.join(fqdir,'shfiles')
for d in [rgdir,shdir]:
    assert op.exists(d)
# dupdir  = op.join(fqdir,'dedup_rg_filtered_indexed_sorted_bamfiles')
gvcfdir  = op.join(shdir,'gvcf_shfiles')
gatkdir = op.join(fqdir,'gatkdir')
vcfdir  = op.join(fqdir,'vcfs')
scheddir = op.join(op.dirname(fqdir),'shfiles/gvcf_shfiles')
for d in [gvcfdir,gatkdir,vcfdir,scheddir]:
    if not op.exists(d):
        os.makedirs(d)
#[os.remove(f) for f in fs(gvcfdir) if f.endswith('.sh')] # don't want to submit the wrong files

    
# create filenames
pool    = op.basename(op.dirname(op.dirname(rgout))) 
samp    = op.basename(rgout).split("---")[1].split("_R1R2")[0].split(".")[1]
dupdir  = op.join(fqdir,'dedup_rg_filtered_indexed_sorted_bamfiles')
dupfile = op.join(dupdir,"%s_rd.bam" % samp)
assert op.exists(dupfile)
rawvcf    = op.join(vcfdir,'raw_%s.g.vcf.gz' % samp)


#get ploidy 
PLOIDY = pickle.load(open(op.join(fqdir,'ploidy.pkl')))
ploidy = int(PLOIDY[pool])
assert type(ploidy) == int
print 'ploidy =',ploidy

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
#SBATCH --time=11:59:00
#SBATCH --nodes=1
#SBATCH --mem=60000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=%s%s%s
#SBATCH --export=all
#SBATCH --output=mark%s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu


source $HOME/.bashrc
module load gatk/4.0.0.0

# call variants
gatk HaplotypeCaller --sample-ploidy %s -R %s --genotyping-mode DISCOVERY -ERC GVCF -I %s -O %s -L %s --minimum-mapping-quality 20

cd $HOME/pipeline
python scheduler.py %s

''' % (s,  pool,  str(tcount).zfill(3),
       str(tcount).zfill(3),
       ploidy,  ref,  dupfile,  rawvcf.replace(".g.vcf.gz","_%s.g.vcf.gz" % s),
       scaff,
       fqdir       
      )
        filE = op.join(gvcfdir,'%s_%s.sh' % (samp,s))
        with open(filE,'w') as o:
            o.write("%s" % text)
        dst = op.join(scheddir,op.basename(filE))
        if not op.exists(dst):
            os.symlink(filE,dst)
#         shfiles.append(filE)
        
else: # non poolseq
    print "this is an individual's file"
    idir = '/scratch/lindb/testdata/intervals/individual'
    ifiles = [f for f in fs(idir) if f.endswith('.list')]
    print 'len(ifiles)=',len(ifiles)
    for scaff in ifiles:
        s = "scaff%s" % scaff.split(".list")[0].split("scaff_")[1]
        text = '''#!/bin/bash
#SBATCH --time=2:59:00
#SBATCH --nodes=1
#SBATCH --mem=20000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=%s%s%s
#SBATCH --export=all
#SBATCH --output=mark%s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu


source $HOME/.bashrc
module load gatk/4.0.0.0

# call variants
gatk HaplotypeCaller --sample-ploidy %s -R %s --genotyping-mode DISCOVERY -ERC GVCF -I %s -O %s -L %s --minimum-mapping-quality 20 

cd $HOME/pipeline
python scheduler.py %s

''' % (s,  pool,  str(tcount).zfill(3),
       str(tcount).zfill(3),
       ploidy,  ref,  dupfile,  rawvcf.replace(".g.vcf.gz","_%s.g.vcf.gz" % s),
       scaff,
       fqdir
      )
        filE = op.join(gvcfdir,'%s_%s.sh' % (samp,s))
        with open(filE,'w') as o:
            o.write("%s" % text)
        dst = op.join(scheddir,op.basename(filE))
        if not op.exists(dst):
            os.symlink(filE,dst)
#         shfiles.append(filE)

# submit to scheduler
pipedir = os.popen('echo $HOME/pipeline').read().replace("\n","")
os.system('python %s %s' % (op.join(pipedir,'scheduler.py'),fqdir))

#the old way I was doing it
# print shfiles
# for s in shfiles:
#     os.chdir(op.dirname(s))
#     os.system('sbatch %s' % s)








