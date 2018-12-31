### usage
# python 03a_combine_and_genotype_by_pool.py /path/to/parentdir-used-in-00_start-pipeline.py-command 
### 

### FIX
# customize time for var levels of ploidy 
# merge with 03b
### 

### imports
import sys
import os
from os import path as op
from os import listdir
import pickle
import numpy as np
def uni(mylist):
    return (np.unique(mylist).tolist())
def ls(DIR):
    return sorted([f for f in listdir(DIR)])
def fs (DIR):
    return sorted([op.join(DIR,f) for f in ls(DIR)])
def createdirs(dirs):
    assert type(dirs) == list
    for d in dirs:
        if not op.exists(d):
            os.makedirs(d)
### 

### args
thisfile, parentdir = sys.argv
if parentdir.endswith("/"):
    parentdir = parentdir[:-1]
###

### reqs
samp2pool = pickle.load(open(op.join(parentdir,'samp2pool.pkl'),'rb'))
poolref   = pickle.load(open(op.join(parentdir,'poolref.pkl'),'rb'))   #key=pool val=/path/to/ref.fa
ploidy    = pickle.load(open(op.join(parentdir,'ploidy.pkl'),'rb'))    #not used yet, but maybe if pooled needs more time
pools = uni(list(samp2pool.values()))
###

# get a list of subdirectory pool dirs created earlier in pipeline
pooldirs = []
for p in pools:
    pooldir = op.join(parentdir,p)
    if op.exists(pooldir):
        pooldirs.append(pooldir)

# get a list of files and partition by those that need to be combined
spp = {}
for d in pooldirs:
    vcfdir = op.join(d,'vcfs')
    # only get vcfs if they're ready for next step (those with .tbi files)
    files = [f.replace(".gz.tbi",".gz") for f in fs(vcfdir) if f.endswith('.gz.tbi')] 
    count = 0
    for f in files:
        count += 1
        splits = op.basename(f).replace("raw_","").split("_")
        sp = splits[0]
        if not sp in spp:
            spp[sp] = {}
        if splits[1].startswith('p'):
#             kind  = 'pooled'    ################################## CHANGE THIS BACK FOR MULTI-POOL GENOTYPING!!!!!!!!!!!!!!!!
            kind2 = op.basename(d) # for single pool genotyping
        else:
            kind = 'unpooled'
#         if not kind in spp[sp]: ################################## CHANGE THIS BACK FOR MULTI-POOL GENOTYPING!!!!!!!!!!!!!!!!
#             spp[sp][kind]  = []
        if not kind2 in spp[sp]:
            spp[sp][kind2] = []
#         spp[sp][kind].append(f) ################################## CHANGE THIS BACK FOR MULTI-POOL GENOTYPING!!!!!!!!!!!!!!!!
        spp[sp][kind2].append(f)
for kind in spp['DF']:
    print(kind,len(spp['DF'][kind]))

# create some dirs
outdir = op.join(parentdir,'snps')
shdir  = op.join(parentdir,'shfiles/select_variants_within_and_across')
createdirs([outdir,shdir])

# make sh files
alreadycreated = [f for f in fs(shdir) if f.endswith('.sh') and 'swp' not in f]
print('len(alreadycreated)=',len(alreadycreated))
shfiles = []
for sp in spp:
    print(sp)
    for kind in spp[sp]:
        writesh = False
        print('\t',kind)
        # get the files that need to be combined
        groups = {}
        pools = []
        for f in spp[sp][kind]:
            POOL = op.basename(op.dirname(op.dirname(f)))
            if not POOL in pools:
                pools.append(POOL)
            scaff = op.basename(f).split("scaff")[1].split(".g.v")[0]
            if scaff not in groups:
                groups[scaff] = []
            groups[scaff].append(f)
        pools = '---'.join([x for x in pools])
        # get ref.fa
        ref = poolref[POOL] # doesn't matter which pool I use, will be same ref (POOL is from iter)
        # get commands
        for scaff in groups:
            cmds = ''''''
            combfile = op.join(outdir,'%s--%s_combined.vcf.gz' % (pools,scaff))
            gfile    = combfile.replace("_combined.vcf.gz","_genotyped.vcf.gz")
            snpfile  = gfile.replace("_genotyped.vcf.gz","_snps.vcf.gz")
            if len(groups[scaff]) > 1:
                # these files need to be combined
                varcmd = '--variant ' + ' --variant '.join([x for x in sorted(groups[scaff])])
                cmd    = '''echo COMBINING GVCFs
gatk CombineGVCFs -R %(ref)s %(varcmd)s -O %(combfile)s 

echo GENOTYPING GVCFs
gatk --java-options "-Xmx4g" GenotypeGVCFs -R %(ref)s -V %(combfile)s -O %(gfile)s 

echo SELECTING VARIANTS
gatk SelectVariants -R %(ref)s -V %(gfile)s --select-type-to-include SNP --restrictAllelesTo BIALLELIC -O %(snpfile)s

''' % locals()
            elif len(groups[scaff]) == 1 and kind != 'pooled': # otherwise if jobs aren't finished, cross-pool sh will skip combining
                # no need to combine, just symlink vcf and tbi file
                f = groups[scaff][0] # skip symlink and just specify original file
#                 tbi = f.replace(".gz",".gz.tbi")
#                 tbilink = combfile.replace(".gz",".gz.tbi")
#                 try:
#                     if not op.exists(combfile):
#                         os.symlink(f,combfile)
#                     if not op.exists(tbilink):
#                         os.symlink(tbi,tbilink)
#                 except:
#                     # this should only happen if the f or tbi file get moved during execution (somehow) and no longer exist
#                     # this probably won't ever happen
#                     continue
#                     os.system('echo continuing')
#                     print('continuing')
                gfile    = combfile.replace("_combined.vcf.gz","_genotyped.vcf.gz")
                snpfile  = gfile.replace("_genotyped.vcf.gz","_snps.vcf.gz")
                cmd = '''echo GENOTYPING GVCFs
gatk --java-options "-Xmx4g" GenotypeGVCFs -R %(ref)s -V %(f)s -O %(gfile)s 

echo SELECTING VARIANTS
gatk SelectVariants -R %(ref)s -V %(gfile)s --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O %(snpfile)s

''' % locals()
            else:
                print(sp,kind,scaff,'has no files')
            cmds = cmds + cmd
            if not cmds == '''''':
                file = op.join(shdir,'%(pools)s--%(scaff)s.sh' % locals())
                if not file in alreadycreated: # this way I can start running the genotyping phase before all files are ready
                    print('\t','\t',file)
                    text = '''#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --ntasks=1
#SBATCH --mem=2000M
#SBATCH --cpus-per-task=1
#SBATCH --job-name=%(pools)s--%(scaff)s
#SBATCH --export=all
#SBATCH --output=%(pools)s--%(scaff)s---%%j.out 

source $HOME/.bashrc
cat $0
echo shfile = %(file)s

# requeue jobs with errors
cd $HOME/pipeline
python genotyping_rescheduler.py %(parentdir)s

# fill up the queue
python genotyping_scheduler.py %(parentdir)s

# genotype current file
module load gatk/4.0.8.1
%(cmds)s

# give time for slurm scheduler to check mem and kill job
sleep 60

# keep running jobs until time runs out
echo getting help from genotyping_helper
cd $HOME/pipeline
python genotyping_helper.py %(parentdir)s

''' % locals()
                    # do not want files created edited by genotyping_rescheduler to change time/mem back to default
                    with open(file,'w') as o:
                        o.write("%s" % text)
                    shfiles.append(file)

# create scheddir queue
scheddir = op.join(parentdir,'shfiles/supervised/select_variants_within_and_across')
if not op.exists(scheddir):
    os.makedirs(scheddir)
    
for sh in shfiles:
    dst = op.join(scheddir,op.basename(sh))
    try:
        os.symlink(sh,dst)
    except:
        print('could not create symlink')
        
# submit to scheduler
# pipedir = os.popen('echo $HOME/pipeline').read().replace("\n","")
# os.system('python %s %s' % (op.join(pipedir,'genotyping_scheduler.py'),parentdir))
    

print(shdir)
