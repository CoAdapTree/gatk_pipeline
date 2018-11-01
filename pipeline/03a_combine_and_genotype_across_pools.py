### usage
# python 03a_combine_and_genotype_by_pool.py /path/to/parentdir/used/in/00_start-pipeline/command 
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
rginfo    = pickle.load(open(op.join(parentdir,'rginfo.pkl'),'rb'))
poolsamps = pickle.load(open(op.join(parentdir,'poolsamps.pkl'),'rb')) #key=pool val=list(sampnames)
poolref   = pickle.load(open(op.join(parentdir,'poolref.pkl'),'rb'))   #key=pool val=/path/to/ref.fa
ploidy    = pickle.load(open(op.join(parentdir,'ploidy.pkl'),'rb'))    #not used yet, but maybe if pooled needs more time
pools = uni(list(samp2pool.values()))
samps = list(rginfo.keys())
###

# get a list of subdirectory pool dirs created earlier in pipeline
pooldirs = []
for p in pools:
    pooldir = op.join(parentdir,p)
    assert op.exists(pooldir)
    pooldirs.append(pooldir)

# get a list of files and partition by those that need to be combined
spp = {}
for d in pooldirs:
    vcfdir = op.join(d,'vcfs')
    files = [f for f in fs(vcfdir) if f.endswith('.gz')]
    count = 0
    for f in files:
        count += 1
        splits = op.basename(f).replace("raw_","").split("_")
        sp = splits[0]
        if not sp in spp:
            spp[sp] = {}
        if splits[1].startswith('p'):
            kind = 'pooled'
        else:
            kind = 'unpooled'
        if not kind in spp[sp]:
            spp[sp][kind] = []
        spp[sp][kind].append(f)

# create some dirs
outdir = op.join(parentdir,'snps')
shdir  = op.join(parentdir,'shfiles/select_variants')
createdirs([outdir,shdir])

# make sh files
shfiles = []
for sp in spp:
    print(sp)
    for kind in spp[sp]:
        if kind == 'pooled':
            TIME = '23:59:59'
        else:
            TIME = '11:59:59'
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
        pool = op.basename(op.dirname(op.dirname(f)))
        ref = poolref[pool]
        # get commands
        for scaff in groups:
            cmds = ''''''
            combfile = op.join(outdir,'%s--%s_combined.vcf.gz' % (pools,scaff))
            gfile    = combfile.replace("_combined.vcf.gz","_genotyped.vcf.gz")
            snpfile  = gfile.replace("_genotyped.vcf.gz","_snps.vcf.gz")
            if len(groups[scaff]) > 1:
                # these files need to be combined
                varcmd = '--variant ' + ' --variant '.join([x for x in sorted(groups[scaff])])
                cmd    = '''gatk CombineGVCFs -R %(ref)s %(varcmd)s -O %(combfile)s
                
gatk --java-options "-Xmx4g" GenotypeGVCFs -R %(ref)s -V %(combfile)s -O %(gfile)s

gatk SelectVariants -R %(ref)s -V %(gfile)s --select-type-to-include SNP -O %(snpfile)s

''' % locals()
            elif len(groups[scaff]) == 1:
                # no need to combine, just symlink vcf and tbi file
                f = groups[scaff][0]
                combfile = op.join(outdir,'%s--%s_combined.vcf.gz' % (pools,scaff))
                tbi = f.replace(".gz",".gz.tbi")
#                 tbilink = op.join(outdir,op.basename(tbi))
                tbilink = combfile.replace(".gz",".gz.tbi")
                try:
                    if not op.exists(combfile):
                        os.symlink(f,combfile)
                    if not op.exists(tbilink):
                        os.symlink(tbi,tbilink)
                except:
                    print("could not create symlink\n%s\n%s\n%s\n%s\n" % (f,combfile,tbi,tbilink)) # in case the tbi file isn't there
                    continue
                gfile    = combfile.replace("_combined.vcf.gz","_genotyped.vcf.gz")
                snpfile  = gfile.replace("_genotyped.vcf.gz","_snps.vcf.gz")
                cmd = '''gatk --java-options "-Xmx4g" GenotypeGVCFs -R %(ref)s -V %(combfile)s -O %(gfile)s
                
gatk SelectVariants -R %(ref)s -V %(gfile)s --select-type-to-include SNP -O %(snpfile)s

''' % locals()
            else:
                print(sp,kind,scaff,'has no files')
            cmds = cmds + cmd
            if not cmds == '''''':
                text = '''#!/bin/bash
#SBATCH --time=%(TIME)s
#SBATCH --nodes=1
#SBATCH --mem=30000M
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=%(pools)s--%(scaff)s
#SBATCH --export=all
#SBATCH --output=%(pools)s--%(scaff)s---%%j.out 

source $HOME/.bashrc
module load gatk/4.0.8.1

%(cmds)s

''' % locals()
                filE = op.join(shdir,'%(pools)s--%(scaff)s.sh' % locals())
                print('\t','\t',filE)
                with open(filE,'w') as o:
                    o.write("%s" % text)
                shfiles.append(filE)
print(shdir)