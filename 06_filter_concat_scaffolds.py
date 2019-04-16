### usage
# python 04_concat_scaffolds.py /path/to/parentdir/used/in/00_start-pipeline/command 
### 

### purpose
# instead of calling snps when combining some of our pools, call snps for indiviual pools
###

### FIX
# I've used 'kit' as a key word that is in all of the library/pool names, this should be changed if not the case
# see creation of combdict
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
poolref = pickle.load(open(op.join(parentdir,'poolref.pkl'),'rb'))
###

### dirs
shdir   = op.join(parentdir,'shfiles/concat')
catdir  = op.join(parentdir,'concatenated_vcfs')
filtdir = op.join(parentdir,'filtered_snps')
createdirs([shdir,catdir,filtdir])
###

# get the snpfiles
snpdir = op.join(parentdir,'snps')
allfiles = fs(snpdir)
snpfiles = [f for f in fs(snpdir) if f.endswith('.gz') and 'snp' in op.basename(f) and f.replace('.gz','.gz.tbi') in allfiles]
os.system('echo "len(snpfiles) = %s"' % str(len(snpfiles)))

# sort snpfiles by combo lib
combdict = {}
for i,snp in enumerate(snpfiles):
    lib = "---".join([x for x in op.basename(snp).split("-") if 'kit' in x])
    if not lib in combdict:
        combdict[lib] = []
    combdict[lib].append(snp)
os.system('echo there are %s keys in combdict' % str(len(combdict.keys())))
for k,kfiles in combdict.items():
    os.system('echo %s %s' % (k,str(len(kfiles))))
# write the sh files
shfiles = []
fcats   = []
for lib in combdict.keys():
    if len(combdict[lib]) in [50289]: # for now this is for Doug fir, will need to change for other spp
        os.system('echo creating sbatch file')
        catout   = op.join(catdir,"%s_concatenated_snps.vcf.gz" % lib)
        filtout  = op.join(filtdir,"%s_filtered_concatenated_snps.vcf.gz" % lib)
        firstlib = lib.split("---")[0]
        ref      = poolref[firstlib]
        # I should have made scaffols be zfill(4) not zfill(3)
        files    = " ".join([snp for snp in sorted(combdict[lib]) if '1000' not in snp])  # I can remove the if-condition once we're done with the pseudo ref
        # (bcftools needs the input files to be sorted)
        files    = files + ' %s ' % [snp for snp in combdict[lib] if '1000' in snp][0] 
        text = '''#!/bin/bash
#SBATCH --time=11:59:59
#SBATCH --mem=15000M
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --job-name=%(lib)s-concat
#SBATCH --output=%(lib)s-concat_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

module load bcftools/1.9

bcftools concat %(files)s -O z -o %(catout)s --threads 32

module load gatk/4.0.8.1

gatk IndexFeatureFile -F %(catout)s

gatk VariantFiltration -R %(ref)s -V %(catout)s -O %(filtout)s --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5" --filter-name "coadaptree_filter"

''' % locals()
        file = op.join(shdir,"%s-concat.sh" % lib)
        if not op.exists(filtout): # so I can run this when files begin to finish from 03_genotype.py
            if not op.exists(filtout.replace(".gz",".gz.tbi")):
                with open(file,'w') as o:
                    o.write("%s" % text)
                shfiles.append(file)
                # some times it's easier to just salloc some resources and run, instead of sbatching 
                with open(file,'r') as o:
                    text = o.read().split("\n")
                lines = []
                for line in text:
                    if ('gatk' in line or 'bcftools' in line) and 'module' not in line:
                        lines.append(line)
                fcat = file.replace(".sh","_catfile.sh")
                fcats.append(fcat)
                with open(fcat,'w') as o:
                    for line in lines:
                        o.write("%s\n" % line)

# os.chdir(shdir)
# for sh in shfiles:
#     os.system('echo %s' % sh)
#     os.system('sbatch %s' % sh)
[print(sh) for sh in shfiles]
[print(f) for f in fcats]
        
