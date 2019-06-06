"""
# usage
# python 06_filter_concat_scaffolds.py /path/to/parentdir
# 

# purpose
# concat all genotyped files for each pool into a vcf file if all files are ready
# filter this vcf file for quality depth, fisher strand, mapping quality, \
MQRankSum, MAF < 0.05, no less than 75% missing data
#

# assumes
# intervals files used to make scatter files in 04_scatter.py will sort properly
# ie in order as they appear of the reference, which is important for concatenating vcfs
#
# GQ for samps should be filtered afterwards and then again for missing data after removing samps
# if I used "--genotypeFilterExpression GQ < 20" I could miss SNPs that have > 75% of samps with GQ >=20
"""

### imports
import sys
import os
from os import path as op
from os import listdir
import pickle
import numpy as np
from coadaptree import fs, createdirs, pklload, get_email_info
### 

### args
thisfile, parentdir = sys.argv
if parentdir.endswith("/"):
    parentdir = parentdir[:-1]
poolref = pklload(op.join(parentdir, 'poolref.pkl'))
email_info = get_email_info(parentdir, 'concat')
###

### dirs
shdir   = op.join(parentdir, 'shfiles/concat')
catdir  = op.join(parentdir, 'concatenated_vcfs')
filtdir = op.join(parentdir, 'filtered_snps')
createdirs([shdir,catdir,filtdir])
###

# get the snpfiles
snpdir = op.join(parentdir, 'snps')
snpfiles = [f.replace('.tbi', '') for f in fs(snpdir) if 'snp' in op.basename(f) and f.endswith('.tbi')]
os.system('echo "len(snpfiles) = %s"' % str(len(snpfiles)))

# sort snpfiles by pool
pools = list(poolref.keys())
combdict = {}
for i,snp in enumerate(snpfiles):
    for p in pools:
        if p in op.basename(snp):
            pool = p
            break
    if not pool in combdict:
        combdict[pool] = []
    combdict[pool].append(snp)
    del pool  # will cause script to error if pool isn't found in snpfile
for pool,poolfiles in combdict.items():
    os.system(f'echo there are {len(poolfiles)} snpfiles for {pool} pool')

# get the expected number of snpfiles per pool
expected = {}
for pool,ref in poolref.items():
    scafdir = op.join(op.dirname(ref), 'intervals')
    scaffiles = [f for f in fs(scafdir) if f.endswith('.list')]
    os.system(f'echo found {len(scaffiles)} intervalfiles for {pool}')
    expected[pool] = len(scaffiles)

# write the sh files
shfiles = []
fcats   = []
for pool,files in combdict.items():
    if len(files) == expected[pool]:
        os.system(f'echo creating sbatch file for {pool}')
        catout = op.join(catdir, f"{pool}_concatenated_snps.vcf.gz")
        filtout = op.join(filtdir, f"{pool}_filtered_concatenated_snps.vcf.gz")
        tbi = filtout.replace(".gz",".gz.tbi")
        file = op.join(shdir,"%s-concat.sh" % pool)
        if op.exists(file) is False: # if sh file hasn't been made before
            nomissing = filtout.replace(".vcf.gz", "_no-missing")
            tablefile = nomissing + "_table.txt"
            ref = poolref[pool]
            files = ' '.join(sorted(files))
            text = f'''#!/bin/bash
#SBATCH --time=11:59:59
#SBATCH --mem=50000M
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --job-name={pool}-concat
#SBATCH --output={pool}-concat_%j.out 
{email_info}

source $HOME/.bashrc
export _JAVA_OPTIONS="-Xms256m -Xmx48g"
echo $0

module load bcftools/1.9
echo "CONCATENATING FILES"
bcftools concat {files} -O z -o {catout} --threads 32
module unload bcftools

module load gatk/4.1.0.0
echo -e "\nFILTERING VARIANTS"
gatk IndexFeatureFile -F {catout}
gatk VariantFiltration -R {ref} -V {catout} -O {filtout} --filter-expression \
"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || AF < 0.05 || AF > 0.95" \
--filter-name "coadaptree_filter"
module unload gatk

module load vcftools/0.1.14
echo -e "\nFILTERING MISSING DATA"
vcftools --gzvcf {filtout} --max-missing 0.75 --recode --recode-INFO-all --out {nomissing}
module unload vcftools

module load gatk/4.1.0.0
echo -e "\nVARIANTS TO TABLE"
gatk VariantsToTable --variant {nomissing}.recode.vcf -F CHROM -F POS -F REF -F ALT -F AF -F DP -F QD \
-F FS -F MQ -F MQRankSum -F ReadPosRankSum -GF AD -GF DP -GF GQ -GF GT -GF SB -O {tablefile}

'''
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

os.chdir(shdir)
for sh in shfiles:
    os.system('echo %s' % sh)
    os.system('sbatch %s' % sh)

print('here are the fcats:')
for f in fcats:
    print(f)