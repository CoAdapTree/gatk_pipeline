"""
# usage
# python 06_filter_concat_scaffolds.py /path/to/parentdir
# 

# purpose
# concat all genotyped files for each pool into a vcf file if all files are ready
# filter this vcf file for:
# - quality depth < 2,
# - fisher strand > 60,
# - mapping quality < 40,
# - MQRankSum < -12.5,
# - MAF < 0.05 (default - user can define when starting pipeline with 00_start)
# - GQ < 20 for individual samps,
# - no less than 75% missing data
#

# assumes
# intervals files used to make scatter files in 04_scatter.py will sort properly
# ie in order as they appear of the reference, which is important for concatenating vcfs
#
"""

### imports
import sys, os, pickle, subprocess
from os import path as op
import numpy as np
from coadaptree import fs, createdirs, pklload, get_email_info
### 

### args
thisfile, parentdir = sys.argv
if parentdir.endswith("/"):
    parentdir = parentdir[:-1]
poolref = pklload(op.join(parentdir, 'poolref.pkl'))
email_info = get_email_info(parentdir, 'concat')
bash_variables = op.join(parentdir, 'bash_variables')
maf = pklload(op.join(parentdir, 'maf.pkl'))
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
    if pool not in combdict:
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
    ref = poolref[pool]
    if len(files) == expected[pool]:
        os.system(f'echo creating sbatch file for {pool}')
        catout = op.join(catdir, f"{pool}_concatenated.vcf.gz")
        snpsout = op.join(catdir, f"{pool}_concatenated_snps.vcf.gz")
        filtout = op.join(filtdir, f"{pool}_filtered_concatenated_snps.vcf.gz")
        tbi = filtout.replace(".gz",".gz.tbi")
        file = op.join(shdir,"%s-concat.sh" % pool)
        if op.exists(file) is False: # if sh file hasn't been made before
            maxmissing = filtout.replace(".vcf.gz", "_max-missing")
            tablefile = maxmissing + "_table.txt"
            tablefile_filtered = tablefile.replace(".txt", "_biallelic-only.txt")
            ref = poolref[pool]
            files = ' '.join(sorted(files))
            text = f'''#!/bin/bash
#SBATCH --time=11:59:59
#SBATCH --mem=50000M
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name={pool}-concat
#SBATCH --output={pool}-concat_%j.out 
{email_info}

source {bash_variables}
export _JAVA_OPTIONS="-Xms256m -Xmx48g"
echo $0

module load bcftools/1.9
echo "CONCATENATING FILES"
bcftools concat {files} -O z -o {catout} --threads 32
module unload bcftools

module load gatk/4.1.0.0
echo -e "\nFILTERING VARIANTS"
gatk IndexFeatureFile -F {catout}
gatk SelectVariants -R {ref} -V {catout} --select-type-to-include SNP -O {snpsout}
gatk VariantFiltration -R {ref} -V {snpsout} -O {filtout} --filter-expression \
"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5" \
--filter-name "coadaptree_filter"
module unload gatk

module load vcftools/0.1.14
echo -e "\nFILTERING MISSING DATA"
vcftools --gzvcf {filtout} --maf {maf} --minGQ 20 --max-missing 0.75 --recode \
--recode-INFO-all --out {maxmissing}
module unload vcftools

module load gatk/4.1.0.0
echo -e "\nVARIANTS TO TABLE"
gatk VariantsToTable --variant {maxmissing}.recode.vcf -F CHROM -F POS -F REF -F ALT -F AF -F DP -F QD \
-F FS -F MQ -F MQRankSum -F ReadPosRankSum -F TYPE -F FILTER -GF AD -GF DP -GF GQ -GF GT -GF SB -O {tablefile} --split-multi-allelic
module unload gatk

echo -e "\nREMOVING MULTIALLELIC, KEEPING noREF SNPs WITH TWO ALT ALLELES"
python $HOME/gatk_pipeline/remove_multiallelic-keep_noREF.py {tablefile} {tablefile_filtered} 

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
    
# balance queue
balance_queue = op.join(os.environ['HOME'], 'gatk_pipeline/balance_queue.py')
subprocess.call([sys.executable, balance_queue, 'genotype', parentdir])
subprocess.call([sys.executable, balance_queue, 'concat', parentdir])
