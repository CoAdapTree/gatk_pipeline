"""
### purpose
# map with bwa, view/sort/index with samtools
###

### usage
# 02_bwa-map_rginfo_mark_realign_lofreq_crisp.py /path/to/ref.fa /path/to/trimmedR1.fastq /path/to/trimmedR2.fastq \
#                                                                                              /sbatch/dir/ sampID
###

### assumes
# outfiles from "bwa index ref.fasta"
#
# export path to lofreq in $HOME/.bashrc
###
"""

import sys, os, balance_queue, subprocess, shutil
from os import path as op
from coadaptree import pklload, pkldump, get_email_info, makedir

# get argument inputs
thisfile, parentdir, samp = sys.argv
pool = pklload(op.join(parentdir, 'samp2pool.pkl'))[samp]
pooldir = op.join(parentdir, pool)
shdir = op.join(pooldir, 'shfiles')
ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
r1r2outs = pklload(op.join(pooldir, 'samp2_r1r2out.pkl'))[samp]

# create dirs 
bwashdir = op.join(shdir, '02_bwa_shfiles')
samdir = op.join(pooldir, '02a_samfiles')
bamdir = op.join(pooldir, '02b_bamfiles')
sortdir = op.join(pooldir, '02c_sorted_bamfiles')
for d in [bwashdir, samdir, bamdir, sortdir]:
    makedir(d)

# get rginfo - THIS CAN STAY EVEN WITH SAMPS SEQUENCED MULTIPLE TIMES - RGID and RGPU are defined with file
rginfo = pklload(op.join(parentdir, 'rginfo.pkl'))
print('pooldir = ', pooldir)
print("RG = ", rginfo[samp])
rglb = rginfo[samp]['rglb']
rgpl = rginfo[samp]['rgpl']
rgsm = rginfo[samp]['rgsm']


def getbwatext(r1out, r2out):
    # bwa: fastq -> sam
    sam = op.basename(r1out).replace("R1_trimmed.fastq.gz", "R1R2_trimmed.sam")
    samfile = op.join(samdir, sam)
    # samtools view: sam -> bam
    bam = op.basename(samfile).replace('.sam', '.bam')
    bamfile = op.join(bamdir, bam)
    # samtools sort: bamfile -> sortfile
    sort = op.basename(bamfile).replace('.bam', '_sorted.bam')
    sortfile = op.join(sortdir, sort)
    flagfile = op.join(sortdir, sort.replace('.bam', '.bam.flagstats'))
    coordfile = op.join(sortdir, sort.replace('.bam', '.bam.coord'))
    
    return (sortfile, f'''# get RGID and RGPU
RGID=$(zcat {r1out} | head -n1 | sed 's/:/_/g' | cut -d "_" -f1,2,3,4)
RGPU=$RGID.{rglb}

# map, sam to bam, sort by coordinate, index
module load bwa/0.7.17
bwa mem -t 32 -M -R "@RG\\tID:$RGID\\tSM:{rgsm}\\tPL:{rgpl}\\tLB:{rglb}\\tPU:$RGPU" \
{ref} {r1out} {r2out} > {samfile}
module unload bwa

module load samtools/1.9
samtools view -@ 32 -q 20 -F 0x0004 -f 0x0002 -Sb {samfile} > {bamfile}
samtools sort -@ 32 {bamfile} > {sortfile}
samtools index {sortfile}
samtools flagstat {sortfile} > {flagfile}
module unload samtools

module load bedtools/2.27.1
bedtools bamtobed -i {sortfile} > {coordfile}

''')


# get bwatext
bwatext = ''''''
sortfiles = []
for r1, r2 in r1r2outs:
    sortfile, text = getbwatext(r1, r2)
    bwatext = bwatext + text
    sortfiles.append(sortfile)
pkldump(sortfiles, op.join(pooldir, '%s_sortfiles.pkl' % samp))

# send it off
email_text = get_email_info(parentdir, '02')
text = f'''#!/bin/bash
#SBATCH --time=23:59:00
#SBATCH --mem=35000M
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --job-name={pool}-{samp}-bwa
#SBATCH --output={pool}-{samp}-bwa_%j.out 
{email_text}

{bwatext}

# mark and build
source $HOME/.bashrc
export PYTHONPATH="${{PYTHONPATH}}:$HOME/pipeline"
export SQUEUE_FORMAT="%.8i %.8u %.12a %.68j %.3t %16S %.10L %.5D %.4C %.6b %.7m %N (%r)"
python $HOME/pipeline/03_mark_build.py {pooldir} {samp}
'''

# create shfile
qsubfile = op.join(bwashdir, f'{pool}-{samp}-bwa.sh')
with open(qsubfile, 'w') as o:
    o.write("%s" % text)

# sbatch file
os.chdir(bwashdir)
print('shdir = ', shdir)
subprocess.call([shutil.which('sbatch'), qsubfile])

balance_queue.main('balance_queue.py', 'bwa')
balance_queue.main('balance_queue.py', 'trim')
