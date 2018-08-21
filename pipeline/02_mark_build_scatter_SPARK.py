###
# execution: 02b_merge_get-coverage.py /path/to/rgout /path/to/fastq.gz-folder/ /path/to/ref.fa <int passed from 01b.py>
###

###
# when using samtools from anaconda, change default_jvm_mem_opts to "-Xms512m -Xmx4g" in <which picard> file
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
rawvcf    = op.join(vcfdir,'raw_%s.vcf' % samp)
snpvcf    = op.join(vcfdir,'snps_%s.vcf' % samp)
indelvcf  = op.join(vcfdir,'indel_%s.vcf' % samp)

#get ploidy and rginfo
PLOIDY = pickle.load(open(op.join(fqdir,'ploidy.pkl')))
ploidy = int(PLOIDY[pool])
print 'ploidy =',ploidy
rginfo = pickle.load(open(op.join(fqdir,'rginfo.pkl')))

# create sh files
shfiles = []
if ploidy > 2: #non poolseq
    print "this is a poolseq file"
    # run it
    text = '''#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --nodes=4
#SBATCH --mem=200000mb
#SBATCH --cpus-per-task=32
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=mark%s
#SBATCH --export=all
#SBATCH --output=mark%s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

source $HOME/.bashrc
module load gatk/4.0.0.0
module load spark/2.2.0
export SPARK_IDENT_STRING=$SLURM_JOBID
export SPARK_WORKER_DIR=$SLURM_TMPDIR

# startup spark
start-master.sh
sleep 1
MASTER_URL=$(grep -Po '(?=spark://).*' $SPARK_LOG_DIR/spark-${SPARK_IDENT_STRING}-org.apache.spark.deploy.master*.out)
NWORKERS=$((SLURM_NTASKS - 1))
SPARK_NO_DAEMONIZE=1 srun -n ${NWORKERS} -N ${NWORKERS} --label --output=$SPARK_LOG_DIR/spark-%j-workers.out start-slave.sh -m ${SLURM_MEM_PER_NODE}M -c ${SLURM_CPUS_PER_TASK} ${MASTER_URL} &
slaves_pid=$!
echo $MASTER_URL

# remove dups
picard MarkDuplicates I=%s O=%s M=%s VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true

# Build bam index for GATK
picard BuildBamIndex I=%s

# call variants
java -jar $EBROOTGATK/gatk-package-4.0.0.0-spark.jar --java-options "-Xmx8g" -T HaplotypeCaller -sample-ploidy %s -R %s -genotyping-mode DISCOVERY -ERC GVCF -I %s -o %s -- --spark-runner SPARK --spark-master $MASTER_URL -num-executors 1 --executor-cores 32 --executor-memory 8000m 

kill $slaves_pid
stop-master.sh
    ''' % (str(tcount).zfill(3), str(tcount).zfill(3),
           rgout, dupfile,  dupstat,
           dupfile,
           ploidy,  ref,  realigned,  rawvcf,
          )
else: # poolseq
#     for i in range(1000):
    print "this is an individual's file"
    for i in range(1):
        j = i + 1 #there is no Scaffold_0
        scaff = "Scaffold_%i"  % j
        text = '''#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --nodes=1
#SBATCH --mem=200000mb
#SBATCH --cpus-per-task=32
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=mark%s
#SBATCH --export=all
#SBATCH --output=mark%s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

source $HOME/.bashrc
module load gatk/4.0.0.0
module load spark/2.2.0
export SPARK_IDENT_STRING=$SLURM_JOBID
export SPARK_WORKER_DIR=$SLURM_TMPDIR

# startup spark
start-master.sh
sleep 1
MASTER_URL=$(grep -Po '(?=spark://).*' $SPARK_LOG_DIR/spark-${SPARK_IDENT_STRING}-org.apache.spark.deploy.master*.out)
NWORKERS=$((SLURM_NTASKS - 1))
SPARK_NO_DAEMONIZE=1 srun -n ${NWORKERS} -N ${NWORKERS} --label --output=$SPARK_LOG_DIR/spark-%%j-workers.out start-slave.sh -m ${SLURM_MEM_PER_NODE}M -c ${SLURM_CPUS_PER_TASK} ${MASTER_URL} &
slaves_pid=$!
echo $MASTER_URL

# remove dups
picard MarkDuplicates I=%s O=%s M=%s VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true

# Build bam index for GATK
picard BuildBamIndex I=%s

# call variants
module load gatk/4.0.0.0
gatk HaplotypeCaller --sample-ploidy %s -R %s -genotyping-mode DISCOVERY -ERC GVCF -I %s -o %s -- --spark-runner SPARK -num-executors 1 --executor-cores 32 --executor-memory 8000m 

kill $slaves_pid
stop-master.sh
    ''' % ("%s_%s" % (scaff,str(tcount).zfill(3)),
           "%s_%s" % (scaff,str(tcount).zfill(3)),
           rgout, dupfile,  dupstat,
           dupfile,
           ploidy,  ref,  realigned,  rawvcf,
          )

print text

# for s in shfiles:
#     os.system('sbatch %s' % s)







