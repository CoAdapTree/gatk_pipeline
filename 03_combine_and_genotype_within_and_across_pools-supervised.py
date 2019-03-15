### usage
# python 03a_combine_and_genotype_by_pool.py /path/to/parentdir-used-in-00_start-pipeline.py-command 
### 

### imports
import sys
import os
from os import path as op
from os import listdir
import pickle
from collections import Counter
def uni(mylist):
    return list(set(mylist))
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

# get a list of files and partition by those that need to be combined (multi-pool/indSeq genotyping)
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
            kind  = 'pooled'          ######################## CHANGE THIS BACK FOR MULTI-POOL GENOTYPING!!!!!!!!!!!!!!!!
            if not kind in spp[sp]:   ######################## CHANGE THIS BACK FOR MULTI-POOL GENOTYPING!!!!!!!!!!!!!!!!
                spp[sp][kind]  = []   ######################## CHANGE THIS BACK FOR MULTI-POOL GENOTYPING!!!!!!!!!!!!!!!!
            spp[sp][kind].append(f)   ######################## CHANGE THIS BACK FOR MULTI-POOL GENOTYPING!!!!!!!!!!!!!!!!
#             kind2 = op.basename(d) # for single pool genotyping
#             if not kind2 in spp[sp]:
#                 spp[sp][kind2] = []
#             spp[sp][kind2].append(f)
#         else:
#             kind = 'unpooled'
#             if not kind in spp[sp]: 
#                 spp[sp][kind]  = []
#             spp[sp][kind].append(f) 
for kind in spp['DF']:
    print(kind,len(spp['DF'][kind]))

# create some dirs
outdir = op.join(parentdir,'snps')
shdir  = op.join(parentdir,'shfiles/select_variants_within_and_across')
createdirs([outdir,shdir])

# get a list of snpfiles that have already been made
snpfiles = [f.replace(".tbi","") for f in fs(outdir) if 'snps' in op.basename(f) and f.endswith('.tbi')]
print('len(snpfiles)=',len(snpfiles))
print(snpfiles[0])


# make sh files
alreadycreated = [f for f in fs(shdir) if f.endswith('.sh') and 'swp' not in f]
print('len(alreadycreated)=',len(alreadycreated))
shfiles = []
newfiles = Counter()
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
            scaff = op.basename(f).split("scatter")[1].split(".g.v")[0]
            if scaff not in groups:
                groups[scaff] = []
            groups[scaff].append(f)
        pools = '---'.join([x for x in sorted(pools)])
        # get ref.fa
        ref = poolref[POOL] # doesn't matter which pool I use, will be same ref (POOL is from iter)
        # get commands
        for scaff in groups:
#             print('len(groups[%s])= '%scaff,len(groups[scaff]))
            cmds = '''''' # I think this might need to be moved about the for loop? it might be working bc ...?
            combfile = op.join(outdir,'%s--%s_combined.vcf.gz' % (pools,scaff))
            gfile    = combfile.replace("_combined.vcf.gz","_genotyped.vcf.gz")
            snpfile  = gfile.replace("_genotyped.vcf.gz","_snps.vcf.gz")
            if snpfile in snpfiles:
                continue
#             if (len(groups[scaff]) == 20 and kind == 'unpooled') or (len(groups[scaff]) == 2 and kind == 'pooled'):
            if len(groups[scaff]) == 2 and kind == 'pooled':
                # these files need to be combined
                varcmd = '--variant ' + ' --variant '.join([x for x in sorted(groups[scaff])])
                cmds   = '''echo COMBINING GVCFs
gatk CombineGVCFs -R %(ref)s %(varcmd)s -O %(combfile)s 

echo GENOTYPING GVCFs
gatk --java-options "-Xmx4g" GenotypeGVCFs -R %(ref)s -V %(combfile)s -O %(gfile)s 

echo SELECTING VARIANTS
gatk SelectVariants -R %(ref)s -V %(gfile)s --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O %(snpfile)s

''' % locals()
            elif len(groups[scaff]) == 1 and kind not in ['unpooled','pooled']: # otherwise if jobs aren't finished, cross-pool sh will skip combining
                # symlink vcf and tbi file
                continue           ####################################################### DELETE THIS! #####################
                f = groups[scaff][0] # skip symlink and just specify original file
#                 gfile    = combfile.replace("_combined.vcf.gz","_genotyped.vcf.gz") # keep this
#                 snpfile  = gfile.replace("_genotyped.vcf.gz","_snps.vcf.gz")
                cmds = '''echo GENOTYPING GVCFs
gatk --java-options "-Xmx4g" GenotypeGVCFs -R %(ref)s -V %(f)s -O %(gfile)s 

echo SELECTING VARIANTS
gatk SelectVariants -R %(ref)s -V %(gfile)s --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O %(snpfile)s

''' % locals()
            else:
                print('\t',sp,kind,scaff,'doesnt have enough files')
                continue
            file = op.join(shdir,'genotype---%(pools)s--%(scaff)s.sh' % locals())
            if not file in alreadycreated: # this way I can start running the genotyping stage before all files are ready
                newfiles[kind] += 1
                if kind == 'unpooled':
                    mem = '15000M'
                    time = '7-0:0:0'
                else:
                    mem = '12000M'
                    time = '11:59:00'
#                     print('\t','\t',file)
                text = '''#!/bin/bash
#SBATCH --time=%(time)s
#SBATCH --ntasks=1
#SBATCH --mem=%(mem)s
#SBATCH --cpus-per-task=1
#SBATCH --job-name=genotype---%(pools)s--%(scaff)s
#SBATCH --export=all
#SBATCH --output=genotype---%(pools)s--%(scaff)s---%%j.out 

source $HOME/.bashrc
cat $0
# next line necessary for rescheduler
echo shfile = %(file)s

# requeue jobs with errors
cd $HOME/gatk_pipeline
python genotyping_rescheduler.py %(parentdir)s

# fill up the queue
python genotyping_scheduler.py %(parentdir)s

# genotype current file
module load gatk/4.0.8.1
%(cmds)s

# keep running jobs until time runs out
echo getting help from genotyping_helper
cd $HOME/gatk_pipeline
python genotyping_helper.py %(parentdir)s %(snpfile)s

# # in case there's time, schedule the next batch
# source $HOME/.bashrc # rescue python env from from evil module
# python 03_combine_and_genotype_within_and_across_pools-supervised.py %(parentdir)s

''' % locals()
                # do not want files created edited by genotyping_rescheduler to change time/mem back to default
                with open(file,'w') as o:
                    o.write("%s" % text)
                shfiles.append(file)

# create scheddir queue
scheddir = op.join(parentdir,'shfiles/supervised/select_variants_within_and_across')
if not op.exists(scheddir):
    os.makedirs(scheddir)
    
for kind in newfiles:
    print(kind,' created ',newfiles[kind],' files')
    
for sh in shfiles:
    dst = op.join(scheddir,op.basename(sh))
    try:
        os.symlink(sh,dst)
    except:
        print('could not create symlink')
        
# # submit to scheduler, balance accounts
# pipedir = os.popen('echo $HOME/gatk_pipeline').read().replace("\n","")
# os.system('python %s %s' % (op.join(pipedir,'genotyping_scheduler.py'),parentdir)) 

print(shdir, len( ls(shdir) ) )
