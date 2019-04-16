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
def pkldump(obj,f):
    with open(f,'wb') as o:
        pickle.dump(obj,o,protocol=pickle.HIGHEST_PROTOCOL)
def pklload(path):
    pkl = pickle.load(open(path,'rb'))
    return pkl
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
print(pooldirs)
exit()
        
        
# get a list of intervals per .list file, or load if already exists
intdir = op.join(parentdir,'intervals')
intpkl = op.join(parentdir,'intervals.pkl')
if not op.exists(intpkl):
#     os.system('echo creating intpkl')
    intervals = {}
    for d in ['pooled','individual']:
        print(d)
        if not d in intervals:
            intervals[d] = {}
        subdir = op.join(intdir,d)
        for i,f in enumerate(fs(subdir)):
            if i % 500 == 0:
                update([(d,i)])
            batch = f.split("_")[-1].replace(".list","")
            if not batch in intervals:
                intervals[batch] = []
            with open(f,'r') as o:
                text = o.read().split("\n")
            intervals[d][batch] = t
    pkldump(intervals,intpkl)
else:
    os.system('echo loading intervals.pkl')
    intervals = pklload(intpkl)
    os.system('echo done loading intervals')
                
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
#             kind = 'unpooled'       ######################## CHANGE THIS TO 'individual' to match intervals dict
#             if not kind in spp[sp]: 
#                 spp[sp][kind]  = []
#             spp[sp][kind].append(f) 
print('len spp keys', len(spp.keys()),spp['DF'].keys(),)
for kind in spp['DF']:
    print(kind,len(spp['DF'][kind]))

# create some dirs
outdir = op.join(parentdir,'snps')
shdir  = op.join(parentdir,'shfiles/select_variants_within_and_across')
createdirs([outdir,shdir])

# get a list of snpfiles that have already been made
snpfiles = [f.replace(".tbi","") for f in fs(outdir) if 'snps' in op.basename(f) and f.endswith('.tbi')]
print('len(snpfiles)=',len(snpfiles))
if len(snpfiles) > 0: 
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
            batch = op.basename(f).split("scatter")[1].split(".g.v")[0]
            if batch not in groups:
                groups[batch] = []
            groups[batch].append(f)
        pools = '---'.join([x for x in sorted(pools)])
        print(pools)
        # get ref.fa
        ref = poolref[POOL] # doesn't matter which pool I use, will be same ref (POOL is from iter)
        # get commands
        for batch in groups:
#             print('len(groups[%s])= '%batch,len(groups[batch]))
            cmds = '''''' # wipe cmds to be safe
            combfile = op.join(outdir,'%s--%s_combined.vcf.gz' % (pools,batch))
            gfile    = combfile.replace("_combined.vcf.gz","_genotyped.vcf.gz")
            snpfile  = gfile.replace("_genotyped.vcf.gz","_snps.vcf.gz")
            tbifile  = snpfile.replace(".vcf.gz",".vcf.gz.tbi")
            if snpfile in snpfiles:
                continue
            if (len(groups[batch]) == 20 and kind == 'unpooled') or (len(groups[batch]) == 2 and kind == 'pooled'):
#             if len(groups[batch]) == 2 and kind == 'pooled':
                # these files need to be combined
                varcmd = '--variant ' + ' --variant '.join([x for x in sorted(groups[batch])])
                
                # first get cmds for each interval, and do combining/genotyping over individual contigs
                combfiles = []
                gfiles = []
                snpfiles = []
                for contig in intervals[kind][batch[1:]]:  # batch[1:] removes "-" at beginning of batch, leaving for now
                    newcombfile = combfile.replace(".vcf.gz","_%s.vcf.gz" % contig)
                    newgfile    = gfile.replace(".vcf.gz","_%s.vcf.gz" % contig)
                    newsnpfile  = snpfile.replace(".vcf.gz","_%s.vcf.gz" % contig)
                    gfiles.append(newgfile)       # keep track of files to delete later
                    combfiles.append(newcombfile) # keep track of files to delete later
                    snpfiles.append(newsnpfile)   # keep track of files to delete later
                    cmd   = '''echo COMBINING GVCFs
gatk CombineGVCFs -R %(ref)s %(varcmd)s -L %(contig)s -O %(newcombfile)s 

echo GENOTYPING GVCFs
gatk --java-options "-Xmx4g" GenotypeGVCFs -R %(ref)s -V %(newcombfile)s -O %(newgfile)s 

echo SELECTING VARIANTS
gatk SelectVariants -R %(ref)s -V %(newgfile)s --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O %(newsnpfile)s

if [ -f %(newsnpfile)s ]; then
    echo 'removing combfile and gfile for %(newsnpfile)s'
    rm %(newgfile)s
    rm %(newcombfile)s
fi

''' % locals()
                    cmds = cmds + cmd
                # next write cmd to combine snpfiles, delete unnecessary files
                delsnp  = ' '.join(snpfiles)
                cmd = '''module load bcftools/1.9
bcftools concat %(snpfiles)s -O z -o %(snpfile)s --threads 1

gatk IndexFeatureFile -F %(snpfile)s

# if all of the contig snpfiles have been combined, delete unneeded files
if [ -f %(tbifile)s ]; then
    echo "tbi file exists, deleting unneeded files"
    rm %(delsnp)s
fi

'''
                cmds = cmds + cmd
                    
            elif len(groups[batch]) == 1 and kind not in ['unpooled','pooled']: # single pool genotyping
                continue           ####################################################### DELETE THIS! #####################
                f = groups[batch][0] # skip symlink and just specify original file
#                 gfile    = combfile.replace("_combined.vcf.gz","_genotyped.vcf.gz") # keep this
#                 snpfile  = gfile.replace("_genotyped.vcf.gz","_snps.vcf.gz")
                cmds = '''echo GENOTYPING GVCFs
gatk --java-options "-Xmx4g" GenotypeGVCFs -R %(ref)s -V %(f)s -O %(gfile)s 

echo SELECTING VARIANTS
gatk SelectVariants -R %(ref)s -V %(gfile)s --select-type-to-include SNP --restrict-alleles-to BIALLELIC -O %(snpfile)s

''' % locals()
            else:
                print('\t',sp,kind,batch,'doesnt have enough files')
                continue
            file = op.join(shdir,'genotype---%(pools)s--%(batch)s.sh' % locals())
            if not file in alreadycreated: # this way I can start running the genotyping stage before all files are ready
#                 print(file)
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
#SBATCH --job-name=genotype---%(pools)s--%(batch)s
#SBATCH --export=all
#SBATCH --output=genotype---%(pools)s--%(batch)s---%%j.out 

source $HOME/.bashrc
cat $0
# next line necessary for rescheduler
echo shfile = %(file)s

# requeue jobs with errors
cd $HOME/gatk_pipeline
python genotyping_rescheduler.py %(parentdir)s

# fill up the queue
python genotyping_scheduler.py %(parentdir)s

# genotype current batch
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
                
                
                
                
                break ############################ FOR TESTING - DELETE LATER ########################################
                
                
                
                

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

