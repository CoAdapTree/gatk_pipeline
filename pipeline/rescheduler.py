###
# usage rescheduler.py /path/to/fq.gzfiles/folder/
###

### imports
import sys
import os
from os import path as op
from os import listdir
import shutil
import pickle
def ls(DIR):
    return sorted([f for f in listdir(DIR)])
def fs (DIR):
    return sorted([op.join(DIR,f) for f in ls(DIR)])
###

### args
thisfile, fqdir = sys.argv
###

if fqdir.endswith("/"):
    fqdir = fqdir[:-1]
dname = op.dirname(fqdir)
print('dname=',dname)
#print dirname
DIR = op.join(dname,'shfiles/gvcf_shfiles') # real deal
#DIR = op.join(op.dname(fqdir),'shfiles/') # practice
#print DIR
os.chdir(DIR)
outs = [f for f in fs(DIR) if f.endswith('out') and not 'checked' in f]
rescheduler = op.join(DIR,'rescheduler.txt')
samp2pool = pickle.load(open(op.join(dname,'samp2pool.pkl'),'rb'))
def vcf2sh(v):
    pooldir  = op.dirname(op.dirname(v))
    gvcfdir  = op.join(pooldir,'shfiles/gvcf_shfiles')
    assert op.exists(gvcfdir)
    bname    = op.basename(v)
    shname   = bname.replace("raw_","").replace(".g.vcf.gz",".sh")
    shfile   = op.join(gvcfdir,shname)
    return shfile
def unlink(linkname):
    try:
        os.unlink(linkname)
        print('unlinked %s' % linkname)
    except:
        print('no symlink to unlink: %s' % linkname)
        pass
def addlink(args):
    trushfile,linkname = args
    print('symlink from: %s' % linkname)
    print('to: %s' % trushfile)
    if not op.exists(linkname):
        os.symlink(trushfile,linkname)
        print('added symlink to queue: %s' % linkname)
    else:
        print('unable to create symlink from %s to %s' % (linkname,trushfile))     


# identify outs that aren't running
sq = os.popen("squeue -u lindb | grep 'R 2'").read().split("\n")
#print(sq)
pids = []
for s in sq:
    if not s == '':
    #     print s
        pid = s.split()[0]
        pids.append(pid)
#print(pids)
runs = []
for out in outs:
    pid = op.basename(out).split(".out")[0].split("_")[1]
    if pid not in pids:
        runs.append(out)
outs = runs

print('running rescheduler.py')
createdrescheduler = False
if len(outs) > 0:
#     print('outs =',outs)
    if not op.exists(rescheduler):
        # reserve the rescheduler
        with open(rescheduler,'w') as o:
            o.write("rescheduler")
        createdrescheduler = True

        # look for errors in outfiles and resubmit the error-causing shfile using more mem or time
        for out in outs:
            #check again to see if job is running
            sq = os.popen("squeue -u lindb | grep 'R 2'").read().split("\n")
            pids = []
            for s in sq:
                if not s == '':
                #     print s
                    pid = s.split()[0]
                    pids.append(pid)
            pid = op.basename(out).split(".out")[0].split("_")[1]
            if pid in pids:
                continue
                        
            print('\nworking on %s' % out)
            with open(out,'r') as OUT:
                o = OUT.readlines()
            # look for mem error
            edited = False
            timelimit  = False
            founderror = False
            cancelled = False
            for line in o[-20:]: # look for an error message
                if 'oom-kill' in line or 'error' in line:
                    print ('found an error')
                    founderror = True
                    break
            if founderror == True:
                for test in o[-20:]: # look for a time error 
                    if 'time limit' in test.lower() or 'cancelled' in test.lower():
                        timelimit = True
                        if 'cancelled' in test.lower():
                            cancelled = True
                        break
                if timelimit == True:
                    # look for time error
                    edited = True
                    # time error could be caused by the original sh command running out of time or ...
                    ## ... gvcf_helper.py ran out of time
                    print('due to time limit')
                    # determine if its original or gvcf_helper
                    helped = False
                    for line in o[::-1]:
                        if 'getting help from gvcf_helper' in line and 'echo' not in line:
                            # we dont necessarily need to increase time if last call was due to gvcf_helper.py
                            helped = True
                            print('helped by gvcf_helper =',helped)
                            break
                    if helped == True or cancelled == True: # if the job ended on a call from gvcf_helper.py or was cancelled
                        # no need to change time this first time
                        print('leaving orginal time as-is')
                        print('cancelled =',cancelled)
                        for line in o[::-1]:
                            if line.startswith('gatk HaplotypeCaller'):
                                vcf = line.split()[-5]
                                trushfile = vcf2sh(vcf)
                                linkname = op.join(DIR,op.basename(trushfile))
                                
                                # add job back to the queue 
                                addlink((trushfile,linkname))

                                break
                    else:
                        for line in o[::-1]:
                            # this was the call from the original sh file
                            if line.startswith('gatk HaplotypeCaller'):
                                print ('adjusting time of original sh file')
                                vcf = line.split()[-5]
                                trushfile = vcf2sh(vcf)

                                print('linked to %s' % trushfile)
                                sh = open(trushfile).read()
                                if '00:00:05' in sh: # for debugging/testing
                                    text = sh.replace('00:00:05','02:59:00')
                                elif '02:59:00' in sh:
                                    text = sh.replace('02:59:00','11:59:00')
                                    print('extending time to 11:59:00')
                                elif '11:59:00' in sh:
                                    text = sh.replace('11:59:00','23:59:00')
                                    print('extending time to 23:59:00')
                                elif '23:59:00' in sh:
                                    text = sh.replace('23:59:00','7-00:00:00')
                                    print('extending time to 7 days')
                                with open(trushfile,'w') as o:
                                    o.write("%s" % text)

                                # add job back to the queue  
                                linkname = op.join(DIR,op.basename(trushfile))
                                addlink((trushfile,linkname))
                                
                                break
                else: # there's a mem oerror
                    edited = True
                    # at t=0, all sh files have mem==8000M, so if gvcf_helper.py caused mem error, the last call needs more mem
                    ## and I dont have to figure out if its the original's call or gvcf_helper's call
                    print('due to mem limit')
                    # find the last job and resubmit with more mem
                    for line in o[::-1]:
                        if line.startswith('gatk HaplotypeCaller'):
                            print('adjusting mem limit of: %s' % trushfile)
                            vcf = line.split()[-5]
                            trushfile = vcf2sh(vcf)
                            print('linked to %s' % trushfile)
                            sh = open(trushfile).read()
                            if '8000M' in sh:
                                text = sh.replace('8000M','20000M')
                                print('increasing mem to 20G')
                            elif '20000M' in sh:
                                text = sh.replace('20000M','30000M')
                                print('increasing mem to 30G')
                            elif '30000M' in sh:
                                text = sh.replace('30000M','50000M')
                                print('increasing mem to 50G')
                            elif '50000M' in sh:
                                text = sh.replace('50000M','100000M')
                                print('increasing mem to 100G')
                            elif '100000M' in sh:
                                text = sh.replace('100000M','120000M')
                                print('increasing mem to 120G')
                            with open(trushfile,'w') as o:
                                o.write("%s" % text)
                            
                            # add job back to the queue  
                            linkname = op.join(DIR,op.basename(trushfile))
                            addlink((trushfile,linkname))

                            break
            else:
                print('no mem or time errors found in %s' % out)

            # move the .out file to a new name so rescheduler doesn't look at it again
            dst = out.replace(".out","_checked.out")
            try:
                shutil.move(out,dst)
                print('moved file: %s' % dst)
            except OSError as e:
                print('could not move outfile to noerror.out: %s' % out)
                pass
    else:
        print('rescheduler was running')
else:
    print('rescheduler found no outfiles to analyze or all outfiles are for jobs currently running')

if createdrescheduler == True:
    try:
        os.remove(rescheduler)
        print('removed rescheduler')
    except:
        print('could not remove rescheduler')
        pass  