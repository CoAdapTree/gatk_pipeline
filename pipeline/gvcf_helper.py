### FIX
# at some point, have the jobID/mem as an input arg instead of script
###

###
# usage: gvcf_helper.py /path/to/fastq.gz/folder/

###
# purpose: to keep running gatk commands until time or memory runs out
###

### imports
import sys
import os
from os import path as op
from os import listdir
import pickle
import shutil
from random import shuffle
def ls(DIR):
    return sorted([f for f in listdir(DIR)])
def fs (DIR):
    return sorted([op.join(DIR,f) for f in ls(DIR)])
###

### args
thisfile, fqdir = sys.argv
###

os.system('source $HOME/.bashrc')
DIR = op.join(op.dirname(fqdir),'shfiles/gvcf_shfiles')
os.chdir(DIR)
workingdir = op.join(DIR,'workingdir')
if not op.exists(workingdir):
    os.makedirs(workingdir)
    
# get job info and current memory/time limits
jobid    = os.popen('echo ${SLURM_JOB_ID}').read().replace("\n","")
jobinfo  = os.popen("sacct -j %s | grep 'lindb'" % jobid).read()
jobmem   = int([x for x in jobinfo.split() if 'mem' in x][0].split(",")[1].split('=')[1].replace("M",""))
timeinfo = os.popen("sacct -j %s --format Timelimit" % jobid).read()
jobtime  = int(timeinfo.split()[-1].split(':')[0])

# get list of remaining gatk calls
shfiles = [f for f in fs(DIR) if f.endswith('.sh')]
shuffle(shfiles) 

# run commands until I run out of time
print('running gvcf_helper.py')
for s in shfiles:
#     print (s)
    if op.exists(s):
        reservation = op.join(workingdir,op.basename(s))
        try:
            shutil.move(s,reservation) # so that other jobs don't rewrite
        except:
            print('could not move shfile %s' % s)
            print('to reservation %s' % reservation)
            continue
        print(reservation)
        with open(reservation,'r') as O:
            o = O.readlines()
#         o = open(reservation,'r').readlines()
        
        # only continue to run jobs that fit in the same memory allocation (dont waste resources if its going to fail)
        mem = int([x for x in o if 'mem' in x][0].split("=")[1].replace("M\n",""))
        if mem > jobmem:
            print ('file exceeds mem limit')
            shutil.move(reservation,s) # put the job back in the queue
            continue
        # only continue to run jobs that might fit in same time allocation
        TIME = int([x for x in o if 'time' in x][0].split("=")[1].split(':')[0])
        if TIME > jobtime:
            print ('file exceeds necessary time')
            shutil.move(reservation,s)
            continue
        
        
        print ('file is ok to proceed')
        for line in o:
            if line.startswith('gatk'):
                cmd = line.replace('\n','')
                print('running cmd:')
                print(cmd)
                os.system('%s' % cmd)
                try:
                    os.system('unlink %s' % reservation)
                    print('unlinking shfile %s' % reservation)
                except OSError as e:
                    print('unable to unlink %s' % reservation)
                    pass
                pipedir = os.popen('echo $HOME/pipeline').read().replace("\n","")
                os.system('python %s %s' % (op.join(pipedir,'scheduler.py'),
                                            fqdir))
                os.system('python %s %s' % (op.join(pipedir,'rescheduler.py'),
                                            fqdir))
                
                break


