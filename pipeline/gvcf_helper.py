### FIX
# at some point, have the jobID/mem as an input arg instead of scripted
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

# os.system('source $HOME/.bashrc')
DIR = op.join(op.dirname(fqdir),'shfiles/gvcf_shfiles')
os.chdir(DIR)
workingdir = op.join(DIR,'workingdir')
if not op.exists(workingdir):
    os.makedirs(workingdir)
    
# get job info and current memory/time limits
jobid    = os.popen('echo ${SLURM_JOB_ID}').read().replace("\n","")
#print('jobid=',jobid)
jobinfo  = os.popen("sacct -j %s | grep 'lindb'" % jobid).read()
#print('jobinfo=',jobinfo)
jobmem   = int([x for x in jobinfo.split() if 'mem' in x][0].split(",")[1].split('=')[1].replace("M",""))
#print('jobmem=',jobmem)
timeinfo = os.popen("sacct -j %s --format Timelimit" % jobid).read()
#print('timeinfo=',timeinfo)
jobtime  = int(timeinfo.split()[-1].split(':')[0])
#print('jobtime=',jobtime)

# get list of remaining gatk calls
shfiles = [f for f in fs(DIR) if f.endswith('.sh')]
shuffle(shfiles) 

# run commands until I run out of time
os.system('echo running gvcf_helper.py')
if len(shfiles) > 0:
    for s in shfiles:
    #     print (s)
        reservation = op.join(workingdir,op.basename(s))
        if op.exists(s):
            try:
                shutil.move(s,reservation) # so that other jobs don't rewrite
            except:
                os.system('echo could not move shfile %s' % s)
                os.system('echo to reservation %s' % reservation)
                continue
            os.system('echo %s' % reservation)
            
            with open(reservation,'r') as O:
                o = O.readlines()

            # only continue to run jobs that fit in the same memory allocation (dont waste resources if its going to fail)
            mem = int([x for x in o if 'mem' in x][0].split("=")[1].replace("M\n",""))
            if mem > jobmem:
                os.system('echo file exceeds mem limit')
                shutil.move(reservation,s) # put the job back in the queue
                continue
            # only continue to run jobs that might fit in same time allocation
            TIME = int([x for x in o if 'time' in x][0].split("=")[1].split(':')[0])
            if TIME > jobtime:
                os.system('echo file exceeds necessary time')
                shutil.move(reservation,s)
                continue


            os.system('echo file is ok to proceed')
            for line in o:
                if line.startswith('gatk'):
                    cmd = line.replace('\n','')
                    os.system('echo running cmd:')
                    os.system('echo %s' % cmd)
                    os.system('%s' % cmd)
                    try:
                        os.unlink(reservation)
                        os.system('echo unlinked shfile %s' % reservation)
                    except:
                        os.system('echo unable to unlink %s' % reservation)
                        pass
                    pipedir = os.popen('echo $HOME/pipeline').read().replace("\n","")
                    os.system('python %s %s' % (op.join(pipedir,'rescheduler.py'),
                                                fqdir))
                    os.system('python %s %s' % (op.join(pipedir,'scheduler.py'),
                                                fqdir))

                    break
else:
    os.system('echo no files to help')
    
    
    
    