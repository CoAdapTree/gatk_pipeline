### FIX
# at some point, have the jobID/mem as an input arg instead of scripted
###

###
# usage: gvcf_helper.py /path/to/fastq.gz/folder/ /path/to/last/.tbi/file.tbi

###
# purpose: to keep running gatk commands until time or memory runs out
###

### imports
import sys
import os
import signal
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
thisfile, fqdir, tbi = sys.argv
###


# kill the job if the previous command ran out of mem
def checktbi(tbi):
    if not op.exists(tbi):
        os.system('echo previous tbi did not finish: %s' % tbi) 
        os.kill(os.getppid(), signal.SIGHUP)
checktbi(tbi)
        
os.system('source $HOME/.bashrc')
DIR = op.join(op.dirname(fqdir),'shfiles/gvcf_shfiles')
os.chdir(DIR)
workingdir = op.join(DIR,'workingdir')
if not op.exists(workingdir):
    os.makedirs(workingdir)
    
# get job info and current memory/time limits
jobid    = os.popen('echo ${SLURM_JOB_ID}').read().replace("\n","")
if not float(jobid) == int(jobid):
    # if I can't retrieve the jobid, just let the job die
    # ensures that there wasn't an error to the slurm request (some times I get timeouts etc as returns)
    exit()
os.system('echo jobid = %s' % jobid)
f = [f for f in fs(DIR) if str(jobid) in f][0]
with open(f,'r') as o:
    text = o.read().split("\n")
for line in text:
    if 'mem-' in line:
        jobmem = int(line.split("=")[-1][:-1])
    if 'time=' in line:
        splits = line.split("=")[-1].split("-")
        if len(splits) > 1:            
            # multi-day job
            jobtime = int(splits[0])*24
        else:
            # xhour job
            jobtime = int(line.split("=")[-1].split(':')[0] )    
#jobinfo  = os.popen("sacct -j %s | grep 'lindb'" % jobid).read()
#print('jobinfo = ',jobinfo)
#jobmem   = int([x for x in jobinfo.split() if 'mem' in x][0].split(",")[1].split('=')[1].replace("M",""))
os.system('echo jobmem = %s' % jobmem)
#timeinfo = os.popen("sacct -j %s --format Timelimit" % jobid).read()
#print('timeinfo = ',timeinfo)
#jobtime  = int(timeinfo.split()[-1].split(':')[0])
os.system('echo jobmem = %s' % jobmem)

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
#             if not mem == jobmem: # don't waste resources on jobs requiring less mem either
            if mem > jobmem: # since I'm not using RAC for jobs > 4G, allow jobs with less mem to run
                os.system('echo file does not match mem limit')
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
                    # make sure the command executed until completion, else kill the job to be swept by rescheduler.py
                    vcf = cmd.split()[-5]
                    tbi = vcf.replace(".gz",".gz.tbi")
                    checktbi(tbi) # in case I exceeded mem
                    
                    pipedir = os.popen('echo $HOME/pipeline').read().replace("\n","")
                    os.system('python %s %s' % (op.join(pipedir,'rescheduler.py'),
                                                fqdir))
                    os.system('python %s %s' % (op.join(pipedir,'scheduler.py'),
                                                fqdir))

                    break
else:
    os.system('echo no files to help')
    
    
    
    