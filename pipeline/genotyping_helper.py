### FIX
# at some point, have the jobID/mem as an input arg instead of scripted
# merge gvcf_helper.py and genotyping_helper.py
###

###
# usage: genotyping_helper.py /path/to/parentdir-used-in-00_start-pipeline.py/

###
# purpose: to keep running gatk GenotypeGVCF commands until time or memory runs out
###

### imports
import os
import sys
import pickle
import shutil
import signal
from os import listdir
from os import path as op
from random import shuffle
def ls(DIR):
    return sorted([f for f in listdir(DIR)])
def fs (DIR):
    return sorted([op.join(DIR,f) for f in ls(DIR)])
###

### args
thisfile, parentdir, outfile = sys.argv
if parentdir.endswith("/"):
    parentdir = parentdir[1:]
###

# balance the queue
pipedir = os.popen('echo $HOME/pipeline').read().replace("\n","")
os.system('python %s genotyp' % op.join(pipedir,'balance_queue.py'))

# make sure the previous outfile was created
def checkoutfile(outfile):
    if not op.exists(outfile):
        os.system('echo previous outfile did not finish: %s' % outfile) 
        os.kill(os.getppid(), signal.SIGHUP)
checkoutfile(outfile)


# os.system('source $HOME/.bashrc')
scheddir = op.join(parentdir,'shfiles/supervised/select_variants_within_and_across')
os.chdir(scheddir)
workingdir = op.join(scheddir,'workingdir')
if not op.exists(workingdir):
    os.makedirs(workingdir)
    
# get job info and current memory/time limits
jobid    = os.popen('echo ${SLURM_JOB_ID}').read().replace("\n","")
if not float(jobid) == int(jobid):
    # if I can't retrieve the jobid, just let the job die
    # ensures that there wasn't an error to the slurm request (some times I get timeouts etc as returns)
    exit()
os.system('echo jobid = %s' % jobid)
f = [f for f in fs(scheddir) if str(jobid) in f]
try:
    assert len(f) == 1
    f = f[0]
except:
    os.system('echo could not find slurm job file, exiting genotyping_helper')
    exit()
with open(f,'r') as o:
    text = o.read().split("\n")
for line in text:
    if 'mem=' in line:
        jobmem = int(line.split("=")[-1][:-1])
    if 'time=' in line:
        splits = line.split("=")[-1].split("-")
        if len(splits) > 1:            
            # multi-day job
            jobtime = int(splits[0])*24 # this will never happen
        else:
            # xhour job
            jobtime = int(line.split("=")[-1].split(':')[0] )    
os.system('echo jobmem = %s' % jobmem)
os.system('echo jobtime = %s' % jobtime)

# get list of remaining gatk calls
shfiles = [f for f in fs(scheddir) if f.endswith('.sh')]
shuffle(shfiles) 

# run commands until I run out of time
os.system('echo running gvcf_helper.py')
badcount = 0
if len(shfiles) > 0:
    for s in shfiles:
    #     print (s)
        if badcount > 500:
            # no need to have the job keep looking for other jobs if most other jobs have different reqs (mem, time)
            os.system('echo exceeded badcount, exiting')
            exit()
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
            if not mem == jobmem: # don't waste resources on jobs requiring less mem either
                os.system('echo file does not match mem limit')
                shutil.move(reservation,s) # put the job back in the queue
                badcount += 1
                continue
            # only continue to run jobs that might fit in same time allocation
            TIME = int([x for x in o if 'time' in x][0].split("=")[1].split(':')[0])
            if TIME > jobtime:
                os.system('echo file exceeds necessary time')
                shutil.move(reservation,s)
                badcount += 1
                continue

            os.system('echo file is ok to proceed')
            os.system('echo reading %s' % op.realpath(reservation))
            gatkcmds = []
            for line in o:
                if line.startswith('gatk'):
                    gatkcmds.append(line)
            for cmd in gatkcmds:
                cmd = cmd.replace('\n','')
                os.system('echo running cmd:')
                os.system('echo %s' % cmd)
                os.system('%s' % cmd)
                # check that outfile was made
                outfile = cmd.split()[-1]
                checkoutfile(outfile)
            try:
                os.unlink(reservation)
                os.system('echo unlinked shfile %s' % reservation)
            except:
                os.system('echo unable to unlink %s' % reservation)
                pass
            
            pipedir = os.popen('echo $HOME/pipeline').read().replace("\n","")
            os.system('python %s %s' % (op.join(pipedir,'genotyping_rescheduler.py'),
                                        parentdir))
            os.system('python %s %s' % (op.join(pipedir,'genotyping_scheduler.py'),
                                        parentdir))
else:
    os.system('echo no files to help')
    
    