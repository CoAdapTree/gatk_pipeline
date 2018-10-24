###
# usage: python scheduler.py /path/to/fastq.gz/folder/
###

### imports
import os
import sys
from os import path as op
from os import listdir
import time
import random
def ls(DIR):
    return sorted([f for f in listdir(DIR)])
def fs (DIR):
    return sorted([op.join(DIR,f) for f in ls(DIR)])
###

### args
thisfile, fqdir = sys.argv
###

### reqs
if fqdir.endswith("/"): #sometimes I run the scheduler from the command line, which appends / which screws up op.dirname()
    fqdir = fqdir[:-1]
DIR       = op.join(op.dirname(fqdir),'shfiles/gvcf_shfiles')
print("DIR=",DIR)
assert op.exists(DIR)
scheduler = op.join(DIR,'scheduler.txt')
os.chdir(DIR)
qthresh   = 700 
###

### defs
print('running scheduler.py')
def sq(command):
    # how many jobs are running
    return int(os.popen(str(command)).read().replace("\n",""))
def delsched(scheduler):
    # stop scheduler
    try:
        os.remove(scheduler)
    except OSError as e:
        pass
def startscheduler(scheduler):
    with open(scheduler,'w') as o:
        # after creating the file, write job id in case i want to cancel process
        jobid = os.popen('echo ${SLURM_JOB_ID}').read().replace("\n","")
        o.write("scheduler id = %s" % jobid)
def sbatchjobs(files):
    for f in files:
        realp = op.realpath(f) # find the file to which the symlink file is linked
        if op.exists(f):
            # print (f)
            try:
                os.unlink(f) # first remove the symlink from the scheddir
                print('unlinked %s' % f)
            except:          # unless gvcf_helper has already done so (shouldnt be the case, but maybe with high qthresh)
                print('unable to unlink symlink %f' % f)
                continue
            os.system('sbatch %s' % realp) # then sbatch the real sh file if & only if the symlink was successfully unlinked    
            
def main(DIR):
    # write a file and reserve scheduling to this call of the scheduler, or pass if another scheduler is running
    startscheduler(scheduler) # reserve right away
    x = sq("squeue -u lindb | grep 'lindb' | wc -l") # number of jobs in the queue
    print ('queue length = ',x)
    if x < qthresh: # if there is room in the queue
        print('scheduler not running')
        print('queue length less than thresh')
        nsbatch = qthresh - x # how many should I submit?
        print ('nsbatch =',nsbatch)
        print (len(fs(DIR)))
        files = [f for f in fs(DIR) if 'scheduler.txt' not in f and '.out' not in f and 'workingdir' not in f][0:nsbatch]
        if len(files) > 0:
            print('submitting %s jobs' % str(len(files)))
            print(files)
            sbatchjobs(files)
        else:
            print('no files to sbatch')
    else:
        print('scheduler was not running, but no room in queue' )
    delsched(scheduler)
###

# main
time.sleep(random.random())  # just in case the very first instances of scheduler.py start at v similar times
if not op.exists(scheduler): # if scheduler isn't running
    main(DIR)
else:
    print('scheduler was running')
