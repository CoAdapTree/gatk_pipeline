###
# usage: python scheduler.py /path/to/fastq.gz/folder/
###

### imports
import os
import sys
from os import path as op
from os import listdir
def ls(DIR):
    return sorted([f for f in listdir(DIR)])
def fs (DIR):
    return sorted([op.join(DIR,f) for f in ls(DIR)])
###

### args
thisfile, fqdir = sys.argv
###

### reqs
if fqdir.endswith("/"):
    fqdir = fqdir[:-1]
DIR       = op.join(op.dirname(fqdir),'shfiles/gvcf_shfiles')
print "DIR=",DIR
assert op.exists(DIR)
scheduler = op.join(DIR,'scheduler.txt')
os.chdir(DIR)
qthresh   = 304
###

### defs
print 'running scheduler.py'
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
        o.write("scheduler")
def sbatchjobs(files):
    for f in files:
        if op.exists(f):
            os.system('sbatch %s' % f)     # maybe there will be dup sbatches? oh well hopefully not
            print f
            try:
                os.system('unlink %s' % f) # remove the symlink from the scheddir
                print 'made it this far'
            except OSError as e:           # unless another scheduler has done so (shouldnt be the case)
                print 'unable to unlink symlink %f' % f
                pass
def main():
    # write a file and reserve scheduling to this call of the scheduler, or pass if another scheduler is running
    startscheduler(scheduler) # reserve right away
    x = sq("squeue -u lindb | grep 'lindb' | wc -l") # number of jobs in the queue
    if x < qthresh: # if there is room in the queue
        print 'scheduler not running'
        print 'queue length less than thresh'
        nsbatch = qthresh - x # how many should I submit?
        files   = [f for f in fs(DIR) if not 'scheduler.txt' in f and '.out' not in f][0:nsbatch]
        if len(files) > 0:
            print 'submitting %s jobs' % str(nsbatch)
            sbatchjobs(files)
        else:
            print 'no files to sbatch'
    else:
        print 'scheduler was not running, but no room in queue' 
    delsched(scheduler)
###

# main
if not op.exists(scheduler): # if scheduler isn't running
    main()
else:
    print 'scheduler was running'