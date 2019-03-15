###
# usage: python genotyping_scheduler.py /path/to/<parentdir>/
###

###
# fix: merge with scheduler.py
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
def luni(mylist):
    return len(set(mylist))
###

### args
thisfile, parentdir = sys.argv
###

### reqs
if parentdir.endswith("/"): #sometimes I run the scheduler from the command line, which appends / which screws up op.dirname()
    parentdir = parentdir[:-1]
scheddir  = op.join(parentdir,'shfiles/supervised/select_variants_within_and_across')
print("scheddir=",scheddir)
assert op.exists(scheddir)
scheduler = op.join(scheddir,'scheduler.txt')
os.chdir(scheddir)
qthresh   = 1
user = os.popen("echo $USER").read().replace("\n","")
###

### defs
print('running scheduler.py')
def checksq(rt):
    exitneeded = False
    if not type(rt) == list:
        os.system('echo "type(sq) != list, exiting rescheduler.py"')
        exitneeded = True
    count = 0
    for s in rt:
        if 'socket' in s.lower():
            os.system('echo "socket in sq return, exiting rescheduler.py"')
            exitneeded = True
        try:
#             print(s)
            assert int(s.split()[0]) == float(s.split()[0])
            count =+ 1
        except:
            os.system('echo "could not assert int == float, %s %s"' % (s[0],s[0]))
            exitneeded = True
    if count == 0 and len(rt) > 0:
        os.system('echo never asserted pid, exiting rescheduler.py')
        exitneeded = True
    if exitneeded == True:
        delsched(globals()['scheduler'])
        exit()
def sq(command):
    # how many jobs are running
    q = [x for x in os.popen(str(command)).read().split("\n") if not x == '']
    checksq(q)
    return len(q)
def delsched(scheduler):
    # stop scheduler
    try:
        os.remove(scheduler)
    except:
        pass
def getpids():
    pids = os.popen('squeue -u lindb -o "%i"').read().split("\n")
    pids = [p for p in pids if not p == '']
    if len(pids) != luni(pids):
        print('len !- luni pids')
        delsched(scheduler)
        exit()
    return pids[1:]
def startscheduler(scheduler):
    with open(scheduler,'w') as o:
        # after creating the file, write job id in case i want to cancel process
        jobid = os.popen('echo ${SLURM_JOB_ID}').read().replace("\n","")
        o.write("scheduler id = %s" % jobid)
    # double check that the scheduler is correct
    with open(scheduler,'r') as o:
        text = o.read()
    if not text.split()[-1] == '=':
        if not text.split()[-1] == jobid:
            os.system('echo another scheduler is in conflict. Allowing other scheduler to proceed. Exiting')
            exit()
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
    x = sq("squeue -u %(user)s | grep genotype " % globals()) # number of genotyping jobs in the queue
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
        print('genotyping_scheduler was not running, but no room in queue' )
    pipedir = os.popen('echo $HOME/gatk_pipeline').read().replace("\n","")
    os.system('python %s geno' % (op.join(pipedir,'balance_queue.py')))
    delsched(scheduler)
def bigbrother(scheduler,DIR):
    # if the scheduler controller has died, remove the scheduler
    with open(scheduler,'r') as o:
        text = o.read().replace("\n","")
    pid = text.split()[-1]
    if not pid == '=':
        pids = getpids()
        if not pid in pids:
            print('controller was not running, so the scheduler was destroyed')
            delsched(scheduler)
            main(DIR)
        else:
            print('controller is running, allowing it to proceed')
###

# main
time.sleep(random.random())  # just in case the very first instances of scheduler.py start at v similar times
if not op.exists(scheduler): # if scheduler isn't running
    main(scheddir)
else:
    print('scheduler was running')
    bigbrother(scheduler,scheddir)
