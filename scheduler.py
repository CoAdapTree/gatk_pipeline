###
# usage: python scheduler.py /path/to/fastq.gz/folder/
###

### imports
import os, sys, time, random, subprocess
from os import path as op
from coadaptree import *
from balance_queue import getsq
###

### args
thisfile, pooldir = sys.argv
###

### reqs
if pooldir.endswith("/"): #sometimes I run the scheduler from the command line, which appends / which screws up op.dirname()
    pooldir = pooldir[:-1]
parentdir = op.dirname(pooldir)
scheddir = op.join(parentdir, 'shfiles/gvcf_shfiles')
print("scheddir=", scheddir)

scheduler = op.join(scheddir, 'scheduler.txt')
os.chdir(scheddir)
cluster = os.environ['CC_CLUSTER']  # which compute canada cluster is this job running on?
qthresh = 950 if cluster == 'cedar' else 950

user = os.environ['USER']
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
            assert int(s.split()[0]) == float(s.split()[0])
            count += 1
        except:
            os.system('echo "could not assert int == float, %s %s"' % (s[0], s[0]))
            exitneeded = True
    if count == 0 and len(rt) > 0:
        os.system('echo never asserted pid, exiting rescheduler.py')
        exitneeded = True
    if exitneeded is True:
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
    except OSError as e:
        pass


def getpids():
    pids = os.popen('squeue -u lindb -h -o "%i"').read().split("\n")
    pids = [p for p in pids if not p == '']
    if len(pids) != len(set(pids)):
        print('len != luni pids')
        delsched(scheduler)
        exit()
    return pids


def startscheduler(scheduler):
    with open(scheduler, 'w') as o:
        # after creating the file, write job id in case i want to cancel process
        jobid = os.environ['SLURM_JOB_ID']
        o.write("scheduler id = %s" % jobid)
    # double check that the scheduler is correct
    with open(scheduler, 'r') as o:
        text = o.read()
    if not text.split()[-1] == '=':
        if not text.split()[-1] == jobid:
            os.system('echo another scheduler is in conflict. Allowing other scheduler to proceed. Exiting')
            exit()


def sbatchjobs(files):
    for f in files:
        realp = op.realpath(f) # find the file to which the symlink file is linked
        if op.exists(f): # as long as the symlink is still there
            # print (f)
            try:
                os.unlink(f) # first try to remove the symlink from the scheddir
                print('unlinked %s' % f)
            except:          # unless gvcf_helper has already done so (shouldnt be the case, but maybe with high qthresh)
                print('unable to unlink symlink %s' % f)
                continue
            os.system('sbatch %s' % realp) # then sbatch the real sh file if & only if the symlink was successfully unlinked    


def main(scheddir):
    # write a file and reserve scheduling to this call of the scheduler, or pass if another scheduler is running
    startscheduler(scheduler) # reserve right away
    x = len(getsq())
    print ('queue length = ', x)
    if x < qthresh: # if there is room in the queue
        print('scheduler not running')
        print('queue length less than thresh')
        nsbatch = qthresh - x # how many should I submit?
        print ('nsbatch =', nsbatch)
        print (len(fs(scheddir)))
        files = [f for f in fs(scheddir) 
                 if 'scheduler.txt' not in f 
                 and '.out' not in f 
                 and 'workingdir' not in f][0:nsbatch]
        if len(files) > 0:
            print('submitting %s jobs' % str(len(files)))
            print(files)
            sbatchjobs(files)
        else:
            print('no files to sbatch')
    else:
        print('scheduler was not running, but no room in queue' )
    balance_queue = op.join(os.environ['HOME'], 'gatk_pipeline/balance_queue.py')
    subprocess.call([sys.executable, balance_queue, 'scatter', parentdir])
    delsched(scheduler)


def bigbrother(scheduler, scheddir):
    # if the scheduler controller has died, remove the scheduler
    with open(scheduler, 'r') as o:
        text = o.read()
    pid = text.split()[-1]
    if not pid == '=':
        pids = getpids()
        if not pid in pids:
            print(f'controller ({pid}) was not running, so the scheduler was destroyed')
            delsched(scheduler)
            main(scheddir)
        else:
            print('controller is running, allowing it to proceed')
###


# main
time.sleep(random.random())  # just in case the very first instances of scheduler.py start at v similar times
if not op.exists(scheduler): # if scheduler isn't running
    main(scheddir)
else:
    print('scheduler was running')
    bigbrother(scheduler, scheddir)

