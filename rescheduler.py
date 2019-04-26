###
# FIX: make more pythonic
###

###
# usage rescheduler.py pooldir
###

### imports
import sys
import shutil
from random import shuffle
from coadaptree import *
from balance_queue import getsq
###

### args
thisfile, pooldir = sys.argv
###

if pooldir.endswith("/"):
    pooldir = pooldir[:-1]
parentdir = op.dirname(pooldir)
os.system('echo parentdir = %s' % parentdir)
#print dirname
scheddir = op.join(parentdir, 'shfiles/gvcf_shfiles') 
#print scheddir
os.chdir(scheddir)
outs = [f for f in fs(scheddir) if f.endswith('out') and 'checked' not in f and 'swp' not in f]
rescheduler = op.join(scheddir, 'rescheduler.txt')
samp2pool = pklload(op.join(parentdir, 'samp2pool.pkl'))


def vcf2sh(v):
    print('v =', v)
    pooldir = op.dirname(op.dirname(v))
    pool = op.basename(pooldir)
    gvcfdir = op.join(pooldir, 'shfiles/04_gvcf_shfiles')
    bname = op.basename(v)
    shname = bname.replace("raw_", "").replace(".g.vcf.gz", ".sh")
    shfile = op.join(gvcfdir,shname)
    return shfile


def unlink(linkname):
    try:
        os.unlink(linkname)
        os.system('echo unlinked %s' % linkname)
    except:
        os.system('echo no symlink to unlink: %s' % linkname)
        pass


def addlink(args):
    trushfile,linkname = args
    os.system('echo symlink from: %s' % linkname)
    os.system('echo to: %s' % trushfile)
    if not op.exists(linkname):
        os.symlink(trushfile,linkname)
        os.system('echo added symlink to queue: %s' % linkname)
    else:
        os.system('echo unable to create symlink from %s to %s' % (linkname,trushfile))     


def delrescheduler(rescheduler,createdrescheduler):
    if createdrescheduler is True:
        try:
            os.remove("%(rescheduler)s" % globals())
            os.system('echo removed rescheduler')
        except:
            os.system('echo could not remove rescheduler')
            pass


def removeworker(scheddir,trushfile):
    # remove worker from workingdir
    workingdir = op.join(scheddir, 'workingdir')
    worker = [f for f in fs(workingdir) if op.basename(trushfile) in f]
    if len(worker) == 1:
        worker = worker[0]
        try:
            os.unlink(worker)
        except:
            os.system('echo could not unlink worker: %s' % worker)


def getpids(sq):
    pids = []
    if len(sq) > 0:
        print('getpids len(sq) > 0')
        pids = []
        for q in sq:
            if not q == '':
                pid = q[0]
                pids.append(pid)
    return pids


def bigbrother(rescheduler):
    print('reschduler = ', rescheduler)
    # if the scheduler controller has died, remove the scheduler
    with open(rescheduler, 'r') as o:
        text = o.read().replace("\n", "")
        os.system('echo %s' % text)
    pid = text.split()[-1]
    if not pid == '=':
        sq = getsq(states=['running'])
        pids = getpids(sq)
        if not pid in pids:
            os.system('echo controller was not running, so the scheduler was destroyed')
            delrescheduler(rescheduler, True)
        else:
            os.system('echo controller is running, allowing it to proceed')


# identify outs that aren't running
createdrescheduler = False
sq = getsq(states=['running'])
print('len(sq) = ', len(sq))
print(sq)
pids = getpids(sq)
#print(pids)
runs = []
for out in outs:
    pid = op.basename(out).split(".out")[0].split("_")[-1]
    if pid not in pids:
        runs.append(out)
outs = runs

os.system('echo running rescheduler.py')
if len(outs) > 0:
    if not op.exists(rescheduler):
        # reserve the rescheduler
        with open(rescheduler, 'w') as o:
            jobid = os.environ['SLURM_JOB_ID']
            o.write("rescheduler id = %s" % jobid)
        createdrescheduler = True
        # double check that the rescheduler is correct
        with open(rescheduler,'r') as o:
            text = o.read()
        if not text.split()[-1] == '=':
            if not text.split()[-1] == jobid:
                os.system('echo another rescheduler is in conflict. Allowing other rescheduler to proceed. Exiting')
                time.sleep(5)
                exit()

        # look for errors in outfiles and resubmit the error-causing shfile using more mem or time
        for out in outs:
            #check again to see if job is running
            sq = getsq(states=['running'])
            pids = getpids(sq)
            pid = op.basename(out).split(".out")[0].split("_")[1]
            if pid in pids:
                continue
            if not op.exists(out):
                continue
            os.system('echo -e \n')
            os.system('echo working on %s' % out)
            with open(out,'r') as OUT:
                o = OUT.readlines()
            # look for mem error
            edited = False
            timelimit  = False
            founderror = False
            cancelled = False
            for line in o[-20:]: # look for an error message
                if 'oom-kill' in line or 'error' in line:
                    if not 'no mem or time errors found in' in line:
                        os.system ('echo found an error')
                        founderror = True
                        break
            if founderror is True:
                for test in o[-20:]: # look for a timeout error 
                    if 'time limit' in test.lower() or 'cancelled' in test.lower():
                        timelimit = True
                        if 'cancelled' in test.lower():
                            cancelled = True
                        break
                if timelimit is True:
                    # look for time error
                    edited = True
                    # time error could be caused by the original sh command running out of time or ...
                    ## ... gvcf_helper.py ran out of time
                    os.system('echo due to time limit')
                    # determine if its original or gvcf_helper
                    helped = False
                    for line in o[::-1]:
                        if 'getting help from gvcf_helper' in line and 'echo' not in line:
                            # we dont necessarily need to increase time if last call was due to gvcf_helper.py
                            helped = True
                            os.system('echo helped by gvcf_helper =%s' % helped)
                            break
                    if helped is True or cancelled is True: # if the job ended on a call from gvcf_helper.py or was cancelled
                        # no need to change time this first time
                        os.system('echo leaving orginal time as-is')
                        os.system('echo cancelled =%s' % cancelled)
                        for line in o[::-1]:
                            if line.startswith('gatk HaplotypeCaller'):
                                vcf = line.split()[-5]
                                trushfile = vcf2sh(vcf)
                                linkname = op.join(scheddir, op.basename(trushfile))
                                
                                # add job back to the queue 
                                addlink((trushfile, linkname))
                                
                                # remove worker from workingdir
                                removeworker(scheddir,trushfile)
                                
                                break
                    else:
                        for line in o[::-1]:
                            # this was the call from the original sh file
                            if line.startswith('gatk HaplotypeCaller'):
                                os.system ('echo adjusting time of original sh file')
                                vcf = line.split()[-5]
                                os.system('echo vcf file = %s' % vcf)
                                trushfile = vcf2sh(vcf)

                                os.system('echo linked to %s' % trushfile)
                                with open(trushfile,'r') as O:
                                    sh = O.read()
#                                 sh = open(trushfile).read()
                                if '00:00:05' in sh: # for debugging/testing
                                    text = sh.replace('00:00:05', '02:59:00')
                                elif '02:59:00' in sh:
                                    text = sh.replace('02:59:00', '11:59:00')
                                    os.system('echo extending time to 11:59:00')
                                elif '11:59:00' in sh:
                                    text = sh.replace('11:59:00', '23:59:00')
                                    os.system('echo extending time to 23:59:00')
                                elif '23:59:00' in sh:
                                    text = sh.replace('23:59:00', '7-00:00:00')
                                    os.system('echo extending time to 7 days')
                                elif '14:30:00' in sh:
                                    text = sh.replace('7-00:00:00', '14-00:00:00')
                                    os.system('echo replacing 14-hour time with 7 days')
                                else:
                                    os.system('echo cound not find replacement')
                                    os.system('cat %s' % trushfile)
                                    break
                                with open(trushfile,'w') as O:
                                    O.write("%s" % text)

                                # add job back to the queue  
                                linkname = op.join(scheddir, op.basename(trushfile))
                                addlink((trushfile, linkname))
                                
                                break
                else: # there's a mem oerror
                    edited = True
                    # at t=0, all sh files have mem==8000M, so if gvcf_helper.py caused mem error, the last call needs more mem
                    ## and I dont have to figure out if its the original's call or gvcf_helper's call
                    os.system('echo due to mem limit')
                    # find the last job and resubmit with more mem
                    for line in o[::-1]:
                        if line.startswith('gatk HaplotypeCaller'):
                            vcf = line.split()[-5]
                            trushfile = vcf2sh(vcf)
                            os.system('echo adjusting mem limit of: %s' % trushfile)
                            os.system('echo linked to %s' % trushfile)
                            with open(trushfile,'r') as O:
                                sh = O.read()
#                             sh = open(trushfile).read()
                            if '4000M' in sh:
                                text = sh.replace("4000M", "12000M")
                                os.system('echo increasing mem to 12G')
                            elif '6500M' in sh:
                                text = sh.replace('6500M', '12000M')
                                os.system('echo increasing mem to 12G')
                            elif '8000M' in sh:
                                text = sh.replace('8000M', '12000M')
                                os.system('echo increasing mem to 12G')
                            elif '12000M' in sh:
                                text = sh.replace("12000M", "20000M")
                                os.system('echo increasing mem to 20G')
                            elif '20000M' in sh:
                                text = sh.replace('20000M', '30000M') # keep it in, i changed last if statment, was 8Gb->20Gb
                                os.system('echo increasing mem to 30G')
                            elif '30000M' in sh:
                                text = sh.replace('30000M', '50000M')
                                os.system('echo increasing mem to 50G')
                            elif '50000M' in sh:
                                text = sh.replace('50000M', '100000M')
                                os.system('echo increasing mem to 100G')
                            elif '100000M' in sh:
                                text = sh.replace('100000M', '120000M')
                                os.system('echo increasing mem to 120G')
                            with open(trushfile, 'w') as O:
                                O.write("%s" % text)
                            
                            # add job back to the queue  
                            linkname = op.join(scheddir, op.basename(trushfile))
                            addlink((trushfile, linkname))
                            
                            # remove worker from workingdir
                            removeworker(scheddir, trushfile)
                            break
            else:
                os.system('echo no mem or time errors found in %s' % out)

            # move the .out file to a new name so rescheduler doesn't look at it again
            dst = out.replace(".out", "_checked.out")
            try:
                shutil.move(out, dst)
                os.system('echo moved file: %s' % dst)
            except OSError as e:
                os.system('echo could not move outfile to noerror.out: %s' % out)
                pass
    else:
        os.system('echo rescheduler was running')
        bigbrother(rescheduler)
else:
    os.system('echo rescheduler found no outfiles to analyze or all outfiles are for jobs currently running')

delrescheduler(rescheduler, createdrescheduler)
if createdrescheduler is not False and op.exists(rescheduler):
    bigbrother(rescheduler)
