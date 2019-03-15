###
# FIX: make more pythonic
###

###
# usage rescheduler.py /path/to/parentdir-used-in-00_start-pipeline.py/
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
thisfile, parentdir = sys.argv
###

if parentdir.endswith("/"):
    parentdir = parentdir[:-1]
print('parentdir=',parentdir)
#print dirname
DIR = op.join(parentdir,'shfiles/supervised/select_variants_within_and_across') 
#print DIR
os.chdir(DIR)
outs = [f for f in fs(DIR) if f.endswith('out') and 'checked' not in f and 'swp' not in f]
rescheduler = op.join(DIR,'rescheduler.txt')
samp2pool = pickle.load(open(op.join(parentdir,'samp2pool.pkl'),'rb'))
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
def delrescheduler(rescheduler,createdrescheduler):
    if createdrescheduler == True:
        try:
            os.remove("%(rescheduler)s" % globals())
            os.system('echo removed rescheduler')
        except:
            os.system('echo could not remove rescheduler')
            pass
def getallpids():
    pids = os.popen('squeue -u lindb -o "%i"').read().split("\n")
    pids = [p for p in pids if not p == '']
    if len(pids) != len(set(pids)):
        print('len !- luni pids')
        delrescheduler(rescheduler,globals()['createdrescheduler'])
        exit()
    return pids[1:]
def checksq(rt):
    exitneeded = False
    if not type(rt) == list:
        os.system('echo "type(sq) != list, exiting rescheduler.py"')
        exitneeded = True
#     if len(rt) == 0:
#         os.system('echo "len(sq) == 0, exiting rescheduler.py"')
#         exitneeded = True
    count = 0
    for s in rt:
        if 'socket' in s.lower():
            os.system('echo "socket in sq return, exiting rescheduler.py"')
            exitneeded = True
        try:
            assert int(s.split()[0]) == float(s.split()[0])
            count += 1
        except:
            os.system('echo "could not assert int == float, %s %s"' % (s[0],s[0]))
            exitneeded = True
    if count == 0 and len(rt) > 0:
        os.system('echo never asserted pid, exiting rescheduler.py')
        exitneeded = True
    if exitneeded == True:
        delrescheduler(rescheduler,globals()['createdrescheduler'])
        exit()
def getsq():
    return [x for x in os.popen("squeue -u lindb | grep 'R 2'").read().split("\n") if not x == '']
def removeworker(DIR,trushfile):
    # remove worker from workingdir
    workingdir = op.join(DIR,'workingdir')
    worker = [f for f in fs(workingdir) if op.basename(trushfile) in f]
    if len(worker) == 1:
        worker = worker[0]
        try:
            os.unlink(worker)
        except:
            os.system('echo could not unlink worker: %s' % worker)
def getpids(sq):
    pids = []
    for s in sq:
        if not s == '':
            pid = s.split()[0]
            pids.append(pid)
    return pids
def bigbrother(rescheduler):
    # if the scheduler controller has died, remove the scheduler
    with open(rescheduler,'r') as o:
        text = o.read().replace("\n","")
    pid = text.split()[-1]
    if not pid == '=':
        pids = getallpids()
        if not pid in pids:
            print('controller was not running, so the scheduler was destroyed')
            delrescheduler(rescheduler,True)
        else:
            print('controller is running, allowing it to proceed')
def checktbis(shfile):
    # if some of the previous files were created, no need in wasting time recreating them
    with open(shfile,'r') as s:
        t = s.read().split("\n")
    for linenum,line in enumerate(t):
        if line.startswith('gatk') and 'IndexFeatureFile' not in line:
            vcf = [x for x in line.split() if x.endswith('.gz')][-1]
            if op.exists("%s.tbi" % vcf):
                t[linenum] = line.replace("gatk","#gatk")
                print('changing')
                print('\t%s' % line)
                print('\t%s' % t[linenum])
            else:
                print('file does not exist: %s.tbi' % vcf)
                print('leaving downstream commands as-is (as-are?)')
                break # no need to keep checking downstream files if upstream files have not been created
    text = '\n'.join(t)
    with open(shfile,'w') as s:
        s.write("%s" % text)


# identify outs that aren't running
sq = getsq()
#print(sq)
checksq(sq)
pids = getpids(sq)
#print(pids)
runs = []
for out in outs:
    pid = op.basename(out).split(".out")[0].split("---")[-1]
    if pid not in pids:
        runs.append(out)
outs = runs

os.system('echo running rescheduler.py')
createdrescheduler = False
if len(outs) > 0:
#     print('outs =',outs)
    if not op.exists(rescheduler):
        # reserve the rescheduler
        with open(rescheduler,'w') as o:
            jobid = os.popen('echo ${SLURM_JOB_ID}').read().replace("\n","")
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
            sq = getsq()
            checksq(sq)
            pids = getpids(sq)
            pid = op.basename(out).split(".out")[0].split("---")[-1] 
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
                if 'oom-kill' in line  or 'error' in line:
                    if not 'no mem or time errors found in' in line:
                        os.system ('echo found an error')
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
                    os.system('echo due to time limit')
                    # determine if its original or gvcf_helper
                    helped = False
                    for line in o[::-1]:
                        if 'getting help from genotyping_helper' in line and 'echo' not in line:
                            # we dont necessarily need to increase time if last call was due to gvcf_helper.py
                            helped = True
                            os.system('echo helped by genotyping_helper =%s' % helped)
                            break
                    if helped == True: # if the job ended on a call from gvcf_helper.py or was cancelled
                        # no need to change time this first time
                        os.system('echo leaving orginal time as-is')
                        os.system('echo cancelled =%s' % cancelled)
                        for line in o[::-1]:
                            if line.startswith('reading') or line.startswith('shfile ='):
                                trushfile = [x for x in line.split() if x.endswith('.sh')][0]
                                os.system ('echo checking tbis of original sh file %s' % trushfile)
                                checktbis(trushfile)
                                # add job back to the queue 
                                linkname = op.join(DIR,op.basename(trushfile))
                                addlink((trushfile,linkname))
                                # remove worker from workingdir
                                removeworker(DIR,trushfile)
                                break  
                    else:
                        for line in o:
                            # this was the call from the original sh file, need to change time
                            if line.startswith('shfile ='):
                                os.system ('echo adjusting time of original sh file')
                                trushfile = [x for x in line.split() if x.endswith('.sh')][0]
                                os.system('echo linked to %s' % trushfile)
                                sh = open(trushfile).read()
                                if '00:00:05' in sh: # for debugging/testing
                                    text = sh.replace('00:00:05','02:59:00')
                                elif '02:59:00' in sh:
                                    text = sh.replace('02:59:00','11:59:00')
                                    os.system('echo extending time to 11:59:00')
                                elif '11:59:00' in sh:
                                    text = sh.replace('11:59:00','23:59:00')
                                    os.system('echo extending time to 23:59:00')
                                elif '23:59:00' in sh:
                                    text = sh.replace('23:59:00','7-00:00:00')
                                    os.system('echo extending time to 7 days')
                                elif '7-00:00:00' in sh:
                                    text = sh.replace('7-00:00:00','14-00:00:00')
                                    os.system('echo extending time to 14 days')
                                else:
                                    os.system('echo cound not find replacement')
                                    os.system('cat %s' % trushfile)
                                    break
                                with open(trushfile,'w') as o:
                                    o.write("%s" % text)
                                # check to see which tbi files were made
                                os.system ('echo checking tbis of original sh file %s' % trushfile)
                                checktbis(trushfile)
                                # add job back to the queue  
                                linkname = op.join(DIR,op.basename(trushfile))
                                addlink((trushfile,linkname))
                                break
                elif cancelled == True:
                    # no need to change time this first time
                    os.system('echo leaving orginal time as-is')
                    os.system('echo cancelled =%s' % cancelled)
                    for line in o[::-1]:
                        if line.startswith('reading') or line.startswith('shfile ='):
                            trushfile = [x for x in line.split() if x.endswith('.sh')][0]
                            os.system ('echo checking tbis of original sh file %s' % trushfile)
                            checktbis(trushfile)
                            # add job back to the queue 
                            linkname = op.join(DIR,op.basename(trushfile))
                            addlink((trushfile,linkname))
                            # remove worker from workingdir
                            removeworker(DIR,trushfile)
                            break  
                else: # there's a mem oerror
                    edited = True
                    # at t=0, all sh files have mem==8000M, so if gvcf_helper.py caused mem error, the last call needs more mem
                    ## and I dont have to figure out if its the original's call or gvcf_helper's call
                    os.system('echo due to mem limit')
                    # find the last job and resubmit with more mem
                    for line in o[::-1]:
                        if line.startswith('reading') or line.startswith('shfile ='):
                            trushfile = [x for x in line.split() if x.endswith('.sh')][0]
                            linkname = op.join(DIR,op.basename(trushfile))
                            os.system('echo adjusting mem limit of: %s' % trushfile)
                            os.system('echo linked to %s' % linkname)
                            with open(trushfile) as trush:
                                sh = trush.read()
                            if '=2000M' in sh:
                                text = sh.replace("=2000M","=12000M")
                                os.system('echo increasing mem to 4G')
                            elif '=4000M' in sh:
                                text = sh.replace('=4000M','=12000M')
                                os.system('echo increasing mem to 6.5G')
                            elif '=6500M' in sh:
                                text = sh.replace('=6500M','=12000M')
                                os.system('echo increasing mem to 12G')
                            elif '=12000M' in sh:
                                text = sh.replace('=12000M','=20000M')
                                os.system('echo increasing mem to 30G')
                            elif '=20000M' in sh:
                                text = sh.replace('=20000M','=30000M')
                                os.system('echo increasing mem to 30G')
                            elif '=30000M' in sh:
                                text = sh.replace('=30000M','=50000M')
                                os.system('echo increasing mem to 50G')
                            elif '=50000M' in sh:
                                text = sh.replace('=50000M','=120000M')
                                os.system('echo increasing mem to 120G')
                            elif '=120000M' in sh:
                                text = sh.replace('=120000M','=200000M')
                                os.system('echo increasing mem to 200G')
                            with open(trushfile,'w') as o:
                                o.write("%s" % text)
                            # check to see which tbi files were made
                            os.system ('echo checking tbis of original sh file %s' % trushfile)
                            checktbis(trushfile)
                            # add job back to the queue  
                            addlink((trushfile,linkname))
                            # try and remove worker from workingdir
                            removeworker(DIR,trushfile)
                            break
            else:
                os.system('echo no mem or time errors found in %s' % out)
            # move the .out file to a new name so rescheduler doesn't look at it again
            dst = out.replace(".out","_checked.out")
            try:
                shutil.move(out,dst)
                os.system('echo moved file: %s' % dst)
            except OSError as e:
                os.system('echo could not move outfile to noerror.out: %s' % out)
                pass
    else:
        os.system('echo rescheduler was running')
        bigbrother(rescheduler)
        
else:
    os.system('echo rescheduler found no outfiles to analyze or all outfiles are for jobs currently running')

delrescheduler(rescheduler,createdrescheduler)
if not createdrescheduler == False and op.exists(rescheduler):
    bigbrother(rescheduler)
