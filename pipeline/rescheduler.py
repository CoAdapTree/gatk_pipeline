###
# usage rescheduler.py /path/to/fq.gzfiles/folder/
###

### imports
import sys
import os
from os import path as op
from os import listdir
import shutil
def ls(DIR):
    return sorted([f for f in listdir(DIR)])
def fs (DIR):
    return sorted([op.join(DIR,f) for f in ls(DIR)])
###

### args
thifile, fqdir = sys.argv
###


dirname = op.dirname(fqdir)
#print dirname
DIR = op.join(dirname,'shfiles/gvcf_shfiles') # real deal
#DIR = op.join(op.dirname(fqdir),'shfiles/') # practice
#print DIR
os.chdir(DIR)
outs = [f for f in fs(DIR) if f.endswith('out') and not 'checked' in f]
rescheduler = op.join(DIR,'rescheduler.txt')

print 'running rescheduler.py'
if len(outs) > 0:
    print 'outs =',outs
    if not op.exists(rescheduler):
        # reserve the rescheduler
        with open(rescheduler,'w') as o:
            o.write("rescheduler")

        # look for errors in outfiles and resubmit the error-causing shfile using more mem or time
        for out in outs:
            edited = False
            print 'working on %s' % out
            o = open(out,'r').readlines()
            # look for mem error
            if 'oom-kill' in o[-1]:
                edited = True
                # at t=0, all sh files have mem==8000M, so if gvcf_helper.py caused mem error, the last call needs more mem
                ## and I dont have to figure out if its the original's call or gvcf_helper's call
                print 'due to mem limit'
                # find the last job and resubmit with more mem
                for line in o[::-1]:
                    if dirname in line and line.endswith('.sh\n'):
                        shfile = line.replace("\n","")
                        print 'linked to %s' % shfile
                        sh = open(shfile).read()
                        if '8000M' in sh:
                            text = sh.replace('8000M','20000M')
                        elif '20000M' in sh:
                            text = sh.replace('20000M','30000M')
                        elif '30000M' in sh:
                            text = sh.replace('30000M','50000M')
                        elif '50000M' in sh:
                            text = sh.replace('50000M','100000M')
                        with open(shfile,'w') as o:
                            o.write("%s" % text)
                        print 'sbatching due to mem limit: %s' % shfile
                        os.system('sbatch %s' % shfile)
                        bname = op.basename(shfile)
                        linkname = op.join(DIR,bname)
                        print 'trying to unlink'
                        try:
                            os.system('unlink %s' % linkname)
                            print 'unlinked %s' % linkname
                        except:
                            # if original file there is no symlink
                            print 'unable to unlink %s' % shfile
                            pass
                        break
            # look for time error
            elif 'TIME LIMIT' in o[-1]:
                edited = True
                # time error could be caused by the original sh command running out of time or gvcf_helper.py ran out of time
                print 'due to time limit'
                # determine if its original or gvcf_helper
                linked = False
                for line in o[::-1]:
                    if 'unlinking' in line:
                        # we dont necessarily need to increase time if last call was due to gvcf_helper.py
                        linked = True
                        break
                if linked == True:
                    # no need to change time this first time
                    print 'leaving orginal time as-is'
                    for line in o[::-1]:
                        if dirname in line and line.endswith('.sh\n'):
                            shfile = line.replace("\n","")
                            print 'sbatching %s' % shfile
                            os.system('sbatch %s' % shfile)
                            # unlink 
                            bname = op.basename(shfile)
                            linkname = op.join(DIR,bname)
                            try:
                                os.system('unlink %s' % linkname)
                                print 'unlinked %s' % linkname
                            except OSError as e:
                                print 'no symlink to unlink: %s' % linkname
                                pass
                else:
                    for line in o[::-1]:
                        if dirname in line and line.endswith('.sh\n'):
                            shfile = line.replace("\n","")
                            print 'linked to %s' % shfile
                            sh = open(shfile).read()
                            if '00:00:05' in sh: # for debugging/testing
                                text = sh.replace('00:00:05','02:59:00')
                            elif '02:59:00' in sh:
                                text = sh.replace('02:59:00','11:59:00')
                            elif '11:59:00' in sh:
                                text = sh.replace('11:59:00','23:59:00')
                            elif '23:59:00' in sh:
                                text = sh.replace('23:59:00','7-00:00:00')
                            with open(shfile,'w') as o:
                                o.write("%s" % text)
                            print 'sbatching due to time limit: %s' % shfile
                            os.system('sbatch %s' % shfile)
                            try:
                                bname = op.basename(shfile)
                                linkname = op.join(DIR,bname)
                                os.system('unlink %s' % linkname)
                                print 'unlinked %s' % linkname
                            except OSError as e:
                                print 'could not unlink: %s' % linkname
                                pass
                            break
            else:
                print 'no mem or time errors found in %s' % out

            # move the .out file to a new name so rescheduler doesn't look at it again
            if edited = True: # files with errors arent running, but files without errors might be running (dont want to move the outfile of a running job)
                dst = out.replace(".out","_checked.out")
                try:
                    shutil.move(out,dst)
                except OSError as e:
                    print 'could not move outfile to noerror.out: %s' % out
                    pass
        try:
            os.remove(rescheduler)
        except OSError as e:
            pass
    else:
        print 'rescheduler was running'
else:
    print 'rescheduler found no outfiles to analyze'

