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

os.system('source $HOME/.bashrc')
DIR = op.join(op.dirname(fqdir),'shfiles/gvcf_shfiles')
os.chdir(DIR)
workingdir = op.join(DIR,'workingdir')
if not op.exists(workingdir):
    os.makedirs(workingdir)

shfiles = [f for f in fs(DIR) if f.endswith('.sh')]
shuffle(shfiles) # so I don't have to code something in to deal with/avoid multiple gvcf_helpers 'helping' the same sh file
# shuffling won't be perfect, esp when len(shfiles) == <a smallish number>, but oh well

# run commands until I run out of time
print 'running gvcf_helper.py'
for s in shfiles:
    if op.exists(s):
#         print s # so that rescheduler can find it (should print to stdout in the sh file's .out file)
        reservation = op.join(workingdir,op.basename(s))
        try:
            shutil.move(s,reservation) # so that other jobs don't rewrite
        except:
            print 'could not move shfile %s' % s
            print 'to reservation %s' % reservation
            continue
        print reservation # so that rescheduler can find it (should print to stdout in the sh file's .out file)
        o = open(reservation,'r').readlines()
        for line in o:
            if line.startswith('gatk'):
                cmd = line.replace('\n','')
                print 'running cmd:'
                print cmd
                os.system('%s' % cmd)
                try:
                    os.system('unlink %s' % reservation)
                    print 'unlinking shfile %s' % reservation
                except OSError as e:
                    print 'unable to unlink %s' % reservation
                    pass
                pipedir = os.popen('echo $HOME/pipeline').read().replace("\n","")
                os.system('python %s %s' % (op.join(pipedir,'scheduler.py'),
                                            fqdir))
                os.system('python %s %s' % (op.join(pipedir,'rescheduler.py'),
                                            fqdir))
                break


