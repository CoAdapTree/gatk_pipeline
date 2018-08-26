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
thifile, fqdir = sys.argv
###

os.system('source $HOME/.bashrc')
DIR = op.join(op.dirname(fqdir),'shfiles/gvcf_shfiles')
os.chdir(DIR)

shfiles = [f for f in fs(DIR) if f.endswith('.sh')]
shuffle(shfiles) # so I don't have to code something in to deal with/avoid multiple gvcf_helpers 'helping' the same sh file
# shuffling won't be perfect, esp when len(shfiles) == <a smallish number>, but oh well

# run commands until I run out of time
print 'running gvcf_helper.py'
for s in shfiles:
    print s # so that rescheduler can find it (should print to stdout in the sh file's .out file)
    o = open(s,'r').readlines()
    for line in o:
        if line.startswith('gatk'):
            cmd = line.replace('\n','')
            print 'running cmd:'
            print cmd
            os.system('%s' % cmd)
            print 'unlinking shfile %s'
            try:
                os.system('unlink %s' % s)
            except OSError as e:
                print 'unable to unlink %s' % s
                pass
            pipedir = os.popen('echo $HOME/pipeline').read().replace("\n","")
            os.system('python %s %s' % (op.join(pipedir,'scheduler.py'),
                                        fqdir))
            os.system('python %s %s' % (op.join(pipedir,'rescheduler.py'),
                                        fqdir))
            break
          

