import sys
import os
import pickle
import random
import string
import math
import shutil
import time
import numpy as np
import pandas as pd
from IPython.display import clear_output
from collections import OrderedDict, Counter
from IPython.display import Markdown, display
from os import path as op
from os import chdir as cd
from os import getcwd as cwd
from os import listdir
from shutil import copy as cp
from shutil import move as mv
from ipyparallel import Client
from decimal import Decimal
from datetime import timedelta
from datetime import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as pl
import matplotlib.dates as mdates

pd.set_option('display.max_columns', 100)
def ls(DIR):
    return sorted([f for f in listdir(DIR)])
def fs(DIR):
    return sorted([op.join(DIR,f) for f in ls(DIR)])
def uni(mylist):
    mylist = list(mylist)
    return np.unique(mylist).tolist()
def luni(mylist):
    return len(uni(mylist))
def suni(mylist):
    return sorted(uni(mylist))
def nrow(df):
    return len(df.index)
def ncol(df):
    return len(df.columns)
def table(lst,exclude=[]):
    c = Counter()
    for x in lst:
        c[x] += 1
    if len(exclude) > 0:
        for ex in exclude:
            if ex in c:
                c.pop(ex)
    return c
def pkldump(obj,f):
    with open(f,'wb') as o:
        pickle.dump(obj,o,protocol=pickle.HIGHEST_PROTOCOL)
def head(df):
    return df.head()
def sbatch(files):
    if not type(files) == list:
        files = [files]
    for f in files:
        os.chdir(op.dirname(f))
        os.system('sbatch %s' % f)
def update(args):
    clear_output(wait=True)
    [print(x) for x in args]
def keys(Dict):
    return list(Dict.keys())
def values(Dict):
    return list(Dict.values())
def setindex(df,colname):
    df.index = df[colname].tolist()
    df.index.names = ['']
    df = df[[c for c in df.columns if not colname in c]]
    return df
def pklload(path):
    pkl = pickle.load(open(path,'rb'))
    return pkl
def gettimestamp(f):
    return time.ctime(os.path.getmtime(f))
def getmostrecent(files,remove=False):
    if not type(files) == list:
        files = [files]
    if len(files) > 1:
        whichout = files[0]
        dt1 = dt.strptime(gettimestamp(whichout)[4:],"%b %d %H:%M:%S %Y")
        for o in files[1:]:
            dt2 = dt.strptime(gettimestamp(o)[4:],"%b %d %H:%M:%S %Y")
            if dt2 > dt1:
                #print(gettimestamp(o),'>',gettimestamp(whichout))
                whichout = o
                dt1 = dt2
        if remove == True:
            [os.remove(f) for f in files if not f == whichout]
        return whichout
    else:
        return files[0]
def formatclock(hrs):
    # format the time
    TIME = dt(1,1,1) + timedelta(hours=hrs)
    # zero out the minutes, add an hour
    if TIME.minute > 0:
        TIME = TIME + timedelta(hours=1) - timedelta(minutes=TIME.minute)
    # round up to 7 days, zero out the hours
    if 3 < (TIME.day-1) < 7 or (3 <= TIME.day-1 and TIME.hour > 0):
        diff = 7 - (TIME.day-1)
        TIME = TIME + timedelta(days=diff) - timedelta(hours=TIME.hour)
    # round up to 3 days, zero out the hours
    if 1 < (TIME.day-1) < 3 or (1 <= TIME.day-1 and TIME.hour > 0):
        diff = 3 - (TIME.day-1)
        TIME = TIME + timedelta(days=diff) - timedelta(hours=TIME.hour)
    # round up to 24 hrs, zero out the hours
    if TIME.day == 1 and 12 < TIME.hour < 24:
        TIME = TIME + timedelta(days=1) - timedelta(hours=TIME.hour)
    # round up to 12 hrs 
    if TIME.day == 1 and 3 < TIME.hour < 12:
        diff = 12 - TIME.hour
        TIME = TIME + timedelta(hours=diff)
    # round up to 3 hrs 
    if TIME.day == 1 and TIME.hour < 3:
        diff = 3 - TIME.hour
        TIME = TIME + timedelta(hours=diff) 
    clock = "%s-%s:00:00"  % (TIME.day -1, str(TIME.hour).zfill(2))
    return clock
def getpids():
    pids = os.popen('squeue -u lindb -o "%i"').read().split("\n")
    pids = [p for p in pids if not p == '']
    if len(pids) != luni(pids):
        print('len !- luni pids')
    return pids[1:]
def getjobs():
    jobs = os.popen('squeue -u lindb -o "%j"').read().split("\n")
    jobs = [j for j in jobs if not j == '']
    if len(jobs) != luni(jobs):
        print('len != luni jobs')
    return jobs[1:]
def getaccounts(pd=False):
    if pd == False:
        accounts = os.popen('squeue -u lindb -o "%a"').read().split("\n")
    else:
        accounts = os.popen('squeue -u lindb -t "pd" -o "%a"').read().split('\n')
    accounts = [a for a in accounts if not a in ['','ACCOUNT']]
    return accounts
def printmd(string):
    string = str(string)
    display(Markdown(string))
