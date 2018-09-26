import sys
import os
import pickle
import random
import string
import math
import shutil
import numpy as np
import pandas as pd
from IPython.display import clear_output
from collections import OrderedDict, Counter
from os import path as op
from os import chdir as cd
from os import getcwd as cwd
from os import makedirs as mkdir
from os import listdir
from shutil import copy as cp
from shutil import move as mv
from ipyparallel import Client
from decimal import Decimal
from datetime import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as pl
import matplotlib.dates as mdates


def ls(DIR):
    return sorted([f for f in listdir(DIR)])
def fs(DIR):
    return (sorted([op.join(DIR,f) for f in ls(DIR)]))
def uni(mylist):
    return (np.unique(mylist).tolist())
def luni(mylist):
    return (len(uni(mylist)))
def suni(mylist):
    return(sorted(uni(mylist)))
def nrow(df):
    return len(df.index)
def ncol(df):
    return len(df.columns)
def table(lst):
    c = Counter()
    for x in lst:
        c[x] += 1
    return(c)
def pkldump(obj,f):
    with open(f,'wb') as o:
        pickle.dump(obj,o,protocol=pickle.HIGHEST_PROTOCOL)
def head(df):
    return df.head()
def sbatch(files):
    for f in files:
        os.chdir(op.dirname(f))
        os.system('sbatch %s' % f)
def update(args):
    clear_output(wait=True)
    print(args)
