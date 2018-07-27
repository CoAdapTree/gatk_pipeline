import numpy as np
import pandas as pd
from collections import OrderedDict, Counter
import os
from os import path as op
from os import listdir as ls
import pickle
import random
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import string
from os import chdir as cd
from os import getcwd as cwd
from os import makedirs as mkdir

import math
import shutil
from shutil import copy as cp
from shutil import move as mv

from matplotlib import pyplot as pl
import matplotlib.dates as mdates

from ipyparallel import Client
from decimal import Decimal
from datetime import datetime as dt
from scipy import stats
interval = stats.norm.interval


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
