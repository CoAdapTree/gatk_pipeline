"""
### purpose
# make sure these modules can be loaded before starting pipline
# this is a list of modules across .py scripts
# also a source for functions common across .py scripts
###

### fix
# export PYTHONPATH in .sh files to import these functions
###
"""

import os
# import sys # commented imports are not used in these funcs below, but are used in other python apps in the pipeline
# import json
# import math
# import time
import pickle
# import random
# import compiler
# import numpy as np
# import pandas as pd
from os import path as op
# from collections import OrderedDict, Counter


def fs(directory):
    return sorted([op.join(directory, f) for f in os.listdir(directory)])


def pkldump(obj, f):
    with open(f, 'wb') as o:
        pickle.dump(obj, o, protocol=pickle.HIGHEST_PROTOCOL)


def pklload(path):
    pkl = pickle.load(open(path, 'rb'))
    return pkl


def get_email_info(parentdir, stage):
    pkl = op.join(parentdir, 'email_opts.pkl')
    if op.exists(pkl):
        email_info = pklload(pkl)
        # email_info = {'email': 'lindb@vcu.edu', 'opts': ['pipeline-finish']}  # for testing

        # make text
        email_text = '''#SBATCH --mail-user=%s''' % email_info['email']
        options = [opt.upper() for opt in email_info['opts'] if opt != 'pipeline-finish']
        # first determine if it's only when the pipeline finishes
        if email_info['opts'] == ['pipeline-finish'] and stage != 'final':
            # if default opt, but it's not the final stage
            return ''
        elif 'pipeline-finish' in email_info['opts'] and stage == 'final':
            email_text = email_text + '\n' + "#SBATCH --mail-type=END"
            if 'END' in options:
                options.remove('END')
        # now for stages earlier than final
        for opt in options:
            email_text = email_text + '\n#SBATCH --mail-type=%s' % opt
        return email_text
    else:
        # no email options
        return ''


def uni(mylist):
    return list(set(mylist))


def luni(mylist):
    return len(uni(mylist))


def makedir(directory):
    if not op.exists(directory):
        os.makedirs(directory)
    return directory


def createdirs(dirs):
    for d in dirs:
        makedir(d)
