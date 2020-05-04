"""
### usage
# usage: 00_start-gatk_pipeline.py -p PARENTDIR [-e EMAIL [-n EMAIL_OPTIONS]]
###

### fix
# dicts with samps as keys assume adaptors, ref, rglb, etc are same for reseqs
# (the second file with samp as key will overwrite the last)
# as of now, that hasn't created an issue
###

### TODO
# assert 0 <= input maf <= 1
###
"""

import os, sys, distutils.spawn, subprocess, shutil, argparse, pandas as pd
import balance_queue
from os import path as op
from coadaptree import fs, pkldump, uni, luni, makedir, askforinput, Bcolors


def create_sh(pooldirs, poolref, parentdir):
    # create sh files
    print(Bcolors.BOLD + '\nwriting sh files' + Bcolors.ENDC)
    for pooldir in pooldirs:
        pool = op.basename(pooldir)
        print(Bcolors.BOLD + '\npool = %s' % pool + Bcolors.ENDC)
        ref = poolref[pool]
        print('\tsending pooldir and ref to 01_trim-fastq.py')
        subprocess.call([shutil.which('python'),
                         op.join(os.environ['HOME'], 'gatk_pipeline/01_trim-fastq.py'),
                         pooldir,
                         ref])
    print("\n")
    balance_queue = op.join(os.environ['HOME'], 'gatk_pipeline/balance_queue.py')
    subprocess.call([sys.executable, balance_queue, 'trim', parentdir])


def get_datafiles(parentdir, f2pool, data):
    # get list of files from datatable, make sure they exist in parentdir, create symlinks in /parentdir/<pool_name>/
    print(Bcolors.BOLD + '\nchecking for existance of fastq files in datatable.txt' + Bcolors.ENDC)
    files = [f for f in fs(parentdir) if 'fastq' in f and 'md5' not in f]
    datafiles = data['file_name_r1'].tolist()
    for x in data['file_name_r2'].tolist():
        datafiles.append(x)
    if len(files) > len(datafiles):
        desc = 'more'
    if len(files) < len(datafiles):
        desc = 'less'
    try:
        print(Bcolors.WARNING +
              'WARN: there are %s fastq files in %s than in datatable.txt' % (desc, parentdir) +
              Bcolors.ENDC)
        print(Bcolors.BOLD + 'Here are the files in %s' % parentdir + Bcolors.ENDC)
        [print(op.basename(x)) for x in files]
        print(Bcolors.BOLD + 'Here are the files in datatable.txt' + Bcolors.ENDC)
        [print(x) for x in datafiles]
        askforinput()

    except NameError:
        pass

    for f in datafiles:
        src = op.join(parentdir, f)
        if not op.exists(src):
            # make sure file in datatable exists
            print("could not find %s in %s\nmake sure file_name in datatable is its basename" % (f, parentdir))
            print("(symlinks in parentdir to fastq files in other dirs works fine, and is the intentional use)")
            sys.exit(1)
        pooldir = op.join(parentdir, f2pool[f])
        dst = op.join(pooldir, f)
        if not op.exists(dst):
            # easy to visualize in cmdline if script is finding correct group of files by ls-ing pooldir
            os.symlink(src, dst)


def make_pooldirs(data, parentdir):
    # make pool dirs
    print(Bcolors.BOLD + "\nmaking pool dirs" + Bcolors.ENDC)
    pools = uni(data['pool_name'].tolist())
    pooldirs = []
    for p in pools:
        DIR = op.join(parentdir, p)
        if op.exists(DIR):
            print("The pooldir already exists, this could overwrite previous data: %s" % DIR)
            askforinput()
        pooldirs.append(makedir(DIR))
        makedir(op.join(DIR, 'shfiles'))
    return pooldirs


def handle_rg_fails(failing, warning, parentdir, data):
    if len(failing) > 0:
        print(Bcolors.FAIL + 'FAIL: The following samples have blank RG info.' + Bcolors.ENDC)
        for fail in failing:
            print(Bcolors.FAIL + "FAIL: %s" % fail + Bcolors.ENDC)
        print('exiting 00_start-pipeline.py')
        exit()
    if len(warning) > 0:
        outputs = []
        for row in data.index:
            samp = data.loc[row, 'sample_name']
            if samp in warning:
                r1 = op.join(parentdir, data.loc[row, 'file_name_r1'])
                outputs.append("\t\t%s\t%s" % (samp, get_rgid(r1)))
        print(Bcolors.WARNING + '\n\n\tWARN: at least one of the samples has a blank RGID in the datatable.\n' +
              '\tWARN: If RGPU is also blank, the pipeline will assign RGPU as: $RGID.$RGLB\n' +
              '\tWARN: The pipeline will automatically assign the following RGIDs.\n' +
              '\n\t\tsample_name\tassigned_RGID' +
              Bcolors.ENDC)
        for output in outputs:
            print(Bcolors.WARNING + output + Bcolors.ENDC)
        askforinput(tab='\t', newline='')


def read_datatable(parentdir):
    # read in the datatable, save info for later
    datatable = op.join(parentdir, 'datatable.txt')
    if not op.exists(datatable):
        print(Bcolors.FAIL + '''FAIL: the datatable is not in the necessary path: %s
FAIL: exiting 00_start-gatk_pipeline.py''' % datatable + Bcolors.ENDC)
        sys.exit(3)
    print(Bcolors.BOLD + 'reading datatable, getting fastq info' + Bcolors.ENDC)
    data = pd.read_csv(datatable, sep='\t')
    rginfo = {}     # key=sampname vals=rginfo
    samp2pool = {}  # key=samp val=pool
    poolref = {}    # key=pool val=ref.fa
    ploidy = {}     # key=pool val=ploidy
    poolsamps = {}  # key=pool val=sampnames
    f2samp = {}     # key=f val=samp
    f2pool = {}     # key=f val=pool
    adaptors = {}   # key=samp val={'r1','r2'} val=adaptor
    warning = []  # whether to print out warning about optional RG info
    failing = []  # whether to print out failing about required RG info
    for row in data.index:
        samp = data.loc[row, 'sample_name']
        adaptors[samp] = {'r1': data.loc[row, 'adaptor_1'],
                          'r2': data.loc[row, 'adaptor_2']}
        pool = data.loc[row, 'pool_name']
        pooldir = op.join(parentdir, pool)
        print('\t{}\tsamp = {}\tpool = {}'.format(row, samp, pool))
        if pool not in poolsamps:
            poolsamps[pool] = []
        if samp not in poolsamps[pool]:
            poolsamps[pool].append(samp)
        if samp in samp2pool:
            if samp2pool[samp] != pool:
                print(Bcolors.FAIL + 'FAIL: there are duplicate sample names with \
different pool assignments: %s' % samp + Bcolors.ENDC)
                print('exiting')
                exit()
        samp2pool[samp] = pool
        df = data[data['pool_name'] == pool].copy()
        if not luni(df['ploidy']) == 1:
            print(Bcolors.WARNING + 
                  "The ploidy values for some elements with pool name '%s' are not the same." % pool +
                  "\n\tHere are the ploidy values: %s" % uni(df['ploidy']) +
                  Bcolors.ENDC)
            askforinput()
        if samp not in ploidy:
            ploidy[samp] = data.loc[row, 'ploidy']
        if pool in poolref:
            if not poolref[pool] == data.loc[row, 'ref']:
                print("ref genome for samples in %s pool seems to have different paths in datatable.txt" % pool)
                sys.exit(1)
        else:
            ref = data.loc[row, 'ref']
            if not op.exists(ref):
                print('ref for %s does not exist in path: %s' % (samp, ref))
                print('exiting 00_start-gatk_pipeline.py')
                exit()
            needed = []
            for suffix in ['.dict', '.amb', '.ann', '.bwt', '.fai', '.pac', '.sa']:
                refext = ref + suffix if suffix != '.dict' else ref.split('.fa')[0] + suffix
                if not op.exists(refext):
                    needed.append(refext)
            if len(needed) > 0:
                print(Bcolors.FAIL + 
                      'FAIL: the following extensions of the reference are needed to continue, \
please create these files' + 
                      Bcolors.ENDC)
                for n in needed:
                    print(Bcolors.FAIL + n + Bcolors.ENDC)
                print('exiting')
                exit()
            printneeded = False
            intdir = op.join(op.dirname(ref), 'intervals')
            if not op.exists(intdir):
                printneeded = True
            elif len([f for f in fs(intdir) if '.list' in f and 'batch_' in f]) == 0:
                printneeded = True
            if printneeded is True:
                print(Bcolors.FAIL + 
                      'FAIL: either the intervals dir doesn not exist or there are not batch_interval.list files\
\nFAIL: intdir should be here: %s\
\nFAIL: interval filenames should be of the form "batch_uniqueIDENTIFIER.list"\
\nFAIL:      these files must begin with "batch" followed by one underscore followed\
\nFAIL:      by "uniqueIDENTIFIER" which can be any string that does not include\
\nFAIL:      an underscore (eg chrXIII or 0013 for batch_0013.list).' % intdir + Bcolors.ENDC)
                exit()
            poolref[pool] = ref
        rginfo[samp] = {}
        for col in ['rglb', 'rgpl', 'rgsm']:  # rg info columns
            rginfo[samp][col] = data.loc[row, col]
        for f in [data.loc[row, 'file_name_r1'], data.loc[row, 'file_name_r2']]:
            if "__" in f:
                print(Bcolors.BOLD + 
                      Bcolors.FAIL + 
                      "FAIL: file names cannot have double underscores, replace __ with _ (single)" + 
                      Bcolors.END)
                exit()
            f2pool[f] = pool
            f2samp[op.join(pooldir, f)] = samp

        # hangle RG info
        rginfo[samp] = {}
        # required RG info
        for col in ['rglb', 'rgpl', 'rgsm']:  # rg info columns
            if not data.loc[row, col] == data.loc[row, col]:
                failing.append('%s\t%s' % (samp, col))
            rginfo[samp][col] = data.loc[row, col]
        # optional RG info
        for col in ['rgid', 'rgpu']:
            if data.loc[row, col] != data.loc[row, col]:
                # if nan
                rginfo[samp][col] = None
                if samp not in warning:
                    warning.append(samp)
            else:
                rginfo[samp][col] = data.loc[row, col]

    # RG info failing/warnings
    handle_rg_fails(failing, warning, parentdir, data)

    pkldump(rginfo, op.join(parentdir, 'rginfo.pkl'))
    pkldump(ploidy, op.join(parentdir, 'ploidy.pkl'))
    pkldump(f2samp, op.join(parentdir, 'f2samp.pkl'))
    pkldump(poolsamps, op.join(parentdir, 'poolsamps.pkl'))
    pkldump(poolref, op.join(parentdir, 'poolref.pkl'))
    pkldump(adaptors, op.join(parentdir, 'adaptors.pkl'))
    pkldump(samp2pool, op.join(parentdir, 'samp2pool.pkl'))
    return data, f2pool, poolref


def check_reqs(parentdir):
    """Check for assumed exports."""
    
    print(Bcolors.BOLD + '\nChecking for exported variables' + Bcolors.ENDC)
    variables = ['SLURM_ACCOUNT', 'SBATCH_ACCOUNT', 'SALLOC_ACCOUNT', 'PYTHONPATH', 'SQUEUE_FORMAT']
    
    # check to see if bash_variables file has been created
    if not op.exists(op.join(parentdir, 'bash_variables')):
        print('\tCould not find bash_variables file in parentdir. Please create this file and add \
in variables from README (eg SLURM_ACCOUNT, SQUEUE_FORMAT, etc). See example in $HOME/gatk-pipeline.')
    else:
        with open(op.join(parentdir, 'bash_variables')) as bv:
            text = bv.read().split("\n")
        needed = []
        for var in variables:
            found = False
            for line in text:
                if var in line:
                    found = True
                    break
            if found is False:
                needed.append(var)
        if len(needed) > 0:
            print(Bcolors.FAIL + '\tFAIL: not all bash variables were found in parentdir/bash_variables file.' + Bcolors.ENDC)
            print(Bcolors.FAIL + '\tFAIL: the following variables must be present' + Bcolors.ENDC)
            for var in needed:
                print(Bcolors.FAIL + '\t%s' % var + Bcolors.ENDC)
            print('exiting pipeline')
    
    # check to see if bash_variables file has been sourced
    for var in variables:
        try:
            print('\t%s = %s' % (var, os.environ[var]))
        except KeyError:
            print(Bcolors.FAIL + '\tCould not find %s in exported vars\n\texport this var in parentdir/bash_variables \
so it can be used later in gatk_pipeline, then source this file before restarting pipeline.' % var + Bcolors.ENDC)
            print('\texiting 00_start-gatk_pipeline.py')
            exit()

    # make sure an environment can be activated (activation assumed to be in parentdir/bash_variables)
    for exe in ['activate']:
        if distutils.spawn.find_executable(exe) is None:
            print(Bcolors.FAIL + '\tcould not find %s in $PATH\nexiting 00_start-gatk_pipeline.py' % exe
                  + Bcolors.ENDC)
            if exe == 'activate':
                print(Bcolors.FAIL + '\t\t(the lack of activate means that the python env is not correctly installed)'
                     + Bcolors.ENDC)
            exit()
    # make sure pipeline can be accessed via $HOME/gatk_pipeline
    if not op.exists(op.join(os.environ['HOME'], 'gatk_pipeline')):
        print('\tcould not find gatk_pipeline via $HOME/gatk_pipeline\n\texiting 00_start-gatk_pipeline.py')
        exit()


def check_pyversion():
    # check python version
    pyversion = float(str(sys.version_info[0]) + '.' + str(sys.version_info[1]))
    if not pyversion >= 3.6:
        text = '''FAIL: You are using python %s. This pipeline was built with python 3.7+.
FAIL: You will need at least python v3.6+.
FAIL: exiting 00_start-gatk_pipeline.py
    ''' % pyversion
        print(Bcolors.BOLD + Bcolors.FAIL + text + Bcolors.ENDC)
        exit()


def get_pars():
    choices = ['all', 'fail', 'begin', 'end', 'pipeline-finish']
    parser = argparse.ArgumentParser(description=print(mytext),
                                     add_help=False,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    requiredNAMED = parser.add_argument_group('required arguments')
    requiredNAMED.add_argument("-p",
                               required=True,
                               default=argparse.SUPPRESS,
                               dest="parentdir",
                               type=str,
                               help="/path/to/directory/with/fastq.gz-files/")
    parser.add_argument("-e",
                        required=False,
                        dest="email",
                        help="the email address you would like to have notifications sent to")
    parser.add_argument("-n",
                        default=None,
                        nargs='+',
                        required=False,
                        dest="email_options",
                        help='''the type(s) of email notifications you would like to receive from the pipeline.\
                        Requires --email-address. These options are used to fill out the #SBATCH flags.
must be one (or multiple) of %s''' % [x for x in choices])
    parser.add_argument("-maf",
                        required=False,
                        default='0.05',
                        dest="maf",
                        help='''At the end of the pipeline, VCF files will be filtered for MAF. If the pipeline is run on a single population/pool, the user can set MAF to 0.0 so as to filter variants based on global allele frequency across populations/pools at a later time.''')
    parser.add_argument('-h', '--help',
                        action='help',
                        default=argparse.SUPPRESS,
                        help='Show this help message and exit.\n')
    args = parser.parse_args()
    if args.parentdir.endswith('/'):
        args.parentdir = args.parentdir[:-1]
    if args.email and args.email_options is None:
        print(Bcolors.FAIL + 'FAIL: --notification-types are required when specifying email' + Bcolors.ENDC)
        print(Bcolors.FAIL + 'FAIL: choices = {%s}\n' % [x for x in choices] + Bcolors.ENDC)
        exit()
    if args.email_options and args.email is None:
        print(Bcolors.FAIL + 'FAIL: specifying --notification-types requires specifying \
--email-address\n' + Bcolors.ENDC)
        exit()
    if args.email_options:
        for choice in args.email_options:
            if not choice.lower() in choices:
                print(Bcolors.FAIL +
                      '''FAIL: There can be multiple options, but they must be from the set:''' +
                      Bcolors.ENDC)
                print(Bcolors.FAIL +
                      '''\t%s\n''' % choices +
                      Bcolors.ENDC)
                exit()
    if args.email:
        if '@' not in args.email:
            print(Bcolors.FAIL + 'FAIL: email address does not have an "@" symbol in it, \
please check input\n' + Bcolors.ENDC)
            exit()
        if 'all' in args.email_options:
            args.email_options = ['all']
        # save email
        epkl = {'email': args.email,
                'opts': args.email_options}
        pkldump(epkl, op.join(args.parentdir, 'email_opts.pkl'))

    pkldump(args.maf, op.join(args.parentdir, 'maf.pkl'))

    return args


def main():
    # parse arguments
    args = get_pars()

    # WARN if version = 3.6, FAIL if < 3.6
    check_pyversion()

    # look for exported vars (should be in /parentdir/bash_variables)
    check_reqs(args.parentdir)
    
    # determine which slurm accounts to use
    balance_queue.get_avail_accounts(args.parentdir, save=True)

    # read in the datatable
    data, f2pool, poolref = read_datatable(args.parentdir)

#     # create bedfiles to parallelize later on
#     create_bedfiles(poolref)

    # create directories for each group of pools to be combined
    pooldirs = make_pooldirs(data, args.parentdir)

    # assign fq files to pooldirs for visualization (good to double check)
    get_datafiles(args.parentdir, f2pool, data)

    # create and sbatch sh files
    create_sh(pooldirs, poolref, args.parentdir)
    
    print(Bcolors.BOLD +
          Bcolors.OKGREEN +
          "\nDone with 00_start-gatk_pipeline.py" +
          Bcolors.ENDC)


if __name__ == '__main__':
    mytext = Bcolors.BOLD + Bcolors.OKGREEN + '''
*****************************************************************************


         ___|               \         |          _   __|
        |      _ \           \    __  |  _     _    |    _|  _ \\  _ \\
        |     (   | __|   /_  \  (    | (   | (  |  |   |    __/  __/
         ___|\___/      _/    _\\\___/_|\__/_|  __/  |  _|  \___|\___|
                                              |
                                              |

                                GATK pipeline

*****************************************************************************


''' + Bcolors.ENDC

    main()
