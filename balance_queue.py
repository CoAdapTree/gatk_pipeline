"""
Distribute priority jobs among accounts.

###
# usage: balance_queue.py phaseOFpipeline
###

###
# purpose: evenly redistributes jobs across available accounts based on priority based on the job name
#          (phaseOFpipeline);
#          helps speed up effective run time
###
"""

import os, shutil, sys, math, subprocess, time


def announceacctlens(accounts, fin):
    """How many priority jobs does each account have?
    
    Positional arguments:
    accounts - dictionary with key = account_name, val = list of jobs (squeue output)
    fin - True if this is the final job announcement, otherwise the first announcement
    """
    print('%s job announcement' % ('final' if fin is True else 'first'))
    if fin is True:
        time.sleep(1)
    for account in accounts:
        print('%s jobs with Priority status on %s' % (str(len(accounts[account])), account))


def checksq(sq):
    """Make sure queue slurm command worked.
    
    Positional arguments:
    sq - list of squeue slurm command jobs, each line is str.split()
       - slurm_job_id is zeroth element of str.split()
    """
    exitneeded = False
    if not isinstance(sq, list):
        print("type(sq) != list, exiting %(thisfile)s" % globals())
        exitneeded = True
    for s in sq:
        if 'socket' in s.lower():
            print("socket in sq return, exiting %(thisfile)s" % globals())
            exitneeded = True
        if not int(s.split()[0]) == float(s.split()[0]):
            print("could not assert int == float, %s" % (s[0]))
            exitneeded = True
    if exitneeded is True:
        print('slurm screwed something up for %(thisfile)s, lame' % globals())
        exit()
    else:
        return sq


def getsq_exit(balancing):
    """Determine if getsq is being used to balance priority jobs.

    Positional arguments:
    balancing - bool: True if using to balance priority jobs, else for other queue queries
    """
    print('no jobs in queue matching query')
    if balancing is True:
        print('exiting balance_queue.py')
        exit()
    else:
        return []


def getsq(grepping=None, states=[], balancing=False):
    """
    Get jobs from squeue slurm command matching crieteria.

    Positional arguments:
    grepping - list of key words to look for in each column of job info
    states - list of states {pending, running} wanted in squeue jobs
    balancing - bool: True if using to balance priority jobs, else for other queue queries

    Returns:
    grepped - list of tuples where tuple elements are line.split() for each line of squeue \
slurm command that matched grepping queries
    """
    if grepping is None:
        grepping = [os.environ['USER']]
    if isinstance(grepping, str):
        # in case I pass a single str instead of a list of strings
        grepping = [grepping]

    # get the queue, without a header
    cmd = [shutil.which('squeue'),
           '-u',
           os.environ['USER'],
           '-h']
    if 'running' in states:
        cmd.extend(['-t', 'RUNNING'])
    elif 'pending' in states:
        cmd.extend(['-t', 'PD'])
    sqout = subprocess.check_output(cmd).decode('utf-8').split('\n')

    sq = [s for s in sqout if not s == '']
    checksq(sq)  # make sure slurm gave me something useful

    # look for the things I want to grep
    grepped = []
    if len(sq) > 0:
        for q in sq:  # for each job in queue
            splits = q.split()
            if 'CG' not in splits:  # grep -v 'CG'
                keepit = 0
                if len(grepping) > 0:  # see if all necessary greps are in the job
                    for grep in grepping:
                        for split in splits:
                            if grep.lower() in split.lower():
                                keepit += 1
                                break
                if keepit == len(grepping) and len(grepping) != 0:
                    grepped.append(tuple(splits))

        if len(grepped) > 0:
            return grepped
    return getsq_exit(balancing)


def adjustjob(acct, jobid):
    """Move job from one account to another."""
    subprocess.Popen([shutil.which('scontrol'),
                      'update',
                      'Account=%s_cpu' % acct,
                      'JobId=%s' % str(jobid)])


def getaccounts(sq, stage):
    """
Count the number of priority jobs assigned to each account.

Positional arguments:
sq - list of squeue slurm command jobs, each line is str.split()
   - slurm_job_id is zeroth element of str.split()
stage - stage of pipeline, used as keyword to filter jobs in queue
    """
    accounts = {}
    for q in sq:
        pid = q[0]
        account = q[2].split("_")[0]
        account = account.split("_")[0]
        if account not in accounts:
            accounts[account] = {}
        accounts[account][pid] = q
#     if len(accounts.keys()) == 3 and stage != 'final': # all accounts have low priority ### use 3 when using RAC
    if len(accounts.keys()) == 2 and stage != 'final':  # all accounts have low priority   ### use 2 when not using RAC
        print('all accounts have low priority, leaving queue as-is')
        announceacctlens(accounts, True)
        exit()
    return accounts


def getbalance(accounts, num):
    """Determine how many jobs should be given from one account to another.

    Positional arguments:
    accounts - dictionary with key = account_name, val = list of jobs (squeue output)
    num - number of accounts to balance among (this needs to be changed to object not number)
    """
    sums = 0
    for account in accounts:
        sums += len(accounts[account].keys())
    bal = math.ceil(sums/num)
    print('bal%i %i= ' % (num, bal))
    return bal


# def checknumaccts(accts, checking, mc):
#     # len(accounts) will never == 2 after pop, since I checked for len(accounts) == 3
#     if len(accts.keys()) == 0:
#         if checking == 'RAC':
#             print('RAC has low priority status, skipping RAC as taker')
#         else:
#             print('moved %s jobs to RAC' % str(mc))
#         exit()
#
#
# def redistribute4g(accounts, bal, rac, mcount=0):
#     if rac in accounts:    # no need to redistribute to rac if rac has low priority
#         accounts.pop(rac)  # drop rac from list to redistribute, exit if nothing to redistribute
#         checknumaccts(accounts, 'rac', '')    # if all jobs are on rac, exit
#         return accounts
#     keys = list(accounts.keys())
#     # print('before loop %s' % keys)
#     for account in keys:
#         # distribute 4G jobs to rac
#         pids = list(accounts[account].keys())
#         mcount = 0
#         for pid in pids:
#             mem = int([m for m in accounts[account][pid] if m.endswith('M')][0].split("M")[0])
#             if mem <= 4000:
#                 # if it can be scheduled on the rac, change the account of the jobid, and remove jobid from list
#                 adjustjob(rac, pid)
#                 accounts[account].pop(pid)
#                 mcount += 1
#                 if mcount == bal:
#                     break
#         print("distributed {} jobs from {} to rac".format(mcount, account))
#         if len(accounts[account].keys()) == 0:
#             accounts.pop(account)
#     checknumaccts(accounts, 'none', mcount)  # if all jobs were redistributed to the rac, exit
#     return accounts


def gettaker(accounts, defs):
    """Determine which job should receive jobs from the one with priority jobs.

    Positional arguments:
    accounts - dictionary with key = account_name, val = list of jobs (squeue output)
    defs - default account names to balance among. (TODO: will 'work' if defs > 2, but won't evenly dist.)
    """
    giver = ''
    keys = list(accounts.keys())
    if len(keys) > 1:
        # if there are at least two accounts, figure out which account has more (assign to final giver if tie)
        maxx = 0
        for acct in keys:
            if len(accounts[acct]) > maxx:
                giver = acct
                maxx = len(accounts[acct])
    else:
        if not len(keys) == 1:
            print('assertion error')
        giver = keys[0]
    # taker = list({'def-saitken', 'def-yeaman'}.symmetric_difference({giver}))[0]
    taker = list(set(defs).symmetric_difference({giver}))[0]  # TODO: will 'work' if defs > 2, but won't evenly dist.
    return giver, taker


def givetotaker(giver, taker, accounts, bal):
    """Give jobs to the account without jobs with priority status.

    Positional arguments:
    giver - account giving jobs to taker
    taker - account receiving jobs from giver
    accounts - dictionary with key = account_name, val = list of jobs (squeue output)
    """
    taken = 0
    pids = list(accounts[giver].keys())
    numtotake = len(pids) - bal
    if bal == 1 and len(pids) == 1:
        numtotake = 1
    printout = 'giver has {} jobs to give. (bal= {}). Giver ({}) is giving {} jobs to taker ({})'.format(len(pids),
                                                                                                         bal,
                                                                                                         giver,
                                                                                                         numtotake,
                                                                                                         taker)
    print("\t %s" % printout)
    if numtotake > 0:
        for pid in pids[::-1]:  # re-assign the newer jobs, hopefully older jobs will eventually run
            adjustjob(taker, pid)
            taken += 1
            if taken == numtotake:
                print("\t redistributed %s jobs from %s to %s" % (str(taken), giver, taker))
                break
    else:
        print("\t giver sees that taker has enough, so giver is not giving")


def get_availaccounts():
    """Query slurm with sshare command to determine accounts available."""
    acctout = subprocess.check_output([shutil.which('sshare'),
                                       '-U',
                                       '--user',
                                       os.environ['USER'],
                                       '--format=Account']).decode('utf-8').split('\n')
    accts = [acct.split()[0].split("_")[0] for acct in acctout if 'def' and 'cpu' in acct]
    defs = [acct for acct in accts if 'def' in acct]
    rac = [acct for acct in accts if 'rrg' in acct]
    if len(defs) == 1:
        # no need to try and balance
        print('there is only one account (%s), no more accounts to f queue.\nexiting balance_queue.py' % defs[0])
        exit()
    if len(rac) == 1:
        rac = rac[0]
    elif len(rac) == 0:
        rac = ''
    return defs, rac


def main(thisfile, phase):
    globals().update({'thisfile': thisfile, 'phase': phase})

    # get accounts available for billing
    defs, rac = get_availaccounts()

    # get the queue
    sq = getsq(grepping=[phase, 'Priority'], balancing=True)

    # get per-account counts of jobs in Priority pending status, exit if all accounts have low priority
    accts = getaccounts(sq, '')
    announceacctlens(accts, False)

    # figure out how many to balance remaining
    # balance = getbalance(accts, 3)

    # redistribute 4G jobs to RAC unless RAC has low priority, exit if all jobs redistributed or no jobs to redistribute
    # accts = redistribute4g(accts, balance, rac)

    # figure out which account to add to
    giver, taker = gettaker(accts, defs)

    # redistribute to taker
    balance = getbalance(accts, 2)
    givetotaker(giver, taker, accts, balance)

    # announce final job counts
    announceacctlens(getaccounts(getsq(grepping=[phase, 'Priority'], balancing=True),
                                 'final'),
                     True)


if __name__ == '__main__':
    # args
    thisfile, phase = sys.argv

    main(thisfile, phase)
