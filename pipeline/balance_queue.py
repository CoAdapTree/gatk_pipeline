###
# usage: balance_queue.py phaseOFpipeline
###

###
# purpose: evenly redistributes jobs accross available accounts based on priority based on the job name (phaseOFpipeline)
#          helps speed up effective time
###

import sys
import math
import os

thisfile, phase = sys.argv

def checksq(rt,phase):
    exitneeded = False
    if not type(rt) == list:
        os.system('echo "type(sq) != list, exiting rescheduler.py"')
        exitneeded = True
    if len(rt) == 0:
        os.system('echo "len(sq) == 0, exiting rescheduler.py"')
        exitneeded = True
    for s in rt:
        if not s == '':
            if 'socket' in s.lower():
                os.system('echo "socket in sq return, exiting rescheduler.py"')
                exitneeded = True
            try:
                assert int(s.split()[0]) == float(s.split()[0])
            except:
                os.system('echo "could not assert int == float, %s %s"' % (s[0],s[0]))
                exitneeded = True
    if exitneeded == True:
        os.system('echo slurm screwed something up for balanace_queue.py %s, lame' % phase)
        exit()
    else:
        return rt
def getsq(phase):
    SQ = os.popen('''squeue -u lindb -t "PD" | grep %s | grep Priority''' % phase).read().split("\n")
    if len(SQ) > 0:
        return checksq(SQ,phase)
    else:
        os.system('echo no jobs in queue to balance')
        exit()
def adjustjob(acct, jobid):
    os.system('scontrol update Account=%s_cpu jobid=%s' % (acct, str(jobid)) )
def getaccounts(SQ):
    accounts = {}
    for q in SQ:
        if not q == '':
            splits = q.split()
            pid = splits[0]
            account = splits[2]
            account = account.split("_")[0]
            if not account in accounts:
                accounts[account] = {}
            accounts[account][pid] = splits
    if len(accounts) == 3: # all accounts have low priority
        os.system('echo all accounts have low priority, leaving queue as-is')
        exit()
    return accounts
def getbalance(accounts):
    sums = 0
    for account in accounts:
        sums += len(accounts[account].keys())
    bal = math.ceil(sums/3)
    return bal
def checknumaccts(accts,checking,mc):
    # len(accounts) will never == 2 after pop, since I checked for len(accounts) == 3
    if len(accts.keys()) == 0:
        if checking == 'RAC':
            os.system('echo all jobs on RAC')
        else:
            os.system('echo moved %s jobs to RAC' % str(mc))
        exit()
def redistribute4G(accounts,bal):
    RAC = 'rrg-yeaman'
    if RAC in accounts:   # no need to redistribute to RAC if RAC has low priority
        accounts.pop(RAC) # drop RAC from list to redistribute, exit if nothing to redistribute
        checknumaccts(accounts,'RAC')    # if all jobs are on RAC, exit
        return accounts
    accts = list(accounts.keys())
    for account in accts: 
        # distribute 4G jobs to RAC
        pids = list(accounts[account].keys())
        mcount = 0
        for pid in pids:
            mem = int([m for m in accounts[account][pid] if m.endswith('M')][0].split("M")[0])
            if mem <= 4000:
                # if it can be scheduled on the RAC, change the account of the jobid, and remove jobid from list
                adjustjob(RAC,pid)
                accounts[account].pop(pid)
                mcount += 1
                if mcount == bal:
                    break
        if len(accounts[account].keys()) == 0:
            accounts.pop(account)
    checknumaccts(accounts,'none',mcount) # if all jobs were redistributed to the RAC, exit
    return accounts
def gettaker(accounts):
    keys = list(accounts.keys())
    if len(keys) == 2:
        # if there are two accounts, figure out which account has more
        maxx = 0
        for acct in keys:
            if len(accounts[acct]) > maxx:
                giver = acct
                maxx = len(accounts[acct])
    else:
        assert len(keys) == 1
        giver = keys[0]
    taker = list(set(['def-saitken','def-yeaman']).symmetric_difference(set([giver])))[0]
    return giver,taker
def givetotaker(giver,taker,accounts,bal):
    taken = 0
    for pids in accounts[giver]:
        for pid in pids:
            adjustjob(taker,pid)
            taken += 1
            if taken == bal:
                os.system('echo redistributed %s jobs from %s to %s' % (str(taken),giver,taker))
                break
def announcefinal(accounts):
    os.system('echo final job announcement')
    for account in accounts:
        os.system('echo %s jobs on %s') % (str(len(accounts[account])),account)
def main(phase):
    # get the queue
    sq = getsq(phase)

    # get per-account counts of jobs in Priority pending status, exit if all accounts have low priority
    accts = getaccounts(sq)
    
    # figure out how many to balance remaining
    balance = getbalance(accts)

    # redistribute 4G jobs to RAC unless RAC has low priority, exit if all jobs redistributed or no jobs to redistribute
    accts = redistribute4G(accts,balance)
    
    # figure out which account to add to
    giver, taker = gettaker(accts)
    
    # redistribute to taker
    givetotaker(giver,taker,accts)
                              
    # announce final job counts
    announcefinal(getaccounts(getsq(phase)))
    
main(phase)
    
    
    
    
    