'''python wrappers for slurm job submission, status and batch control
'''

import os,sys,re,time,random,string,subprocess
from glob import glob

MAX_RETRY = 3

def random_filename(chars=string.hexdigits, length=8, prefix='', suffix='', \
                    verify=True, attempts=10, chosen=[]):
    for attempt in range(attempts):
        filename = ''.join([random.choice(chars) for i in range(length)])
        filename = prefix + filename + suffix
        if not verify or (not filename in chosen and not os.path.exists(filename)):
            return filename


def flatten_list(li):
    if not any([type(el) == list for el in li]):
        return li
    
    flat_li = []
    for this_el in li:
        if type(this_el) == list:
            this_el = flatten_list(this_el)
            flat_li.extend(this_el)
        else:
            flat_li.append(this_el)

    return flat_li

def slurm_script(cmd,jobname,scriptdir, runtime=1440,mem=4096,outdir=None,partition='general',**kwargs):
    '''takes one or more shell executable commands and produces a slurm srun/sbatch script
    if cmd is a list, will produce a multistep script accordingly
    runtime is max run duration IN MINUTES (e.g. 60 = 1h; 1440 = 1d; 10080 = 1w)
    mem is max process ram IN MBYTES (e.g. 4096 = 4G)'''

    if outdir is None:
        outdir = scriptdir

    if not os.path.exists(scriptdir): os.makedirs(scriptdir)
    if not os.path.exists(outdir): os.makedirs(outdir)

    if type(cmd) == list:
        #flatten_list() collapses lists potentially containing lists-of-lists into single-depth lists
        cmd = '\nsrun '.join(flatten_list(cmd))
    cmd = '\nsrun ' + cmd

    outstr = os.path.join(outdir,jobname)+'-%J-%N'

    ss_name = os.path.join(outdir,'%s.slurm.sh' % jobname)

    ss_head = '#!/bin/bash\n#SBATCH -J %s\n#SBATCH -n 1\n#SBATCH -t %s\n#SBATCH -p %s\n#SBATCH --mem-per-cpu=%s\n#SBATCH -o %s.out\n#SBATCH -e %s.err' % (ss_name,runtime,partition,mem,outstr,outstr)
    for k,v in kwargs.items():
        ss_head += '\n#SBATCH --%s=%s' % (k,v)
    ss_body = ss_head+cmd

    fh = open(ss_name,'w')
    fh.write(ss_body+'\n')
    fh.close()

    return ss_name
    
def jobs_submit(cmds,jobname_base,scriptdir, runtime=1440,mem=4096, num_batches=None, batch_size=None,partition='general',outdir=None ,force_source=False,**kwargs):
    '''takes a list of commands (which may themselves be lists; see list handling of cmd in slurm_script())
    writes slurm srun scripts and submits these, returning a jobsdict { <jobid> : <slurmscript> }
    kwargs is passed uninspected and unmolested to slurm_script, key/value pairs will be #SBATCH --<key>=<value> lines in resulting script
    '''

    if isinstance(cmds,str):
        cmds = [cmds]

    #print >> sys.stderr, '\n\n',num_batches,'\n\n'

    if num_batches:
        tot = len(cmds)
        batch_size = tot / num_batches
        print >> sys.stderr, 'batching invoked, %s batches requested (%s jobs per batch)' % (num_batches,batch_size)
        
    if batch_size:
        cmds = [cmds[i:i+batch_size] for i in xrange(0,len(cmds),batch_size)]

    jobsdict = {}
    chosen = []
    
    print >> sys.stderr,'Adding jobs'
    for cmd in cmds:
        if force_source:
            cmd = ['source ~/.bashrc',cmd]
            
        print >> sys.stderr,'.',
        
        time.sleep(0.1) # try not to stompy SLURM

        jobname = random_filename(prefix=jobname_base+'-',chosen=chosen)
        chosen.append(jobname)
        
        ss_name = slurm_script(cmd,jobname,scriptdir, runtime,mem,outdir,partition,**kwargs)

        match = None
        while match is None:
            ret = subprocess.Popen('sbatch %s' % (ss_name),shell=True,stdout=subprocess.PIPE).stdout.read()
            match = re.search('Submitted batch job (\d+)',ret)

        jid = match.groups()[0]

        jobsdict[jid] = ss_name

    return jobsdict
                                                   
def get_jobs_status(jobids=None,toplevel=True):
    '''returns status of the jobs indicated (jobsdict or list of job ids) or all jobs if no jobids supplied.
    Set toplevel=False for job step data
    '''
    if jobids is None:
        sacct_return = subprocess.Popen('sacct -p -l',shell=True,stdout=subprocess.PIPE).stdout.readlines()
    else:
        if type(jobids) == dict:
            qjobs = jobids.keys()
        else:
            qjobs = jobids        
        sacct_return = subprocess.Popen('sacct -j %s -p -l' % (','.join(qjobs),),shell=True,stdout=subprocess.PIPE).stdout.readlines()

    jobs_status = {}
    for el in sacct_return[1:]:
            d =  dict(zip(sacct_return[0].strip().split('|'),el.strip().split('|')))
            if not '.' in d['JobID'] or not toplevel: 
                jobs_status[d['JobID']] = d

    return jobs_status

def previous_submissions(scriptdir,runsafe_script,check_done=True):
    '''
    indended only for run_safe submissions (run_safe.py "runsafe_script" donefile)
    returns number of submission attempts as read from slurm.sh scripts in scriptdir
    '''
    if check_done:
        if os.path.exists(runsafe_script.rsplit('.',1)[0]+'.done'):
            return []
    prev_sub = subprocess.Popen('grep -l "%s" %s' % (runsafe_script,os.path.join(scriptdir,'*.slurm.sh')),shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE).stdout.readlines()
    return prev_sub
    

def get_status_counts(jobids=None):
    '''returns the counts of all jobs by status category
    '''
    from collections import defaultdict
    
    jobs_status = get_jobs_status(jobids)
    
    status_counts = defaultdict(int)
    for jd in jobs_status.values():
        status_counts[jd['State'].split()[0]] += 1

    return dict(status_counts)

def last_job_output(jobid,outdir,stdstream='out'):
    globstr = os.path.join(outdir,'*%s*.%s' % (jobid,stdstream))
    cand = glob(os.path.join(outdir,'*%s*.%s' % (jobid,stdstream)))
    if len(cand) == 1:
        try:
            return open(cand[0]).readlines()[-1]
        except:
            print >> sys.stderr, 'no output for %s' % cand[0]
            return ''
    else:
        print >> sys.stderr, 'unique output stream for %s not found: %s' % (globstr,cand)
        return None

def wait_for_jobs(jobsdict,restart_partition='general',sleeptime = 20,restart_z=None,restart_stragglers_after=0.75,kill_if_all_ssusp=False):
    '''loops checking status until no jobs are waiting or running / all are finished.
    wait/run states:

    CF  CONFIGURING     Job has been allocated resources, but are waiting for them to become ready for use (e.g. booting).
    CG  COMPLETING      Job is in the process of completing. Some processes on some nodes may still be active.
    PD  PENDING         Job is awaiting resource allocation.    
    R   RUNNING         Job currently has an allocation.
    RS  RESIZING        Job is about to change size.    
    S   SUSPENDED       Job has an allocation, but execution has been suspended.

    done states:

    CA  CANCELLED       Job was explicitly cancelled by the user or system administrator.  The job may or may not have been initiated.
    CD  COMPLETED       Job has terminated all processes on all nodes.
    F   FAILED          Job terminated with non-zero exit code or other failure condition.
    NF  NODE_FAIL       Job terminated due to failure of one or more allocated nodes.
    PR  PREEMPTED       Job terminated due to preemption.
    TO  TIMEOUT         Job terminated upon reaching its time limit.
    
    '''

    run_status = ['CONFIGURING','COMPLETING','PENDING','RUNNING','RESIZING','SUSPENDED']
    done_status = ['CANCELLED','COMPLETED','FAILED','NODE_FAIL','PREEMPTED','TIMEOUT']

    import time,datetime
    print >> sys.stderr, 'running %s jobs' % len(jobsdict)
    status = get_status_counts(jobsdict)
    t = time.time()
    maxllen = 0

    while any([k in run_status for k in status.keys()]):
        time.sleep(sleeptime)
        status = get_status_counts(jobsdict)
        pctdone = sum([status.get(rs,0) for rs in done_status]) / float(sum(status.values()))
        
        #CHECK SUSPENDED; RESTART STRAGGLERS, ETC

        outl = '\r%s %s (%3d%%)' % (str(datetime.timedelta(seconds=int(time.time() - t))),status.__repr__(),pctdone*100)
        if len(outl) < maxllen:
            pad = maxllen - len(outl)
            outl += ' '*pad
        else:
            maxllen = len(outl)
            
        sys.stderr.write(outl)
        sys.stderr.flush()

    print >> sys.stderr, '\ncompleted iteration in',str(datetime.timedelta(seconds=int(time.time() - t)))
            
def run_until_done(to_run_dict,jobname_base,scriptdir, runtime,mem, num_batches, partition='general' ,force_source=False,MAX_RETRY=MAX_RETRY,**kwargs):
    '''given to-run dictionary as populated by run_safe.add_cmd (see run_safe.py in py_util) and scheduling parameters
    submits jobs that have not yet completed per run_safe .done files until all jobs finish or until identical job lists are submitted MAX_RETRY times
    see jobs_submit and wait_for_jobs in this module for more details
    kwargs go to jobs_submit; see jobs_submit and slurm_script for handling of additional arguments
    '''
    from run_safe import unfinished_cmds
    cmds = unfinished_cmds(to_run_dict)
    
    retries = 0
    last_cmds = []
    while len(cmds) > 0:
        print >> sys.stderr, '%s: %s cmds to run in %s batches on queue %s, logs in %s' % (jobname_base,len(cmds),num_batches,partition,scriptdir)
        #code to halt execution on recurrent errors
        if set(last_cmds) == set(cmds):
            if retries > MAX_RETRY:
                errstr = 'maximum number of retry attempts (%s) exceeded with identical jobs lists.  Check logs (%s) for recurrent errors' % (MAX_RETRY,scriptdir)
                raise IOError, errstr
            else:
                retries += 1
        last_cmds = cmds
        
        jobsdict = jobs_submit(cmds,jobname_base,scriptdir, runtime,mem, num_batches,partition=partition ,force_source=force_source, **kwargs)
        time.sleep(20)
        wait_for_jobs(jobsdict,restart_partition=partition,sleeptime = 20)
        time.sleep(20)
        cmds = unfinished_cmds(to_run_dict)
    print >> sys.stderr, 'DONE\n'
                                                                                                                                            
