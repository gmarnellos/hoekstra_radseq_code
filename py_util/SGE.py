import os, sys, subprocess, re

import random, string
def random_filename(chars=string.hexdigits, length=8, prefix='', suffix='', \
                        verify=True, attempts=10):
    for attempt in range(attempts):
        filename = ''.join([random.choice(chars) for i in range(length)])
        filename = prefix + filename + suffix
        if not verify or not os.path.exists(filename):
            return filename

class Qopen ():
    '''analogous to subprocess.Popen
    when created, takes an argument list, as well as stdin, stdout, stderr and cwd.

    immediately submits to qsub, capturing qprocess id to return from poll()
    '''

    qsub = 'qsub'

    def __init__(self,args,stdin=None,stdout=None,stderr=None,cwd=None,name=None):
        '''basically does all the work:
        kicks off the process, gets qsub pid and listens for qsub error
        '''
        
        if cwd:
            self.cwd = cwd
        else:
            self.cwd = os.getcwd()

        if name:
            self.name = name
        else:
            self.name = random_filename(prefix='qsub-')

        (self.stdin, self.stdout, self.stderr) = (stdin,stdout,stderr)
        
        self.runargs = args
        self.qsubargs = self.compose_arguments()
        self.pid, self.err = subprocess.Popen(self.qsubargs,stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=self.cwd).communicate()
        self.pid = self.pid.strip()
        
        if self.err:
            sys.stderr.write('%s\n[%s]\n%s\n' % (' '.join(args), self.pid, self.err))

    def compose_arguments(self):
        args = [self.qsub]
        args.extend(['-d',self.cwd])
        args.append(self.write_sh())
        return args

    def write_sh(self):
        if self.stdin:
            if isinstance(self.stdin,file):
                self.stdin = self.stdin.name
            self.runargs.extend(['<', self.stdin])
        if self.stdout:
            if isinstance(self.stdout,file):
                self.stdout = self.stdout.name
            self.runargs.extend(['>', self.stdout])
        if self.stderr:
            if isinstance(self.stderr,file):
                self.stderr = self.stderr.name
            self.runargs.extend(['2>', self.stderr])

        #now write
        sh = os.path.join(self.cwd,self.name+'.sh')
        fh = open(sh,'w')
        fh.write(' '.join(self.runargs))
        fh.close()
        os.chmod(sh,0755)
        return sh

    def poll(self):
        '''to mimic subprocess.Popen, basically checks to see if the job is _NOT_ on the queue
        in which case returns 0 (not none)
        '''

        if subprocess.Popen(['qstat',self.pid],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0]:
            return None
        else:
            return 0

    def wait(self):
        raise NotImplementedError, 'qsub operations resist urgency...'

        
