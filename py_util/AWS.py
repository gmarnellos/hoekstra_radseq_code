'''
AWS.py

tools and objects for creating, managing and processing on amazon web services
(EC2 and S3) based "clusters"

Includes objects and methods intended for use on mpi ipengines, 
requires ipython 0.9 or higher
'''
from paths import paths
from IPython.kernel import client
import os, sys

def get_keys_from_file(filename=os.path.join(paths['home'],'.aws/aws_keys')):
    '''returns aws_keys dict from 'keyname:keyval,keyname:keyval' string
    in file filename'''

    return dict([i.split(':') for i in 
                 open(filename).read().split(',')])

def push_engine_furl(target_ip,
                     username=os.environ['USER'],
                     furl='~/.ipython/ipcontroller-engine.furl',
                     key_file='~/.ssh/gsg-keypair'):
    '''scp's the furl file pointing to the current local instance ipcontroller to target_ip
    '''
    

    scp = 'scp %s -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no %s %s%s:%s 2> /dev/null' % \
        ((key_file and '-i '+key_file or ''),
         furl,
         (username and username+'@' or ''),
         target_ip,furl)

    #print scp
    return os.system(scp)

def push_file(target_ip,sourcefile,targetfile,
              username=os.environ['USER'],
              key_file='~/.ssh/gsg-keypair'):
    '''scp's a file to target_ip
    '''

    scp = 'scp -r %s -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no %s %s%s:%s 2> /dev/null' % \
        ((key_file and '-i '+key_file or ''),
         sourcefile,
         (username and username+'@' or ''),
         target_ip,targetfile)

    return os.system(scp)

def pull_file(target_ip,sourcefile,targetfile,
              username=os.environ['USER'],
              key_file='~/.ssh/gsg-keypair'):
    '''scp's a file to target_ip
    '''

    scp = 'scp -r %s -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no %s%s:%s %s 2> /dev/null' % \
        ((key_file and '-i '+key_file or ''),
         
         (username and username+'@' or ''),
         target_ip,sourcefile,targetfile)

    return os.system(scp)
                              
def engines_import(module, mec=None):
    '''imports module (or more than one, separated by commas) on all engines
    '''
    if not mec:
        mec = client.MultiEngineClient()
    return mec.execute('import %s' % module)
    
def launch_engine(target_ip, username=os.environ['USER'], mec=None, key_file='~/.ssh/gsg-keypair'):
    '''launches an engine on target_ip
    '''
    if not mec:
        mec = client.MultiEngineClient()

#    ssh = r'ssh %s -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no %s%s ipengine 2\> /dev/null \> /dev/null \& 2> /dev/null' % \    
    ssh = r'ssh %s -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no %s%s ipengine \& &' % \
        ((key_file and '-i '+key_file or ''), 
         (username and username+'@' or ''), target_ip)
    return os.system(ssh)

def launch_bootstrap_engine(target_ip,mec=None,**kwargs):
    '''launches an engine on target_ip ONLY IF no engine is currently attached from there
    '''
    if not mec:
        mec = client.MultiEngineClient()
    try:
        ids_on_ip = engine_ids_by_ip(mec)[target_ip]
    except KeyError:
        pass
    else:
        mec.kill(targets=ids_on_ip)
    
    return launch_engine(target_ip,mec=mec,**kwargs)

def write_clusterconfig(engines_dict,controller_host=None,
                        base_clusterconfig=os.path.join(paths['home'],'py_util/base_clusterconfig-ec2.py'),
                        final_clusterconfig=os.path.join(paths['home'],'.ipython/full_clusterconfig-ec2.py')):
    '''takes an engine description dict, and a base configuration filename 
    (i.e. clusterconfig.py _without_ engines = {}), writes a full config to final_clusterconfig
    '''

    if controller_host is None:
        import net
        controller_host = net.get_ip_address('eth0')

    fh = open(final_clusterconfig,'w')
    fh.write( open(base_clusterconfig).read() )
    fh.write( 'controller = %s\n' % ({'host':controller_host}, ) )
    fh.write( 'engines = %s\n' % (engines_dict, ) )
    fh.close()


def write_runclusterconfig(host_list,controller_host=None,mec=None,eng_per_proc=1,**kwargs):
    '''takes a list of hostnames and number of engines per processor, 
    makes an engine description dict and calls write_clusterconfig
    '''
    if not mec:
        mec = client.MultiEngineClient()

    engines_import('os',mec)
        
    mec.execute("n = int(os.environ['NUMPROCS'])")
    engines_dict = dict(zip(host_list,[i*eng_per_proc for i in mec.gather("n")]))

    write_clusterconfig(engines_dict,controller_host,**kwargs)


def write_bstrapclusterconfig(host_list,controller_host=None,**kwargs):
    '''takes a list of hostnames, writes a clusterconfig for one engine per host
    (calls write_clusterconfig)
    '''
    engines_dict = dict.fromkeys(host_list,1)
    
    write_clusterconfig(engines_dict,controller_host,**kwargs)

    
def engine_ids_by_ip(mec=None,interface='eth0'):
    if not mec:
        mec = client.MultiEngineClient()
    if not mec.get_ids():
        return {}
    engines_import('net',mec)
    mec.execute("ip = net.get_ip_address('%s')" % interface)
    ids_by_ip = {}
    for id in mec.get_ids():
        [ip] = mec.pull('ip',[id])
        try:
            ids_by_ip[ip].append(id)
        except KeyError:
            ids_by_ip[ip] = [id]
    return ids_by_ip

def fill_instance_engines(mec=None,eng_per_proc=2,**kwargs):
    if not mec:
        mec = client.MultiEngineClient()

    engines_import('os',mec)

    ids_by_ip = engine_ids_by_ip(mec)
    for ip,ids in ids_by_ip.items():
        mec.execute("p = int(os.environ['NUMPROCS'])")
        [procs] = mec.pull('p',ids[0])
        for i in range( (procs*eng_per_proc) - len(ids) ):
            launch_engine(ip,mec=mec,**kwargs)

def scatter_and_run(queue,sleeptime=60,mec=None,verbose=False):
    '''not tested, probably not finished :)
    '''
    from time import sleep
    if not mec:
        mec = client.MultiEngineClient()
    mec.scatter('q',queue)

    mec.execute('res = AWS.run_queue(q)',block=False)
    
    while any([i[1]['pending'] != 'None' for i in mec.queue_status()]):
        done = 0
        for i in mec.get_ids():
            if verbose: print >> sys.stderr, 'Engine %s:' % i,
            [status] = mec.pull(['on','tot'],[i])
            if verbose: print >> sys.stderr, '[%s/%s]' % tuple(status)
            done += status[0]-1

        if verbose: print >> sys.stderr,'total progress: %s of %s' % (done, len(queue))
        sleep(sleeptime)

    return mec.gather('res')


def get_instance_list(ec2_conn):
    '''given a valid boto ec2 connection object
    returns a list of all current instance objects
    '''
    instancelist = []
    for r in ec2_conn.get_all_instances():
        instancelist.extend(r.instances)
    return instancelist

def engines_get_from_bucket(bucket,mec=None,filestring='*',dir_root='/mnt'):
    if not mec:
        mec = client.MultiEngineClient()

    engines_import('boto',mec)

    raise NotImplementedError, 'unfinished!'

def engines_get_from_bucket_by_s3fs(bucket,mec=None,filestring='*',mount_root='/mnt/s3fs',target_root='/mnt',gz=False):
    '''uses s3fs to mount bucket at mount_root/<bucket>, 
    creates directory target_root/<bucket>,
    does cp -R mount_root/<bucket>/filestring target_root/<bucket>
    if gz does gunzip target_root/<bucket>/*.gz

    NB: should be run with bootstrap engines only (i.e. 1 per host)
    '''

    if not mec:
        mec = client.MultiEngineClient()

    engines_import('os',mec)

    target = os.path.join(target_root,bucket)
    mount = os.path.join(mount_root,bucket)

    makeTdirs = "os.makedirs('%s')" % target
    print makeTdirs
    makeMdirs = "os.makedirs('%s')" % mount
    print makeMdirs


    mountbucket = "os.system('s3fs %s %s')" % (bucket,mount)
    print mountbucket
    

    cpsource = os.path.join(mount,filestring)
    copy = "os.system('cp -R %s %s')" % (cpsource,target)
    print copy

    mec.execute(makeTdirs)
    mec.execute(makeMdirs)
    mec.execute(mountbucket)
    mec.execute(copy)

    if gz:
        zipfiles = os.path.join(target,'*.gz')
        unzip = "os.system('gunzip %s')" % zipfiles
        print unzip

        mec.execute(unzip)


def put_to_bucket(filestring,bucket,s3_conn):
    try:
        b = s3_conn.get_bucket(bucket)
    except S3ResponseError:
        b = s3_conn.create_bucket(bucket)
    filebase = os.path.basename(filestring)
    k = b.new_key(filebase)
    k.set_contents_from_filename(filestring)

def put_all_to_bucket(bucket,s3_conn,directory):
    try:
        b = s3_conn.get_bucket(bucket)
    except S3ResponseError:
        b = s3_conn.create_bucket(bucket)

    for f in os.listdir(directory):
        filebase = os.path.basename(f)
        k = b.new_key(filebase)
        k.set_contents_from_filename(filestring)

def get_from_bucket(keyname,bucket,s3_conn,basedir):
    if not os.path.exists(basedir): os.makedirs(basedir)
    b = s3_conn.get_bucket(bucket)
    targetfile = os.path.join(basedir,bucket,keyname)
    k = b.get_key(keyname)
    k.get_contents_to_filename(targetfile)

def get_all_from_bucket(bucket,s3_conn,basedir):
    if not os.path.exists(basedir): os.makedirs(basedir)
    b = s3_conn.get_bucket(bucket)
    for k in b.get_all_keys():
        targetfile = os.path.join(basedir,bucket,k.name)
        k.get_contents_to_filename(targetfile)

def engines_put_to_bucket_by_s3fs(filestring,mount_root='/mnt/s3fs',mec=None,s3_conn=None,bucket=None,gz=True):
    '''transfers files in filestring (can be glob) to bucket
    if bucket is none, attempts to transfer to a bucket named the last dir in the path of filestring
    creates bucket if it doesn't exist
    '''
    import boto
    
    if not mec:
        mec = client.MultiEngineClient()

    if not bucket:
        bucket = filestring.split('/')[-2]

    print 'bucket:',bucket

    if not s3_conn:
        aws_keys = get_keys_from_file()
        s3_conn = boto.connect_s3(**aws_keys)

    if not bucket in [b.name for b in s3_conn.get_all_buckets()]:
        print 'create bucket:',bucket
        s3_conn.create_bucket(bucket)

    engines_import('os',mec)

    if gz:
        zip = "os.system('gzip %s')" % filestring
        print zip
        filestring = filestring+'.gz'

        mec.execute(zip)

    mount = os.path.join(mount_root,bucket)
    makeMdirs = "os.makedirs('%s')" % mount
    print makeMdirs


    mountbucket = "os.system('s3fs %s %s')" % (bucket,mount)
    print mountbucket
    

    copy = "os.system('cp %s %s')" % (filestring,mount)
    print copy

    mec.execute(makeMdirs)
    mec.execute(mountbucket)
    mec.execute(copy)

    

####ENGINE TOOLS

def run_queue(q):
    '''runs a list-of-lists queue, 
    where every item in each inner list is executed in list order (ascending index)
    for now, items in queues are controllers (or anything else with a .run() method that has a "now" option
    (see Util.Controller for an example)
    '''
    
    res = []
    tot = len(q)
    print 'total jobs:',tot
    for on,subq in enumerate(q):
        res.append([])
        for j,job in enumerate(subq):
            print 'running',job['name']
            res[on].append(job.run(1))

    return res
