#!/usr/bin/env python
'''does no sorting; just merges bams and respects headers
usage:
merge_bams_with_headers_samtools.py out.bam in1.bam in2.bam in3.bam ...
'''
import os,sys,hashlib
from subprocess import Popen,PIPE
tmp_root = '/tmp'

def add_header_line(line,header):
    if line in header:
        return None
    if not any([l.startswith(line.split('\t')[0]) for l in header]) or header[-1].startswith(line.split('\t')[0]):
        header.append(line)
        return -1
    ins_idx = list(reversed([l.startswith(line.split('\t')[0]) for l in header])).index(True)
    header.insert(-ins_idx,line)
    return -ins_idx

def merge_header_from_bams(bams):
    header = []
    for b in bams:
        for l in Popen('samtools view -H  %s' % b,shell=True,stdout=PIPE).stdout:
            ins_idx = add_header_line(l,header)
    return header

if __name__ == "__main__":
    outbam = sys.argv[1]
    bams = sys.argv[2:]
    header_basename = '%s.txt' % hashlib.md5(' '.join(map( os.path.abspath, map(os.path.expanduser,bams)))).hexdigest()
    header_dir = os.path.join(tmp_root,os.environ['USER'],'bamheaders')
    header_file = os.path.join( header_dir,header_basename)

    try:
        os.makedirs(header_dir)
    except:
        pass

    header = merge_header_from_bams(bams)

    fh = open(header_file,'w')
    fh.writelines(header)
    fh.close()

    cmd = 'samtools merge -h %s %s %s' % (header_file, outbam, ' '.join(bams)) 
    print >> sys.stderr, cmd
    ret = os.system(cmd)
    print >> sys.stderr, 'DONE'

