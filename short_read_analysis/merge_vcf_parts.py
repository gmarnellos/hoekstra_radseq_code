#!/usr/bin/env python

'''given as arguments vcf parts (vcf files containing sections of a single reference) and the reference from which they are derived
outputs (stdout) a single vcf

all parts must have identical field order (#CHROM line must be identical)

NOW ONLY WORKS ON GZIPPED FILES
'''

from short_read_analysis import map_reads_by_indiv_stampy
from subprocess import Popen
import os, sys, gzip

def smartopen(filename,*args,**kwargs):
    '''opens with open unless file ends in .gz, then use gzip.open

    in theory should transparently allow reading of files regardless of compression
    '''
    if filename.endswith('.gz'):
        return gzip.open(filename,*args,**kwargs)
    else:
        return open(filename,*args,**kwargs)

def head_and_start(vcf):
    fh = smartopen(vcf)
    while 1:
        l = fh.readline()
        if l == '':
            return None,None
        if l.startswith('#CHROM'):
            start = fh.readline()
            fh.close()
            return l,start

def order_vcf_parts(vcf_parts,ref):
    vcf_part_start = {}
    header = head_and_start(vcf_parts[0])[0]

    print >> sys.stderr, 'get start position for %s vcf parts' % (len(vcf_parts))

    for vcf in vcf_parts:
        head,start = head_and_start(vcf)
        if start:
            key = tuple(start.split()[:2])
            if vcf_part_start.has_key(key):
                raise ValueError, '%s starts at %s, but this start is already present in %s' % (vcf,key,vcf_part_start[key])
            if head != header:
                raise ValueError, 'header for %s does not match, expect:\n%s\nthis was\n%s' % (vcf,header,head)
            vcf_part_start[key] = vcf

    order = []

    print >> sys.stderr, 'get contig order from reference %s' % (ref)

    for l in smartopen(ref):
        if l.startswith('>'):
            order.append(l[1:].strip().split()[0])

    print >> sys.stderr, 'get sort order'
    ordered_vcf_parts = [vcf_part_start[key] for key in sorted(vcf_part_start.keys(), key = lambda x: (order.index(x[0]), int(x[1])) )]
    print >> sys.stderr, 'order parts done.'
    return ordered_vcf_parts

if __name__ == "__main__":

    import argparse

    ds =  ' [%(default)s]'
    #create command line parser
    parser = argparse.ArgumentParser(description='merge parts into single VCF')
    parser.add_argument('-v','--vcf_parts',nargs='+',help='VCF parts (specify any number of files)')
    parser.add_argument('-r','--reference_fasta',help='reference fasta for VCF')
    parser.add_argument('-o','--output',help='output filename for VCF')
    
    opts = parser.parse_args()

    ordered_parts = order_vcf_parts(opts.vcf_parts,opts.reference_fasta)

    print >> sys.stderr, 'merge %s vcfs, write output to %s' % (len(ordered_parts),opts.output)

    cmd = 'cp %s %s' % (ordered_parts[0],opts.output)
    ret = os.system(cmd)
    if ret != 0:
        raise OSError, 'cp failed'
    
    for i,f in enumerate(ordered_parts[1:],1):
        if i%10==0:
            print >> sys.stderr, '\r\t %s / %s' % (i,len(ordered_parts)),
        cmd = 'zcat %s | grep -v "^#" | gzip -c - >> %s' % (f,opts.output)
        ret = os.system(cmd)
        if ret != 0:
            errstr = 'vcf %s failed merge' % f
            raise OSError, errstr
    
    #outfh = smartopen(opts.output,'w')

    #for l in smartopen(ordered_parts[0]):
    #    outfh.write(l)
    #for i,f in enumerate(ordered_parts[1:],1):
    #    if i%10==0:
    #        print >> sys.stderr, '\r\t %s / %s' % (i,len(ordered_parts)),
    #    for l in smartopen(f):
    #        if not l.startswith('#'):
    #            outfh.write(l)
    #outfh.close()
    print >> sys.stderr, '\tdone.'
    

