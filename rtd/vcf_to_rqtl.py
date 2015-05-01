#!/usr/bin/env python
'''invocation:
vcf_to_rqtl.py path_to_vcf.vcf "P" Q G
where:
P = comma-separated pair of prefixes identifying cross founder species/strains
    for instance, if founding parents of strain "BW" are BW01, BW02, BW03
    and founding parents of strain "PO" are PO01,PO02,PO03,
    this paramter would be "BW,PO"
Q = GATK "QD" score threshold for inclusion
G = individual genotype call quality threshold

'''
# functionality from short_read now copied here
#from short_read_analysis import variant_detection,extract_genotypes_from_mclgr

import os,sys,re,numpy
from collections import defaultdict

def hwe(indiv_gt):
    '''given an indiv_gt dict as in vcf_data

    returns the observed and expected classes for hardy-weinberg equilibrium'''

    alleles = set(reduce(lambda x,y:x+y,[re.split('[/|]',v['GT']) for v in indiv_gt.values()]))
    if len(alleles) == 1:
        return 0
    if len(alleles) > 2:
        return None

    Aid, Bid = sorted(list(alleles))

    A = 0
    a = 0

    AA = 0
    Aa = 0
    aa = 0

    for v in indiv_gt.values():
        if set(re.split('[/|]',v['GT'])) == set([Aid]):
            AA += 1
            A += 2
        elif set(re.split('[/|]',v['GT'])) == set([Aid,Bid]):
            Aa += 1
            A += 1
            a += 1
        elif set(re.split('[/|]',v['GT'])) == set([Bid]):
            aa += 1
            a += 2
    try:
        p = A/float(A+a)
    except ZeroDivisionError:
        print >> sys.stderr, indiv_gt
    q = 1-p

    tot = AA+Aa+aa

    expAA = p**2 * tot
    expAa = 2*p*q*tot
    expaa = q**2 *tot

    X2 = sum([(o-e)**2/e for o,e in zip([AA,Aa,aa],[expAA,expAa,expaa]) if e != 0])

    return X2
        
def maf(indiv_gt,alt_as_minor=False):
    '''returns the counts of the minor allele for an indiv_gt dict as in vcf_data above'''

    A = 0
    a = 0

    for v in indiv_gt.values():
        A += v['GT'].count('0')
        a += v['GT'].count('1')

    if alt_as_minor:
        return a
    else:
        return min(A,a)

def fract_het(indiv_gt):
    gts = [v['GT'] for v in indiv_gt.values()]

    return (gts.count('0/1') + gts.count('0|1') + gts.count('1|0'))/float(len(gts))

def load_vcf_header(vcfh):
    '''given an open VCF file handle (that has not yet advanced past the header line)
    returns the headers, number of expected elements, and position of FORMAT field'''

    while 1:
        line = vcfh.readline()
        if line.startswith('#CHROM'):
            headers = line[1:].split()
            exp_elements = len(line.split())
            FORMAT = headers.index('FORMAT')
            return (headers, exp_elements, FORMAT)
    raise IOError, 'header not found; try zeroing vcfh position'

def load_vcf_line(vcfh,headers,exp_elements,FORMAT,indiv_gt_phred_cut=None,store_indiv_prefix=None,drop_indiv=None,biallelic_sites=True,skip_noGQ=True):
    line = vcfh.readline()
    if len(line) == 0:
        return 'EOF'
    
    fields = line.split()
    if len(fields) != exp_elements:
        print >>sys.stderr, 'unexpected length, line %s (exp %s obs %s)' % (i,exp_elements,len(fields))
        return None

    sd = dict(zip(headers[:FORMAT],fields[:FORMAT]))
    key = (sd['CHROM'],sd['POS'])
    try:
        infostr = sd.pop('INFO')
        sd.update(dict([el.split('=') for el in infostr.split(';') if '=' in el]))
    except KeyError:
        pass

    if biallelic_sites is True:
        if ',' in sd['ALT']:
            return None

    sd['indiv_gt'] = {}
    for ind,gt in zip(headers[FORMAT+1:],fields[FORMAT+1:]):
        if store_indiv_prefix is not None and not ind.startswith(store_indiv_prefix): continue
        #drop_indiv HERE
        if drop_indiv is not None and ind in drop_indiv: continue
        if ':' in gt:
            try:
                this_gt = dict(zip(fields[FORMAT].split(':'),gt.split(':')))
                if indiv_gt_phred_cut is None or float(this_gt['GQ']) >= indiv_gt_phred_cut:
                    sd['indiv_gt'][ind] = this_gt
            except:
                if skip_noGQ:
                    pass
                else:
                    print >> sys.stderr, 'parse failed for genotype string:\n\n',this_gt,'\n'
                    raise

    if len(sd['indiv_gt']) > 0: 
        sd['hwe'] = hwe(sd['indiv_gt'])
        sd['mac'] = maf(sd['indiv_gt'])
        sd['maf'] = sd['mac'] / (2. * len(sd['indiv_gt']))
        if sd['hwe'] is None:
            sd['mafbyhwe'] = None
        else:
            sd['mafbyhwe'] = sd['maf'] / numpy.log1p(sd['hwe'])
        sd['fh'] = fract_het(sd['indiv_gt'])
        sd['totcov'] = sum([int(gtd['DP']) for gtd in sd['indiv_gt'].values()])
        sd['numind'] = len(sd['indiv_gt'])
    else:
        return None

    return sd

def load_valid_vcf_line(vcfh,headers,exp_elements,FORMAT,indiv_gt_phred_cut=None,store_indiv_prefix=None,drop_indiv=None,biallelic_sites=True,skip_noGQ=True):
    while 1:
        sd = load_vcf_line(vcfh,headers,exp_elements,FORMAT,indiv_gt_phred_cut,store_indiv_prefix,drop_indiv,biallelic_sites,skip_noGQ)
        if sd == 'EOF':
            return None
        elif sd is not None:
            return sd

def load_vcf(vcf,cutoff_fn=None,ding_on=100000,store_only=None,indiv_gt_phred_cut=None,store_indiv_prefix=None,drop_indiv=None,biallelic_sites=True,skip_noGQ=True):
    '''populates and returns a site:properties dict from vcf file

    if store_only is set, must be a list of fields to retain

    if cutoff_fn is set, each sd (dict of fully parsed line, subject to store_only filter) is passed to cutoff_fn
    must be callable, and truth value of the return determines retention of that site
    '''

    vcf_data = {}
    
    i = 0

    vcfh = open(vcf)
    headers, exp_elements, FORMAT = load_vcf_header(vcfh)
    
    while 1:
        
        if i % ding_on == 0: print >> sys.stderr, 'reading',i
        i += 1

        sd = load_vcf_line(vcfh,headers,exp_elements,FORMAT,indiv_gt_phred_cut,store_indiv_prefix,drop_indiv,biallelic_sites,skip_noGQ)
        if sd is None:
            continue
        elif sd == 'EOF':
            break
        key = (sd['CHROM'],sd['POS'])
        if cutoff_fn is None or cutoff_fn(sd):
            if store_only is not None:
                keep_sd = {}
                for k in sd:
                    if k in store_only:
                        keep_sd[k] = sd[k]
                sd = keep_sd
            if len(sd) > 0:
                vcf_data[key] = sd

    return vcf_data

def genotypes_from_vcf_obj(vcf,min_indiv=40):
    pm = {}
    gt = defaultdict(dict)
    for k,v in vcf.items():
        if len(v['indiv_gt']) < min_indiv or v['mac'] == 0: continue
        sname = '%s.%s' % k
        pm[sname] = [v['REF'],v['ALT']]
        for ind,gtdict in v['indiv_gt'].items():
            gt[ind][sname] = [int(i) for i in gtdict['GT'].split('/')]

    return pm,gt

def genotypes_by_parent(pm,gt,parents,hybrids=None,remove_targets=None):
    '''Given pm and gt per tabulate_genotypes,
    and dictionary parents of the form
    {'A': ['BW1','BW2']...}
    + optional list hybrids
    returns loci and genotypes suitable for output_cross_radtag_genotypes
    + a list of target individuals to remove (e.g. parents, problematic indivs)'''

    print >> sys.stderr, '%s snps to evalute' % (len(pm))

    fixedsnps = {}
    for site in pm.keys():
        gts = defaultdict(list)
        for pk in parents.keys():
            for indiv in parents[pk]:
                try:
                    gts[pk].extend(gt[indiv].get(site,[]))
                except:
                    print indiv,site
                    raise
        gts['A'] = list(set(gts['A']))
        gts['B'] = list(set(gts['B']))
        if len(gts['A']) == 1 and len(gts['B']) == 1 and gts['A'][0] != gts['B'][0]:
            gtlookup = {}
            for pk in gts.keys():
                gtlookup[gts[pk][0]] = pk
            fixedsnps[site] = gtlookup
    print >> sys.stderr, '%s fixed snps between parents' % (len(fixedsnps))
    
    if hybrids:
        f1het = []
        for site in fixedsnps.keys():
            #print >>sys.stderr,[gt[polarized_loci,polarized_geno = extract_genotypes_from_mclgr.genotypes_by_parent(pm,gt,parents,hybrids=hybrids,remove_targets=reduce(lambda x,y: x+y,parents.values()) + hybrids)h][site] for h in hybrids if gt[h][site]]
            #print >>sys.stderr,[(len(gt[h][site]), set([fixedsnps[site][gt[h][site][i]] for i in range(2)])) for h in hybrids if gt[h][site]]
            if all([len(gt[h].get(site,[]))==2 and set(['A','B'])==set([fixedsnps[site][gt[h].get(site,[])[i]] for i in range(2)]) for h in hybrids if gt[h].get(site,[])]):
                f1het.append(site)
        loci = f1het
        print >> sys.stderr, '%s parent-fixed snps het in hybrids' % (len(f1het))        
    else:
        loci = fixedsnps.keys()


    genotypes = {}
    for indiv in gt.keys():
        genotypes[indiv] = {}
        for site in fixedsnps.keys():
            g = ''.join(sorted([fixedsnps[site].get(allele,'') for allele in gt[indiv].get(site,[])]))
            if len(g) == 2:
                genotypes[indiv][site] = g
    if remove_targets is not None:
        for t in remove_targets:
            del genotypes[t]

    return loci, genotypes

def output_cross_radtag_genotypes(loci,genotypes,filename,lg0='X'):
    '''Given list loci and dictionary genotype per genotypes_by_parent, writes file <filename>
    suitable for RQTL

    overloads 20101202:
    - if loci is a dict per maploci from load_cross_radtag_genotypes below, sort by map position in output
    - if filename is not string, use as filehandle (permits passing sys.stdout, for instance)
    '''

    def sortkey(x):
        if x == '':
            return 0
        else:
            return x

    if isinstance(loci,list):
        locnames = loci
        lgs = ['1' for i in range(len(loci))]
        mps = [str(i+1) for i in range(len(loci))]
    elif isinstance(loci,dict):
        for k,v in loci.items():
            if v[0] == 0:
                loci[k] = (lg0,v[1])
        locnames,lgs,mps = zip(*[(loc,str(lg),str(mp)) for loc,(lg,mp) in sorted(loci.items(),key=lambda x:[sortkey(v) for v in x[1]])])

    mID_lookup = dict([(m,str(i)) for i,m in enumerate(sorted(genotypes.keys()))])

    if isinstance(filename,str):
        fh = open(filename ,'w')
        #open(filename+'.mIDlookup','w').write('\n'.join(['%s\t%s' % (i,m) for m,i in sorted(mID_lookup.items())]))
    else:
        fh = filename
        #open(filename.name+'.mIDlookup','w').write('\n'.join(['%s\t%s' % (i,m) for m,i in sorted(mID_lookup.items())]))
        
    fh.write('ID,')
    fh.write(','.join(['%sr' % l for l in locnames]))
    fh.write('\n')
    fh.write(',')
    fh.write(','.join(lgs))
    fh.write('\n')
    fh.write(',')
    fh.write(','.join(mps))
    fh.write('\n')

    out_geno = {}
    
    for mID in genotypes.keys():
        #fh.write(mID_lookup[mID]+',')
        fh.write(mID+',')
        out_geno[mID] = dict([(mkr,genotypes[mID][mkr]) for mkr in locnames if genotypes[mID].has_key(mkr)])
        fh.write(','.join([genotypes[mID].get(mkr,'-') for mkr in locnames]))
        fh.write('\n')

    fh.close()

    return out_geno,mID_lookup


#other output methods
def write_structure_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):
    '''given a parsed vcf data structure per load_vcf and an outfile name

    prints a structure-formatted genotype file.

    if keys_to_write is supplied, only vcf_data items corresponding to those keys will be written
    '''

    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile,'w')
    ofh.write('\t'.join(['%s.%s' % (c,p) for c,p in keys_to_write]))
    ofh.write('\n')

    for ind in indiv_to_write:
        ofh.write(ind)
        for k in keys_to_write:
            try:
                ofh.write('\t'+(vcf_data[k]['indiv_gt'][ind]['GT'].replace('/',' ')))
            except KeyError:
                ofh.write('\t-9 -9')
        ofh.write('\n')

    ofh.close()

        
def write_spagedi_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):
    '''generates output intended for SPAGeDi

    currently treats all individuals as originating from a single population;
    this will need to be elaborated upon
    '''

    from short_read_analysis import preprocess_radtag_lane
    lookup = dict([(l['sampleid'],l['population']) for l in preprocess_radtag_lane.get_table_as_dict('DB_library_data') if l.has_key('population')])

    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile,'w')
    #write header
    ofh.write('%s\t1\t0\t%s\t1\t2\n0\nInd\tPop\t%s\n' % \
              (len(indiv_to_write),len(keys_to_write), '\t'.join(['%s.%s' % (c,p) for c,p in keys_to_write])))

    #write genotypes
    for ind in indiv_to_write:
        ofh.write('%s\t%s' % (ind,lookup.get(ind,'pop1')))
        for k in keys_to_write:
            try:
                gt = '/'.join([str(int(i)+1) for i in vcf_data[k]['indiv_gt'][ind]['GT'].split('/')])
                ofh.write('\t'+gt)
            except KeyError:
                ofh.write('\t0/0')
        ofh.write('\n')

    ofh.write('END\n')
    ofh.close()
    

def write_tassel_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):
    '''given vcf_data per load_vcf above and an outfile name

    writes tassel input format'''

    xd = { '0':'REF','1':'ALT' }

    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile,'w')

    ofh.write('%s\t%s:%s\n' % (len(indiv_to_write),len(keys_to_write),'2'))
    if len(set([c for c,p in keys_to_write])) == 1:
        ofh.write('%s\n' % '\t'.join([p for c,p in keys_to_write]))
    else:
        ofh.write('%s\n' % '\t'.join(['%s.%s' % (c,p) for c,p in keys_to_write]))

    for ind in indiv_to_write:
        ofh.write('%s' % ind)
        for k in keys_to_write:
            try:
                gt = ':'.join([vcf_data[k][xd[i]] for i in vcf_data[k]['indiv_gt'][ind]['GT'].split('/')])
            except:
                gt = '?:?'
            ofh.write('\t' + gt)
        ofh.write('\n')

    ofh.close()
    
def write_plink_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):

    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile,'w')

    idx = 0
    for ind in indiv_to_write:
        idx += 1
        ofh.write('FAM%s\t%s\t0\t0\t1\t0' % (idx,ind))
        for k in keys_to_write:
            try:
                gt = ' '.join([str(int(i)+1) for i in vcf_data[k]['indiv_gt'][ind]['GT'].split('/')])
            except:
                gt = '0 0'
            ofh.write('\t' + gt)
        ofh.write('\n')

    ofh.close()

    mapout = os.path.splitext(outfile)[0] + '.map'
    ofh = open(mapout,'w')
    infout = open(outfile + '.info','w')
    
    chrom_translation = dict([(k,i+1) for i,k in enumerate(set([k[0] for k in keys_to_write]))])
    open(mapout+'.xlat','w').write('\n'.join(['%s\t%s' % (k,v) for k, v in chrom_translation.items()]))

    for k in keys_to_write:
        ofh.write('%s\t%s.%s\t%s\t%s\n' % (chrom_translation[k[0]], k[0], k[1], 0, k[1]))
        infout.write('%s\t%s\n' % (k[1],k[1]))

    ofh.close()
    infout.close()


def write_flapjack_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):

    xd = { '0':'REF','1':'ALT' }


    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile,'w')

    if len(set([c for c,p in keys_to_write])) == 1:
        ofh.write('\t%s\n' % '\t'.join([p for c,p in keys_to_write]))
    else:
        #ofh.write('\t%s\n' % '\t'.join(['%s.%s' % (c,p) for c,p in keys_to_write]))
        raise NotImplementedError, 'currently only supports single-chromosome output'

    for ind in indiv_to_write:
        ofh.write('%s' % ind)
        for k in keys_to_write:
            try:
                gt = '/'.join([vcf_data[k][xd[i]] for i in vcf_data[k]['indiv_gt'][ind]['GT'].split('/')])
            except:
                gt = '?/?'
            ofh.write('\t' + gt)
        ofh.write('\n')

    ofh.close()
    
    mapout = os.path.splitext(outfile)[0] + '.fj.map'
    ofh = open(mapout,'w')
    
    for k in keys_to_write:
        ofh.write('%s\t%s\t%s\n' % (k[1], 1, k[1]))

    ofh.close()

def write_peas_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):
    '''given vcf_data per load_vcf above and an outfile name

    writes peas input format'''

    xd = { '0':'REF','1':'ALT' }

    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile,'w')

    ofh.write('SNPID\t%s\n' % '\t'.join('%s.%s' % (k[0],k[1]) for k in  keys_to_write))
    ofh.write('Chrom\t%s\n' % '\t'.join([k[0] for k in keys_to_write]))
    ofh.write('Position\t%s\n' % '\t'.join([k[1] for k in keys_to_write]))
    ofh.write('AlleleState\t%s\n' % '\t'.join(['%s/%s' % (vcf_data[k]['REF'],vcf_data[k]['ALT']) for k in keys_to_write]))
    ofh.write('Strand\t%s\n' % '\t'.join(['+'] * len(keys_to_write)))
    
    for ind in indiv_to_write:
        ofh.write('%s' % ind)
        for k in keys_to_write:
            try:
                gt = ''.join([vcf_data[k][xd[i]] for i in vcf_data[k]['indiv_gt'][ind]['GT'].split('/')])
            except:
                gt = '??'
            ofh.write('\t' + gt)
        ofh.write('\n')

    ofh.close()

def write_bimbam_genotypes(vcf_data, outfile, keys_to_write = None, indiv_to_write = None):
    '''does not suppport phased or multiple chromosome input vcfs
    '''

    xd = { '0':'REF','1':'ALT' }


    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    ofh = open(outfile,'w')

    if len(set([c for c,p in keys_to_write])) == 1:
        ofh.write('%s\n%s\n' % (len(indiv_to_write),len(keys_to_write)))
        ofh.write('IND,%s\n' % ','.join([ind for ind in indiv_to_write]))
    else:
        #ofh.write('\t%s\n' % '\t'.join(['%s.%s' % (c,p) for c,p in keys_to_write]))
        raise NotImplementedError, 'currently only supports single-chromosome output'

    for k in keys_to_write:
        ofh.write('%s' % k[1])
        for ind in indiv_to_write:
            try:
                gt = ''.join([vcf_data[k][xd[i]] for i in vcf_data[k]['indiv_gt'][ind]['GT'].split('/')])
            except:
                gt = '??'
            ofh.write(',' + gt)
        ofh.write('\n')

    ofh.close()
    
    mapout = os.path.splitext(outfile)[0] + '.bb.map.txt'
    ofh = open(mapout,'w')
    
    for k in keys_to_write:
        ofh.write('%s,%s\n' % (k[1], k[1]))

    ofh.close()

def write_fastPHASE_genotypes(vcf_data,outbase, keys_to_write = None, indiv_to_write = None):
    '''
    given vcf per load_vcf, writes one output file PER CHROM, name as:
    <outbase>_<CHROM>.inp

    it is recommended that outbase include path information to a directory in which these will be written, i.e.:

    outbase="beachmouse_fastPHASE/beaches_l55_k4_QD8_GQ12"

    will produce files:

    beachmouse_fastPHASE/beaches_l55_k4_QD8_GQ12_agouti_bac.inp
    beachmouse_fastPHASE/beaches_l55_k4_QD8_GQ12_mc1r_bac.inp
    ...etc
    
    directories will be created as necessary.
    '''

    
    xd = { '0':'REF','1':'ALT' }


    if keys_to_write is None:
        keys_to_write = vcf_data.keys()

    keys_to_write.sort(key = lambda x: (x[0],int(x[1])))

    if indiv_to_write is None:
        indiv_to_write = set()
        for k in keys_to_write:
            v = vcf_data[k]
            indiv_to_write = indiv_to_write.union(set(v['indiv_gt'].keys()))
        indiv_to_write = sorted(list(indiv_to_write))

    outroot = os.path.dirname(outbase)

    try:
        os.makedirs(outroot)
    except:
        pass


    this_chrom = None
    this_len = 0
    poslist = []
    for k in keys_to_write:
        if k[0] != this_chrom:
            if this_chrom is not None:
                outfile = '%s_%s.inp' % (outbase,this_chrom)
                posfile = '%s_%s.pos' % (outbase,this_chrom)
                open(posfile,'w').write(('\t'.join(poslist)) + '\n')
                ofh = open(outfile,'w')
                ofh.write('%s\n%s\n' % (len(indiv_to_write), this_len))
                for ind in indiv_to_write:
                    ofh.write('# %s\n' % ind)
                    ofh.write('%s\n%s\n' % tuple([''.join(li) for li in Util.dezip(chrom_data[ind])]))
                ofh.close()
            chrom_data = defaultdict(list)
            this_chrom = k[0]
            this_len = 0
            poslist = []
        for ind in indiv_to_write:
            try:
                chrom_data[ind].append(vcf_data[k]['indiv_gt'][ind]['GT'].split('/'))
            except:
                chrom_data[ind].append(['?','?'])
        this_len += 1
        poslist.append(k[1])

    #flush:
    outfile = '%s_%s.inp' % (outbase,this_chrom)
    posfile = '%s_%s.pos' % (outbase,this_chrom)
    open(posfile,'w').write(('\t'.join(poslist)) + '\n')
    ofh = open(outfile,'w')
    ofh.write('%s\n%s\n' % (len(indiv_to_write), this_len))
    for ind in indiv_to_write:
        ofh.write('# %s\n' % ind)
        ofh.write('%s\n%s\n' % tuple([''.join(li) for li in Util.dezip(chrom_data[ind])]))
    ofh.close()



if __name__ == "__main__":

    #parent_str = 'Ep,Ti'
    #qd = 6
    #gq = 20
    min_indiv = 2
    fh = 0.7
    site_before = numpy.inf #polymorphism must occur before this base in a fragment
    #chi2crit = 30
    
    #vcfn,qd,gq,chi2crit = sys.argv[1:]
    vcfn,parent_str,qd,gq = sys.argv[1:]
    
    
    outbase = os.path.splitext(vcfn)[0]

    cut_fn = lambda sd: sd.has_key('QD') and float(sd['QD']) >= float(qd) and len(sd['indiv_gt']) >= min_indiv and sd['fh'] < fh


    print >> sys.stderr, 'loading vcf',vcfn
    vcf = load_vcf(vcfn,cutoff_fn=cut_fn,indiv_gt_phred_cut=float(gq))
    print >> sys.stderr, '%s sites loaded' % len(vcf)

    print >> sys.stderr, 'convert to pm/gt matrices'
    pm,gt = genotypes_from_vcf_obj(vcf,min_indiv=min_indiv)
    print >> sys.stderr, 'length pm: %s length gt: %s' % (len(pm),len(gt))

    parents_prefixes = dict(zip(['A', 'B'],parent_str.split(',')))
    parents = dict([(l,[k for k in gt.keys() if k.startswith(p)]) for l,p in parents_prefixes.items()])

    polarized_loci,polarized_geno = genotypes_by_parent(dict([(k,v) for k,v in pm.items() if int(k.split('.')[1]) < site_before]),gt,parents,remove_targets=reduce(lambda x,y: x+y,parents.values()))

    #chi2-free output:
    ret = output_cross_radtag_genotypes(polarized_loci,polarized_geno,'%s_QD%s-GQ%s_%sbp.csv' % (outbase,qd,gq,site_before))
    

    """ #ditch chi2
    print >> sys.stderr, 'filter X linked, chi2 critical %s' % chi2crit
    xsites,autsites = extract_genotypes_from_mclgr.filter_Xlinked_loci(polarized_loci, polarized_geno,float(chi2crit))
    print >> sys.stderr, '%s X linked, %s autosomal' % (len(xsites),len(autsites))

    print >> sys.stderr, 'write output'
    ret = extract_genotypes_from_mclgr.output_cross_radtag_genotypes(xsites,polarized_geno,'%s_QD%s-GQ%s_%sbp_Xchi%s.csv' % (outbase,qd,gq,site_before,chi2crit))
    ret = extract_genotypes_from_mclgr.output_cross_radtag_genotypes(autsites,polarized_geno,'%s_QD%s-GQ%s_%sbp_autchi%s.csv' % (outbase,qd,gq,site_before,chi2crit))
    print >> sys.stderr, 'wrote:'
    print >> sys.stderr, '%s_QD%s-GQ%s_%sbp_Xchi%s.csv' % (outbase,qd,gq,site_before,chi2crit)
    print >> sys.stderr, '%s_QD%s-GQ%s_%sbp_autchi%s.csv' % (outbase,qd,gq,site_before,chi2crit)
    print >> sys.stderr, 'done'
    """
