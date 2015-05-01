from Seq import Sequence, Fasta
from Util import Controller
from paths import paths
from subprocess import Popen, PIPE
from copy import deepcopy

import os, sys, re, subprocess, cPickle, gzip
import numpy

import Trees, Util, Seq, GFF


##########
# store and manipulate alignments
##########

class AlnSeq (Sequence):
    '''sequence subclass (string)

    '-' chars read as gaps'''

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-','-':'-'}

    def __init__(self, seq):
        Sequence.__init__(self,seq)

    def get_ungapped(self):
        '''returns a Seq.Sequence object with '-' chars stripped'''
        ungapped = self.replace('-','')
#        ungapped = ''.join([str(c) for c in list(self) if c != '-'])
        return Sequence(ungapped)

    def aln_len(self):
        '''returns the length of the alignment'''
        return len(self)

    def seq_len(self):
        '''returns the length of the ungapped sequence'''
        return len(self.get_ungapped())

    def seq_coord(self,col,nearest=False):
        '''returns the sequence position (base) at the supplied column (col)

        since gff starts at 1, sequence starts at 1
        since sequence starts at 1, alignments start at 1 (sorry, Stuart)'''

        if self[col] == '-':
            if nearest:
                up = len(self[col-1::-1]) - len(self[col-1::-1].lstrip('-'))
                down = len(self[col+1:]) - len(self[col+1:].lstrip('-'))
                if (up<down and self[col-1::-1].lstrip('-')) or len(self[col+1:].lstrip('-')) == 0:
                    return self.seq_coord(col-1-up)
                else:
                    return self.seq_coord(col+1+down)
            else:
                return None
        else:
            return AlnSeq(self[:col]).seq_len()
        

    def aln_coord(self,base):
        '''returns alignment column of the given base

        same rules as seq_coord (sequences start at 1)'''
        
        col = 0
        b = 0
        while b <= base:
            if self[col] != '-': b += 1
            col += 1
        return col-1
            
class AlnFasta (Fasta):
	'''store and operate on alignments stored in Fasta format

	'''

	def __init__(self,filename=None, seqs_to_load=None, seqobj=AlnSeq):
		Fasta.__init__(self, filename, seqs_to_load, seqobj)

	def load(self, filename, seqs_to_load=None,seqobj=AlnSeq):
		Fasta.load(self, filename, seqs_to_load, seqobj)

	def seq_len(self,include_seq=None):
		'''return lengths of each sequence as dict.  Calls seq_len() of each member sequence

		set list "include_seq" to return only some seqs'''
		include_seq = include_seq and include_seq or [k for k in self.keys() if k[0:2] != "__"]
		return dict([(k, v.seq_len()) for k, v in self.data.items() if k in include_seq])


	def match_cols(self,match_if_gapped=False,drop_n=False):
		'''returns a list of perfect match columns'''
		pmc = []
		for col in range(len(self.values()[0])):
			if match_if_gapped:
				bases = [v[col] for v in self.values() if v[col] != '-']
			else:
				bases = [v[col] for v in self.values()]
                        if drop_n:
                            bases = [b for b in bases if b.upper() != 'N']
                        #print bases, set(bases)
			if len(set(bases)) <= 1:
				pmc.append(col)
		return pmc
		
	def nonmatch_cols(self,match_if_gapped=False,drop_n=False):
		return sorted(list(set(range(len(self.values()[0]))) - set(self.match_cols(match_if_gapped,drop_n))))
		
	def extract_polymorphism_table(self,drop_ends=False,snp_only=False,filename=None,snp_flank=5,match_if_gapped=False,drop_n=False):
		nmc = self.nonmatch_cols(match_if_gapped,drop_n)
                if len(nmc) == 1:
                	poly_cols = [(nmc[0],nmc[0]+1)]
                else:
                	poly_cols = Util.get_consecutive_value_boundaries(nmc)
		fixed_cols = self.match_cols(match_if_gapped)
		#print poly_cols
		if drop_ends:
			poly_cols = poly_cols[1:-1]
			
		table = {}
		pos_in = []
		for id in self.keys():
			table[id] = []
		if snp_only:
			for s,e in poly_cols:
				if e == s+1:
					pre_cols = range(s-snp_flank,s)
					post_cols = range(e,e+snp_flank)
					pre_check = [(pre in fixed_cols) for pre in pre_cols]
					post_check = [(post in fixed_cols) for post in post_cols]
					#print >> sys.stderr, pre_cols, post_cols, pre_check, post_check
					if all(pre_check) and all(post_check):
						if any([self[id][s]=='-' for id in self.seq_names()]) and not match_if_gapped:
							continue
						for id in self.seq_names():
							if self[id][s] == '-':
								table[id].append(None)
							else:
								table[id].append(self[id][s])
						pos_in.append(str(s))
				#table[id] = [self[id][s:e] for s,e in poly_cols \
				#			if e == s+1 and \
				#			all([(pre in fixed_cols) for pre in range(s-snp_flank,s)]) and \
				#			all([(post in fixed_cols) for post in range(e,e+snp_flank)])]
		else:
			for id in self.keys():
				table[id] = [self[id][s:e] for s,e in poly_cols]

		if snp_only:
			pos = list(pos_in)
		else:
			pos = [str(s) for s,e in poly_cols]
			
		for i in xrange(len(pos)):
			pass
		
		if filename:
			out = open(filename,'w')
			out.write('id\t%s\n' % '\t'.join(pos))
			for k,v in table.items():
				out.write('%s\t%s\n' % (k,'\t'.join([str(i) for i in v])))
			out.close()
		else:
			return (pos,table)

class AlnRefMaf ():
    '''Store and operate on REFERENCED alignments stored in maf format

    uses Seq.Fasta (and therefore Seq.Sequence) to store full, unaligned sequences,
    and Util.AssociativeArray to store alignment column info.

    currently limited to operation with REFERENCED alignments like those produced by Aln.RoastController
    when run with the ref_sp flag.  Will treat the first sequence encountered as reference, unless specified.

    currently only supports object population at construction (no piece-wise building after instance creation)
    this could change, but is not necessary

    fastas, if specified, must be a dict of species_name:fasta_file
    (easier to use fasta_dir; constructor will look for both non-extension fastas (e.g. blastz input)
    or anything with .fa[sta]?
    
    if neither is present, looks in the same dir as maf, and the parent dir

    minsize is the shortest block that will be populated
    '''

    sequences = Fasta()

    def __init__(self,maf,fastas=None,fasta_dir=None,ref_sp=None,minsize=10,populate=True,load_seqs=True,seqrange=None,init_sequences=None):


        self.filename = maf
        self.ref_sp = ref_sp or get_reference(maf)
        self.assoc_arrays = {}
        self.seqrange = seqrange

        if init_sequences:
            self.sequences.update(init_sequences)

        if seqrange:
            self.offset = seqrange[1]
        else:
            self.offset = 0

        to_load = {}
        for (seqname,length),match_dict in describe_maf(maf,seqrange).items():
            ref_seqname = '%s.%s' % fasta_seqname_from_maf_seqname(seqname)
            matches = ['%s.%s' % fasta_seqname_from_maf_seqname(m[0]) for m in match_dict.keys()]
            self.assoc_arrays[ref_seqname] = Util.AssociativeArray(matches, length, dtype=int)

            matches.append(ref_seqname)
            for maf_seqname in matches:
                sp,fasta_seqname = maf_seqname.split('.',1)
                try:
                    to_load[sp].append(fasta_seqname)
                except KeyError:
                    to_load[sp] = [fasta_seqname]

        if load_seqs:
            AlnRefMaf.sequences.update(self.load_sequences(to_load,fastas,fasta_dir))

        if populate: 
            self.populate_assoc_arrays(minsize)

    def load_sequences(self,to_load,fastas,fasta_dir):

        def first_valid_file(filelist,sp):
            for f in filelist:
                if os.path.basename(f) == sp:
                    return f
            for f in filelist:
                if re.search(sp+'\.fa[sta]*(\.gz)?$',os.path.basename(f)):
                    return f
        
        if fastas is None:
            fastas = {}
            if fasta_dir:
                files = [os.path.join(fasta_dir,f) for f in os.listdir(fasta_dir)]
            else:
                thisdir = os.path.dirname(self.filename)
                updir = os.path.dirname(thisdir)
                files = []
                for di in [thisdir,updir]:
                    files.extend([os.path.join(di,f) for f in os.listdir(di)])
            
            for sp in to_load.keys():
                fastas[sp] = first_valid_file(files,sp)

        loaded_seqs = Fasta()
        for sp,contigs in to_load.items():
            print >>sys.stderr,sp, fastas[sp], '\n\ttargets:',contigs
            if set([sp+'.'+c for c in contigs]) <= set(AlnRefMaf.sequences.keys()):
                print >>sys.stderr,'\tall loaded',AlnRefMaf.sequences.keys()
                continue
            
            sp_seqs = Fasta(fastas[sp])
            #print  >>sys.stderr,'from:\n',sp_seqs.seq_len()
            #print len(contigs), len(these_seqs.seq_len())

            for seqname,seq in sp_seqs.items():
                if sp_seqs.filename.endswith(sp):
                    clean_seqname = seqname.split(':')[1]
                else:
                    clean_seqname = seqname

                if clean_seqname in contigs:
                    new_seqname = sp+'.'+clean_seqname
                    #print  >>sys.stderr,'take',seqname,'as',new_seqname
                    loaded_seqs[new_seqname] = seq

        return loaded_seqs
        
    def populate_assoc_arrays(self,minsize):
        '''populates each associative array with sequence coords in all non-reference species indexed on the reference

        e.g. traverses each non-reference sequence in the maf, assigning the sequence coordinate of each non-reference base
        to the array element at the index of the reference base to which it aligns
        '''

        keyfollows = False
        includeblock = False
        key = ()
        
        for l in Util.smartopen(self.filename):
            if l.startswith('a'):
                keyfollows = True
                includeblock = False

            match = re.search('^s\s(?P<seqid>.+?)\s+(?P<start>\d+)\s+(?P<len>\d+)\s+(?P<strand>[+-])\s+(?P<total>\d+)\s+(?P<seq>.+?)$',l)
            if match:
                if keyfollows:
                    keyfollows = False
                    key = match.groupdict()
                    key['seqid'] = '.'.join(fasta_seqname_from_maf_seqname(key['seqid']))
                    key['start'],key['len'],key['total'] = int(key['start']),int(key['len']),int(key['total'])
                    if self.seqrange:
                        if (key['seqid'] == self.seqrange[0] and key['start'] >= self.seqrange[1] and 
                            key['start']+key['len'] <= self.seqrange[2]and key['len'] >= minsize):
                            includeblock = True
                        else:
                            includeblock = False
                    elif key['len'] >= minsize:
                        includeblock = True
                    else:
                        includeblock = False
                elif includeblock:
                    comp = match.groupdict()
                    comp['seqid'] = '.'.join(fasta_seqname_from_maf_seqname(comp['seqid']))
                    comp['start'],comp['len'],comp['total'] = int(comp['start']),int(comp['len']),int(comp['total'])
                    if self.seqrange:
                        idx = key['start']-self.seqrange[1]
                    else:
                        idx = key['start']
                    step = int(comp['strand']+'1')
                    val = step * comp['start']
                    for i, ref_base in enumerate(key['seq']):
                        if ref_base != '-':
                            if comp['seq'][i] != '-':
                                self.assoc_arrays[key['seqid']][comp['seqid']][idx] = val
                                val += step
                            idx += 1

    def get_seqs_in_reference_range(self,start,end,seqname=None,abs_coords=False,include_ref=True,returnobj_filename=None,size_factor = 2):
        if abs_coords:
            start,end = start-self.offset,end-self.offset

        if seqname is None:
            seqs_present = self.assoc_arrays.keys()
            if len(seqs_present) != 1:
                raise ValueError, 'seqname can only be omitted for AlnRefMaf objects with one seqid'
            seqname = seqs_present[0]

        seqs_retr = Fasta()
        if returnobj_filename is None:
            seqs_retr.filename = '/tmp/%s_%s_%s-all_aln_seqs.fa' % (seqname,start,end)

        if include_ref:
            seqs_retr[seqname] = self.sequences[seqname][start+self.offset:end+self.offset]

        for k in self.assoc_arrays[seqname].keys():
            try:
                bounds = self.assoc_arrays[seqname][k][start:end][self.assoc_arrays[seqname][k][start:end].nonzero()][[0,-1]]
                #print k, bounds, abs(bounds[1]-bounds[0])

                if (end-start) / size_factor < abs(bounds[1] - bounds[0]) < (end-start) * size_factor: 
                    if bounds[0] < 0:
                        hit = self.sequences[k][bounds[1]:bounds[0]].rc()
                    else:
                        hit = self.sequences[k][bounds[0]:bounds[1]]
                    seqs_retr[k] = hit
            except:
                pass

        return seqs_retr

    def ref_sp_seqnames(self):
        return self.assoc_arrays.keys()


##Useful helper functions:

def get_AlnRefMaf_for_region(region,ref_sp='dm3',aln_dir=paths['UCSC_drosophila_align'],fasta_dir=paths['UCSC_drosophila_sequence'],**kwargs):

    if os.path.exists(os.path.join(aln_dir,region['seqid']+'.maf')):
        maf = os.path.join(aln_dir,region['seqid']+'.maf')        
    elif os.path.exists(os.path.join(aln_dir,region['seqid']+'.maf.gz')):
        maf = os.path.join(aln_dir,region['seqid']+'.maf.gz')
    else:
        raise OSError, 'no file '+os.path.join(aln_dir,region['seqid']+'.maf[.gz]')

    seqrange = (ref_sp+'.'+region['seqid'], region['start'], region['end'])
    print >> sys.stderr, maf,fasta_dir, seqrange
    return AlnRefMaf(maf=maf,fasta_dir=fasta_dir,seqrange=seqrange,**kwargs)

def get_reference(maf):
    '''returns the reference species for alignment stored in self.filename
    
    looks for # roast.v3 E=ceratitis_capitata 
    '''
    for l in Util.smartopen(maf):
        match = re.search(r'#\sroast.+?\sE=(.+?)\s',l)
        if match:
            return match.groups()[0]

def fasta_seqname_from_maf_seqname(maf_seqname):
    '''given a seqname from maf sp.escaped-seqname, returns (sp,unescaped-seqname)
    '''

    def unescape_maf_periods(string):
        '''minihack to get back period-having names
        '''
        idx = string.rfind('__')
        while idx != -1:
            string = string[:idx]+'.'+string[idx+2:]
            idx = string.rfind('__')
        return string

    return (maf_seqname.split('.',1)[0], unescape_maf_periods(maf_seqname.split('.',1)[1]))
        
        

def describe_maf(maf,seqrange=None):
    '''gets the names, lengths and grouping of all sequences in a maf alignment

    returns a dict of grouped sequences and their lengths, keyed on first (seqname,seqlen) in each block,
    with values as dicts of match sequences, keyed on (seqname,seqlen), 
    with values as total number of bases anywhere aligned to outer key seq
    e.g.

    {('dmel.X_2281862_2298462___', 16601): {('dana.CH902617__1', 251596): 16,
                                            ('dana.CH902625__1', 166543): 214,
                                            ('dana.CH902632__1', 369476): 384,
                                            ('dere.CH954177__1', 166403): 423,
                                            ('dere.CH954181__1', 167147): 7,
                                            ('dere.CH954183__1', 1025277): 15125,
                                              ...
                                             }
                                          }

    if range is specified as 3-tuple (seqid,start,end), only matches in that range will be returned
    '''
    description = {}
    keyfollows = False
    if seqrange:
        includeblock = False
        key = (seqrange[0],seqrange[2]-seqrange[1])
        description[key] = {}
    else:
        includeblock = True
        key = ()

    for l in Util.smartopen(maf):
        if l.startswith('a'):
            keyfollows = True

        match = re.search(r'^s\s(.+?)\s+(\d+)\s+(\d+)\s+[+-]\s+(\d+)',l)
        if match:
            if keyfollows:
                keyfollows = False
                if seqrange:
                    start,end = int(match.groups()[1]), int(match.groups()[1])+int(match.groups()[2])
                    if match.groups()[0] == seqrange[0] and start >= seqrange[1] and end <= seqrange[2]:
                        includeblock = True
                    else:
                        includeblock = False
                else:
                    key = (match.groups()[0],int(match.groups()[3]))
                    if not description.has_key(key): description[key] = {}
            elif includeblock:
                subkey = (match.groups()[0],int(match.groups()[3]))
                matchlen = int(match.groups()[2])
                try:
                    description[key][subkey] += matchlen
                except KeyError:
                    description[key][subkey] = matchlen

    return description

def project_maf(maf,fasta,sp_order,outfile=None):
    '''given a maf file, a species order sp_order with reference as [0], and sequence in fasta
    and optional outfile (if none generates by appending .fa)

    returns tuple (projected_maf,ordered_maf,fasta)'''

    if outfile:
        pass
    else:
        outfile = maf+'.fa'

    basefile = outfile.rsplit('.',1)[0]
    projected_maf = basefile+'_proj-'+sp_order[0]+'.maf'
    ordered_maf = basefile+'_projord-'+sp_order[0]+'.maf'
    
    maf_project = '%s/maf_project %s %s > %s' % (paths['tba_dir'],maf,sp_order[0],projected_maf)
    maf_order = '%s/maf_order %s %s all > %s' % (paths['tba_dir'],projected_maf,' '.join(sp_order),ordered_maf)
    maf2fasta = '%s/maf2fasta %s %s fasta iupac2n > %s' % (paths['tba_dir'],fasta,ordered_maf,outfile)
    os.system(maf_project)
    os.system(maf_order)
    os.system(maf2fasta)

    return (projected_maf,ordered_maf,outfile)

def get_best_match_to_query(filename,qseqnum='2',source=None,sourceparse=lambda x:x,qfasta=None,tfasta=None,return_n=1):
	'''extracts best single match to each query sequence from a tabular file such as lastz .general table
	
	returns GFF formatted intervals from the target sequence
	
	if source is not None, sets ['source'] in GFF region to sourceparse(hit['source'])
	
	returns <return_n> best matches
	'''
	
	scores = Util.parse_tabular(filename)
	
	#break into lists by query
	qname = 'name'+qseqnum
	if qseqnum == '2':
		tseqnum = '1'
	elif qseqnum == '1':
		tseqnum = '2'
	else:
		raise ValueError, 'qseqnum really should be 1 or 2 (is %s)' % qseqnum
		
	qseqnames =  set([h[qname] for h in scores])
	
	by_qseq = {}
	for q in qseqnames:
		by_qseq[q] = [h for h in scores if h[qname] == q]
		
	best_matches = GFF.File()
	for query,hits in by_qseq.items():
		hits.sort(key=lambda x:float(x['score']),reverse=True)
		for h in hits[:return_n]:
			r = GFF.Region(fields=
				{
				'seqid':h['name'+tseqnum],
				'start':h['start'+tseqnum],
				'end':h['end'+tseqnum],
				'score':h['score'],
				'strand':h['strand'+tseqnum],
				'attribute_qcovPct':h['covPct'],
				'attribute_qcoverage':h['coverage'],
				'attribute_qidPct':h['idPct'],
				'attribute_qidentity':h['identity'],
				'attribute_qname':h['name'+qseqnum],
				'attribute_qstart':h['start'+qseqnum],
				'attribute_qend':h['end'+qseqnum],
				'attribute_qstrand':h['strand'+qseqnum],
				})
			if source is not None:
				try:
					r['source'] = sourceparse(h[source])
				except:
					pass
			if qfasta is not None:
				r['attribute_qfasta'] = os.path.abspath(qfasta)
			if tfasta is not None:
				r['attribute_tfasta'] = os.path.abspath(tfasta)
			best_matches.append(r)
	
	return best_matches

def get_zcut_match_to_query(filename,qseqnum='2',zcut=1,gff=True):
    '''extracts and sums score information for all sequences in a tabular file (e.g. lastz .general table)
	intended for use with fragmentary output (e.g. potentially many "correct" matches)

    if gff=False, returns dict:


    { qname1 : ('tname1:start1-end1',score1)
      qname2 : ('tname2:start2-end2',score2)
    }

    if gff=True, returns list of gff regions, coords in target genome
    IMPLEMENT GFF!
    

    AND scores (original tabular data from .general or other table)

    N.B. - qseqnum determines which sequence (1 or 2) should be treated as 'query'
    RETURNS ONLY REGIONS SCORING 
    '''

    if qseqnum == '2':
        tseqnum = '1'
    else:
        tseqnum = '2'

    qname = 'name'+qseqnum
    qstart = 'start'+qseqnum
    qend = 'end'+qseqnum
    tname = 'name'+tseqnum
    tstart = 'start'+tseqnum
    tend = 'end'+tseqnum

    scores = Util.parse_tabular(filename)

    '''NO LONGER PRE-SPLIT
    sum_by_strand = {}.fromkeys(set([d['name'+qseqnum] for d in scores]),None)
    for k in sum_by_strand.keys():
        sum_by_strand[k] = {}
        for d in scores:
            if d['name'+qseqnum] == k:
                sum_by_strand[k][d['name'+tseqnum]] = {'+':0,'-':0}
    '''

    #sort scores by qname, then tname, then tstart, then tend
    scores.sort(key=lambda x: (x['name'+qseqnum],x['name'+tseqnum],x['start'+tseqnum],x['end'+tseqnum]))


    last = scores[0]
    best_score = float(last['score'])
    best_tname = last[tname]
    best_start = int(last[tstart])
    matches = {best_tname:{best_start:best_score}}

    if gff:
        best_matches = GFF.File()
    else:
        best_matches = {}
    qlens = {}

    for score_dict in scores[1:]:
        if score_dict[qname] == last[qname]:
            if matches.has_key(score_dict[tname]):
                for k in matches[score_dict[tname]].keys():
                    if k + int(score_dict['size'+qseqnum]) >= int(score_dict[tstart]): #currently just requires overlap, not contains
                        matches[score_dict[tname]][k] += float(score_dict['score'])
                matches[score_dict[tname]][int(score_dict[tstart])] = float(score_dict['score'])
            else:
                matches[score_dict[tname]] = {int(score_dict[tstart]):float(score_dict['score'])}
        else:
            #remove overlaps
            for k,v in matches.items():
                if len(v) > 1:
                    drops = set([])
                    v_sorted = sorted(v.items(),key=lambda x:x[1],reverse=True)
                    for i,(v_start,v_score) in enumerate(v_sorted):
                        for next_start,next_score in v_sorted[i+1:]:
                            if (v_start + int(last['size'+qseqnum]) > next_start):
                                drops.add(next_start)
                    for d in drops:
                        v.pop(d)

            flatmatches = sorted([((k,st,st+int(score_dict['size'+qseqnum])),sc) \
                                  for k, v in matches.items() for st,sc in v.items()],key=lambda x:x[1],reverse=True)
            if max(Util.zscore([i[1] for i in flatmatches])) >= zcut:
                if gff:
                    reg = GFF.Region(fields=
                                     {'seqid':flatmatches[0][0][0],
                                      'start':flatmatches[0][0][1],
                                      'end':flatmatches[0][0][2],
                                      'score':flatmatches[0][1],
                                      'attribute_Name':last[qname]} )
                    best_matches.append(reg)
                else:
                    best_matches[last[qname]] = '%s:%s-%s' % flatmatches[0][0]
            #set new 'last'
            last = score_dict
            best_score = float(last['score'])
            best_tname = last[tname]
            best_start = int(last[tstart])
            matches = {best_tname:{best_start:best_score}}

    #repeat last processing step
    for k,v in matches.items():
        if len(v) > 1:
            drops = set([])
            v_sorted = sorted(v.items(),key=lambda x:x[1],reverse=True)
            for i,(v_start,v_score) in enumerate(v_sorted):
                for next_start,next_score in v_sorted[i+1:]:
                    if (v_start + int(last['size'+qseqnum]) > next_start):
                        drops.add(next_start)
            for d in drops:
                v.pop(d)

    flatmatches = sorted([((k,st,st+int(score_dict['size'+qseqnum])),sc) \
                          for k, v in matches.items() for st,sc in v.items()],key=lambda x:x[1],reverse=True)
    if max(Util.zscore([i[1] for i in flatmatches])) >= zcut:
        if gff:
            reg = GFF.Region(fields=
                             {'seqid':flatmatches[0][0][0],
                              'start':flatmatches[0][0][1],
                              'end':flatmatches[0][0][2],
                              'score':flatmatches[0][1],
                              'attribute_Name':last[qname]} )
            best_matches.append(reg)
        else:
            best_matches[last[qname]] = '%s:%s-%s' % flatmatches[0][0]
 
    return best_matches,scores

def load_wig(filename,seqlen=None):
    '''loads WIG format files (like phastCons output) into a list, padding w/ 0s as necessary
    '''
    li = []
    if seqlen is None and 'maf' in filename:
        try_fa = filename[:filename.find('maf')]+'fa'
        if os.path.exists(try_fa):
            seqlen = Seq.Fasta(try_fa).seq_len().values()[0]

    for l in open(filename):
        match = re.search('start=(\d+)',l)
        if match:
            target = int(match.groups()[0])
            while len(li) < target:
                li.append(0.0)
        else:
            li.append(float(l))

    if seqlen:
        while len(li) < seqlen+1:
            li.append(0.0)
        
    return li

def split_maf_and_fasta(maf,fasta=None,outdir=None,clear_empties=1):
    '''assumes maf and fasta in TBA/roast format (e.g. fasta has sp:scaff:otherstuff and maf has sp.scaff)

    deposits .fa and .maf files where sequences are sp named and files are scaff named in dir outdir'''

    if outdir is None:
        outdir = os.path.join(os.path.dirname(maf),'single_contigs')

    try:
        os.makedirs(outdir)
    except:
        pass
    
    if fasta:
        all_fa = Fasta(fasta)
        all_seqs = []
        for seq in all_fa.seq_names():
            (sp,seqname) = seq.split(':')[0:2]
            fa = Fasta()
            fa[sp] = all_fa[seq]
            fa.write_to_file(os.path.join(outdir,seqname+'.fa'))
            all_seqs.append(seqname)
    else:
        all_seqs = list(set(re.findall(r'a\sscore.+?\ns\s\w+\.(\w+)',open(maf).read())))

    li = re.split(r'(a\sscore=[\.\d]+\n)',open(maf).read())
    header = li[0]
    mafs = {}
    for seq in all_seqs:
        mafs[seq] = deepcopy([header])
    for i in range(1,len(li),2):
        match = re.search(r'^s\s.+?\.(.+?)\s',li[i+1])
        if match:
            seqname = match.groups()[0]
            mafs[seqname].extend(li[i:i+2])

    for k,vl in mafs.items():
        open(os.path.join(outdir,k+'.maf'),'w').writelines(vl)
        if clear_empties and len(vl) == 1:
            os.unlink(os.path.join(outdir,k+'.maf'))
            os.unlink(os.path.join(outdir,k+'.fa'))
            del all_seqs[all_seqs.index(k)]

    if fasta:
        return zip([os.path.join(outdir,s+'.fa') for s in all_seqs],[os.path.join(outdir,s+'.maf') for s in all_seqs])
    else:
        return [os.path.join(outdir,s+'.maf') for s in all_seqs]

def split_maf_to_segments(maf,size=30000,outdir=None):
    '''takes a maf file, splits into multiple mafs spanning roughly size (in bases)
    '''
    
    if outdir is None:
        outdir = os.path.dirname(maf)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    li = re.split(r'(a\sscore=[-\.\d]+\n)',open(maf).read())
    header = li[0]
    maf_seg = [header]
    curr_start = None
    curr_end = 1
    seqid = os.path.splitext(os.path.basename(maf))[0]
    for i in range(1,len(li),2):
        match = re.search(r'^s\s.+?\..+?\s+(\d+)\s+(\d+)',li[i+1])
        if match:
            (start,length) = match.groups()
            if curr_start is None: curr_start = int(start)
            maf_seg.extend(li[i:i+2])
            curr_end = int(start)+int(length)
        if curr_end > curr_start + size or i+1 >= len(li):
            fname = os.path.join(outdir,'%s_%s_%s.maf' % (seqid,curr_start,curr_end))
            open(fname,'w').write(''.join(maf_seg))
            curr_start = None
            maf_seg = [header]

    
    

##########
# Generate alignments
##########

class AlignmentController (Controller):
    '''generic class with useful tools for controlling apps that generate alignments

    e.g. Lagan and TBA'''

    def __init__(self,aln=None,tree=paths["PD_tree"],newickstr=None,
                 tempdir=None,stdin=None,stderr=None,stdout=None,executor=Popen,use_defIOE=0,name=None,cleanup=0,args=None):
        '''just handles alignment specific items--IO and tempfiles from ancestor

        '''
        Controller.__init__(self,tempdir,stdin,stderr,stdout,executor,use_defIOE,name,cleanup,args)
        self["aln"] = aln
        self["tree"] = tree
        self["newickstr"] = newickstr

    def __setitem__(self, key, item):
        Controller.__setitem__(self,key,item)
        if key == "newickstr" and isinstance(item,str):
            self.data["tree"] = Trees.Tree(item)
        if key == "tree" and isinstance(item,str):
            self["newickstr"] = open(item).read()

    def fasta_filename(self):
        try:
            return self["fasta"].filename 
        except:
            return None

    def treestr(self,sp_list=None,ref_sp=None):
        '''takes list of species for a subtree or nothing

        returns newick tree of subset of species given, 
        or subset in fasta headers, or in aln fasta headers (if None)
        if ref_sp specified, will call Trees.redorder_newick() to set ref_sp first
        '''

        if sp_list is None:
            if "fasta" in self.data.keys():
                if isinstance(self["fasta"], str):
                    self["fasta"] = Fasta(self["fasta"])
                sp_list = self["fasta"].seq_names()
            elif "aln" in self.data.keys():
                if self['aln'] is not None:
                    if isinstance(self['aln'], str):
                        if self['aln'].endswith('.maf'):
                            match = re.search(r'#\sroast.+?\s(\(.+?\))\s/',open(self['aln']).read())
                            if match:
                                sp_list = re.findall(r'\w+',match.groups()[0])
                            else:
                                sp_list = self['tree'].get_taxa()
                        else:
                            self['aln'] = AlnFasta(self['aln'])
                            sp_list = self["aln"].seq_names()
            #elif "fastas" in self.data.keys():
            #    sp_list = [os.path.basename(f.filename).split('.')[0] for f in self["fastas"]]
        subtree = self["tree"].get_subset_tree(sp_list)
        subtree_newick = subtree.to_string(plain_newick=1).strip(';')
        if ref_sp:
            subtree_newick = Trees.reorder_newick(subtree_newick,ref_sp)
        return subtree_newick

    def space_treestr(self,sp_list=None,ref_sp=None):
        '''makes wacko lagan tree format
        '''
        return self.treestr(sp_list,ref_sp).replace(',',' ')

class LaganController (AlignmentController):
    '''run lagan.

    lagan requires fasta, tree (filenames or Seq.Fasta and Trees.Tree) and aln (filename)

    tree will be subtree'd from PD_tree if not provided'''

    command = paths["mlagan"]

    controller_type = 'mlagan'
    
    def __init__(self,fasta=None,aln=None,tree=paths["PD_tree"],newickstr=None,
                 tempdir=None,stdin=None,stderr=None,stdout=None,executor=Popen,use_defIOE=0,name=None,cleanup=0):
        '''just handles lagan specific items--alignment, IO and tempfiles from ancestor

        '''
        AlignmentController.__init__(self,aln,tree,newickstr,tempdir,stdin,stderr,stdout,executor,use_defIOE,name,cleanup)

        self["fasta"] = fasta
        self["tailed"] = 0 #ticker to ensure A-tailing only happens once

    def __setitem__(self, key, item):
        AlignmentController.__setitem__(self,key,item)
        if key == "fasta" and isinstance(item,str):
            self.data["fasta"] = Fasta(item)
            
    def lagan_treestr(self):
        '''makes wacko lagan tree format
        '''
        return AlignmentController.space_treestr(self)

    def prepare_fastas(self,flush_memory=0):
        '''adds trailing 'A' and dumps to individual files in self["tempdir"]

        if flush_memory=1, fasta sequences are removed from self["fasta"]
        can be reloaded at any time with self["fasta"].load(self["fasta"].filename)'''

        if self["fasta"] and not self["tailed"]:
            a = Sequence('A')
            for seq in self["fasta"].seq_names():
                self["fasta"][seq] += a
            self["tailed"] = 1
        seq_files = self["fasta"].split_multifasta(self["tempdir"],basenames_only=1,ext='')
        #add files created this way to list of files generated by this object
        self['files'].extend([os.path.join(self['tempdir'],f) for f in seq_files])
        #add all .anchors files (pairwise seq name combos) to list of files generated
        seq_files.sort()
        for i in range(0,len(seq_files)):
            for j in range(i+1,len(seq_files)):
                self['files'].append(os.path.join(self['tempdir'],'%s%s.anchors' % (seq_files[i],seq_files[j])))
                
        if flush_memory:
            for seq in self["fasta"].seq_names():
                self["fasta"][seq] = ''
        return seq_files

    def compose_arguments(self):
        '''build the args list for subprocess.Popen (or whatever)

        if aln is still undefined at this point, output becomes tempdir/fasta_filename.aln
        if no fasta_filename, becomes "alignment.aln"'''

        os.environ["LAGAN_DIR"] = paths["LAGAN_DIR"]

        if self["aln"] is None:
            if self.fasta_filename():
                alnfile = os.path.basename(self.fasta_filename())+'.aln'
            else:
                alnfile = self["name"]+'.aln'
            self["aln"] = os.path.join(self["tempdir"], alnfile)

        args = [self.command]
        args.extend(self.prepare_fastas())
        args.extend(['-tree',"'%s'" % self.lagan_treestr()])
        args.extend(['-out',self["aln"]])
        args.append("-verbose")

        return args

    def handle_results(self):
        '''to be run upon completion of alignment; removes 'A' tails

        to do: clean up tempdirs, etc'''

        new_aln = AlnFasta(self["aln"])
        for seq in new_aln.seq_names():
            new_aln[seq] = new_aln[seq][:-1]
        new_aln.metadata.update(self['fasta'].metadata)   
        new_aln.write_to_file()
        Controller.handle_results(self)
        return new_aln

class AllBZController (AlignmentController):
    '''run all_bz from TBA package to generate all pairwise blastz alignments

    '''

    controller_type = 'AllBZ'

    command = os.path.join(os.path.dirname(paths['TBA']),'all_bz')

    def __init__(self,fastas=None,blastz_specs=None,tree=paths["PD_tree"],newickstr=None,ref_sp=None,preblat=None,
                 tempdir=None,stdin=None,stderr=None,stdout=None,executor=Popen,use_defIOE=0,name=None,cleanup=0):

        AlignmentController.__init__(self,None,tree,newickstr,tempdir,stdin,stderr,stdout,executor,use_defIOE,name,cleanup)

        if blastz_specs is None: blastz_specs = os.path.join(os.path.dirname(self.command),'blastz.specs')
        self.data["fastas"] = fastas
        self['blastz_specs'] = blastz_specs
        self['species'] = []
        self['mafs'] = []
        self['ref_sp'] = ref_sp
	if self['ref_sp']:
		self['ref_sp_idx'] = [self.sp_from_filename(i) for i in self['fastas']].index(self['ref_sp'])
        self['preblat'] = preblat

    def __setitem__(self, key, item):
        '''for now, direct assignment of a list of strings to self['fastas']
        will cause those strings to be treated as fasta filenames, and those files will be loaded as Fasta objects

        since this is a bit odd, self['fastas'] is best set through the __init__ constructor'''
        
        AlignmentController.__setitem__(self,key,item)
        if key == "fastas":
            if isinstance(item, list):
                for i in range(len(item)):
                    if isinstance(item[i],str):
                        if not os.path.exists(item[i]):
                            raise IOError, "file %s doesn't exist!" % item[i] 
                        item[i] = Fasta(item[i])
                    elif isinstance(item[i],Fasta):
                        pass
                    else:
                        raise TypeError, "items in fastas list must be either Seq.Fasta objects or strings"
                self.data['fastas'] = item
            else:
                raise TypeError, "fastas member variable must be a list"
        if key == 'blastz_specs':
            if isinstance(item,dict) or isinstance(item,list):
                self.prepare_blastz_specs()

    def sp_from_filename(self,filen):
	if isinstance(filen,str):  #if not string, maybe a Seq.Fasta instance, try to get a .filename
		filename = filen
	else:
		filename = filen.filename
	filename = filename.strip('.gz')
        if re.search(r'(\.fa|\.fasta|\.fas)$',filename,re.I):
            return os.path.basename(filename).split('.')[-2]
        else:
            #assume files are already formatted, so filename is species
            return os.path.basename(filename)

    def prepare_blastz_specs(self):
        '''assumes self['blastz_specs'] is a blastz parameters hash, which looks like:

        {
        GROUP : ([member species], {OTHER_GROUP : "param string"})
        }

        or a list of parameter tuples to apply to the whole set of species:

        [(param, value)]'''
        specs = self['blastz_specs']
        self.data['blastz_specs'] = os.path.join(self['tempdir'],'blastz'+self['name']+'.specs')
        bz_specfile = open(self.data['blastz_specs'],'w')

        if isinstance(specs,dict):
            for group,(species,params) in specs.items():
                bz_specfile.write('#define %s %s\n' % (group,' '.join(species)))
            for group,(species,params) in specs.items():
                for othergroup,paramstr in params.items():
                    bz_specfile.write('%s : %s\n\t%s\n' % (group,othergroup,paramstr))
        elif isinstance(specs,list):
            filenames = []
            for fa in self['fastas']:
                filenames.append(self.sp_from_filename(fa))

            print '#define ALL %s\n' % (' '.join(filenames), )
            print 'ALL : ALL\n\t%s\n' % (' '.join(['%s=%s' % (k,v) for (k,v) in specs]))
            
            bz_specfile.write('#define ALL %s\n' % (' '.join(filenames), ))
            bz_specfile.write('ALL : ALL\n\t%s\n' % (' '.join(['%s=%s' % (k,v) for (k,v) in specs])))

        bz_specfile.close()

    def prepare_fastas(self,flush_memory=0):
        '''copies fasta files to self["tempdir"], adds required header info

        if flush_memory=1, fasta sequences are removed from self["fasta"]
        can be reloaded at any time with self["fastas"]['species_name'].load(self["fasta"].filename)'''

        for (i,fasta) in enumerate(self['fastas']):
            #if a fasta item is still a string, check to see if it needs to be loaded
            # i.e. if filename ends in a fasta extension, likely has not been formatted for blastz
            if isinstance(fasta,str):
                sp = self.sp_from_filename(fasta)
                if re.search(r'(\.fa|\.fasta|\.fas)(.gz)?$',fasta,re.I):
			self['fastas'][i] = Fasta(fasta)

			outfile = os.path.join(self['tempdir'],sp)

			if self['ref_sp'] and self['preblat'] and sp != self['ref_sp']:
				self['fastas'][i] = find_subset_matching_reference(self['fastas'][self['ref_sp_idx']],self['fastas'][i],tempdir=self['tempdir'])
                        
			format_fasta_instance_for_bz(self['fastas'][i],outfile)

                elif os.path.dirname(fasta) != self['tempdir']:
			#if a filename is given without an extension, assume it's bz formatted
			#symlink it in the working directory
			os.system('ln -s %s %s' % (fasta, self['tempdir']))
            else:
                sp = self.sp_from_filename(fasta.filename)
                outfile = os.path.join(self['tempdir'],sp)
                
		if self['ref_sp'] and self['preblat'] and sp != self['ref_sp']:
			self['fastas'][i] = find_subset_matching_reference(self['fastas'][self['ref_sp_idx']],self['fastas'][i],tempdir=self['tempdir'])
                
                format_fasta_instance_for_bz(self['fastas'][i],outfile)

            #self['files'].append(outfile) cannot remove these files--required by TBA
            self['species'].append(sp)

        if flush_memory:
            try:
                for seq in self["fasta"].seq_names():
                    self["fasta"][seq] = ''
            except:
                pass

    def calc_output_sp_pairs(self):
        '''returns a list of output species pairs in the order that they will be created

        requires that self['sp_order'] exist'''

        if len(self['sp_order']) < 2:
            raise KeyError, "too few sp in self['sp_order']"
        
        sp_order = self['sp_order']
        sp_pairs = []
        for i in range(len(sp_order)):
            for j in range(i+1,len(sp_order)):
                sp_pairs.append('%s.%s' % (sp_order[i],sp_order[j]))
        return sp_pairs

    def space_treestr(self,sp_list=None,ref_sp=None):
        '''special instance of the method captures the tree order, which all_bz uses in file naming

        call to set self['sp_order'] - TO DO: pull sp_order calc out to another sub, to call independently'''
        treestr = AlignmentController.space_treestr(self,sp_list,ref_sp)
        self['sp_order'] = [s.strip('()') for s in treestr.split()]
        return treestr
                

    def compose_arguments(self,flush_memory=0):

        self.prepare_fastas(flush_memory)

        sargs = [self.command]
        sargs.append('-')
        if self['ref_sp']:
            sargs.extend(['F='+self['ref_sp']])

        sargs.append('%s' % (self.space_treestr(self['species'],self['ref_sp']),))
        sargs.append(self['blastz_specs'])

        bz_out = Popen(sargs,stdout=PIPE)
        scriptname = os.path.join(self['tempdir'],self['name']+'-bz.sh')
        bz_script = open(scriptname,'w')
        bz_script.write('cd %s\n' % self['tempdir'])
        if self['ref_sp']:
            for line in bz_out.stdout:
                if re.search(self['ref_sp'],line):
                    bz_script.write(line)
        else:
            for line in bz_out.stdout:
                bz_script.write(line)
        bz_script.close()

        os.chmod(scriptname,0777)
        args = ['bash',scriptname]
        for n in self.calc_output_sp_pairs():
            if self['ref_sp'] is None or re.search(self['ref_sp'],n):
                print self['ref_sp'],n
                self['files'].append(os.path.join(self['tempdir'],n+'.orig.maf'))
                self['mafs'].append(os.path.join(self['tempdir'],n+'.sing.maf'))
                
        self['files'].append(scriptname)

        return args

    def handle_results(self):
        '''runs cleanup if called for, a la ancestor

        currently returns the sing.maf files, but might return the list of commands to run,
        if a future revision enables running with the '-' flag to push to cluster'''
        
        AlignmentController.handle_results(self)
        return self['mafs']


#blastz/roast/TBA helpers

def format_fasta_instance_for_bz(fasta_obj,outfile):
    '''helper function for use in properly formatting fasta headers to conform to blastz input requirements

    returns the modified fasta'''
    new = Fasta()
    for seq,length in fasta_obj.seq_len().items():
        new['%s:%s:1:+:%s' % (os.path.basename(outfile),seq,length)] = fasta_obj[seq]
        
    if outfile: 
        new.write_to_file(outfile,metadata=0)

    return new

def prep_roast_from_prereq(item,Dispatcher_obj):
    '''takes a triplet dict from a Util.Dispatcher queue (controller/runner/prereqs)
    in which the prereqs set has a single item, the name of the AllBZ controller
    responsible for pre-roast alignments
    
    '''
    if item['controller'].controller_type == 'roast' and len(item['prereqs']) == 1:
        #set roast starting point (mafs) with the list of files in the dispatcher results dict
        # .returns stored under the key matching the current item's prereq (AllBZ instance)
        item['controller']['mafs'] = Dispatcher_obj.returns[list(item['prereqs'])[0]]
        

class RoastController (AlignmentController):
    '''run roast on a directory of multifasta files;
    assumes one file per species, named (.\.)?species.fa
    or Fasta objects of same (or any combination thereof)
    in list argument 'fastas'.  Otherwise, must be provided with mafs, e.g.:
    ref_sp.other_sp1.sing.maf
    ref_sp.other_sp2.sing.maf
    etc--see AllBZController to generate these (using the ref_sp mode of AllBZController will save time)

    handle resulting output in a variety of ways

    TO-DO: extend AlignmentController treestr() method to try getting sp_from_mafs()

    '''
    
    controller_type = 'roast'

    command = paths[controller_type]

    def __init__(self,fastas=None,mafs=None,aln=None,blastz_specs=None,ref_sp=None,tree=paths["PD_tree"],newickstr=None,
                 tempdir=None,stdin=None,stderr=None,stdout=None,executor=Popen,use_defIOE=0,name=None,cleanup=0):
        '''just handles roast specific items--alignment, IO and tempfiles from ancestor

        '''
        AlignmentController.__init__(self,aln,tree,newickstr,tempdir,stdin,stderr,stdout,executor,use_defIOE,name,cleanup)

        self["fastas"] = fastas
        self['blastz_specs'] = blastz_specs
        self['mafs'] = mafs
        self['aln'] = aln
        self['ref_sp'] = ref_sp
        
    def __setitem__(self, key, item):
        AlignmentController.__setitem__(self,key,item)
        

    def sp_from_filename(self,filen):
	if isinstance(filen,str):  #if not string, maybe a Seq.Fasta instance, try to get a .filename
		filename = filen
	else:
		filename = filen.filename
	filename = filename.strip('.gz')

        if re.search(r'(\.fa|\.fasta|\.fas)$',filename,re.I):
            return os.path.basename(filename).split('.')[-2]
        else:
            #assume files are already formatted, so filename is species
            return os.path.basename(filename)

    def sp_from_mafs(self):
        d = {}
        for m in self['mafs']:
            sp = os.path.basename(m).split('.')[:2]
            d.update({sp[0]:'',sp[1]:''})
        return d.keys()

    def compose_arguments(self):
	'''workaround nasty bug that tries to double-quote in the tree string, 
	by concatenating args, adding stdin and stdout as appropriate, then pushing that to a .sh script
	'''
        if not self['ref_sp']:
            raise ValueError, 'ref_sp must be set before calling compose arguments'
        args = [self.command]
        args.append('E='+self['ref_sp'])
        if not self['mafs']:
            sys.stderr.write('mafs not supplied\n')
            if self['fastas']:
                #sys.stderr.write('fastas: %s\n' % (' '.join(self['fastas']),))
                bzer = AllBZController(fastas=self['fastas'],
                                       blastz_specs=self['blastz_specs'],ref_sp=self['ref_sp'],
                                       tempdir=self['tempdir'],use_defIOE=1,cleanup=self['cleanup'],executor=self['executor'])
                self['mafs'] = bzer.run(1)
            else:
                raise SystemError, 'either fastas or mafs must be set'
        args.append('"%s"' % (self.space_treestr(self.sp_from_mafs(),self['ref_sp']), ) )
        args.append(self['tempdir'])
        if self['aln'] is None:
            self['aln'] = os.path.join(self['tempdir'],'tba'+self['name']+'.maf')
        args.append(self['aln'])
	if self['stdout']: args.extend(['>',self['stdout'].name])
	if self['stderr']: args.extend(['2>',self['stderr'].name])

	script = os.path.join(self['tempdir'],self['name']+'-roast.sh')
	open(script,'w').write('cd %s\n%s\n' % (self['tempdir'],' '.join(args)))
	
	args = ['bash',script]
        return args

    def handle_results(self):
        indiv_contigs = split_maf_and_fasta(self['aln'],os.path.join(self['tempdir'],self['ref_sp']),os.path.join(self['tempdir'],'single_contigs'))
        sp_order = [self['ref_sp']]+[s for s in self.sp_from_mafs() if s != self['ref_sp']]
        single_contig_results = [project_maf(maf,fa,sp_order) for (fa,maf) in indiv_contigs]
        AlignmentController.handle_results(self)
        return single_contig_results

def find_subset_matching_reference(ref_fasta_obj, subj_fasta_obj, pad_multiplier=5, blat_args={'minIdentity':'0.75'},
                                   tempdir=None, outfile=None, message_sock=None,logfile=None):
	'''helper method runs blat to find regions in subject matching the reference query.
	pads these regions by (pad_multiplier X length of the query seq)
	if outfile in None, returns a Fasta object of these match regions
        else writes resulting object to outfile
	'''
        if logfile: message_sock = open(logfile,'a')
        if message_sock: message_sock.write('finding subset regions for %s in %s (output to %s %s)\n' % 
                                            (ref_fasta_obj.filename,subj_fasta_obj.filename,tempdir,outfile))
        if logfile: message_sock.close()

	if not os.path.exists(ref_fasta_obj.filename): 
            pathdir = os.path.dirname(ref_fasta_obj.filename)
            if not os.path.exists(pathdir): os.makedirs(pathdir)
            ref_fasta_obj.write_to_file()
	if not os.path.exists(subj_fasta_obj.filename): 
            pathdir = os.path.dirname(subj_fasta_obj.filename)
            if not os.path.exists(pathdir): os.makedirs(pathdir)
            subj_fasta_obj.write_to_file()
	
        if logfile: message_sock = open(logfile,'a')
        if message_sock: message_sock.write('run blat controller...')
        if logfile: message_sock.close()

	blater = BlatController(ref_fasta_obj.filename,subj_fasta_obj.filename,blat_args,tempdir=tempdir,logfile=logfile)	
	(Qgff,Sgff) = blater.run(1)

        if logfile: message_sock = open(logfile,'a')
        if message_sock: message_sock.write('done\n')
        if logfile: message_sock.close()

        qlens = ref_fasta_obj.seq_len()
        slens = subj_fasta_obj.seq_len()
	for i, reg in enumerate(Sgff):
		
		qlen = qlens[Qgff[i]['seqid']]
		slen = slens[Sgff[i]['seqid']]
		pad = qlen * pad_multiplier

		Sgff[i]['start'] -= pad
		if Sgff[i]['start'] < 1: Sgff[i]['start'] = 1

		Sgff[i]['end'] += pad
		if Sgff[i]['end'] > slen: Sgff[i]['end'] = slen
	
	Sgff.collapse_overlapping_regions()
        result_fasta = subj_fasta_obj.substr_from_gff(Sgff.data,name_key=None)
        if outfile:
            if not os.path.exists(os.path.dirname(outfile)): os.makedirs(os.path.dirname(outfile))
            result_fasta.write_to_file(outfile)
        else:
            return result_fasta

def unwrap_if_gz(fasta):
    import os
    if isinstance(fasta,str):
        if fasta.endswith('.gz'):
            unzipped = fasta.rstrip('.gz')
            if not os.path.exists(unzipped):
                os.system('gunzip -c %s > %s' % (fasta,unzipped))
            return unzipped
        else:
            return fasta
    else:
        if fasta.filename.endswith('.gz'):
            fasta.filename = fasta.filename.rstrip('.gz')
        return fasta

		
class BlatController(Controller):
    '''takes either fasta instances or fasta files as query and subject,
    generates .blat file in tempdir named based on instance name

    handle_results to get GFF.File lists of match regions in query and subject
    '''
    
    controller_type = 'blat'

    command = paths[controller_type]

    def __init__(self,query=None,subject=None,blatargs={},tempdir=None,stdin=None,stderr=None,stdout=None,
                 executor=Popen,use_defIOE=0,name=None,cleanup=0,args=None,logfile=None):

        Controller.__init__(self,tempdir,stdin,stderr,stdout,executor,use_defIOE,name,cleanup,args,logfile)

        self['query'] = unwrap_if_gz(query)
        self['subject'] = unwrap_if_gz(subject)
        self['blatargs'] = blatargs
        self['logfile'] = logfile

    def unwrap_if_gz(self,fasta):
        import os
        if isinstance(fasta,str):
            if fasta.endswith('.gz'):
                unzipped = fasta.rstrip('.gz')
                if not os.path.exists(unzipped):
                    os.system('gunzip -c %s > %s' % (fasta,unzipped))
                return unzipped
            else:
                return fasta
        else:
            if fasta.filename.endswith('.gz'):
                fasta.filename = fasta.filename.rstrip('.gz')
            return fasta

    def compose_arguments(self):
        if self['logfile']: open(self['logfile'],'a').write('starting compose_arguments\n')
        from Seq import write_fastas_to_files

        args = [self.command]
        blatargs_list = ['-%s=%s' % (k,v) for k,v in self['blatargs'].items()]
        args.extend(blatargs_list)
        infiles = write_fastas_to_files([self['subject'],self['query']])
        args.extend(infiles)
        outfile = os.path.join(self['tempdir'],self['name']+'.blat')
        self['outfile'] = outfile
        args.append(outfile)

        if self['logfile']: open(self['logfile'],'a').write('finished compose_arguments\n')
        return(args)

    def handle_results(self):
        if self['logfile']: open(self['logfile'],'a').write('starting handle_results\n')
        import GFF
        fh = open(self['outfile'])
        lines = fh.readlines()[5:]
        Qgff = GFF.File()
        Sgff = GFF.File()
        for line in lines:
            Qreg = GFF.Region()
            Sreg = GFF.Region()

            fields = line.split('\t')
            [Qlen,Slen] = [float(fields[10]), float(fields[14])]
            [Qhitlen, Shitlen] = [ (int(fields[12]) - int(fields[11])) , (int(fields[16]) - int(fields[15])) ]
            [Qreg['seqid'], Qreg['start'], Qreg['end'], Qreg['strand'], Qreg['score'], Qreg['attribute_coverage']] = \
                [fields[9], fields[11], fields[12], '+', float(fields[0])/Qhitlen, Qhitlen / Qlen]
            [Sreg['seqid'], Sreg['start'], Sreg['end'], Sreg['strand'], Sreg['score'], Sreg['attribute_coverage']] = \
                [fields[13], fields[15], fields[16], fields[8], float(fields[0])/Shitlen, Shitlen / Slen]
            Qgff.append(Qreg)
            Sgff.append(Sreg)

        Controller.handle_results(self)
        if self['logfile']: open(self['logfile'],'a').write('finished handle_results\n')
        return (Qgff, Sgff)


class Bl2seqController(Controller):
    '''Runs bl2seq on two (Seq.Fasta,seqid) tuples 'query' and 'subject' with blast args 'blastargs' (dict)
    

    returns qgff, sgff tuple (a la BlatController)
    '''

    controller_type = 'bl2seq'

    command = paths[controller_type]

    def __init__(self,query=None,subject=None,blastargs={'-p':'blastn'},tempdir=None,stdin=None,stderr=None,stdout=subprocess.PIPE,
                 executor=Popen,use_defIOE=0,name=None,cleanup=0,args=None,logfile=None):

        Controller.__init__(self,tempdir,stdin,stderr,stdout,executor,use_defIOE,name,cleanup,args,logfile)

        self['query'] = query
        self['subject'] = subject
        self['blastargs'] = blastargs
        

    def prepare_fastas(self):
        q = os.path.join(self['tempdir'],self['name']+'-q.fa')
        s = os.path.join(self['tempdir'],self['name']+'-s.fa')
        self['files'].extend([q,s])
        self['query'][0].write_to_file(q,ids=[self['query'][1]])
        self['subject'][0].write_to_file(s,ids=[self['subject'][1]])
        return (q,s)

    def compose_arguments(self):
        if self['logfile']: open(self['logfile'],'a').write('starting compose_arguments\n')
        from Seq import write_fastas_to_files

        args = [self.command]
        blastargs_list = reduce(lambda x,y:x+y, list(self['blastargs'].items()))
        args.extend(blastargs_list)
        args.extend(['-D','1'])
        infiles = self.prepare_fastas()
        args.extend(reduce(lambda x,y:x+y, zip(['-i','-j'],infiles)))
        #outfile = os.path.join(self['tempdir'],self['name']+'.blat')
        #self['outfile'] = outfile
        #args.append(outfile)

        if self['logfile']: open(self['logfile'],'a').write('finished compose_arguments\n')

        return(args)
                    
    def handle_results(self):
        if self['logfile']: open(self['logfile'],'a').write('starting handle_results\n')
        import GFF
        lines = self['out'][0].split('\n')
        Qgff = GFF.File()
        Sgff = GFF.File()
        i = 0
        for line in lines:

            fields = line.split('\t')
            if line.startswith('#') or len(fields) != 12: continue            

            Qreg = GFF.Region()
            Sreg = GFF.Region()


            
            if fields:
                [Qreg['seqid'], Qreg['start'], Qreg['end'], Qreg['strand'], Qreg['score'],Qreg['attribute_e'],Qreg['attribute_ID']] = \
                    [fields[0], fields[6], fields[7], '+', float(fields[11]),fields[10],str(i)]

                if fields[8] < fields[9]:
                    s_start,s_end,s_strand = fields[8], fields[9], '+'
                else:
                    s_start,s_end,s_strand = fields[9], fields[8], '-'
                [Sreg['seqid'], Sreg['start'], Sreg['end'], Sreg['strand'], Sreg['score'],Sreg['attribute_e'],Sreg['attribute_ID']] = \
                    [fields[1], s_start, s_end, s_strand, float(fields[11]),fields[10],str(i)]
                Qgff.append(Qreg)
                Sgff.append(Sreg)
                i+=1
        Controller.handle_results(self)
        if self['logfile']: open(self['logfile'],'a').write('finished handle_results\n')
        return (Qgff, Sgff)
    
class MummerController(Controller):
    '''Runs mummer on two Seq.Fasta instances 'query' and 'subject' with blast args 'blastargs' (dict)
    
    returns qgff, sgff tuple (a la BlatController)
    '''

    controller_type = 'mummer'

    command = 'mummer'

    def __init__(self,query=None,subject=None,mumargs={'-l':'20'},tempdir=None,stdin=None,stderr=None,stdout=subprocess.PIPE,
                 executor=Popen,use_defIOE=0,name=None,cleanup=0,args=None,logfile=None):

        Controller.__init__(self,tempdir,stdin,stderr,stdout,executor,use_defIOE,name,cleanup,args,logfile)

        self['query'] = unwrap_if_gz(query)
        self['subject'] = unwrap_if_gz(subject)
        self['mumargs'] = mumargs
        self['baseargs'] = ['-mum', '-b', '-F', '-c']

    def compose_arguments(self):
        if self['logfile']: open(self['logfile'],'a').write('starting compose_arguments\n')
        from Seq import write_fastas_to_files

        args = [self.command]
        mumargs_list = reduce(lambda x,y:x+y, list(self['mumargs'].items()))
        args.extend(mumargs_list)
        args.extend(self['baseargs'])
        infiles = write_fastas_to_files([self['subject'],self['query']])
        args.extend(infiles)

        if self['logfile']: open(self['logfile'],'a').write('finished compose_arguments\n')

        return(args)
                    
    def handle_results(self):
        if self['logfile']: open(self['logfile'],'a').write('starting handle_results\n')
        import GFF
        lines = self['out'][0].split('\n')
        Qgff = GFF.File()
        Sgff = GFF.File()
        i = 0
        Qseqid = ''
        Qstrand = ''
        for line in lines:
#            print >>sys.stderr,line                    
            if line.startswith('#'): 
                continue

            if line.startswith('>'):
                Qseqid = line[1:].strip()
                if Qseqid.endswith('Reverse'):
                    Qseqid = Qseqid[:(-1 * len('Reverse'))].strip()
                    Qstrand = '-'
                else:
                    Qstrand = '+'
#                print >>sys.stderr,Qstrand
                continue

            fields = line.split()
            if len(fields) == 4:


                fields[1:] = [int(field) for field in fields[1:]]
                
                Qreg = GFF.Region()
                Sreg = GFF.Region()

                if Qstrand == '-':
                    Qstart = fields[2] - fields[3]
                    Qend = fields[2]
                else:
                    Qstart = fields[2]
                    Qend = fields[2] + fields[3]

                [Qreg['seqid'], Qreg['start'], Qreg['end'], Qreg['strand'], Qreg['score'],Qreg['attribute_ID']] = \
                    [Qseqid, Qstart, Qend, Qstrand, fields[3], str(i)]

                [Sreg['seqid'], Sreg['start'], Sreg['end'], Sreg['strand'], Sreg['score'],Sreg['attribute_ID']] = \
                    [fields[0], fields[1], fields[1]+fields[3], '+', fields[3],str(i)]
#                print >>sys.stderr,Qseqid, Qstart, Qend, Qstrand, fields[3], str(i),fields[0], fields[1], fields[1]+fields[3], '+', fields[3],str(i)
                Qgff.append(Qreg)
                Sgff.append(Sreg)
                i+=1

        Controller.handle_results(self)
        if self['logfile']: open(self['logfile'],'a').write('finished handle_results\n')
        return (Qgff, Sgff)



class RefListAlnPipelineMetacontroller():
    '''chains together blat, blastz and tba/roast to generate alignments for a list of single-sequence-containing Seq.Fasta objects
    and a list of filenames of comparative species' sequence

    NB: uses ref_seq_list object.filename to determine where to put results, 
    so this should be set for each member of the list before calling

    wraps functionality in a .run() method for compatiblity with other tools built for Controller objects
    '''

    def __init__(self, ref_seq_list, comp_sp_files_list, blastz_specs=None):
        self.ref_seq_list = ref_seq_list
        self.comp_sp_files_list = comp_sp_files_list
        self.blastz_specs = blastz_specs
        base = os.path.basename(ref_seq_list[0].filename)
        base = base.strip('.gz')
        self.ref_sp = os.path.splitext(base)[0]

    def run(self,null=None):
        from copy import deepcopy

        fastas_di = {}
        for comp_sp in self.comp_sp_files_list:
            print >> sys.stderr,'loading sequence',comp_sp,
            if os.path.exists(comp_sp+'.pkl'):
                comp_seq = cPickle.load(open(comp_sp+'.pkl','rb'))
            else:
                comp_seq = Fasta(comp_sp)
            print >> sys.stderr,'done'
            for i,ref_seq in enumerate(self.ref_seq_list):
                try: #hack to get around "TypeError: 'NoneType' object is unsubscriptable" --should get fixed for real once engine debugging is set up
                    tempdir = os.path.dirname(ref_seq.filename)
                    name = ref_seq.seq_names()[0]
                    ss_file = os.path.join(tempdir, os.path.basename(comp_seq.filename))
                    ss_file_noext = os.path.splitext(ss_file.strip('.gz'))[0]
                    print >> sys.stderr,'dir:',tempdir,'name:',name,'ss_file:',ss_file,'ss_file_noext:',ss_file_noext
                    if os.path.exists(ss_file_noext) or os.path.exists(ss_file):
                        print >> sys.stderr,'genome-wide anchors for', name,'in',comp_sp,'present, skipping'
                    else:
                        print >> sys.stderr,'finding genome-wide anchors for', name,'in',comp_sp,'[%s/%s]' % (i+1,len(self.ref_seq_list)),
                        subseq = find_subset_matching_reference(ref_seq,comp_seq,tempdir=tempdir,blat_args={'minIdentity':'0.85'})
                        subseq.filename = ss_file
                        subseq.write_to_file()
                        print >> sys.stderr,'done'
                    try:
                        fastas_di[name].append(ss_file)
                    except KeyError:
                        fastas_di[name] = [ss_file]
                except TypeError:
                    pass

        results = []

        for i, ref_seq in enumerate(self.ref_seq_list):
            tempdir = os.path.dirname(ref_seq.filename)
            if not os.path.exists(os.path.join(tempdir,'tba*.maf')):
                try: #hack to get around "TypeError: 'NoneType' object is unsubscriptable" --should get fixed for real once engine debugging is set up
                    name = ref_seq.seq_names()[0]
                    print >> sys.stderr,'computing alignment for', name,'[%s/%s]' % (i+1,len(self.ref_seq_list)),
                    mafs = AllBZController(fastas=deepcopy([ref_seq]+fastas_di[name]),
                                           tempdir=tempdir,
                                           ref_sp=self.ref_sp,
                                           blastz_specs=self.blastz_specs,
                                           name='AllBZ-'+name,
                                           use_defIOE=1).run(now=1)
                    alns = RoastController(mafs=mafs,
                                           tempdir=tempdir,
                                           ref_sp=self.ref_sp,
                                           name='Roast-'+name,
                                           use_defIOE=1).run(now=1)
                    results.append(alns)
                    print >> sys.stderr,'done'
                except TypeError:
                    pass

        return results
                

##########
# beyond primary sequence homology
#  (e.g. mapping, etc.)
##########

def blast(query,subject,**kwargs):
    allQgff = GFF.File()
    allSgff = GFF.File()

    for seqid in query.seq_names():
        (Qgff,Sgff) = Bl2seqController((query,seqid),(subject,subject.seq_names()[0]),stderr=open('/dev/null','w'),**kwargs).run(1)
        allQgff.extend(Qgff)
        allSgff.extend(Sgff)
        
    return (allQgff,allSgff)

def mum(query,subject,**kwargs):
    allQgff = GFF.File()
    allSgff = GFF.File()

    if len(subject.seq_names()) != 1:
        raise ValueError, 'mapping with mummer is only confirmed to work with a single subject seq'

    if query.filename is None:
        query.filename = Util.random_filename(prefix='/n/home/brantp/tmp/',suffix='-mumquery.fa')
    if subject.filename is None:
        subject.filename = Util.random_filename(prefix='/n/home/brantp/tmp/',suffix='-mumsubject.fa')
        wipesubj=True
    else:
        wipesubj=False

	print >>sys.stderr, 'q: %s s: %s' % (query.filename,subject.filename)

    for seqid in query.seq_names():
        thisq = Fasta()
        thisq.filename = '/tmp/%s-%s.fa' % (os.path.basename(query.filename),seqid)
        thisq[seqid] = query[seqid]
        (Qgff,Sgff) = MummerController(thisq,subject,stderr=open('/dev/null','w'),**kwargs).run(1)
        allQgff.extend(Qgff)
        allSgff.extend(Sgff)

    if wipesubj:
        print >> sys.stderr, 'removing %s' % subject.filename
        os.unlink(subject.filename)
        subject.filename=None

    return (allQgff,allSgff)

def mumblast(query,subject,**kwargs):
    '''calls both blast and mum, 

    passes blastkwargs to blast, mumkwargs to mum
    NOTE: if either is unspecified, sets mapping-specific defaults defined below
    '''

    default_blastargs = {'blastargs': {'-e': '10', '-p': 'blastn'}}
    default_mumargs = {'mumargs': {'-l': '10'}}

    allQgff = GFF.File()
    allSgff = GFF.File()
    
    (Qgff,Sgff) = blast(query,subject,**kwargs.get('blastkwargs',default_blastargs))
    
    allQgff.extend(Qgff)
    allSgff.extend(Sgff)
    
    (Qgff,Sgff) = mum(query,subject,**kwargs.get('mumkwargs',default_mumargs))

    allQgff.extend(Qgff)
    allSgff.extend(Sgff)
    
    return (allQgff,allSgff)

def map_by_aggregate_short_matches(query,subject,Qseqid=None,Sseqid=None,
                                   winsize=1000,winstep=None,
                                   method=mumblast,method_args={},
                                   scorefunct=lambda li: sum([r['score'] for r in li]), 
                                   return_seqlens=False, return_reverse=False, iplot_style=False):

    '''map from AlnRefMaf query to Fasta subject by short sequences matches aggregated in windows

    default winstep is winsize/2

    currently has no notion of precise sequence coords of non-reference orthologous sequence matches
    (i.e. matches from the Dpse aligned sequence to a given Dmel interval will be scored in that interval,
    but not to a specific base position)

    if multiple sequence Fasta or multi-region AlnRefMaf is supplied, qseqid or sseqid must be provided

    method can be any valid function that takes a query Fasta (which can have many sequences) and a single seq subject Fasta
    it must return a Qgff,Sgff tuple pair, a la BlatController

    scorefunct takes the subject GFF list and returns a numerical score (default: sum of region scores)
    '''
    
    def get_AlnRefMaf_seqid(inst):
        if len(inst.assoc_arrays.keys()) == 1:
            return inst.assoc_arrays.keys()[0]
        else:
            raise ValueError, 'instance has more than one region'
        
    def get_Fasta_seqid(inst):
        if len(inst.keys()) == 1:
            return inst.keys()[0]
        else:
            raise ValueError, 'instance has more than one region'

    Qmap = GFF.File()
    Smap = GFF.File()

    if winstep is None:
        winstep = winsize/2

    if isinstance(query,AlnRefMaf) and isinstance(subject,Fasta):
        if Qseqid is None: Qseqid = get_AlnRefMaf_seqid(query)
        if Sseqid is None: Sseqid = get_Fasta_seqid(subject)
        
        if len(subject.seq_names()) == 1:
            one_seq_subject = subject
        else:
            one_seq_subject = Fasta()
            one_seq_subject.filename = '/tmp/%s-%s.fa' % (os.path.basename(subject.filename),Sseqid)
            one_seq_subject[Sseqid] = subject[Sseqid]

        Qsize = len(query.assoc_arrays[Qseqid].array[0])
        Ssize = len(one_seq_subject[Sseqid])
        
        print >>sys.stderr, 'map from %s (%sbp) to %s (%sbp) in %sbp windows, step %sbp' % (Qseqid,Qsize,Sseqid,Ssize,winsize,winstep)
        
        pair_idx = 1
        for Qwinstart in xrange(0,Qsize-winsize,winstep):
            Qwinend = Qwinstart + winsize

            (Qgff,Sgff) = method(query.get_seqs_in_reference_range(Qwinstart,Qwinend,seqname=Qseqid),one_seq_subject,**method_args)

            for Swinstart in xrange(0,Ssize-winsize,winstep):
                Swinend = Swinstart + winsize

                Swinhits = [ r for r in Sgff if Swinstart <= r['start'] and r['end'] <= Swinend ]

                if Swinhits:
                    new_score = scorefunct(Swinhits)
                    Qwinreg = GFF.Region(fields={'seqid':Qseqid,'start':Qwinstart,'end':Qwinend,'score':new_score,'attribute_ID':str(pair_idx)})
                    Swinreg = GFF.Region(fields={'seqid':Sseqid,'start':Swinstart,'end':Swinend,'score':new_score,'attribute_ID':str(pair_idx)})
                    #print >>sys.stderr,'adding',Qwinreg
                    #print >>sys.stderr,'adding',Swinreg
                    Qmap.append(Qwinreg)
                    Smap.append(Swinreg)
                    pair_idx += 1

        if iplot_style: return_seqlens = return_reverse = True
        if return_seqlens:
            if return_reverse:
                return (Smap,Qmap),(Ssize,Qsize)
            else:
                return (Qmap,Smap),(Qsize,Ssize)
        else:
            if return_reverse:
                return (Smap,Qmap)
            else:
                return (Qmap,Smap)
            
            
            


##########
# evaluate and display
# (some items in this section 
#    may move to iplot)
##########

class PhastconsController (AlignmentController):
    '''class to run phastcons and associated programs

    intended use as "pipeline" controller for msa_view, phyloFit and phastCons'''

    controller_type = 'phastCons'


    def __init__(self,aln=None,rho='0.3',gff_out=None,bg_aln=None,phy_out=None,
                 tree=paths["PD_tree"],newickstr=None,
                 tempdir=None,stdin=None,stderr=None,stdout=None,
                 executor=Popen,use_defIOE=0,name=None,cleanup=0,args=None):

        if rho:
            rho = str(rho)

        if stdout is None:
            if aln:
                stdout = os.path.join(tempdir,os.path.basename(aln)+'_pcons-rho'+rho+'.scores')
        if gff_out is None:
            if aln:
                gff_out = os.path.join(tempdir,os.path.basename(aln)+'_pcons-rho'+rho+'.gff')

        AlignmentController.__init__(self,aln,tree,newickstr,
                                     tempdir,stdin,stderr,stdout,
                                     executor,use_defIOE,name,cleanup,args)
        self['rho'] = rho
        self['gff_out'] = gff_out

        if phy_out:
            self['phy_out'] = phy_out
        elif aln and aln.endswith('.maf'):
            self['phy_out'] = aln
            self['msa-format'] = 'MAF'
        else:
            self['phy_out'] = os.path.join(self['tempdir'],self['name']+'align.phylip')
            self['msa-format'] = 'PHYLIP'
        if bg_aln:
            self['bg_aln'] = bg_aln
        else:
            self['bg_aln'] = self['phy_out']
        
        self['mod-root'] = os.path.join(tempdir,os.path.basename(self['bg_aln']))
        self['mod'] = self['mod-root']+'.mod'



    def run_sub_app(self,args,stdin=None,stdout=None,stderr=None,now=1,cwd=None):

        print >> sys.stderr, "calling %s for stdo %s stde %s\n" % (self["executor"],stdout,stderr)

        if cwd is None: 
            cwd = self["tempdir"]

        if stderr is None:
            stderr = os.path.join(self['tempdir'],
                                  self['name']+'-'+args[0].rsplit('/',1)[-1]+'.stderr')

        if isinstance(stdout,str):
            stdout = open(stdout,'w')
        if isinstance(stderr,str):
            stderr = open(stderr,'a')
        
        print >> sys.stderr, "running %s in %s" % (' '.join(args),cwd)
        if now:
            self["executor"](args,
                             stdin=stdin,
                             stdout=stdout,
                             stderr=stderr,
                             cwd=cwd).wait()
        else:
            return self["executor"](args,
                                    stdin=stdin,
                                    stdout=stdout,
                                    stderr=stderr,
                                    cwd=cwd)
        

    def msa_view(self,now=1,cwd=None):
        
        if isinstance(self['aln'],str):
            alignment = self['aln']
        elif isinstance(self['aln'],AlnFasta):
            alignment = self['aln'].filename
        args = ['msa_view',alignment,'--in-format', 'FASTA', '--out-format', 'PHYLIP']

        return self.run_sub_app(args=args,stdout=self['phy_out'],now=now)

    def phyloFit(self,now=1,cwd=None):

        if not os.path.exists(self['bg_aln']):
            self.msa_view(now=1)

        if self['msa-format'] == 'MAF':
            #need to get rid of all '.' after the sp-sequence separator
            lines = open(self['bg_aln']).readlines()
            fh = open(self['bg_aln'],'w')
            for l in lines:
                match = re.search(r'^s\s.+?\.(.+?)\s',l)
                if match:
                    seq = match.groups()[0]
                    l = re.sub(seq,seq.replace('.','__'),l)
                fh.write(l)
            fh.close()

        args = ['phyloFit',self['bg_aln'],'--tree',self.treestr(),
                '--subst-mod=HKY85+Gap', '--msa-format', self['msa-format'], '--out-root',self['mod-root']]
        
        return self.run_sub_app(args=args,now=now)
        
    def compose_arguments(self):
        if not os.path.exists(self['bg_aln']+'.mod'):
            self.phyloFit(now=1)

        args = ['phastCons',self['phy_out'],'--msa-format',self['msa-format']]
        if self['gff_out']:
            args.extend(['--most-conserved',self['gff_out']])
        args.extend(['--score', '--refidx', '1', '--rho', self['rho'], self['mod']])

        #expected-length = 23.8; target-coverage = 0.393        
        #args.extend(['--score', '--refidx', '1', '--target-coverage','0.393','--expected-length','23.8', '--rho', self['rho'],self['bg_aln']+'.mod'])
        return args

def phastcons_list(fa_list,rho='0.3',tempdir=None):
    for f in fa_list:
        pconer = PhastconsController(aln=f,rho=rho,tempdir=tempdir)
        run = pconer.run(1)

def pairwise_id(seq1,seq2, cut=0):
    ids = []
    for i,v in enumerate(seq1):
        if v == seq2[i]: 
            ids.append(1)
        else:
            ids.append(0)
    
    id = sum(ids)/float(len(ids))
    if id > cut:
        return id
    else:
        return cut


def id_matrix(seq1,seq2, window=20, cut=0.6, rev=False, revcomp_seq2=False,normed=False):
    if revcomp_seq2:
        seq2 = seq2.rc()
        j_pos = 'len(m[i])-j-1'
    else:
        j_pos = 'j'

    if rev:
        m = numpy.ones([len(seq1),len(seq2)],float)
    else:
        m = numpy.zeros([len(seq1),len(seq2)],float)
    for i in range(0,len(m)-window):
        r = []
        for j in range(0,len(m[i])-window):
            v = pairwise_id(seq1[i:i+window],seq2[j:j+window])
            if v <= cut:
                v = 0
            else:
                if normed:
                    normv = (v-cut)/(1-cut)
                    v = normv
            if rev:
                v = 1-v

            m[i][eval(j_pos)] = v

    return m

def gff_matrix(gff1,gff2,seq_lens,rev=False,norm_lower=0.1,min_topscore=30):
    '''generate a matrix representing matches between 2 sequences in a pair of matched-rank gffs

    assumes only 1 sequence per GFF list
    '''
    scores = numpy.array([r['score'] for r in gff1])
    print >> sys.stderr, scores,scores.min(),max(scores.max(),min_topscore)
    #scores = scores - scores.min()
    
    newmax = max(scores.max(),min_topscore)
    shiftdown = scores.min() - (norm_lower*newmax)
    print >> sys.stderr,scores,newmax
    scores = (scores - shiftdown) / newmax
    #scores = Util.normalize(scores,norm_lower)
    
    if rev:
        scores = 1 - scores

    n1,n2 = gff1[0]['seqid'], gff2[0]['seqid']

    if rev:
        f = numpy.ones([seq_lens[n1],seq_lens[n2]],float)
        r = numpy.ones([seq_lens[n1],seq_lens[n2]],float)
    else:
        f = numpy.zeros([seq_lens[n1],seq_lens[n2]],float)
        r = numpy.zeros([seq_lens[n1],seq_lens[n2]],float)

    for i,reg in enumerate(gff1):
        if gff2[i]['strand'] == '-':
            for subi in range(len(reg)):
                try:
                    r[reg['start']+subi][gff2[i]['start']-subi] = scores[i]
                except IndexError:
                    print sys.stderr >> len(r),reg['start']+subi,len(gff2),i,len(scores)
        else:
            for subi in range(len(reg)):
                try:
                    r[reg['start']+subi][gff2[i]['start']+subi] = scores[i]
                except IndexError:
                    print >> sys.stderr, len(r),reg['start']+subi,subi,len(gff2),i,len(scores)

    return (f,r)
                

def bothstrands_id_matrix(seq1,seq2, window=20, cut=0.6, rev=False,normed=False):
    from copy import deepcopy
    f = id_matrix(seq1,seq2,rev=rev,cut=cut,window=window)
    r = id_matrix(seq1,deepcopy(seq2),rev=rev,cut=cut,window=window,revcomp_seq2=True)
    return (f,r)

#DOTPLOT MOVED TO iplot


if __name__ == "__main__":
    import os, time,glob
    testfile = os.path.join(os.path.expanduser('~'),"py_util/unit_test_data/test.aln")

    if os.path.exists(testfile):
        test_aln = AlnFasta(testfile)
        print test_aln.seq_len()
    else:
        print "no file at %s" % (testfile,)

    print "-----test mlagan-----\n"
    import Seq
    fasta_file = os.path.join(paths["py_testdata"],"seq.fasta")
    fasta_obj = Seq.Fasta(fasta_file)
    from subprocess import PIPE
    import SGE
    laganer = LaganController(fasta_file,use_defIOE=1,cleanup=0,executor=SGE.Qopen)
    print laganer.lagan_treestr()
    print laganer["fasta"].seq_len()
    print str(laganer.prepare_fastas())
    print ' '.join(laganer.compose_arguments())
    runner = laganer.run()
    while runner.poll() != 0:
        print '.'
        time.sleep(10)
    align = laganer.handle_results()
    print align.seq_len()
    #print align

    '''
    print "-----test TBA-----\n\non:\n"
    fasta_files = glob.glob(os.path.join(paths["py_testdata"],"eve.*.fa"))
    print fasta_files
    '''
    
    '''
    print "-----all_bz then TBA-----\n
    bzer = AllBZController(fasta_files,use_defIOE=1,cleanup=0)
    #    print bzer
    bzer['blastz_specs'] = [('C',1)]
    print bzer['blastz_specs']
    print bzer['tempdir']
    mafs = bzer.run(1)
    print bzer['files']
    print bzer.calc_output_sp_pairs()

    tbaer = TBAController(mafs=mafs,use_defIOE=1,cleanup=0)
    print tbaer
    tbaer.run(1)

    del bzer
    del tbaer
    '''

    '''
    print "-----TBA with fasta files-----\n"
    tbaer = TBAController(fastas=fasta_files,use_defIOE=1,cleanup=1)
    print tbaer
    tbaer.run(1)
    '''
