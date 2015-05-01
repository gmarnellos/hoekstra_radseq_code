#!/usr/bin/env python

'''
give fastq file(s) (1 or 4 line format)
preprocess_radtag_lane.py /path/to/flowcell/s_N_1_sequence.txt [s_N_2_sequence.txt]

for which individual data is present in <LIBRARY_DATA> gdoc spreadsheet (see config.py)

generates tabular uniqued read data:
ID nreads sequence mean_qual comma,delim,indiv comma,delim,depths_by_ind \
                                  comma,delim,unique,read2,seq comma,delim,unique,read2,counts

NOTES ON FILENAMES:
flowcell is extracted from the name of the directory containing the fastq file(s)
Thus, the path supplied to the fastq file(s) in running preprocess_radtag_lane.py
must include at least the containing folder.  The name of the containing  folder
must match the flowcell designation in <LIBRARY_DATA> gdoc spreadsheet 

'''

import os, sys, re, numpy, gzip
import gdata.spreadsheet.service

from collections import defaultdict
from subprocess import Popen, PIPE
from editdist import distance
from glob import glob
from config import EMAIL,PASS,SOURCE,LIBRARY_DATA,ADAPTER_DATA,RTDROOT

def dezip(values_in):
    '''opposite of zip(), i.e. 
    >>> dezip([("a",1),("b",2),("c",3)])
    (["a","b","c"],[1,2,3])'''

    if isinstance(values_in,tuple):
        values_in = [values_in]
	
    lol = []
    for i in range(len(values_in[0])):
        lol.append([])
    for l in values_in:
        for it,li in zip(l,lol):
	    li.append(it)
    return tuple(lol)


def smartopen(filename,*args,**kwargs):
    '''opens with open unless file ends in .gz, then use gzip.open

    in theory should transparently allow reading of files regardless of compression
    '''
    if filename.endswith('.gz'):
        return gzip.open(filename,*args,**kwargs)
    else:
        return open(filename,*args,**kwargs)


def get_read_count(filename,lnum=None,use_cache=True):

    if use_cache:
        rcc = filename+'.rc.cache'
        try:
            filesize,rc = open(rcc).readline().strip().split()
            if float(filesize) == os.path.getsize(filename):
                print >> sys.stderr, 'read count from cached value: %s' % rc
                return int(rc)
        except:
            pass

    if lnum is None:
        if smartopen(filename).read(1) == '@':
	    lnum = 4
	else:
	    lnum = 1
	
    if filename.endswith('.gz'):    
        print >> sys.stderr, 'getting read count for compressed file',filename,'...',
        rc = int(Popen('gunzip -c %s | wc -l' % filename,shell=True,stdout=PIPE).stdout.read().split()[0]) / lnum
        print >> sys.stderr, rc
    else:
        print >> sys.stderr, 'getting read count for file',filename,'...',
        rc = int(Popen('wc -l %s' % filename,shell=True,stdout=PIPE).stdout.read().split()[0]) / lnum
        print >> sys.stderr, rc

    if use_cache:
        open(rcc,'w').write('%s\t%s\n' % (os.path.getsize(filename),rc))
        
    return rc

def get_baseQ(qstr):
    q = [ord(c) for c in qstr]
    if any([i<66 for i in q]):
        return 33
    elif any([i>74 for i in q]):
        return 64
    else:
        return None

def fq_path_from_db_dict(db_dict,index_lookup='DB_multiplex_indices'):
    if type(index_lookup) == str:
        index_lookup =  get_table_as_dict(index_lookup,suppress_fc_check=True)
    idxlookup = dict([(d['idx'],d['seq']) for d in index_lookup])
    idxnumlookup = dict([(d['seq'],d['idx']) for d in index_lookup])

    cand = glob(db_dict['datapath']+'/*')
    r1 = None
    r2 = None
    for c in cand:
        if os.path.basename(c).startswith('Sample') and 'lane%s' % db_dict['lane'] in os.path.basename(c).lower() and os.path.basename(c).endswith('.R1.fastq.gz'):
            if db_dict.has_key('index'):
                if db_dict['index'].isdigit():
                    alt_idx = idxlookup
                    idx = db_dict['index']
                else:
                    if '-' in db_dict['index']:
                        idx = db_dict['index'].split('-')[0]
                    else:
                        idx = db_dict['index']
                    alt_idx = idxnumlookup
                if os.path.basename(c).endswith('%s.R1.fastq.gz' % idx) or os.path.basename(c).endswith('%s.R1.fastq.gz' % alt_idx[idx]):
                    r1 = c
                    if os.path.exists(c.replace('.R1.','.R2.')):
                        r2 = c.replace('.R1.','.R2.')
                    break
            else:
                if os.path.basename(c).endswith('None.R1.fastq.gz') or os.path.basename(c).endswith('noidx.R1.fastq.gz'):
                    r1 = c
                    if os.path.exists(c.replace('.R1.','.R2.')):
                        r2 = c.replace('.R1.','.R2.')
                    break
        if os.path.basename(c).startswith('s_%s_1_sequence' % db_dict['lane']) and os.path.basename(c).endswith('.txt.gz'):
            if db_dict.has_key('index'):
                if db_dict['index'].isdigit():
                    alt_idx = idxlookup
                    idx = db_dict['index']
                else:
                    if '-' in db_dict['index']:
                        idx = db_dict['index'].split('-')[0]
                    else:
                        idx = db_dict['index']
                    alt_idx = idxnumlookup
                if os.path.basename(c).endswith('index%s.txt.gz' % idx) or os.path.basename(c).endswith('index%s.txt.gz' % alt_idx[idx]):
                    r1 = c
                    if os.path.exists(c.replace('_1_sequence','_2_sequence')):
                        r2 = c.replace('_1_sequence','_2_sequence')
                    break
            else:
                r1 = c
                if os.path.exists(c.replace('_1_sequence','_2_sequence')):
                    r2 = c.replace('_1_sequence','_2_sequence')
                break
    
    return r1,r2

def get_legacy_to_DB_lookup(table_dict,mouseDB = 'Hoekstra lab mouse database',mouse_sheet='Mice',prefix_for_missing='BC',remove_suffixes=['_pilot']):
    #db_td = get_table_as_dict(mouseDB,target_worksheet=mouse_sheet,suppress_fc_check=True)
    db_td = no_net_get_table_as_dict(mouseDB)
    failures = []
    transtable = {}
    #identity pass-through for "DB" records
    for sid in set([d['sampleid'] for d in table_dict if d['idtype'] == 'DB']):
        transtable[sid] = sid
    #find translation for "legacy" records
    for sid in set([d['sampleid'] for d in table_dict if d['idtype'] == 'legacy' and not d.has_key('altid')] + [d['altid'] for d in table_dict if d.has_key('altid') and d['altidtype'] == 'legacy']): #ADD ALTID TO SET
        for suff in remove_suffixes:
            if sid.endswith(suff):
                oldsid=sid
                sid = sid[:-len(suff)]
                print >> sys.stderr, 'found suffix %s in %s; now %s' % (suff,oldsid,sid)
                break
        legid = sid
        candidate_records = [rec for rec in db_td if rec.get('legacyaka','') == legid]
        if len(candidate_records) != 1 and sid[0].isdigit():
            legid = prefix_for_missing + sid
            candidate_records = [rec for rec in db_td if rec.get('legacyaka','') == legid]
        if len(candidate_records) == 1:
            transtable[sid] = candidate_records[0]['id']
        else:
            #print legid,candidate_records
            failures.append((sid,candidate_records))

    return transtable,failures
                                                                           

def create_empty_table(table_name):
    try:
        key, gd_client = get_spreadsheet_key(table_name)
        print >> sys.stderr, 'table %s exists, skip' % table_name
    except:
        import gdata.docs.data
        import gdata.docs.client
        client = gdata.docs.client.DocsClient(source=SOURCE)
        client.ssl = True  # Force all API requests through HTTPS
        client.http_client.debug = False  # Set to True for debugging HTTP requests
        
        client.ClientLogin(EMAIL,PASS,client.source)
        
        new_spreadsheet = client.Create(gdata.docs.data.SPREADSHEET_LABEL, table_name , writers_can_invite=False)
        print >> sys.stderr, 'Spreadsheet "%s" created' % new_spreadsheet.title.text


def get_spreadsheet_key(target_sheet,gd_client=None):
    '''returns the key string for a spreadsheet given its name'''

    if gd_client is None:
        gd_client = gdata.spreadsheet.service.SpreadsheetsService()
        gd_client.email = EMAIL
        gd_client.password = PASS
        gd_client.source = SOURCE

    gd_client.ProgrammaticLogin()

    feed = gd_client.GetSpreadsheetsFeed()
    key = [entry.id.text.rsplit('/', 1)[1] for entry in feed.entry if entry.title.text == target_sheet][0]

    return key,gd_client

def get_worksheet_key(target_worksheet,sheet_key,gd_client):
    feed = gd_client.GetWorksheetsFeed(sheet_key)
    key = [entry.id.text.rsplit('/', 1)[1] for entry in feed.entry if entry.title.text == target_worksheet][0]

    return key,gd_client

def get_all_worksheet_names(sheet_key,gd_client):
    feed = gd_client.GetWorksheetsFeed(sheet_key)
    names = [entry.title.text for entry in feed.entry]
    return names,gd_client

def get_table_as_dict(target_sheet,sq=None,gd_client=None,suppress_fc_check=False,target_worksheet=None):
    key,gd_client = get_spreadsheet_key(target_sheet,gd_client)
    if sq is not None:
        filter_el = sq.split(' and ')
        filter_pairs = [s.replace('"','').split('=') for s in filter_el]
        filter_fn_str = r'lambda d: all([d.get(k,"") == v for k,v in %s])' % filter_pairs
        print >> sys.stderr, 'structured queries no longer work (thanks google!)\n' +\
              'replaced by list filter:\n%s' % filter_fn_str
        filter_fn = eval(filter_fn_str)

    if target_worksheet is None:
        names = get_all_worksheet_names(key,gd_client)[0]
        recs = []
        for sheet_name in names:
            print >> sys.stderr, 'add %s' % sheet_name
            recs.extend(get_table_as_dict(target_sheet,sq=sq,gd_client=gd_client,suppress_fc_check=suppress_fc_check,target_worksheet=sheet_name))
        return recs
    else:
        wskey,gd_client = get_worksheet_key(target_worksheet,key,gd_client)
    #if sq is not None:
    #    q = gdata.spreadsheet.service.ListQuery()
    #    q.sq = sq
    #    feed = gd_client.GetListFeed(key,wskey,query=q)
    #else:
    #    feed = gd_client.GetListFeed(key,wskey)
    feed = gd_client.GetListFeed(key,wskey)
    
    recs = []
    for entry in feed.entry:
        #d = []
	#for el in entry.content.text.split(','):
	#print entry.content.ToString()
        try:
	    #recs.append(dict(re.findall('([\d\w_-]+?):\s(.+?)(?:(?:,\s)|$)',entry.content.text)))
            l = [m.strip(' ,:') for m in re.split('([^\s]+?:\s)',entry.content.text) if m]
            d = dict(zip(l[::2],l[1::2]))
            recs.append(d)
	    if not suppress_fc_check and not all([k in recs[-1].keys() for k in ['flowcell','lane','pool']]):
	        print >> sys.stderr, 'missing keys:', dict(re.findall('(.+?):\s(.+?)(?:(?:,\s)|$)',entry.content.text))
	        print >> sys.stderr, 'line was:\n',entry.content.text
	except:
	    print >> sys.stderr, 'invalid:', entry.content.text#.split(',')
    if sq is not None:
        print >> sys.stderr, 'loaded %s records' % len(recs)
        recs = filter(filter_fn,recs)
        print >> sys.stderr, 'return %s records' % len(recs)
    return recs

def no_net_get_table_as_dict(target_sheet,host='heroint3'):
    print >> sys.stderr, 'retrieve %s via %s ...' % (target_sheet,host),
    from subprocess import Popen,PIPE
    cmd = 'ssh %s "python -c \\"from rtd.preprocess_radtag_lane import get_table_as_dict; td=get_table_as_dict(\'%s\',suppress_fc_check=True); print td.__repr__()\\"" 2> /dev/null | tail -n 1' % (host,target_sheet)
    #cmd = 'ssh %s "python -c \\"from rtd.preprocess_radtag_lane import get_table_as_dict; td=get_table_as_dict(\'%s\',suppress_fc_check=True); print td.__repr__()\\"" 2> /dev/null' % (host,target_sheet)
    #return cmd
    td_str = Popen(cmd,shell=True,stdout=PIPE).stdout.read()
    #return td_str
    td = eval(td_str)
    print >> sys.stderr, '%s records' % len(td)
    return td

def get_adapter_index_lookup(verbose=True):
    '''returns dict of dicts:
    { <adaptersversion> : { <well> : <index_seq> } }
    '''
    
    key,gd_client = get_spreadsheet_key(ADAPTER_DATA)
    feed = gd_client.GetListFeed(key)

    d = []
    for entry in feed.entry:
        tl = [[st.strip() for st in el.split(':')] for el in entry.content.text.split(',')]
        try:
            d.append(dict(tl))
        except:
            print >> sys.stderr, 'tuple list did not parse:\n\t%s' % tl

    idxlookup = defaultdict(dict)
    for el in d:
        idxlookup[el['adaptersversion']].update({el['well']:el['idx']})
    if verbose==True:
        print >> sys.stderr, 'loaded adapter lookup from %s lines in %s' % (len(d),ADAPTER_DATA) 
    return idxlookup

def get_individual_data_for_lane(filename=None,idxlookup=None,fc=None,lane=None,index=None,transtable=None,get_transtable_from_mouseDB=False):
    '''given a fastq file, treats the directory immediately above as the flowcell ID, returns dict:
    { <sequence_index_tag> : **ROW_FROM_LIBRARY_DATA } 
    '''
    if idxlookup is None:
        idxlookup = get_adapter_index_lookup()

    if fc is None and lane is None:
        fc = os.path.basename(os.path.dirname(filename))
        lane = os.path.basename(filename)[2]
        #print >> sys.stderr, fc,lane,idxlookup
        fbase = os.path.basename(filename)
        #
        #key,gd_client = get_spreadsheet_key(LIBRARY_DATA)
        #
        #q = gdata.spreadsheet.service.ListQuery()
        #
        if 'index' in fbase:
        #    q.sq = 'flowcell="%s" and lane="%s" and index="%s"' % (fc,lane,fbase.split('index')[-1].split('.')[0])
            sq = 'flowcell="%s" and lane="%s" and index="%s"' % (fc,lane,fbase.split('index')[-1].split('.')[0])
        else:
    	#    q.sq = 'flowcell="%s" and lane="%s"' % (fc,lane)
            sq = 'flowcell="%s" and lane="%s"' % (fc,lane)
    else:
        #key,gd_client = get_spreadsheet_key(LIBRARY_DATA)
        #
        #q = gdata.spreadsheet.service.ListQuery()
        if index is not None:
        #    q.sq = 'flowcell="%s" and lane="%s" and index="%s"' % (fc,lane,index)
            sq = 'flowcell="%s" and lane="%s" and index="%s"' % (fc,lane,index)
        else:
        #    q.sq = 'flowcell="%s" and lane="%s"' % (fc,lane)
            sq = 'flowcell="%s" and lane="%s"' % (fc,lane)



    #feed = gd_client.GetListFeed(key,query=q)
    #
    #recs = []
    #
    #for entry in feed.entry:
    #    tl = [[st.strip() for st in el.split(':')] for el in entry.content.text.split(',')]
    #    try:
    #        recs.append(dict(tl))
    #    except:
    #        print >> sys.stderr, 'tuple list did not parse:\n\t%s' % tl

    #recs = [dict([[st.strip() for st in el.split(':')] for el in entry.content.text.split(',')]) for entry in feed.entry]

    #try with table_as_dict
    recs = get_table_as_dict(LIBRARY_DATA,sq=sq)

    if len(recs) == 0:
        raise ValueError, 'no records returned'
    print >> sys.stderr, "%s records found for %s" % (len(recs), sq)
    adaptersversions = set([r['adaptersversion'] for r in recs])
    print >> sys.stderr, "adapters used: %s" % adaptersversions
    idxs = reduce(lambda x,y: x+y, [idxlookup[adver].values() for adver in adaptersversions])
    idxlens = set([len(idx) for idx in idxs])
    if len(idxlens) != 1:
        raise ValueError, 'non-uniform index lengths %s for %s' % (idxlens,filename)

    try:
        sampleids = [r['sampleid'] for r in recs]
    except KeyError:
        try: #permit backup sample ID use
            sampleids = [r['sampleid2'] for r in recs]
        except KeyError:
            print >> sys.stderr, 'not all samples have ID:'
            for d in recs:
                print >> sys.stderr, d.get('sampleid','MISSING'),d
            raise
    wells = [r['adapter'] for r in recs]

    if len(set(sampleids)) != len(sampleids):
        raise ValueError, '%s sampleids, %s unique' % (len(sampleids),len(set(sampleids)))
    if len(set(wells)) != len(wells):
        raise ValueError, '%s wells, %s unique' % (len(wells),len(set(wells)))

    if transtable is None and get_transtable_from_mouseDB:
        transtable,failures = get_legacy_to_DB_lookup(recs)

    indiv_data = {}
    for r in recs:
        if transtable is not None:
            if transtable.has_key(r['sampleid']):
                r['sampleid'] = transtable[r['sampleid']]
                indiv_data[idxlookup[r['adaptersversion']][r['adapter']]] = r
            else:
                print >> sys.stderr, '%s missing from transtable; check database lookup' % r['sampleid']
        else:
            indiv_data[idxlookup[r['adaptersversion']][r['adapter']]] = r

    return indiv_data

def match_index(t,idx_d,idx_len=None,mismatch_allowed=1):
    '''given an index read sequence a dictionary of form {"<read_sequence>":"<index_number>" ...}
    returns <index_number> if the best match is the only index within mismatch_allowed

    DOES NOT check to make sure all indices are idx_len (even if left to get idx_len, i.e. no idx_len supplied)
    '''

    if idx_len is None:
        idx_len = list(set([len(k) for k in idx_d]))[0]
	
    tagdist = sorted([(distance(t_this,t[:idx_len]),t_this) for t_this in idx_d.keys()])
    if tagdist[0][0] <= mismatch_allowed and tagdist[1][0] > mismatch_allowed:
        return idx_d[tagdist[0][1]]
    else:
        return None

def sam_line_to_fastq(samline,idx_field=None,idx_d=None,idx_len=None):
    '''converts a sam line to 4-line fastq string.
    If idx_field is present, matches strings in this field with match_index using idx_d (required) and idx_len (which can be none)
    If idx_field is present, return is (string,idx) else return is string
    
    TO DO:
    CURRENTLY ASSUMES READ 1 (should figure this out from flag)
    '''

    qname,flag,rname,pos,mapq,cigar,mrnm,mpos,tlen,seq,qual,opt = samline.split(None,11)

    opts = dict([(s.split(':')[0],s.split(':')[-1]) for s in opt.split()])

    if idx_field is not None and idx_d is not None and idx_field in opts:
        if idx_len is None:
            idx_len = list(set([len(k) for k in idx_d]))[0]
        fqstr = '@%s#%s/1\n%s\n+\n%s\n' % (qname,opts[idx_field][:idx_len],seq,qual)
        return (fqstr,match_index(opts[idx_field][:idx_len],idx_d,idx_len))
    else:
        fqstr = '@%s#0/1\n%s\n+\n%s\n' % (qname,seq,qual)
        return fqstr
    

def assign_read_to_indiv(line,indiv_data,mismatch_allowed=1, \
		indiv_reads_out_pattern=None,fhdict=None,passfh=None,read2_has_idx=None, \
		trim_Q2=False,min_readlen=None,lnum=4,output_lnum=4,baseQ_in=None,baseQ_out=None):
    '''given a fastq line (actually a list of [read_name,seq,qual_str]), and an indiv_data object (see get_individual_data_for_lane)

    assigns the read to an individual based on the index tag, strips the index sequence and quality positions,
    converts quality to list of integers, and returns the sampleid, sequence and quality
    
    if a pattern is specified for output (like "/path/to/per-indiv-data/%s_s_1_1_sequence.txt")
    will also generate per-individual fastqs.
    
    using a single fhdict and passfh is highly recommended (i.e. creating beforehand and passing as arguments),
    but will be generated if absent.

	FUTURE PLANS:
    if min_readlen is set, will "pass" reads shorter than min_readlen
    if trim_Q2 is True will remove all terminal quality 2 bases.
    If this reduces a read to less then min_readlen good bases, sends to pass
    
    
    returns indiv,read,qual
    
    Paired-Ends (PE) HANDLING:
    if line and indiv_reads_out_pattern are 2-tuples, treats reads as paired-end.
    This requires that read2_has_idx be either True or False
    if False, both reads handled per the index bases of line[0]
    if True, both reads assesssed for index bases, if they DO NOT DISAGREE both reads handled per consensus

    fhdict keys for PE (line is 2-tuple) are 2-tuples (<indiv>,<readnum>) i.e. (BW001,1)

    if passfh supplied, must also be 2-tuple

    returns indiv, (read1, read2), (q1, q2)
    '''

    idxlen = len(indiv_data.keys()[0])

    if isinstance(line, tuple) and len(line) == 2:
        if (isinstance(indiv_reads_out_pattern, tuple) and len(indiv_reads_out_pattern) == 2) or indiv_reads_out_pattern is None:
            if read2_has_idx is not None:
                if indiv_reads_out_pattern is not None:
                    if fhdict is None:
                        fhdict = {}
                    if passfh is None:
                        passfh = [smartopen(p % 'pass','w') for p in indiv_reads_out_pattern]

                indiv = None
                heads = [l[0] for l in line]
                ss = [l[1] for l in line]
                qstrs = [l[2] for l in line]
                
                if baseQ_in is None:
                	bqs = list(set([get_baseQ(qs) for qs in qstrs if get_baseQ(qs) is not None]))
                	if len(bqs) == 1:
                		baseQ_in = bqs[0]
                	else:
                		raise ValueError,'bqs: %s' % bqs
                if baseQ_out is None:
                	baseQ_out = baseQ_in

                if len(set([h.split()[0][:-1] for h in heads])) != 1:
                    raise ValueError, 'read headers not identical prior to last character; %s' % heads 

                if read2_has_idx: #check that indices are concordant
                    ts = [s[:idxlen] for s in ss]
                    tqs = [qstr[:idxlen] for qstr in qstrs]
                    tagdists = [sorted([(distance(t_this,t),t_this) for t_this in indiv_data.keys()]) for t in ts]
                    try:
                        indiv_cand = [indiv_data[tagdist[0][1]]['sampleid'] for tagdist in tagdists \
                                      if tagdist[0][0] <= mismatch_allowed and tagdist[1][0] > mismatch_allowed]
                    except:
                        indiv_cand = [indiv_data[tagdist[0][1]]['sampleid2'] for tagdist in tagdists \
                                      if tagdist[0][0] <= mismatch_allowed and tagdist[1][0] > mismatch_allowed]
                    if len(set(indiv_cand)) == 1:
                        indiv = indiv_cand[0]
                        read = [s[idxlen:] for s in ss]
                        qual = [[ord(c)-baseQ_in for c in qstr[idxlen:]] for qstr in qstrs]
                    
                else: #dump both reads per the first
                    t = ss[0][:idxlen] #tag from read1
                    ts = [t]*2 # hack for getting tag into both reads, below
                    tqs = [qstrs[0][:idxlen]]*2
                    tagdist = sorted([(distance(t_this,t),t_this) for t_this in indiv_data.keys()])
                    if tagdist[0][0] <= mismatch_allowed and tagdist[1][0] > mismatch_allowed:
                        indiv = indiv_data[tagdist[0][1]]['sampleid']
                        read = [ss[0][idxlen:],ss[1]]
                        qual = [[ord(c)-baseQ_in for c in qstrs[0][idxlen:]],[ord(c)-baseQ_in for c in qstrs[1]]]

                if indiv is None:
                    read = ss
                    qual = [[ord(c)-baseQ_in for c in qstr] for qstr in qstrs]
                    if passfh is not None:
                        for id,s,q,fh in zip(heads,read,qual,passfh):
                            fh.write(as_fq_line(id,s,q,baseQ_out,output_lnum,))
                else:
                    if indiv_reads_out_pattern is not None:
                        for h,t,tq,s,q,rn,pat in zip(heads,ts,tqs,read,qual,[1,2],indiv_reads_out_pattern):
                            newhead = '%s %s:%s' % (h,t,tq)
                            try:
                                fhdict[(indiv,rn)].write(as_fq_line(newhead,s,q,baseQ_out,output_lnum))
                            except KeyError:
                                fhdict[(indiv,rn)] = smartopen(pat % indiv,'w')
                                fhdict[(indiv,rn)].write(as_fq_line(newhead,s,q,baseQ_out,output_lnum))
                                
                qual = [numpy.array(q,dtype=int) for q in qual]
                
            else:
                raise ValueError, 'read2_has_idx cannot be None for PE reads'
        else:
            raise ValueError, 'PE handling invoked, but indiv_out_pattern does not match; must be 2-tuple or None, is: %s' % indiv_reads_out_pattern
    else:

        if indiv_reads_out_pattern is not None:
            if fhdict is None:
                fhdict = {}
            if passfh is None:
                passfh = smartopen(indiv_reads_out_pattern % 'pass','w')

        head,s,qstr = line

        if baseQ_in is None:
            if get_baseQ(qstr) is None:
                raise ValueError,'could not determine qual base (33 or 64): %s' % qstr
            else:
                baseQ_in = get_baseQ(qstr)
        if baseQ_out is None:
            baseQ_out = baseQ_in


        t = s[:idxlen]
    
        tagdist = sorted([(distance(t_this,t),t_this) for t_this in indiv_data.keys()])
        if tagdist[0][0] <= mismatch_allowed and tagdist[1][0] > mismatch_allowed:
            indiv = indiv_data[tagdist[0][1]]['sampleid']
            read = s[idxlen:]
            qual = [ord(c)-baseQ_in for c in qstr[idxlen:]]
            if indiv_reads_out_pattern is not None:
                newhead = '%s:%s:%s' % (head,t,qstr[:idxlen])
                try:
                    fhdict[indiv].write(as_fq_line(newhead,read,qual,baseQ_out,output_lnum))
                except KeyError:
                    fhdict[indiv] = smartopen(indiv_reads_out_pattern % indiv,'w')
                    fhdict[indiv].write(as_fq_line(newhead,read,qual,baseQ_out,output_lnum))
        else:
            indiv = None
            read = s
            qual = [ord(c)-baseQ_in for c in qstr]
            if passfh is not None:
                passfh.write(as_fq_line(head,s,qual,baseQ_out,output_lnum))

        qual = numpy.array(qual,dtype=int)
    return indiv,read,qual

def get_fastq_properties(fq):
    fh = smartopen(fq)
    chr1 = fh.readline()[0]
    if chr1 == '@':
        lnum = 4
        print >> sys.stderr, 'using 4-line fastq'
    else:
        lnum = 1
        print >> sys.stderr, 'using 1-line fastq'

    baseQ = None
    qfh = smartopen(fq)
    while baseQ is None:
        baseQ = get_baseQ(next_read_from_fh(qfh,lnum)[2])
    qfh.close()
    print >> sys.stderr, 'using quality encoding base-%s' % baseQ
    return lnum,baseQ

def next_read_from_fh(fh,lnum=None):

    if lnum is None:
        if smartopen(fh.name).read(1) == '@':
	    lnum = 4
	else:
	    lnum = 1

    if lnum == 1:
        return fh.readline().strip().rsplit(':',2)
    elif lnum == 4:
        rl = [fh.readline().strip() for i in range(lnum)]
        return [ rl[0][1:], rl[1], rl[3] ]

def as_fq4_lines(id,s,q,baseQ=None):
    if baseQ is None:
        return '\n'.join(['@'+id] + [s,'+',q+'\n'])
    else:
        return '\n'.join(['@'+id] + [s,'+',(''.join([chr(n+baseQ) for n in q]))+'\n'])

def as_fq1_line(id,s,q,baseQ):
    return ':'.join([id] + [s,(''.join([chr(n+baseQ) for n in q]))+'\n'])

def as_fq_line(id,s,q,baseQ,lnum):
    if lnum == 4:
        return as_fq4_lines(id,s,q,baseQ)
    elif lnum == 1:
        return as_fq1_line(id,s,q,baseQ)

def store_SR(all_quality,indiv,read,qual):
    store_read1(all_quality,indiv,read,qual)

def store_PE(all_quality,indiv,read,qual):
    store_read1(all_quality,indiv,read[0],qual[0])
    store_read2(all_quality,indiv,read[0],read[1])

def store_read1(all_quality,indiv,read,qual):
    try:
        all_quality[read]['count'][indiv] += 1
        all_quality[read]['sum_quality'] += qual
    except KeyError:
        all_quality[read]['count'] = defaultdict(int)
        all_quality[read]['count'][indiv] += 1
        all_quality[read]['sum_quality'] = qual 


def store_read2(all_quality,indiv,read1,read2):
    try:
        all_quality[read1]['read2'][read2] += 1
    except KeyError:
        all_quality[read1]['read2'] = defaultdict(int)
        all_quality[read1]['read2'][read2] += 1


def write_uniqued(all_quality,outfile,baseQ):
    ofh = smartopen(outfile,'w')

    for seq in all_quality.keys():
        aqd = all_quality[seq]
        ind,indcount = dezip(sorted([(k,v) for k,v in aqd['count'].items()],reverse=True,key = lambda x:x[1]))
        c = sum(indcount)
        q = ''.join([chr(i+baseQ) for i in aqd['sum_quality'] / c])
        if aqd.has_key('read2') and any([v > 1 for v in aqd['read2'].values()]):
            r2,r2count = dezip(sorted([(k,v) for k,v in aqd['read2'].items() if v > 1],reverse=True,key = lambda x:x[1]))
        else:
            r2,r2count = '.','.'
        
        line =  '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (seq,c,q, ','.join(ind), ','.join([str(i) for i in indcount]), ','.join(r2), ','.join([str(i) for i in r2count]))
        ofh.write(line)

    ofh.close()


if __name__ == '__main__':

    import argparse

    ds =  ' [%(default)s]'
    #create command line parser
    parser = argparse.ArgumentParser(description='performs index parsing and unique sequence tabulation')

    parser.add_argument('-w','--write_reads_by_indiv',action='store_true',help='enables creation of new .fastq files, one for each individual in reads_by_individual folder'+ds)
    parser.add_argument('-u','--no_uniqued',action='store_true',help='disables creation of .uniqued file'+ds)
    parser.add_argument('-r2idx','--read2_has_idx',action='store_true',help='if specified, individual index (barcode) is also present in read 2'+ds)
    parser.add_argument('-mc','--set_mincycles',default=0,type=int,help='truncate reads to this length (if not 0)'+ds)
    parser.add_argument('-s','--cutsite',default='AATTC',help='sequence left behind by restriction enzyme at read1 end of library NOT NECESSARILY FULL R.E. SITE'+ds)
    parser.add_argument('-iq','--base_Q_in',default=None,type=int,help='integer offset for quality scores IN INPUT.  Usually 33 for "sanger" style (newer illumina runs, input for BWA) or 64 for illumina/solexa (older illumina). If None, ascertain from data'+ds)
    parser.add_argument('-oq','--base_Q_out',default=33,type=int,help='integer offset for quality scores IN OUTPUT.  Usually 33 for "sanger" style (newer illumina runs, input for BWA) or 64 for illumina/solexa (older illumina). value None will output according to input'+ds)
    parser.add_argument('-ol','--output_lnum',default='4',choices=['1','4'],type=int,help='number of lines per record in fastq output if -w is specified (either older 1-line, or newer 4-line)'+ds)
    parser.add_argument('-fc','--flowcell',default=None,type=str,help='flowcell name (if None, derive from sequence infile path)'+ds)
    parser.add_argument('-l','--lane',default=None,type=str,help='lane (if None, derive from sequence infile name)'+ds)
    parser.add_argument('-idx','--index',default=None,type=str,help='multiplex index (if None, no index applied)'+ds)

    parser.add_argument('-suf','--suffix',default=None,type=str,help='suffix for .uniqued file (permits processing split files)'+ds)

    parser.add_argument('--force_db_id',action='store_true',help='force mouse database ids for individuals \n(replacing legacy sampleid from DB_library_data)'+ds)

    parser.add_argument('-e','--estimate_error',action='store_true',help='invokes clustering to estimate error rate after completion of preprocessing'+ds)
    parser.add_argument('-ee','--est_err_engine',default='local',choices=['local','parallel','lsf'],help='use this engine for parallelizing error estimate (rtd_run -pe argument)'+ds)
    parser.add_argument('-ec','--est_err_cores',default=1,type=int,help='parallelize error estimate run over this number of cores (serial if less than 2) REQUIRES GNU PARALLEL'+ds)
    parser.add_argument('-ep','--est_err_parts',default=4,type=int,help='number of query files to split error estimate simliarity calculation into (see rtd_run -np argument)'+ds)
    parser.add_argument('-et','--est_err_threads',default=4,type=int,help='number of MCL expansion threads in clustering (see rtd_run -te argument)'+ds)
    parser.add_argument('-er','--est_err_radius',default=2,type=int,help='MCL radius argument (-I) for error estimate clustering'+ds)

    parser.add_argument('infiles',nargs='+',help='either 1 or 2 fastq files corresponding to reads from a single lane, and optionally read 2 sequences for that lane')

    opts = parser.parse_args()

    write_reads_by_indiv = opts.write_reads_by_indiv

    set_mincycles = opts.set_mincycles
    read2_has_idx = opts.read2_has_idx
    nticks = 20
    cutsite = opts.cutsite
    baseQ_in = opts.base_Q_in
    baseQ_out = opts.base_Q_out

    if len(opts.infiles) == 1:
        r1 = opts.infiles[0]
        fh = smartopen(r1)
        chr1 = fh.read(1)
        fh.seek(0)
    elif len(opts.infiles) == 2:
        r1,r2 = opts.infiles[0:2]
        fh = (smartopen(r1),smartopen(r2))
        chr1 = fh[0].read(1)
        fh[0].seek(0)
    else:
        raise ValueError, 'either one or two input fastq files must be specified, got %s' % opts.infiles
    # here forward, if isinstance(fh,tuple) then we're PE

    if chr1 == '@':
        lnum = 4
        print >> sys.stderr, 'using 4-line fastq'
    else:
        lnum = 1
        print >> sys.stderr, 'using 1-line fastq'

    qfh = smartopen(r1)
    while baseQ_in is None:
        baseQ_in = get_baseQ(next_read_from_fh(qfh,lnum)[2])
    qfh.close()
    
    if baseQ_out is None:
        baseQ_out = baseQ_in
	
    print >> sys.stderr, 'read qualities base %s, write qualities base %s' % (baseQ_in, baseQ_out)

    idxlookup = get_adapter_index_lookup()
    indivs = []

    if opts.flowcell is None:        
        fc = os.path.basename(os.path.dirname(r1))
    else:
        fc = opts.flowcell

    if opts.lane is None:
        lane = os.path.basename(r1)[2]
    else:
        lane = opts.lane

    nreads = get_read_count(r1,lnum)

    #index info append
    if opts.index is None or opts.index == 'None':
        index = None
        idxstr = ''
    else:
        index = opts.index
        idxstr = '_index%s' % index

    if opts.suffix is not None:
        idxstr = idxstr+'_'+opts.suffix
    
    if isinstance(fh,tuple):
        nreads2 = get_read_count(r2,lnum)
        if nreads != nreads2:
            raise ValueError, '%s read count: %s; %s read count: %s' % (r1, nreads, r2, nreads2)

    indiv_data = get_individual_data_for_lane(idxlookup=idxlookup,fc=fc,lane=lane,index=index,get_transtable_from_mouseDB=opts.force_db_id)

    adaptersversions = set([r['adaptersversion'] for r in indiv_data.values()])
    idxs = reduce(lambda x,y: x+y, [idxlookup[adver].values() for adver in adaptersversions])    
    idxlen = len(indiv_data.keys()[0])
    line = next_read_from_fh(smartopen(r1),lnum)
    readlen = len(line[-2]) - idxlen
    print >> sys.stderr, '%s\n\t%s bp reads, %s / %s %s bp idxs used' % (r1,readlen,len(indiv_data),len(idxs),idxlen)


    if write_reads_by_indiv:
	outroot = os.path.dirname(r1)
        indiv_reads_out_base = os.path.join(outroot,'reads_by_individual/%s_lane%s%s/' % (fc,lane,idxstr))
        try:
            os.makedirs(indiv_reads_out_base)
	    os.system('chmod g+w '+indiv_reads_out_base.rstrip('/'))
	    print >> sys.stderr, indiv_reads_out_base,'created for individual output'
        except OSError:
            print >> sys.stderr, indiv_reads_out_base,'exists, using'
            
    else:
        indiv_reads_out_base = None

    #prep filehandles for faster processing, see assign_reads_to_indiv docstring
    if indiv_reads_out_base is not None:
        fhdict = {}	    
        if isinstance(fh,tuple):
            indiv_reads_out_pattern = tuple([os.path.join(indiv_reads_out_base,'%s_'+os.path.basename(r1)),os.path.join(indiv_reads_out_base,'%s_'+os.path.basename(r2))])
	    passfh = tuple([smartopen(indiv_reads_out_pattern[0] % 'pass', 'w'),smartopen(indiv_reads_out_pattern[1] % 'pass', 'w')])
	else:
            indiv_reads_out_pattern = os.path.join(indiv_reads_out_base,'%s_'+os.path.basename(r1))
            passfh = smartopen(indiv_reads_out_pattern % 'pass', 'w')
    else:
        indiv_reads_out_pattern = None
        fhdict = None
        passfh = None

    all_quality = defaultdict(dict)
    tickon = nreads/nticks
    if tickon < 1:
    	tickon = 1
    print >> sys.stderr, '\tloading'

    if isinstance(fh,tuple): #PE
        outfile = os.path.join(os.path.dirname(r1), '%s_lane%s_PE%s.uniqued.gz' % (fc,lane,idxstr))
        for i in xrange(nreads):
            if i%tickon==0: print >> sys.stderr, '\t\t%s / %s' % (i,nreads)

            l = tuple([next_read_from_fh(h,lnum) for h in fh])
            indiv,read,qual = assign_read_to_indiv(l,indiv_data,indiv_reads_out_pattern=indiv_reads_out_pattern,\
            										fhdict=fhdict,passfh=passfh,read2_has_idx=read2_has_idx,\
            										baseQ_in=baseQ_in,baseQ_out=baseQ_out,lnum=lnum,output_lnum=opts.output_lnum)
            if indiv is not None and not opts.no_uniqued:
                #store PE
                store_PE(all_quality,indiv,read,qual)
    else: #SR
        outfile = os.path.join(os.path.dirname(r1), '%s_lane%s_SR%s.uniqued.gz' % (fc,lane,idxstr))
        for i in xrange(nreads):
            if i%tickon==0: print >> sys.stderr, '\t\t%s / %s' % (i,nreads)

            l = next_read_from_fh(fh,lnum)
            
            indiv,read,qual = assign_read_to_indiv(l,indiv_data,indiv_reads_out_pattern=indiv_reads_out_pattern,\
            										fhdict=fhdict,passfh=passfh,baseQ_in=baseQ_in,baseQ_out=baseQ_out,lnum=lnum,output_lnum=opts.output_lnum)
            if indiv is not None and not opts.no_uniqued:
                #store SR
                store_SR(all_quality,indiv,read,qual)
                

    if fhdict is not None:
        for v in fhdict.values():
            v.close()
    if passfh is not None:
        if isinstance(passfh,tuple):
            for pfh in passfh:
	        pfh.close()
	else:
            passfh.close()


    if opts.no_uniqued:
        print >> sys.stderr, '.uniqued output disabled, skip postprocessing on .uniqued'
    else:
        write_uniqued(all_quality,outfile,baseQ_out)
        print >> sys.stderr, 'output written to',outfile

        # summarize stats buggy; pool_lookup throws errors and call to Util isn't portable
        # disable until fixed
        #print >> sys.stderr, 'generate preprocess summary (summarize_sequencing_stats.py)'
        #ret = os.system(os.path.join(RTDROOT,'summarize_sequencing_stats.py %s > %s.stats' % (outfile,outfile)))

        if opts.estimate_error:
            os.system(os.path.join(RTDROOT,'estimate_error_by_clustering.py %s %s %s %s %s %s %s' % (outfile, opts.cutsite, opts.est_err_engine, opts.est_err_cores, opts.est_err_parts, opts.est_err_threads, opts.est_err_radius)))
    
    #done
