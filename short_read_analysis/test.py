#!/usr/bin/env python

import os,sys
import run_safe
from rtd import preprocess_radtag_lane,overlap_preprocess
#import rtd.overlap_preprocess

recs = preprocess_radtag_lane.get_table_as_dict('DB_adapt_trim_seqs',suppress_fc_check=True)

print >> sys.stderr, recs

#def get_adaptseq(table='DB_adapt_trim_seqs'):
#	return dict([(d['adapterstype'],{'r1':d['r1'],'r2':d['r2']}) for d in preprocess_radtag_lane.get_table_as_dict(table,suppress_fc_check=True)])
def get_adaptseq(table='DB_adapt_trim_seqs'):
       return dict([(d['adapterstype'],{'r1':d['r1'],'r2':d['r2']}) for d in recs ])

adaptseq = get_adaptseq()
#print >> sys.stderr, 'use adapterstype: %s\nadaptA: %s\nadaptB: %s' % (adapterstype,adaptA,adaptB)
print >> sys.stderr, adaptseq
