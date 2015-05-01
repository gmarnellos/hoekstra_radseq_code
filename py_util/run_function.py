#!/usr/bin/env python

'''
run_function Aln get_best_match_to_query "{'filename':'some_lastz_output.general'}" outfile
given a module and a function name, plus a string that can be eval'd to either a list or dict to supply as *args or **kwargs

prints whatever is returned to outfile
'''
import sys

mod,fnc,arg,out = sys.argv[1:5]

exec('import %s' % mod)
exec('fnc_obj = %s.%s' % (mod,fnc))

ofh = open(out,'w')

arg_obj = eval(arg)

if isinstance(arg_obj,list):
	ofh.write('%s' % fnc_obj(*arg_obj))
elif isinstance(arg_obj,dict):
	ofh.write('%s' % fnc_obj(**arg_obj))
else:
	raise TypeError, 'argument string %s must evaluate to dict or list (is %s)' % (arg,type(arg_obj))