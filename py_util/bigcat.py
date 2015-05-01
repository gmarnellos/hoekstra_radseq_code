#!/usr/bin/env python
'''given a string to glob and  target file, cats all matching files to target

e.g.
bigcat "lastz-masked-chain-general_mus/pairs/*/*.general" lastz-masked-chain-general_mus/out.general

(QUOTES REQUIRED)
'''

import sys
from Util import cat
from glob import glob

def bigcat(globstring,target):
	sourcelist = glob(globstring)
	cat(sourcelist,target)

if __name__ == '__main__':
	globstring,target = sys.argv[1:3]
	bigcat(globstring,target)