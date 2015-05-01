import copy
import urllib
import re,os,sys
import numpy
import Util


from UserDict import UserDict
from UserList import UserList

class Region (UserDict):
    '''Store and operate on GFF regions (single line)

    setting self["__line"] populates object
    setting self["__attrib_str"] populates attributes

    getting __line, __attrib_str rebuilds strings from current values'''

    FIELDTYPES = { 'seqid':str,
                   'source':str,
                   'type':str,
                   'start':int,
                   'end':int,
                   'strand':str,
                   'score':(float,int,str),
                   'frame':(int,str)
                   }
    
    def __init__(self, line=None, version=3, fields=None):
        UserDict.__init__(self)
        self["__version"] = version
        self["strand"] = '.'
        self["score"] = '.'
        self["source"] = '.'
        self["type"] = '.'
        self["phase"] = '.'
        self["attributes"] = {}
        self["__line"] = line
        if fields: self.update(fields)
#        if line: self.__parse(line)


    def parse_gff3_attrib(self, attrib_str):
        "parses gff3 formatted attribute strings.  DOES NOT CLEAR PREVIOUS KEYS!"
#        self["attributes"]={}
        try:
            for (k,v) in [(kv.strip().split('=')) for kv in attrib_str.strip(' ;').split(';')]: 
                (k,v) = (urllib.unquote(k),urllib.unquote(v))
                self["attributes"][k]=v
        except ValueError:
            print >> sys.stderr, 'failed to parse',attrib_str

    def compose_gff3_attrib(self):
        "compose gff3 formatted attribute string. Don't call without checking for attributes!"
        self.data["__attrib_str"] = ';'.join(["%s=%s" % (urllib.quote(k),urllib.quote(v)) for (k,v) in self["attributes"].iteritems()])
        
    def parse_gff2_attrib(self, attrib_str):
        "parses gff2 formatted attribute strings.  DOES NOT CLEAR PREVIOUS KEYS!"
#        self["attributes"]={}
        for (k,v) in [(kv.split(' ')) for kv in attrib_str.split(';')]: 
            (k,v) = (urllib.unquote(k),urllib.unquote(v.strip('"')))
            self["attributes"][k]=v

    def compose_gff2_attrib(self, attrib_str):
        "compose gff2 formatted attribute string. Don't call without checking for attributes!"
        self.data["__attrib_str"] = ';'.join(['%s "%s"' % (urllib.quote(k),urllib.quote(v)) for (k,v) in self["attributes"].iteritems()])


    def __parse(self, line):
        "parse a gff formatted line, populate self. This should be re-written with regex"
#        print "parsing: %s" % line
	if line.startswith('#'):
		return None
        fields = line.split("\t")
        if len(fields) == 8:
            (self["seqid"], self["source"], self["type"], self["start"], self["end"],
             self["score"], self["strand"], self["phase"] ) = fields[0:8]
        elif len(fields) == 9:
            (self["seqid"], self["source"], self["type"], self["start"], self["end"],
             self["score"], self["strand"], self["phase"], self["__attrib_str"] ) = fields[0:9]
        else:
            print "line:\n%s\n\thas %s fields! (gff has 8 or 9)" % (line, len(fields))

    def __compose_line(self):
        "compose a gff line in current version format"
        try:
            self.data["__line"] = "\t".join([self["seqid"], self["source"], self["type"], 
                                             str(self["start"]), str(self["end"]), str(self["score"]), 
                                             self["strand"], self["phase"] ])
            if "attributes" in self.keys():
                self.data["__line"] += "\t%s" % self["__attrib_str"]
        except:
            self.data["__line"] = None

    def __setitem__(self, key, item):
        if key == "__line" and item:
            item = item.strip("\n")
            self.__parse(item)
        if key == "__attrib_str" and item:
            getattr(self, "parse_gff%s_attrib" % self["__version"])(item)
        if (key == "start" or key == "end") and item:
            item = int(item)
        if key == "score" and item:
            item = item == '.' and item or float(item)
        if re.search(r'^attribute_',key):
            self.data['attributes'][key[10:]] = item
            del key
            del item
        try:
            UserDict.__setitem__(self, key, item)
        except NameError:
            pass


    def __getitem__(self, key):
        if isinstance(key,str):
            if key == "__line":
                self.__compose_line()
            if key == "__attrib_str":
                getattr(self, "compose_gff%s_attrib" % self["__version"])()
            if re.search(r'^attribute_',key):
                return self.data['attributes'][key[10:]]
        return UserDict.__getitem__(self,key)

    def update(self,d):
        '''replace update so special assignments work'''
        for k,v in d.items():
            self[k]=v

    def __repr__(self):
        if self["__line"] is not None:
            return self["__line"]
        else:
            return self.full_repr()

    def full_repr(self):
        return UserDict.__repr__(self)

    def __len__(self):
        try:
            return int(self["end"]) - int(self["start"]) + 1
        except:
            return None
	
	
    def as_ncbi_table(self,attribs):
		'''returns a string suitable for NCBI tbl2asn
		
		attribs is a dict where keys are allowed NCBI 'qualifiers' and values are the corresponding GFF attribute field
		'''
		if self['strand'] == '+':
			start = self['start']
			end = self['end']
		else:
			start = self['end']
			end = self['start']
		tbl = '%d\t%d\t%s\n' % (start,end,self['type'])
		for k,v in attribs.items():
			if self['attributes'].get(v):
				tbl += '\t'*3
				tbl += '%s\t%s\n' % (k,self['attributes'][v])
				
		return tbl

    def rc(self, seqlen):
        '''given the length of the sequence 'seqid', 

        return a region copy on the opposite strand at reversed coordinates'''

        rc_region = copy.deepcopy(self)

        if self["strand"] == '+':
            rc_region["strand"] = '-'
        elif self["strand"] == '-':
            rc_region["strand"] = '+'
        
        oldstart = self["start"]
        oldend = self["end"]
        rc_region["start"] = seqlen - oldend + 1
        rc_region["end"] = seqlen - oldstart + 1
        
        return rc_region

    def coords(self):
        '''returns seqid,start,end,strand tuple'''
        return (self['seqid'],self['start'],self['end'],self['strand'])

    def coords_str(self,delim='_'):
        return delim.join([str(i) for i in self.coords()])

    def contains(self,region):
        '''given another region, returns True if that region is entirely within self, False otherwise'''

        if self["seqid"] == region["seqid"] \
                and self["start"] <= region["start"] \
                and self["end"] >= region["end"]:
            return True
        else:
            return False

    def overlaps(self,region):
        '''given another region, return the number of bases overlapped if that region overlaps self
        
        if non-overlapping, return False'''
        
        overlap = 0

        if self["seqid"] == region["seqid"]:
            if self.contains(region):
                #region contained entirely in self
                overlap = len(region)
            elif region.contains(self):
                #self contained entirely within region
                overlap = len(self)
            elif self["start"] <= region["start"] and self["end"] >= region["start"]:
                #region starts (but does not end) in self
                overlap = (self["end"] - region["start"]) + 1 # +1 since end-start=0; overlap of 1 base
            elif region["start"] <= self["start"] and region["end"] >= self["start"]:
                #self starts (but does not end) in region
                overlap = (region["end"] - self["start"]) + 1

        return overlap and overlap or False

    def to_relative(self,new_origin):
        region = copy.deepcopy(self)
        (region['start'],region['end']) = (region['start'] - new_origin,region['end'] - new_origin)
        return region

    def to_absolute(self,rel_origin):
        region = copy.deepcopy(self)
        (region['start'],region['end']) = (region['start'] + rel_origin,region['end'] + rel_origin)
        return region


class File (UserList):
    '''
    contains GFF.Region objects; reads and writes GFF files

    init argument is filename; if specified loads from gff (otherise call load)

    printing or assigning from (e.g. __repr__) produces tabular GFF
    '''

    def __init__(self, filename=None, version=3):
        UserList.__init__(self)
        if filename: self.load(filename, version)

    def load (self, filename, version):
        "load regions from GFF file"
        import gzip
        if filename.endswith('.gz'):
            opener = gzip.open
        else:
            opener = open

        try:
            fsock = opener(filename, "r")
            try:
                for line in fsock.readlines():
                    if line[0] == "#":
                        pass
                    else:
                        self.append(Region(line,version=version))
            finally:
                self.filename = filename
                self.version = version
                fsock.close()
        except IOError:
            pass
    
    def __repr__(self):
        return "".join([str(region)[-1] == "\n" and str(region) or str(region)+"\n" for region in self])
        
    def rc_regions(self,id_len_dict):
        '''revcomps regions with seqid matching keys in the id_len_dict arg 
        (reguires lengths as keys, not currently validated)'''
        for i in range(len(self)):
            if self[i]["seqid"] in id_len_dict.keys():
                self[i] = self[i].rc(id_len_dict[self[i]["seqid"]])

    def write_to_file(self, outfile=None, set_filename=0):
        '''output current fasta to specified file; if no file is specified, self.filename is used 

        if set_filename is specified, changes the filename of this instance'''
        
        try: 
            fsock = open(outfile or self.filename, "w")
            fsock.write(self.__repr__())
            fsock.close()
            if set_filename:
                self.filename = outfile
        except IOError:
            pass

    def as_ncbi_table(self,attribs):
		'''uses the as_ncbi_table method of Region to return a tbl2asn-compatible table
		'''
		
		seqids = set([r['seqid'] for r in self])
		
		tbl = ''
		
		for s in seqids:
			tbl += '>Feature %s\n' % s
			tbl += ''.join([r.as_ncbi_table(attribs) for r in self if r['seqid'] == s])
			
		return tbl

    def sort(self):
        '''replaces the list.sort() method to specify seqid,start,end as the sort tuple'''
        UserList.sort(self,key=lambda x: x.coords())

    def regions_containing(self,region):
        '''takes a region object

        returns a list of regions that satisfy region.contains()'''
        subset = File()
        subset.extend([thisreg for thisreg in self if region.contains(thisreg)])
        return subset

    def regions_contained_by(self,region):
        print >> sys.stderr, 'deprecated, use regions_containing'
        return self.regions_containing(region)

    def regions_overlapping(self,region):
	    '''takes a region object
	
	    returns regions overlapping (e.g. region.overlaps() == True)'''
	    subset = File()
	    subset.extend([thisreg for thisreg in self if region.overlaps(thisreg)])
	    return subset
	
    def regions_overlapped_by(self,region):
        print >> sys.stderr, 'deprecated, use regions_containing'
        return self.regions_overlapping(region)

    def collapse_overlapping_regions(self):
        '''at present, simply removes regions 'under' the leftmost of a set of overlapping regions

        this is clearly not the perfect behavior--something more like a 'merge attributes' process would be preferred'''
        self.sort()
        i = 0
        while i < len(self)-1:
            if self[i].overlaps(self[i+1]):
                self[i]['end'] = max(self[i]['end'],self[i+1]['end'])
                del self[i+1]
            else:
                i += 1

    def unannotated_regions(self,fields={'type':'intergenic_region'}):
        '''bit of a specialty method--returns _un_annotated regions
        i.e. intervals between regions in self

        known bug: doesn't include trailing (i.e. last annot to end of seq) region
        to add: set-ify attribute keys for flanking annotations and make 'a-b' values for keys present in both'''

        unannotated = File()
        self.collapse_overlapping_regions()

        #hack to seed the upstream unannotated interval
        #last_region = {'seqid':self[0]['seqid'],'end':0}
        last_region = {'seqid' : None}
        for region in self:
            if last_region['seqid'] == region['seqid']:
                new_reg = Region()
                new_reg.update(fields)
                new_reg['start'] = last_region['end']+1
                new_reg['end'] = region['start']-1
                new_reg['seqid'] = region['seqid']
                new_reg['attribute_Name'] = new_reg.coords_str()
                if len(new_reg):
					unannotated.append(new_reg)
            last_region = region
        return unannotated
                
    def to_relative(self,new_origin):
        new_gff = File()
        for r in self:
            new_gff.append(r.to_relative(new_origin))
        return new_gff

    def to_absolute(self,rel_origin):
        new_gff = File()
        for r in self:
            new_gff.append(r.to_absolute(rel_origin))
        return new_gff

    def get_seqids(self):
        return sorted(list(set([r['seqid'] for r in self])))

    def to_coverage_arrays(self,dtype=bool,endpoints=None):
        '''coverts a collection of regions into a dict of coverage counts (if dtype is numeric) or True/False (if dtype is bool)
        
        endpoints is a dict of seqid:length pairs for use in defining array sizes'''

        ca = {}
        
        seqids = self.get_seqids()
        for s in seqids:
            if endpoints is not None and s in endpoints.keys():
                al = endpoints[s]
            else:
                al = max([r['end'] for r in self if r['seqid'] == s]) + 1
            ca[s] = numpy.zeros(al,dtype=dtype)
			
        for r in self:
            ca[r['seqid']][r['start']:r['end']+1] += 1
			
        return ca

def get_multifile_coverage_arrays(gff_list,dtype=bool,endpoints_list=None):
	'''given a list of GFF.File objects, returns a coverage array dict where arrays are y = longest array and x = number of files,
	see GFF.File.to_coverage_arrays for details
	e.g.
	
	{
	'agouti':array([[True,True,False,True],[False,True,False,True]])
	}
	
	where two GFFs contain regions with seqid 'agouti' where the longest sequence (or last end) is 4.  inner arrays ordered per input GFF list order
	'''
	
	#get individual coverage arrays
	indiv_cas = []
	if endpoints_list is not None and len(endpoints_list) == len(gff_list):
		for i,g in enumerate(gff_list):
			indiv_cas.append(g.to_coverage_arrays(dtype=dtype,endpoints_list=endpoints_list[i]))
	else:
		for g in gff_list:
			indiv_cas.append(g.to_coverage_arrays(dtype=dtype))

	#init merged coverage array
	merge_cas = {}.fromkeys(set(Util.flatten_list([a.keys() for a in indiv_cas])),None)
	
	for k in merge_cas.keys():
		slen = max([len(ca.get(k,[])) for ca in indiv_cas])
		merge_cas[k] = numpy.zeros((len(indiv_cas),slen),dtype=dtype)
		
	for i,ca in enumerate(indiv_cas):
		for seqid,ar in ca.items():
			merge_cas[seqid][i][:len(ar)] = ar
		
		
	return merge_cas

def partition_by_mapping_gffs(gff_list,required_sources=[],min_sources=2,min_length=50):
	'''given a list of GFF.File objects, partitions each seqid into segments covered by at least min_sources
	and requiring that required_sources be among them.  Regions must be at least min_length in reference
	
	returns a gff of such regions
	'''
	
	if any([len(set([r['source'] for r in g])) != 1 for g in gff_list]):
		raise ValueError, 'more than one source per GFF! %s' % [set([r['source'] for r in g]) for g in gff_list]
	
	mca = get_multifile_coverage_arrays(gff_list)
	sources = [g[0]['source'] for g in gff_list]
	required_idx = []
	for source in required_sources:
		required_idx.append(sources.index(source))
	
	include_bases = {}
	partition = File()
	
	for seqid, ar in mca.items():
		include_bases[seqid] = []
		for col in xrange(ar.shape[1]):
			sl = ar[:,col]
			#print sl, len(numpy.nonzero(sl)), [sl[i] for i in required_idx]
			if len(numpy.nonzero(sl)[0]) >= min_sources and all([sl[i] for i in required_idx]):
				include_bases[seqid].append(col)
		bounds = Util.get_consecutive_value_boundaries(include_bases[seqid])
		for start,end in bounds:
			if end-start > min_length:
				r = Region(fields={ \
							'seqid':seqid,
							'start':start,
							'end':end,
							})
				partition.append(r)
	
	return partition

from SQLite_DB import DB
class FileDB (DB):
    '''nascent sqlite database format for gff files

    intended use:
    gffdb = GFF.FileDB(gff_filename.gff3)
    all_rows = gffdb.select_as_gff()
    genes_on_2R = gffdb.select_as_gff('seqid = "?" and type = "?"',('2R','gene'))

    attribute_grab = gffdb.select_as_gff('__line like "?"', ('%%Name=%s%%' % (genename,)) )
    in_region = gffdb.select_regions_contained_by_as_gff(region)

    gffdb.cursor.execute('select avg(end-start) from gff')
    mean_len = gffdb.cursor.next()[0] #since matches are field tuples, even with only one query field

    gffdb.cursor.execute('select * from gff limit 100')
    for row in gffdb.cursor:
        print row['start'],row['end'],row['type']

    etc...
    
    see SQLite_DB.py for base class'''
    
    schema = [("id", "INTEGER PRIMARY KEY AUTOINCREMENT"),
              ("seqid","TEXT NOT NULL"),
              ("source", "TEXT NOT NULL"),
              ("type", "TEXT NOT NULL"),
              ("start", "INTEGER NOT NULL"),
              ("end", "INTEGER NOT NULL"),
              ("score", "REAL NOT NULL"),
              ("strand", "TEXT NOT NULL"),
              ("phase", "TEXT NOT NULL"),
              ("__line", "TEXT NOT NULL")
              ]
    
    table = "gff"

    def __init__(self,filename,create=False):
        self.filename = filename
        self.DBfile = self.calcDBfilename(filename)
        DB.__init__(self,self.DBfile,create)
        if create:
            try:
                self.populate_from_file()
            except IOError:
                pass

    def calcDBfilename(self,gff_filename):
        '''calculates the path/.filename.DB database flatfile name'''
        filename = gff_filename.rstrip('.gz')            
        return os.path.join(os.path.dirname(filename),'.'+os.path.basename(filename)+'.DB')

    def insertGFFRegion(self,GFFRegion,commit=1):
        '''makes an insert tuple for current schema from GFF.Region object,
        calls self.insert from ancestor'''

        vals = [GFFRegion[field[0]] for field in self.schema[1:]]
        self.insert(vals,commit)

    def insertGFFlist(self,GFFlist,commit=1):
        for region in GFFlist:
            self.insertGFFRegion(region,0)
        if commit: self.connection.commit()

    def populate_from_file(self):
        '''current implementation opens and reads into DB line-by-line'''

        buffersize = 10000000

        if os.path.exists(self.filename):
            import gzip
            if self.filename.endswith('.gz'):
                opener = gzip.open
            else:
                opener = open

            fh = opener(self.filename)
            lines = fh.readlines(buffersize)
            while lines:
                regions = []
                for line in lines:
                    region = Region(line)
                    if region['__line']:
                        regions.append(region)
                self.insertGFFlist(regions,0)
                lines = fh.readlines(buffersize)
            fh.close()
            self.connection.commit()
        else:
            raise OSError, "File %s doesn't exist" % self.filename
            
        '''too memory intensive--need to go line-by-line
        all_regions = File(self.filename)
        self.insertGFFlist(all_regions)
        '''

    def select_as_gff(self,where=None,select_args=None):
        '''given a SQL where clause, gives back results as a GFF.File
        if select_args is give, assumed SQL '?' value insertion'''
        select = 'SELECT __line FROM %s' % self.table
        if where:
            select = select + ' WHERE %s' % where
        if select_args:
            self.cursor.execute(select, select_args)
        else:
            self.cursor.execute(select)
        hits = self.cursor.fetchall()
        regions = File()
        regions.extend([Region(hit[0]) for hit in hits])
        return regions

    def select_regions_contained_by_as_gff(self,Region,minscore=0,to_relative=False):
        where = "seqid = ? and start > ? and end < ? and score > ? order by start"
        hits = self.select_as_gff(where,Region.coords()[:3]+(minscore,))
        if to_relative:
            hits = hits.to_relative(Region['start'])
        return hits

    def get_upstream_annotation(self,Region,type=None):
        if type is None: #use the type of the specified region
            type = Region['type']

        where = "seqid = ? and type = ? and end < ? order by end desc limit 1"
        try:
            return self.select_as_gff(where,(Region['seqid'],type,Region['start']))[0]
        except IndexError:
            return None

    def get_downstream_annotation(self,Region,type=None):
        if type is None: #use the type of the specified region
            type = Region['type']

        where = "seqid = ? and type = ? and start > ? order by start limit 1"
        try:
            return self.select_as_gff(where,(Region['seqid'],type,Region['end']))[0]
        except IndexError:
            return None

    def get_flanking_annotations(self,Region,type=None):
        if type is None: #use the type of the specified region
            type = Region['type']

        flanks = File()
        flanks.append(self.get_upstream_annotation(Region,type))
        flanks.append(self.get_downstream_annotation(Region,type))
        return flanks

    def get_neighborhood_region(self,feat_string,type=None):
        '''searches for a feature string match (e.g. %Name=eve;%) in __line

        returns a new region object spanning the interval from the previous
        to the next annotations of type 'type' (or the type of the annotation found)
        '''

        (annot,) = self.select_as_gff('__line like "%s"' % feat_string)
        (neighborhood,) = self.get_flanking_annotations(annot,type).unannotated_regions()
        return neighborhood
        

    def close(self):
        self.connection.close()
    
class TabularAnnotation ():
    '''a generic tabular object intended to provide a basis for parsers from tabular formats to GFF

    '''
    def_schema = {'fields' : {0:('seqid',r'.'),1:('start',r'^\d+$'),2:('end',r'^\d+$')} }
    def_delim = None

    def __init__(self,source=None,schema=None,delim=None):
        if schema is not None:
            self.schema = schema
        else:
            self.schema = TabularAnnotation.def_schema
        if delim is not None:
            self.delim = delim
        else:
            self.delim = TabularAnnotation.def_delim
        if isinstance(source,str):
            self.source = open(source)
        elif isinstance (source, (file, list)):
            self.source = source

    def convert(self, source=None):
        '''do the work to convert the list of records in self.schema format into a GFF.File object

        '''
        if source is None:
            if self.source and self.source is not None:
                source = self.source
            else:
                raise TypeError, 'no valid source'

        gff_records = File()
        
        for record in source:
            regions = self.parse(record)
            if regions is not None: gff_records.extend(regions)
        return gff_records

    def populate_valid_fields_dict(self,record,schema_fields_dict):
        fields = {}
        for key,(field,regex) in schema_fields_dict.items():
            try:
                if isinstance(record[key], str) and re.search(regex,record[key].strip()):
                    fields[field] = record[key].strip()
                elif isinstance(record[key],Region.FIELDTYPES[field]):
                    fields[field] = record[key]
            except (IndexError, KeyError):
                pass
        return fields

    def get_subfeature_start_stop_pairs(self,record):
        (startcol,endcol,delim) = (self.schema['subfeatures']['starts'],
                                   self.schema['subfeatures']['ends'],
                                   self.schema['subfeatures']['delim'])
        coords = zip([pos for pos in record[startcol].split(delim) if re.search(r'^\d+$',pos)],
                     [pos for pos in record[endcol].split(delim) if re.search(r'^\d+$',pos)])
        return coords

    def parse(self,record):
        '''basic parse method returns a region built from self.schema

        can be overwritten to perform more complex parse operations,
        but should always return a list of GFF.Region objects'''
        if isinstance(record, str):
            record = record.split(self.delim)
        regions = []
        region = Region()
        region.update(self.populate_valid_fields_dict(record,self.schema['fields']))
        if 'static' in self.schema.keys():
            region.update(self.schema['static'])

        subregions = []
        if 'subfeatures' in self.schema.keys():
            '''will pick out and populate subfeatures from columns specified

            needs at least starts and ends; looks for parent_attr to populate attribute_parent.
            all other keys as fields or attribute_fields'''
                      
            
            fields = self.populate_valid_fields_dict(record,self.schema['subfeatures']['subfeature_fields'])
            if 'subfeature_static' in self.schema['subfeatures'].keys():
                fields.update(self.schema['subfeatures']['subfeature_static'])

            coords = self.get_subfeature_start_stop_pairs(record)

            for (fields['start'],fields['end']) in coords:
                subregion = Region()
                subregion.update(fields)
                subregions.append(subregion)

        if region['__line'] is not None: 
            regions = [region]
        if subregions:
            regions.extend(subregions)
        if regions:
            return regions

'''
tabschemas are schema hashes for use with the TabularAnnotation converter

They can be passed to the object constructor, or written to the member variable self.schema
any time before convert() is called
'''
tabschemas = { 
    'UCSC_dmel_genes+exons' : 
    {'fields' : {1:('attribute_Name',r'.'),
                 2:('seqid',r'.'),
                 3:('strand',r'^[\+-]$'),
                 4:('start',r'^\d+$'),
                 5:('end',r'^\d+$')}, 
     'static' : {'type':'gene',
                 'source':'UCSC'},
     'subfeatures' : {'starts' : 9, 
                      'ends' : 10, 
                      'delim' : ',',
                      'subfeature_fields' : {1:('attribute_Parent',r'.'),
                                             2:('seqid',r'.'),
                                             3:('strand',r'^[\+-]$') }, 
                      'subfeature_static': {'type':'exon','source':'UCSC'} } },

    'UCSC_dmel_genes' : 
    {'fields' : {1:('attribute_Name',r'.'),
                 2:('seqid',r'.'),
                 3:('strand',r'^[\+-]$'),
                 4:('start',r'^\d+$'),
                 5:('end',r'^\d+$')}, 
     'static' : {'type':'gene',
                 'source':'UCSC'}, },
    
    'UCSC_human_genes+exons' :  
    {'fields' : {0:('attribute_Name',r'.'),
                 1:('seqid',r'.'),
                 2:('strand',r'^[\+-]$'),
                 3:('start',r'^\d+$'),
                 4:('end',r'^\d+$') }, 
     'static' : {'type':'gene',
                 'source':'UCSC' },
     'subfeatures' : {'starts' : 8, 
                      'ends' : 9, 
                      'delim' : ',',
                      'subfeature_fields' : {0:('attribute_Parent',r'.'),
                                             1:('seqid',r'.'),
                                             2:('strand',r'^[\+-]$') }, 
                      'subfeature_static': {'type':'exon','source':'UCSC'} } },

    'UCSC_human_genes' :  
    {'fields' : {0:('attribute_Name',r'.'),
                 1:('seqid',r'.'),
                 2:('strand',r'^[\+-]$'),
                 3:('start',r'^\d+$'),
                 4:('end',r'^\d+$') }, 
     'static' : {'type':'gene',
                 'source':'UCSC' }, },

    'UCSC_dmel_pcons' :
    {'fields' : {1:('seqid',r'.'),
                 2:('start',r'^\d+$'),
                 3:('end',r'^\d+$'),
                 4:('attribute_lod',r'.'),
                 5:('attribute_score',r'^\d+$') }, },
    }

                
if __name__ == "__main__":
    
    print "------------------test region object---------------"
    print "building piecemeal:"
    gff_region = Region()
    print gff_region
    gff_region["seqid"] = '2R'
    gff_region["start"] = '50'
    gff_region["end"] = '55'
    print gff_region["__line"]
    gff_region["attributes"]["ID"] = '1'
    print gff_region["__line"]
    del gff_region

    print "assigning __line in an existing instance"
    gff_region = Region()
    gff_line = "FNYB.scaffold_0\tBLAT_homology\ttransferred_annot\t25753\t27044\t0.95\t+\t."
    print gff_line
    gff_region["__line"] = gff_line
    print gff_region
    print "length: %s" % len(gff_region)
    gff_region["attributes"]["species"] = "Bactrocera_dorsalis"
    print gff_region["__line"]
    print "-----\n%s\n-----\n" % gff_region
    gff_line = "FNYB.scaffold_0\tBLAT_homology\ttransferred_annot\t25753\t27044\t0.95\t+\t.\tName=Guelph"
    gff_region["__line"] = gff_line
    print gff_region["__line"]
    print "-----\n%s\n-----\n" % gff_region
    gff_region["__attrib_str"] = "Name=gelfling"
    print "-----\n%s\n-----\n" % gff_region


    print "constructing with __line"
    gff_line = "FNYB.scaffold_0\tBLAT_homology\ttransferred_annot\t25753\t27044\t0.95\t+\t."
    print gff_line
    gff_region = Region(gff_line)
    print "-----\n%s\n-----\n" % gff_region
    gff_line = "FNYB.scaffold_0\tBLAT_homology\ttransferred_annot\t25753\t27044\t0.95\t+\t.\tName=Guelph"
    print gff_line
    gff_region = Region(gff_line)
    print "-----\n%s\n-----\n" % gff_region
    print "-----\n%s\n-----\n" % gff_region.full_repr()
    


    print "------------------test file object---------------"
    test_gff_file = "/home/brant/py_util/unit_test_data/test.gff3"
    gff = File(test_gff_file)
    print gff
    gff.write_to_file("/home/brant/temp/test.gff3")
    del gff
    print "print?"
    try: 
        print gff
    except: 
        print "nope"
    gff = File("/home/brant/temp/test.gff3")
    try:
        print gff
    except:
        print "file write/load test failed!"
