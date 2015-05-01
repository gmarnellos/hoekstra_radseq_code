'''
collected tools for database-like interactions with google spreadsheets
'''

import os, sys, re, numpy, gzip
import gdata.spreadsheet.service

from collections import defaultdict

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
        try:
            gd_client.source = SOURCE
        except NameError:
            pass

    gd_client.ProgrammaticLogin()

    feed = gd_client.GetSpreadsheetsFeed()
    key = [entry.id.text.rsplit('/', 1)[1] for entry in feed.entry if entry.title.text == target_sheet][0]

    return key,gd_client

def get_table_as_dict(target_sheet,sq=None,gd_client=None,suppress_fc_check=False):

    key,gd_client = get_spreadsheet_key(target_sheet,gd_client)
    if sq is not None:
        q = gdata.spreadsheet.service.ListQuery()
        q.sq = sq
        feed = gd_client.GetListFeed(key,query=q)
    else:
        feed = gd_client.GetListFeed(key)

    recs = []
    for entry in feed.entry:
        #d = []
	#for el in entry.content.text.split(','):
	    
        try:
            recs.append(dict(re.findall('(.+?):\s(.+?)(?:(?:,\s)|$)',entry.content.text)))
            if not suppress_fc_check and not all([k in recs[-1].keys() for k in ['flowcell','lane','pool']]):
                print >> sys.stderr, 'missing keys:', dict(re.findall('(.+?):\s(.+?)(?:(?:,\s)|$)',entry.content.text))
                print >> sys.stderr, 'line was:\n',entry.content.text
        except:
            print >> sys.stderr, 'invalid:', entry.content.text#.split(',')

    return recs

def update_row(target_sheet,pkey_dict,upd_dict,gd_client=None):

    key,gd_client = get_spreadsheet_key(target_sheet,gd_client)
    feed = gd_client.GetListFeed(key)
    recs = [dict([[st.strip() for st in el.split(':')] for el in entry.content.text.split(',')]) for entry in feed.entry]
    hit_idx = []
    for i,rec in enumerate(recs):
        if all([v == rec[k] for k,v in pkey_dict.items() if k in rec]):
            hit_idx.append(i)

    if len(hit_idx) != 1:
        raise ValueError, 'invalid number of records match primary key: %s' % hit_idx

    #otherwise, run the update
    idx = hit_idx[0]
    recs[idx].update(upd_dict)
    el = gd_client.UpdateRow(feed.entry[idx],recs[idx])
