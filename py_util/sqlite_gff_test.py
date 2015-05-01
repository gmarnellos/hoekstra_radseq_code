#gff sqlite action

import GFF, os, sqlite3

def InsertGFFRegion(curobj, vals):
    curobj.execute('INSERT INTO gff VALUES (null,?,?,?,?,?,?,?,?,?,?)', vals)
                   

gff_filename = r"G:\AllBrantsStuff\python\ephinaroun\sqlite\dmel-all-r4.3.filtered.gff"
DB_filename = os.path.join(os.path.dirname(gff_filename),'.'+os.path.basename(gff_filename)+'.DB')

gff = GFF.File(gff_filename)
connection = sqlite3.connect(DB_filename)
cursor = connection.cursor()

try:
    cursor.execute('drop table gff')
    connection.commit()
except:
    pass
cursor.execute('''CREATE TABLE gff (
				id INTEGER PRIMARY KEY AUTOINCREMENT,
				sequence_name TEXT NOT NULL,
				source TEXT NOT NULL,
				type TEXT NOT NULL,
				start INTEGER NOT NULL,
				end INTEGER NOT NULL,
				score REAL NOT NULL,
				strand TEXT NOT NULL,
				phase TEXT NOT NULL,
				gff_string TEXT NOT NULL,
				attribute_Name TEXT
			       )''')
connection.commit()

        

for region in gff:
    vals = (region["seqid"],
            region["source"],
            region["type"],
            region["start"],
            region["end"],
            region["score"],
            region["strand"],
            region["phase"],
            region["__line"],
            "Name" in region["attributes"].keys() and region["attributes"]["Name"])
    InsertGFFRegion(cursor,vals)

connection.commit()

cursor.execute('select * from gff')
print cursor.fetchall()
