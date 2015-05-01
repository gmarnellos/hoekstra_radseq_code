import os, sqlite3

class DB ():
    '''generic SQLite file object

    '''
    schema = [("id", "INTEGER PRIMARY KEY AUTOINCREMENT")
              
              ]

    table = "table"

    def __init__(self,DBfile,create=False):
        self.DBfile = DBfile
        #print DBfile
        if create:
            try:
                os.makedirs(os.path.dirname(DBfile))
            except OSError:
                pass             
            self.prepare_connection()
            self.reset_table()
        elif os.path.exists(DBfile):
            self.prepare_connection()
        else:
            raise OSError, 'file does not exist--use create=True'
        
        
    def prepare_connection(self):
        self.connection = sqlite3.connect(self.DBfile)
        self.connection.row_factory = sqlite3.Row
        self.cursor = self.connection.cursor()


    def create_table(self):
        '''performs SQL create'''

        create_str = r'CREATE TABLE %s (' % self.table
        create_str = create_str + ', '.join(['%s %s' % (k,v) for k,v in self.schema]) + ')'
        
        self.cursor.execute(create_str)
        self.connection.commit()


    def reset_table(self):
        '''use to clear any existing data and create table'''
        try:
            self.cursor.execute('drop table '+self.table)
            self.connection.commit()
        except:
            pass
        self.create_table()

    def insert(self,vals,commit=1):
        valstr = ',?'*len(vals)
        insert = 'INSERT INTO %s VALUES (null %s)' % (self.table,valstr)
        self.cursor.execute(insert, vals)
        if commit: self.connection.commit()
