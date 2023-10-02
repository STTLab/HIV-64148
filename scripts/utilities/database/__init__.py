import sqlite3

class SQLDatabase(object):
    def __init__(self, db_name):
        self.db = self.connect(db_name)

    def connect(self, db_name):
        pass

    def query(self):
        pass
