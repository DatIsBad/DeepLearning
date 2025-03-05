import sqlite3

class Database:
    def __init__(self) -> None:
        self.db_path = 'DeepLearning\DATABASE\R_Enzime.db'
        self.schema_file = 'DeepLearning\DATABASE\create.sql'

        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()


    def create_database(self):
        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()

        with open("DeepLearning\DATABASE\create.sql", "r") as file:
            sql_script = file.read()

        self.cursor.executescript(sql_script)
        print("Database schema created successfully!")

        self.connection.commit()
        self.connection.close()


    def reset_database(self):
        # Connect to the database
        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()

        # Step 1: Drop all existing tables
        self.cursor.execute("PRAGMA foreign_keys = OFF;")  # Disable foreign key constraints temporarily
        self.cursor.execute("BEGIN TRANSACTION;")
        tables = self.cursor.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()

        for table in tables:
            print(f"Dropping table: {table[0]}")
            self.cursor.execute(f"DROP TABLE IF EXISTS {table[0]};")

        self.cursor.execute("PRAGMA foreign_keys = ON;")  # Re-enable foreign key constraints

        # Step 2: Recreate tables using the schema file (optional)
        with open(self.schema_file, "r") as file:
            sql_script = file.read()
            self.cursor.executescript(sql_script)

        # Commit changes and close the connection
        self.connection.commit()
        self.connection.close()
        print("Database reset successfully!")


    def insert_data(self, data:list):
        self.connection = sqlite3.connect(self.db_path) 
        self.cursor = self.connection.cursor()


        insert_query = '''
        INSERT INTO samples (enzyme, sample, orf, rec_sequence, size, fragment, filename, line) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        '''
        list_of_tuples = []
        for sublist in data:
            list_of_tuples.append( (sublist[0], sublist[1], sublist[2], ", ".join(sublist[3]), sublist[4], sublist[5], sublist[6], sublist[7]) )
        

        self.cursor.executemany(insert_query, list_of_tuples)

        # Commit the transaction
        self.connection.commit()

        # Close the connection
        self.connection.close()

        
    def dropAll(self):
        # Connect to the database
        self.connection = sqlite3.connect(self.db_path) 
        self.cursor = self.connection.cursor()

        # Fetch all table names
        self.cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = self.cursor.fetchall()

        # Drop each table
        for table in tables:
            table_name = table[0]
            self.cursor.execute(f"DROP TABLE IF EXISTS {table_name}")
            print(f"Table {table_name} dropped.")

        # Commit and close
        self.connection.commit()
        self.connection.close()

        print("All tables dropped successfully!")




