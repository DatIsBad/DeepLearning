import sqlite3


class Database:
    def __init__(self) -> None:
        self.db_path = "DATABASE\R_Enzime.db"
        self.schema_file = "DATABASE\create.sql"
        self.data_path = "DATA"

        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()

    def create_database(self):
        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()

        with open(self.schema_file, "r") as file:
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
        self.cursor.execute(
            "PRAGMA foreign_keys = OFF;"
        )  # Disable foreign key constraints temporarily
        self.cursor.execute("BEGIN TRANSACTION;")
        tables = self.cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table';"
        ).fetchall()

        for table in tables:
            print(f"Dropping table: {table[0]}")
            self.cursor.execute(f"DROP TABLE IF EXISTS {table[0]};")

        self.cursor.execute(
            "PRAGMA foreign_keys = ON;"
        )  # Re-enable foreign key constraints

        # Step 2: Recreate tables using the schema file (optional)
        with open(self.schema_file, "r") as file:
            sql_script = file.read()
            self.cursor.executescript(sql_script)

        # Commit changes and close the connection
        self.connection.commit()
        self.connection.close()
        print("Database reset successfully!")

    def insert_data(self, data: list):
        file_names = self.select_file_names()

        for tuple in file_names:
            for name in tuple:
                if name == data[0][6]:
                    print("File you want to insert into database was already inserted.")
                    return

        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()

        insert_query = """
        INSERT INTO samples (enzyme, sample, orf, rec_sequence, size, fragment, filename, line) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """
        list_of_tuples = []
        for sublist in data:
            list_of_tuples.append(
                (
                    sublist[0],
                    sublist[1],
                    sublist[2],
                    ", ".join(sublist[3]),
                    sublist[4],
                    sublist[5],
                    sublist[6],
                    sublist[7],
                )
            )

        self.cursor.executemany(insert_query, list_of_tuples)
        print("Insert data was sucessful")

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

    def select_all_samples_data(self):
        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()

        self.cursor.execute("SELECT * from samples")  # select all data from samples
        data = self.cursor.fetchall()  # fetch data from executed query

        # Close the connection
        self.connection.close()

        return data

    def select_samples_data(self, select=None, enzyme=None, id=None, filename=None):
        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()
        columns_str = ""
        if select:
            columns_str = ", ".join(select)  # convert list to string

        else:
            columns_str = "*"

        # query base
        query = f"SELECT {columns_str} from samples WHERE 1=1"
        param = []

        # append conditions dynamically
        if enzyme:
            query += " AND enzyme = ?"
            param.append(enzyme)

        if id:
            query += " AND id = ?"
            param.append(id)

        if filename:
            query += " AND filename = ?"
            param.append(filename)

        # execute
        self.cursor.execute(query, param)
        data = self.cursor.fetchall()

        # Close the connection
        self.connection.close()

        return data

    def select_file_names(self):
        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()

        self.cursor.execute(
            "SELECT filename from samples GROUP BY filename"
        )  # select all data from samples
        fileNames = self.cursor.fetchall()  # fetch data from executed query

        # Close the connection
        self.connection.close()

        return fileNames

    def erase_data(self, filename):
        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()

        insert_query = """
        DELETE FROM samples WHERE filename = ?
        """

        self.cursor.execute(insert_query, (filename,))

        print("Insert data was sucessful")

        self.connection.commit()
        self.connection.close()

    def get_size_range(self):
        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()
        self.cursor.execute("SELECT MIN(size), MAX(size) FROM samples")
        result = self.cursor.fetchone()
        self.connection.close()
        return result if result else (0, 2000)

    def get_distinct_enzymes(self):
        self.connection = sqlite3.connect(self.db_path)
        self.cursor = self.connection.cursor()
        self.cursor.execute("SELECT DISTINCT enzyme FROM samples")
        enzymes = [row[0] for row in self.cursor.fetchall()]
        self.connection.close()
        return enzymes
