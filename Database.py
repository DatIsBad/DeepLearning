import sqlite3
import os

class Database:
    def __init__(self) -> None:
        self.db_path = "DATABASE\R_Enzime.db"
        self.schema_file = "DATABASE\create.sql"
        self.data_path = "DATA"

        if os.path.exists(self.db_path):
            print("Database already exists.")
        else:
            self.create_database()

    def create_database(self):
        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()

        with open(self.schema_file, "r") as file:
            sql_script = file.read()

        cursor.executescript(sql_script)
        print("Database schema created successfully!")

        connection.commit()
        connection.close()

    def reset_database(self):
        # Connect to the database
        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()

        # Step 1: Drop all existing tables
        cursor.execute(
            "PRAGMA foreign_keys = OFF;"
        )  # Disable foreign key constraints temporarily
        cursor.execute("BEGIN TRANSACTION;")
        tables = cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name != 'sqlite_sequence';"
        ).fetchall()

        for table in tables:
            print(f"Dropping table: {table[0]}")
            cursor.execute(f"DROP TABLE IF EXISTS {table[0]};")

        cursor.execute(
            "PRAGMA foreign_keys = ON;"
        )  # Re-enable foreign key constraints

        # Step 2: Recreate tables using the schema file (optional)
        with open(self.schema_file, "r") as file:
            sql_script = file.read()
            cursor.executescript(sql_script)

        # Commit changes and close the connection
        connection.commit()
        connection.close()
        print("Database reset successfully!")

    def insert_data(self, data: list):
        file_names = self.select_file_names()

        for tuple in file_names:
            for name in tuple:
                if name == data[0][6]:
                    print("File you want to insert into database was already inserted.")
                    return

        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()

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

        cursor.executemany(insert_query, list_of_tuples)
        print("Insert data was sucessful")

        # Commit the transaction
        connection.commit()

        # Close the connection
        connection.close()

    def dropAll(self):
        # Connect to the database
        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()

        # Fetch all table names
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = cursor.fetchall()

        # Drop each table
        for table in tables:
            table_name = table[0]
            cursor.execute(f"DROP TABLE IF EXISTS {table_name}")
            print(f"Table {table_name} dropped.")

        # Commit and close
        connection.commit()
        connection.close()

        print("All tables dropped successfully!")

    def select_all_samples_data(self):
        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()

        cursor.execute("SELECT * from samples")  # select all data from samples
        data = cursor.fetchall()  # fetch data from executed query

        # Close the connection
        connection.close()

        return data

    def select_samples_data(self, select=None, enzyme=None, id=None, filename=None):
        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()
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
        cursor.execute(query, param)
        data = cursor.fetchall()

        # Close the connection
        connection.close()

        return data

    def select_file_names(self):
        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()

        cursor.execute(
            "SELECT filename from samples GROUP BY filename"
        )  # select all data from samples
        fileNames = cursor.fetchall()  # fetch data from executed query

        # Close the connection
        connection.close()

        return fileNames

    def erase_data(self, filename):
        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()

        insert_query = """
        DELETE FROM samples WHERE filename = ?
        """

        cursor.execute(insert_query, (filename,))

        print("Insert data was sucessful")

        connection.commit()
        connection.close()

    def get_size_range(self):
        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()
        cursor.execute("SELECT MIN(size), MAX(size) FROM samples")
        result = cursor.fetchone()
        connection.close()
        return result if result else (0, 2000)

    def get_distinct_enzymes(self):
        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()
        cursor.execute("SELECT DISTINCT enzyme FROM samples")
        enzymes = [row[0] for row in cursor.fetchall()]
        connection.close()
        return enzymes


    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------

    def insert_alignment(self, enzyme_a_id, enzyme_b_id, algorithm, score, aligned_seq_a, aligned_seq_b):
        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()
        
        insert_query = """
        INSERT INTO alignments (enzyme_a_id, enzyme_b_id, algorithm, score, aligned_seq_a, aligned_seq_b) 
        VALUES (?, ?, ?, ?, ?, ?)
        """

        id_1, id_2 = sorted((enzyme_a_id, enzyme_b_id))

        cursor.execute(insert_query, (id_1, id_2, algorithm, score))

        # Commit the transaction
        connection.commit()

        # Close the connection
        connection.close()


    # ------------------------------------------
    def insert_alignments_bulk(self, alignments_data):
        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()
        
        insert_query = """
        INSERT INTO alignments (enzyme_a_id, enzyme_b_id, algorithm, score, aligned_seq_a, aligned_seq_b) 
        VALUES (?, ?, ?, ?, ?, ?)
        """

        normalized_data = [
            (*sorted((a, b)), algo, score, aligned_seq_a, aligned_seq_b)
            for (a, b, algo, score, aligned_seq_a, aligned_seq_b) in alignments_data
        ]

        cursor.executemany(insert_query, normalized_data)

        # Commit the transaction
        connection.commit()

        # Close the connection
        connection.close()


    # ------------------------------------------
    def select_existing_alignment_combination(self, id_pairs, algorithm):
        # Normalizuj kombinace: (menší_id, větší_id)
        normalized_input = set(tuple(sorted((a, b))) for a, b in id_pairs)

        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()

        query = "SELECT enzyme_a_id, enzyme_b_id FROM alignments WHERE algorithm = ?"

        # Normalizuj i výsledky z databáze
        existing = set(tuple(sorted((row[0], row[1]))) for row in cursor.fetchall())

        cursor.execute(query, algorithm)

        # Close the connection
        connection.close()

        return existing
    

    # ------------------------------------------
    # Metoda která nám vrací hodnotu score daných dvojic daného algoritmu
    def get_alignment_scores(self, id_pairs, algorithm):
        normalized_pairs = [tuple(sorted((a, b))) for a, b in id_pairs]

        if not normalized_pairs:
            return {}

        # Sestavení WHERE podmínky s mnoha OR 
        conditions = " OR ".join(
            "(enzyme_a_id = ? AND enzyme_b_id = ?)" for _ in normalized_pairs
        )

        query = f"""
            SELECT enzyme_a_id, enzyme_b_id, score
            FROM alignments
            WHERE algorithm = ? AND ({conditions})
        """

        # Příprava parametrů
        parameters = [algorithm]
        for pair in normalized_pairs:
            parameters.extend(pair)

        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()

        cursor.execute(query,parameters)
        results = cursor.fetchall()

        # Close the connection
        connection.close()

        # Výstup jako slovník {(id1, id2): score}
        return {
            tuple(sorted((row[0], row[1]))): row[2]
            for row in results
        }
    

    # ------------------------------------------
    # Metoda která nám vrací hodnotu score daných dvojic daného algoritmu
    def get_alignment_data(self, id_pairs, algorithm):
        normalized_pairs = [tuple(sorted((a, b))) for a, b in id_pairs]

        if not normalized_pairs:
            return {}

        # Sestavení WHERE podmínky s mnoha OR 
        conditions = " OR ".join(
            "(enzyme_a_id = ? AND enzyme_b_id = ?)" for _ in normalized_pairs
        )

        query = f"""
            SELECT enzyme_a_id, enzyme_b_id, score, aligned_seq_a, aligned_seq_b
            FROM alignments
            WHERE algorithm = ? AND ({conditions})
        """

        # Příprava parametrů
        parameters = [algorithm]
        for pair in normalized_pairs:
            parameters.extend(pair)

        connection = sqlite3.connect(self.db_path)
        cursor = connection.cursor()

        cursor.execute(query,parameters)
        results = cursor.fetchall()

        # Close the connection
        connection.close()

        # Výstup jako slovník {(id1, id2): {score: ?; aligned_seq_a: ?; aligned_seq_b: ?}}
        return {
            tuple(sorted((row[0], row[1]))): {
                "score": row[2],
                "aligned_seq_a": row[3],
                "aligned_seq_b": row[4]
            }
            for row in results
        }   

