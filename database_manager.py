from Database import Database

class DatabaseManager:
    def __init__(self):
        self.db = Database()

    def create_database(self):
        self.db.create_database()

    def insert_data(self, data):
        self.db.insert_data(data)

    def get_all_samples(self):
        return self.db.select_all_samples_data()

    def get_samples_by_filename(self, filename):
        return self.db.select_samples_data(filename=filename)

    def get_unique_filenames(self):
        return [f[0] for f in self.db.select_file_names()]

    def delete_by_filename(self, filename):
        self.db.erase_data(filename)
        self.get_all_samples()

    def reset_database(self):
        self.db.reset_database()