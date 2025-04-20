from Database import Database

class DatabaseManager:
    def __init__(self):
        self.db = Database()

    # ------------------- SAMPLES -------------------

    def create_database(self):
        self.db.create_database()

    def reset_database(self):
        self.db.reset_database()

    def insert_data(self, data):
        self.db.insert_data(data)

    def get_all_samples(self):
        return self.db.select_all_samples_data()

    def get_samples_by_filename(self, filename):
        return self.db.select_samples_data(filename=filename)

    def get_samples_by_enzyme(self, enzyme):
        return self.db.select_samples_data(enzyme=enzyme)

    def get_samples_by_id(self, enzyme_id):
        return self.db.select_samples_data(id=enzyme_id)

    def get_unique_filenames(self):
        return [f[0] for f in self.db.select_file_names()]

    def delete_by_filename(self, filename):
        self.db.erase_data(filename)

    def get_size_range(self):
        return self.db.get_size_range()

    def get_distinct_enzymes(self):
        return self.db.get_distinct_enzymes()

    # ------------------- ALIGNMENTS -------------------

    def insert_alignment(self, enzyme_a_id, enzyme_b_id, algorithm, score):
        self.db.insert_alignment(enzyme_a_id, enzyme_b_id, algorithm, score)

    def insert_alignments_bulk(self, alignments_data):
        self.db.insert_alignments_bulk(alignments_data)

    def get_scores_for_pairs(self, id_pairs, algorithm):
        return self.db.get_alignment_scores(id_pairs, algorithm)


    def get_alignment_data(self, id_pairs, algorithm):
        return self.db.get_alignment_data(id_pairs, algorithm)