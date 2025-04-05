class FilterManager:
    def __init__(self, db_manager):
        self.db = db_manager

    def filter_data(self, enzyme=None, only_fragment=False, filename=None):
        if filename:
            data = self.db.get_samples_by_filename(filename)
        else:
            data = self.db.get_all_samples()

        filtered = []
        for row in data:
            match = True
            if enzyme and enzyme.lower() not in row[1].lower():
                match = False
            if only_fragment and not row[6]:
                match = False
            if match:
                filtered.append(row)

        return filtered