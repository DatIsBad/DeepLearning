import json
from pathlib import Path

class GroupManager:
    def __init__(self, path="groups.json"):
        self.groups = {}
        self.path = Path(path)
        self.load_groups()

    def add_to_group(self, group_name, ids):
        if group_name not in self.groups:
            self.groups[group_name] = set()
        self.groups[group_name].update(ids)

    def get_group_ids(self, group_name):
        return list(self.groups.get(group_name, []))

    def get_all_group_names(self):
        return list(self.groups.keys())

    def save_groups(self):
        with open(self.path, "w", encoding="utf-8") as f:
            json.dump({k: list(v) for k, v in self.groups.items()}, f, indent=4)

    def load_groups(self):
        if self.path.exists():
            with open(self.path, "r", encoding="utf-8") as f:
                loaded = json.load(f)
                self.groups = {k: set(v) for k, v in loaded.items()}