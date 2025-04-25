import tkinter as tk
from tkinter import ttk
import json
from processPrediction import ProcessPrediction
from database_manager import DatabaseManager
from group_manager import GroupManager

class PredictionTab:
    def __init__(self, parent, db: DatabaseManager, group_manager: GroupManager):
        self.db = db
        self.frame = ttk.Frame(parent)
        self.groups = group_manager

        self.predictor = ProcessPrediction(self.db)

        self._build_interface()

    def _build_interface(self):
        label = ttk.Label(self.frame, text="Režim predikce:")
        label.grid(row=0, column=0, sticky="w")

        self.mode_var = tk.StringVar(value="Vlastní sekvence")
        self.mode_menu = ttk.Combobox(self.frame, textvariable=self.mode_var, values=["Vlastní sekvence", "Skupina enzymů", "Soubor"], state="readonly")
        self.mode_menu.grid(row=0, column=1, sticky="ew")
        self.mode_menu.bind("<<ComboboxSelected>>", self._update_input_visibility)

        self.seq_input = tk.Text(self.frame, height=4, width=60)
        self.seq_input.grid(row=1, column=0, columnspan=2, sticky="ew")

        group_names = self.groups.get_all_group_names()

        self.group_var = tk.StringVar()
        self.group_menu = ttk.Combobox(self.frame, textvariable=self.group_var, values=group_names, state="readonly")
        self.group_menu.grid(row=2, column=0, columnspan=2, sticky="ew")

        files = self.db.get_unique_filenames()
        self.file_var = tk.StringVar()
        self.file_menu = ttk.Combobox(self.frame, textvariable=self.file_var, values=files, state="readonly")
        self.file_menu.grid(row=3, column=0, columnspan=2, sticky="ew")

        self.predict_button = ttk.Button(self.frame, text="Predikuj", command=self._run_prediction)
        self.predict_button.grid(row=4, column=0, columnspan=2, pady=5)

        self.output_label = ttk.Label(self.frame, text="Zatím nebyla provedena žádná predikce.", wraplength=600, justify="left")
        self.output_label.grid(row=5, column=0, columnspan=2, sticky="w")

        self._update_input_visibility()

    def _update_input_visibility(self, *_):
        mode = self.mode_var.get()
        self.seq_input.grid_remove()
        self.group_menu.grid_remove()
        self.file_menu.grid_remove()

        if mode == "Vlastní sekvence":
            self.seq_input.grid()
        elif mode == "Skupina enzymů":
            self.group_menu.grid()
        elif mode == "Soubor":
            self.file_menu.grid()

    def _run_prediction(self):
        self.predictor.load_model()
        mode = self.mode_var.get()

        if mode == "Vlastní sekvence":
            seq = self.seq_input.get("1.0", "end").strip()
            if seq:
                result = self.predictor.predict(seq)
                self.output_label.config(text=f"Predikovaný motiv: {result}")
            else:
                self.output_label.config(text="Zadejte sekvenci.")

        elif mode == "Skupina enzymů":
            try:
                group = self.group_var.get()
                group_items = self.groups.get_group_ids(group)
                results = [f"{list(item)[1]}: {self.predictor.predict(list(item)[9])}" for item in group_items if isinstance(item, tuple) and len(item) > 9]
                self.output_label.config(text="\n".join(results))
            except:
                self.output_label.config(text="Chyba při načítání skupiny.")

        elif mode == "Soubor":
            fname = self.file_var.get()
            if fname:
                all_seqs = self.db.get_samples_by_filename(fname)
                results = []
                for row in all_seqs:
                    seq_text = row[1]
                    result = self.predictor.predict(seq_text)
                    results.append(f"{row[1]}: {result}")
                self.output_label.config(text="\n".join(results))
