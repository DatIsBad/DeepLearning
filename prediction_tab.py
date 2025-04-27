import tkinter as tk
from tkinter import ttk
import os
from processPrediction import ProcessPrediction

class PredictionTab:
    def __init__(self, parent, db, group_manager):
        self.db = db
        self.groups = group_manager
        self.frame = ttk.Frame(parent)
        self.model_var = tk.StringVar()
        self.mode_var = tk.StringVar(value="Vlastní sekvence")
        self.predictor = None
        self._build_interface()

    def _build_interface(self):
        top_frame = ttk.LabelFrame(self.frame, text="Výběr modelu a režimu")
        top_frame.pack(fill="x", padx=10, pady=10)
        ttk.Label(top_frame, text="Model:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.model_menu = ttk.Combobox(top_frame, textvariable=self.model_var, values=self._get_model_folders(), state="readonly")
        self.model_menu.grid(row=0, column=1, sticky="ew", padx=5, pady=5)
        ttk.Label(top_frame, text="Režim:").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.mode_menu = ttk.Combobox(top_frame, textvariable=self.mode_var, values=["Vlastní sekvence", "Skupina enzymů", "Soubor"], state="readonly")
        self.mode_menu.grid(row=1, column=1, sticky="ew", padx=5, pady=5)
        self.mode_menu.bind("<<ComboboxSelected>>", self._update_input_visibility)
        top_frame.columnconfigure(1, weight=1)

        input_frame = ttk.LabelFrame(self.frame, text="Vstupní data")
        input_frame.pack(fill="x", padx=10, pady=10)
        self.seq_input = tk.Text(input_frame, height=4, width=80)
        self.seq_input.pack(fill="x", padx=5, pady=5)
        self.group_var = tk.StringVar()
        self.group_menu = ttk.Combobox(input_frame, textvariable=self.group_var, values=self.groups.get_all_group_names(), state="readonly")
        self.group_menu.pack(fill="x", padx=5, pady=5)
        self.file_var = tk.StringVar()
        self.file_menu = ttk.Combobox(input_frame, textvariable=self.file_var, values=self.db.get_unique_filenames(), state="readonly")
        self.file_menu.pack(fill="x", padx=5, pady=5)
        self._update_input_visibility()

        bottom_frame = ttk.Frame(self.frame)
        bottom_frame.pack(fill="x", padx=10, pady=10)
        self.predict_button = ttk.Button(bottom_frame, text="Predikuj", command=self._run_prediction)
        self.predict_button.pack(fill="x", pady=5)
        self.output_label = ttk.Label(self.frame, text="Výsledky se objeví zde.", justify="left", wraplength=800)
        self.output_label.pack(fill="both", padx=10, pady=10)

    def _get_model_folders(self):
        model_dir = "MODEL"
        if not os.path.exists(model_dir):
            return []
        return [name for name in os.listdir(model_dir) if os.path.isdir(os.path.join(model_dir, name))]

    def _update_input_visibility(self, *_):
        mode = self.mode_var.get()
        self.seq_input.pack_forget()
        self.group_menu.pack_forget()
        self.file_menu.pack_forget()
        if mode == "Vlastní sekvence":
            self.seq_input.pack(fill="x", padx=5, pady=5)
        elif mode == "Skupina enzymů":
            self.group_menu.pack(fill="x", padx=5, pady=5)
        elif mode == "Soubor":
            self.file_menu.pack(fill="x", padx=5, pady=5)

    def _run_prediction(self):
        # Získat vybraný název modelové složky
        model_folder = self.model_var.get()
        if not model_folder:
            self.output_label.config(text="⚠️ Vyber model.")
            return

        # Sestavit plnou cestu k model.pth uvnitř složky modelu
        model_path = os.path.join("MODEL", model_folder, "model.pth")

        # Vytvořit instanci ProcessPrediction s nastaveným modelem
        self.predictor = ProcessPrediction(self.db, model_path=model_path)
        self.predictor.load_model()

        # Získat režim vstupu
        mode = self.mode_var.get()

        if mode == "Vlastní sekvence":
            seq = self.seq_input.get("1.0", "end").strip()
            if seq:
                result = self.predictor.predict(seq)
                self.output_label.config(text=f"Predikovaný motiv: {result}")
            else:
                self.output_label.config(text="⚠️ Zadejte sekvenci.")

        elif mode == "Skupina enzymů":
            group = self.group_var.get()
            if not group:
                self.output_label.config(text="⚠️ Vyber skupinu.")
                return
            results = []
            ids = self.groups.get_group_ids(group)
            for enzyme_id in ids:
                if isinstance(enzyme_id, tuple):
                    enzyme_id = enzyme_id[0]  # vezme první prvek z tuple
                samples = self.db.get_samples_by_id(enzyme_id)
                for sample in samples:
                    seq = sample[3]  # Ověř, že 4. sloupec je sekvence!
                    result = self.predictor.predict(seq)
                    results.append(f"{sample[1]}: {result}")
            self.output_label.config(text="\n".join(results))

        elif mode == "Soubor":
            filename = self.file_var.get()
            if not filename:
                self.output_label.config(text="⚠️ Vyber soubor.")
                return
            samples = self.db.get_samples_by_filename(filename)
            results = []
            for sample in samples:
                seq = sample[3]  # Opět, 4. sloupec je sekvence
                result = self.predictor.predict(seq)
                results.append(f"{sample[1]}: {result}")
            self.output_label.config(text="\n".join(results))
