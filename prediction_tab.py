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
        
        # Model selection + Train button
        ttk.Label(top_frame, text="Model:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        
        model_train_frame = ttk.Frame(top_frame)
        model_train_frame.grid(row=0, column=1, columnspan=2, sticky="ew", padx=5, pady=5)
        model_train_frame.columnconfigure(0, weight=3)
        model_train_frame.columnconfigure(1, weight=1)

        self.model_menu = ttk.Combobox(model_train_frame, textvariable=self.model_var, values=self._get_model_folders(), state="readonly")
        self.model_menu.grid(row=0, column=0, sticky="ew", padx=(0, 5))

        self.train_button = ttk.Button(model_train_frame, text="Vytrénovat model", command=self._train_model)
        self.train_button.grid(row=0, column=1, sticky="ew")

        # Mode selection
        ttk.Label(top_frame, text="Režim:").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.mode_menu = ttk.Combobox(top_frame, textvariable=self.mode_var, values=["Vlastní sekvence", "Skupina enzymů", "Soubor"], state="readonly")
        self.mode_menu.grid(row=1, column=1, sticky="ew", padx=5, pady=5, columnspan=2)

        self.mode_menu.bind("<<ComboboxSelected>>", self._update_input_visibility)
        top_frame.columnconfigure(1, weight=1)
        top_frame.columnconfigure(2, weight=1)

        # Input data section
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

        # Predict button
        bottom_frame = ttk.Frame(self.frame)
        bottom_frame.pack(fill="x", padx=10, pady=10)

        self.predict_button = ttk.Button(bottom_frame, text="Predikuj", command=self._run_prediction)
        self.predict_button.pack(fill="x", pady=5)

        # Scrollable output for prediction results
        output_frame = ttk.LabelFrame(self.frame, text="Výsledky predikce / Tréninku")
        output_frame.pack(fill="both", expand=True, padx=10, pady=10)

        self.output_text = tk.Text(output_frame, wrap="word")
        self.output_text.pack(side="left", fill="both", expand=True)

        scrollbar = ttk.Scrollbar(output_frame, command=self.output_text.yview)
        scrollbar.pack(side="right", fill="y")

        self.output_text.config(yscrollcommand=scrollbar.set)


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
        model_folder = self.model_var.get()
        if not model_folder:
            self.output_text.delete("1.0", "end")
            self.output_text.insert("1.0", "⚠️ Vyber model.")
            return

        model_path = os.path.join("MODEL", model_folder, "model.pth")
        self.predictor = ProcessPrediction(self.db, model_path=model_path)
        self.predictor.load_model()

        mode = self.mode_var.get()
        results = []

        if mode == "Vlastní sekvence":
            seq = self.seq_input.get("1.0", "end").strip()
            if seq:
                result = self.predictor.predict(seq)
                results.append(f"Predikovaný motiv: {result}")
            else:
                results.append("⚠️ Zadejte sekvenci.")

        elif mode == "Skupina enzymů":
            group = self.group_var.get()
            if not group:
                results.append("⚠️ Vyber skupinu.")
            else:
                ids = self.groups.get_group_ids(group)
                for enzyme_id in ids:
                    if isinstance(enzyme_id, tuple):
                        enzyme_id = enzyme_id[0]
                    samples = self.db.get_samples_by_id(enzyme_id)
                    for sample in samples:
                        seq = sample[3]
                        result = self.predictor.predict(seq)
                        results.append(f"{sample[1]}: {result}")

        elif mode == "Soubor":
            filename = self.file_var.get()
            if not filename:
                results.append("⚠️ Vyber soubor.")
            else:
                samples = self.db.get_samples_by_filename(filename)
                for sample in samples:
                    seq = sample[3]
                    result = self.predictor.predict(seq)
                    results.append(f"{sample[1]}: {result}")

        self.output_text.delete("1.0", "end")
        self.output_text.insert("1.0", "\n".join(results))

    def _train_model(self):
        model_folder = self.model_var.get()
        if not model_folder:
            self.output_text.delete("1.0", "end")
            self.output_text.insert("1.0", "⚠️ Vyber model pro trénink.")
            return

        model_path = os.path.join("MODEL", model_folder, "model.pth")
        self.predictor = ProcessPrediction(self.db, model_path=model_path)

        # Chytrý výběr top_n podle velikosti databáze
        total_sequences = len(self.db.get_all_samples())
        if total_sequences > 500:
            top_n = 35
        elif total_sequences > 300:
            top_n = 20
        elif total_sequences > 150:
            top_n = 10
        else:
            top_n = 5

        self.output_text.delete("1.0", "end")
        self.output_text.insert("1.0", f"Trénuji model '{model_folder}' s top_n={top_n}...\n")

        try:
            self.predictor.auto_train(model_name=model_folder, top_n=top_n)
            self.output_text.insert("end", "✅ Trénink dokončen a model uložen.")
        except Exception as e:
            self.output_text.insert("end", f"❌ Chyba při tréninku: {str(e)}")
