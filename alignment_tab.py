import tkinter as tk
from tkinter import ttk

class AlignmentTab:
    def __init__(self, parent, db_manager, similarity_module):
        self.notebook = parent
        self.db = db_manager
        self.sim = similarity_module

        self.mode = tk.StringVar(value="file_vs_file")
        self.algorithm = tk.StringVar(value="needleman_wunsch")

        self._build_interface()

    def _build_interface(self):
        self.frame = ttk.Frame(self.notebook)
        self.frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Výběr režimu
        mode_frame = ttk.LabelFrame(self.frame, text="Režim porovnání")
        mode_frame.pack(fill="x", pady=5)

        ttk.Radiobutton(mode_frame, text="Soubor vs. Soubor", variable=self.mode, value="file_vs_file", command=self._on_mode_change).pack(side="left", padx=5)
        ttk.Radiobutton(mode_frame, text="Soubor vs. Enzym", variable=self.mode, value="file_vs_enzyme", command=self._on_mode_change).pack(side="left", padx=5)
        ttk.Radiobutton(mode_frame, text="Enzym vs. Enzym", variable=self.mode, value="enzyme_vs_enzyme", command=self._on_mode_change).pack(side="left", padx=5)

        # Výběr vstupních sekvencí
        input_frame = ttk.LabelFrame(self.frame, text="Vstupy")
        input_frame.pack(fill="x", pady=5)

        self.combo_file1 = ttk.Combobox(input_frame, state="readonly")
        self.combo_file2 = ttk.Combobox(input_frame, state="readonly")
        self.combo_enzyme = ttk.Combobox(input_frame, state="readonly")
        self.combo_enzyme1 = ttk.Combobox(input_frame, state="readonly")
        self.combo_enzyme2 = ttk.Combobox(input_frame, state="readonly")

        self.combo_file1.pack(padx=5, pady=2, fill="x")
        self.combo_file2.pack(padx=5, pady=2, fill="x")
        self.combo_enzyme.pack(padx=5, pady=2, fill="x")
        self.combo_enzyme1.pack(padx=5, pady=2, fill="x")
        self.combo_enzyme2.pack(padx=5, pady=2, fill="x")

        # Výběr algoritmu
        algo_frame = ttk.LabelFrame(self.frame, text="Algoritmus")
        algo_frame.pack(fill="x", pady=5)
        ttk.Radiobutton(algo_frame, text="Needleman-Wunsch", variable=self.algorithm, value="needleman_wunsch").pack(side="left", padx=5)
        ttk.Radiobutton(algo_frame, text="Smith-Waterman", variable=self.algorithm, value="smith_waterman").pack(side="left", padx=5)

        # Tlačítko zarovnat
        ttk.Button(self.frame, text="Zarovnat", command=self._on_align).pack(pady=10)

        # Výsledky (zatím placeholder)
        self.results_frame = ttk.Frame(self.frame)
        self.results_frame.pack(fill="both", expand=True)

        self._on_mode_change()

    def _on_mode_change(self):
        mode = self.mode.get()
        self.combo_file1.pack_forget()
        self.combo_file2.pack_forget()
        self.combo_enzyme.pack_forget()
        self.combo_enzyme1.pack_forget()
        self.combo_enzyme2.pack_forget()

        if mode == "file_vs_file":
            self.combo_file1.pack(padx=5, pady=2, fill="x")
            self.combo_file2.pack(padx=5, pady=2, fill="x")
        elif mode == "file_vs_enzyme":
            self.combo_file1.pack(padx=5, pady=2, fill="x")
            self.combo_enzyme.pack(padx=5, pady=2, fill="x")
        elif mode == "enzyme_vs_enzyme":
            self.combo_enzyme1.pack(padx=5, pady=2, fill="x")
            self.combo_enzyme2.pack(padx=5, pady=2, fill="x")

    def _on_align(self):
        # TODO: Přidat logiku výpočtu podle zvoleného režimu a algoritmu
        for widget in self.results_frame.winfo_children():
            widget.destroy()

        label = ttk.Label(self.results_frame, text="Zde budou výsledky zarovnání.")
        label.pack()
