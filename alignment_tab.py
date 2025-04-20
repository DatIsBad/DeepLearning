import tkinter as tk
import time
from tkinter import ttk
from processProperties import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
from ProcessFiles import fetch_sequence


class AlignmentTab:
    def __init__(self, parent, db_manager, similarity_module, group_manager):
        self.db = db_manager
        self.sim = similarity_module
        self.groups = group_manager

        self.frame = ttk.Frame(parent)
        self.frame.pack(fill="both", expand=True, padx=10, pady=10)

        self.input_type_a = tk.StringVar(value="Soubor")
        self.input_type_b = tk.StringVar(value="Soubor")
        
        # Hlavní rozdělení na ovládací část (controls_frame) a výslednou tabulku (results_frame)
        self.controls_frame = ttk.Frame(self.frame)
        self.controls_frame.pack(fill="x")

        self.results_frame = ttk.Frame(self.frame)
        self.results_frame.pack(fill="both", expand=True)
        
        self._tree = None
        self._seqs_a = []
        self._seqs_b = []

        self._build_interface()


    # Metoda pro vytvoření ui tabu Zarovnání
    def _build_interface(self):
        # -------------------- Sekce-1 (Vstupy) --------------------
        input_frame = ttk.LabelFrame(self.controls_frame, text="Zarovnat: A proti B")
        input_frame.pack(fill="x", pady=5)

        # Vstup A
        # Soubor / Enzym / Skupina / Vlastní sekvence
        row_a = ttk.Frame(input_frame)
        row_a.pack(fill="x", pady=2)
        ttk.Label(row_a, text="Vstup A:").pack(side="left", padx=5)
        self.combo_type_a = ttk.Combobox(
            row_a,
            textvariable=self.input_type_a,
            state="readonly",
            values=["Soubor", "Enzym", "Skupina", "Vlastní sekvence"],
        )
        self.combo_type_a.pack(side="left")
        self.combo_type_a.bind(
            "<<ComboboxSelected>>", lambda e: self._update_input_widget("A")
        )
        self.input_widget_a = ttk.Combobox(row_a, state="readonly")
        self.input_widget_a.pack(side="left", padx=10)

        # Vstup B 
        # Soubor / Enzym / Skupina / Vlastní sekvence
        row_b = ttk.Frame(input_frame)
        row_b.pack(fill="x", pady=2)
        ttk.Label(row_b, text="Vstup B:").pack(side="left", padx=5)
        self.combo_type_b = ttk.Combobox(
            row_b,
            textvariable=self.input_type_b,
            state="readonly",
            values=["Soubor", "Enzym", "Skupina", "Vlastní sekvence"],
        )
        self.combo_type_b.pack(side="left")
        self.combo_type_b.bind(
            "<<ComboboxSelected>>", lambda e: self._update_input_widget("B")
        )
        self.input_widget_b = ttk.Combobox(row_b, state="readonly")
        self.input_widget_b.pack(side="left", padx=10)

        # Algoritmy
        # Výběr algoritmu pro výpočet zarovnání
        algo_frame = ttk.LabelFrame(self.controls_frame, text="Algoritmus")
        algo_frame.pack(fill="x", pady=5)
        self.algo_var = tk.StringVar(value="needleman_wunsch")
        ttk.Radiobutton(
            algo_frame,
            text="Needleman-Wunsch",
            variable=self.algo_var,
            value="needleman_wunsch",
        ).pack(side="left", padx=5)
        ttk.Radiobutton(
            algo_frame,
            text="Smith-Waterman",
            variable=self.algo_var,
            value="smith_waterman",
        ).pack(side="left", padx=5)

        # Tlačítka
        # tlačítka která spustí algoritmus pro zarovnání
        button_frame = ttk.Frame(self.controls_frame)
        button_frame.pack(pady=10)
        ttk.Button(button_frame, text="Zarovnat", command=self._on_align).pack(
            side="left", padx=5
        )
        self.cancel_requested = False
        ttk.Button(button_frame, text="Zrušit", command=self._cancel_alignment).pack(
            side="left", padx=5
        )

        # Progressbar
        # bar který uživateli nějak vyzobrazí jak dlouho nám bude trvat než se všechno dopočítá
        self.progress_var = tk.DoubleVar(value=0)
        self.progress_bar = ttk.Progressbar(
            self.controls_frame, variable=self.progress_var, maximum=100
        )
        self.progress_bar.pack(fill="x", padx=10, pady=5)

        self.progress_label = ttk.Label(self.controls_frame, text="")
        self.progress_label.pack()

        # Scrollovací rám s tabulkou
        self.table_container = ttk.Frame(self.results_frame)
        self.table_container.pack(fill="both", expand=True)

        self.tree_scroll_y = ttk.Scrollbar(self.table_container, orient="vertical")
        self.tree_scroll_y.pack(side="right", fill="y")
        self.tree_scroll_x = ttk.Scrollbar(self.table_container, orient="horizontal")
        self.tree_scroll_x.pack(side="bottom", fill="x")

        self.table_frame = ttk.Frame(self.table_container)
        self.table_frame.pack(side="left", fill="both", expand=True)

        self._update_input_widget("A")
        self._update_input_widget("B")

    def _update_input_widget(self, side):
        var = self.input_type_a if side == "A" else self.input_type_b
        current_frame = self.input_widget_a if side == "A" else self.input_widget_b
        parent = current_frame.master
        current_frame.destroy()

        if var.get() == "Vlastní sekvence":
            new_widget = tk.Text(parent, height=4, width=40)
        else:
            new_widget = ttk.Combobox(parent, state="readonly")
            if var.get() == "Soubor":
                new_widget["values"] = self.db.get_unique_filenames()
            elif var.get() == "Enzym":
                enzymes = sorted(set(row[1] for row in self.db.get_all_samples()))
                new_widget["values"] = enzymes
            elif var.get() == "Skupina":
                new_widget["values"] = self.groups.get_all_group_names()

        new_widget.pack(side="left", padx=10)
        if side == "A":
            self.input_widget_a = new_widget
        else:
            self.input_widget_b = new_widget

    def _cancel_alignment(self):
        self.cancel_requested = True
        self.alignment_completed = False

    def check_for_fragment(self, enzyme_name):
        for row in self.db.get_all_samples():
            if row[1] == enzyme_name:
                return bool(row[6])
        return False


    # Metoda která spustí zarovnání na základě vstupů A a B
    def _on_align(self):
        # --------------------- 
        # Událost kliknutí na buňku v tabulce výsledků – zobrazí detailní zarovnání
        def on_select(event):
            # získání x a y v tabulce, pro zjištění které zarovnání bylo zvoleno
            item = self._tree.identify_row(event.y)
            col = self._tree.identify_column(event.x)

            if not item or not col:
                return
            
            # odečíst 1 v "col_idx" kvůli názvu sloupce, 
            # ale tabulka má navíc i první sloupec s názvy enzymů – tedy posunout o 1 ještě dále
            row_idx = self._tree.index(item)
            col_idx = (
                int(col.replace("#", "")) - 1
            )  
            col_idx -= 1
            if col_idx < 0:
                return
            name_a, seq_a = self._seqs_a[row_idx]

            # ochrana proti kliknutí mimo rozsah
            if col_idx >= len(self._seqs_b):
                return  
            name_b, seq_b = self._seqs_b[col_idx]

            # Zarovnání se vypočítá znovu a data se pošlou do metody _show_alignment_detail
            # Zde se pak podrobněji vyzobrazí data o zarovnání
            if self.algo_var.get() == "needleman_wunsch":
                aligned1, aligned2, score = self.sim.needleman_wunsh_alignment(
                    seq_a, seq_b
                )
                self._show_alignment_detail(
                    name_a, name_b, seq_a, seq_b, aligned1, aligned2, score
                )
            else:
                alignment = self.sim.smith_waterman_alignment(seq_a, seq_b)
                score = alignment[0][2]
                aligned1 = []
                aligned2 = []
                for item in alignment:
                    aligned1.append(item[0])
                    aligned2.append(item[1])

                self._show_alignment_detail(
                    name_a, name_b, seq_a, seq_b, aligned1, aligned2, score
                )

        # --------------------- 
        # Metoda pro získání sekvencí podle typu vstupu
   
        def fetch_sequences_with_ids(widget, source_type):
            if source_type == "Vlastní sekvence":
                return [("Vlastní_A", widget.get("1.0", "end").strip().replace("\n", "").replace(" ", ""), None)]

            elif source_type in ("Soubor", "Enzym"):
                identifier = widget.get()
                return [
                    (row[1], fetch_sequence(row[7], row[8]), row[0])  # (name, sequence, id)
                    for row in self.db.get_all_samples()
                    if (source_type == "Soubor" and row[7] == identifier) or
                    (source_type == "Enzym" and row[1] == identifier)
                ]

            elif source_type == "Skupina":
                group_name = widget.get()
                group_data = self.groups.get_group_ids(group_name)
                return [
                    (row[1], row[9], row[0])
                    for row in group_data
                    if isinstance(row, tuple) and len(row) >= 10
                ]

            else:
                return [("Vlastní_B", widget.get().strip().replace("\n", "").replace(" ", ""), None)]
            

         # -----------------------------
        
        # Hlavní logika
        type_a = self.input_type_a.get()
        type_b = self.input_type_b.get()
        algorithm = self.algo_var.get()

        seqs_a = fetch_sequences_with_ids(self.input_widget_a, type_a)
        seqs_b = fetch_sequences_with_ids(self.input_widget_b, type_b)
        self._seqs_a = [(n, s) for n, s, _ in seqs_a]
        self._seqs_b = [(n, s) for n, s, _ in seqs_b]

        if not seqs_a or not seqs_b:
            return

        is_single = len(seqs_a) == 1 and len(seqs_b) == 1
        if is_single:
            name_a, seq_a, _ = seqs_a[0]
            name_b, seq_b, _ = seqs_b[0]
            if algorithm == "needleman_wunsch":
                aligned1, aligned2, score = self.sim.needleman_wunsh_alignment(seq_a, seq_b)
            else:
                alignment = self.sim.smith_waterman_alignment(seq_a, seq_b)
                score = alignment[0][2]
                aligned1 = [a[0] for a in alignment]
                aligned2 = [a[1] for a in alignment]
            self._show_alignment_detail(name_a, name_b, seq_a, seq_b, aligned1, aligned2, score)
            return

        # -----------------------------
        # Příprava kombinací pro výpočet
        id_pairs = []
        id_pair_to_seq = {}
        for name_a, sa, id_a in seqs_a:
            for name_b, sb, id_b in seqs_b:
                if id_a is not None and id_b is not None:
                    try:
                        key = tuple(sorted((int(id_a), int(id_b))))
                        id_pairs.append(key)
                        id_pair_to_seq[key] = (sa, sb)
                    except:
                        continue

        known_scores = {}
        chunk_size = 400
        for i in range(0, len(id_pairs), chunk_size):
            chunk = id_pairs[i:i+chunk_size]
            partial = self.db.get_scores_for_pairs(chunk, algorithm)
            known_scores.update(partial)
        to_compute = [pair for pair in id_pairs if pair not in known_scores]

        # -----------------------------
        # Tabulka výstupu
        self.table_frame.pack_forget()
        self.table_frame.pack(fill="both", expand=True)
        for widget in self.table_frame.winfo_children():
            widget.destroy()

        columns = ["enzym"] + [name for name, _, _ in seqs_b]
        self._tree = ttk.Treeview(self.table_frame, columns=columns, show="headings")
        for col in columns:
            self._tree.heading(col, text=col)
        self._tree.pack(fill="both", expand=True)
        self._tree.bind("<ButtonRelease-1>", on_select)
        
        # Vytvoření tabulky, pokud ještě neexistuje
        columns = ["enzym"] + [name for name, _, _ in seqs_b]
        if not self._tree:
            self._tree = ttk.Treeview(self.table_frame, columns=columns, show="headings")
            for col in columns:
                self._tree.heading(col, text=col)
            self._tree.pack(fill="both", expand=True)
            self._tree.bind("<ButtonRelease-1>", on_select)
        else:
            self._tree.delete(*self._tree.get_children())
            self._tree["columns"] = columns
            for col in columns:
                self._tree.heading(col, text=col)


        # -----------------------------
        # Výpočet chybějících skóre + progressbar
        total = len(seqs_a) * len(seqs_b)
        count = 0
        start_time = time.time()
        to_insert = []

        for name_a, sa, id_a in seqs_a:
            row_scores = []

            if self.cancel_requested:
                break

            for name_b, sb, id_b in seqs_b:
                if self.cancel_requested:
                    break

                score = None
                if id_a is not None and id_b is not None:
                    try:
                        key = tuple(sorted((int(id_a), int(id_b))))
                        if key in known_scores:
                            score = known_scores[key]
                        else:
                            if algorithm == "needleman_wunsch":
                                seq1, seq2, score = self.sim.needleman_wunsh_alignment(sa, sb)
                            else:
                                alignment = self.sim.smith_waterman_alignment(sa, sb)
                                seq1, seq2, score = alignment[0][0], alignment[0][1], alignment[0][2] if alignment else 0
                            known_scores[key] = score
                            to_insert.append((*key, algorithm, score, seq1, seq2))
                    except:
                        score = 0
                else:
                    if algorithm == "needleman_wunsch":
                        _, _, score = self.sim.needleman_wunsh_alignment(sa, sb)
                    else:
                        alignment = self.sim.smith_waterman_alignment(sa, sb)
                        score = alignment[0][2] if alignment else 0

                row_scores.append(score)
                count += 1
                self.progress_var.set((count / total) * 100)

                elapsed = time.time() - start_time
                if count > 0:
                    est_total = elapsed / count * total
                    remaining = max(0, est_total - elapsed)
                    self.progress_label.config(
                        text=f"Zpracováno: {count} / {total} | Odhadovaný čas: {int(remaining)} s"
                    )
                self.progress_label.update()
                self.progress_bar.update()

            self._tree.insert("", "end", values=[name_a] + row_scores)

        self.progress_var.set(0)
        self.progress_label.config(text="Zarovnání dokončeno.")
        self.cancel_requested = False

        if to_insert:
            self.db.insert_alignments_bulk(to_insert)



        

    # Metoda pro zobrazení detailů zarovnání dvou enzymů
    # Nové okno zobrazí základní statistiky (shody, neshody, mezery, identita), distribuci AK,
    # barevně zvýrazněné zarovnání a volitelnou heatmapu matic.
    def _show_alignment_detail(
        self, name_a, name_b, seq_a, seq_b, aligned1, aligned2, score
    ):
        # Výpočet základních statistik zarovnání
        matches = sum(1 for x, y in zip(aligned1, aligned2) if x == y)
        mismatches = sum(
            1 for x, y in zip(aligned1, aligned2) if x != y and x != "-" and y != "-"
        )
        gaps = aligned1.count("-") + aligned2.count("-")
        identity = (
            matches / max(len(aligned1), len(aligned2)) * 100
            if max(len(aligned1), len(aligned2)) > 0
            else 0
        )

        # Vytvoření nového okna pro zobrazení detailu zarovnání
        top = tk.Toplevel(self.frame)
        top.title(f"Detail zarovnání: {name_a} vs {name_b}")
        top.geometry("800x500")

        info_frame = ttk.Frame(top)
        info_frame.pack(fill="x")

        # Výpis textových statistik do horní části okna
        stats_text = tk.Text(info_frame, height=10, wrap="word")
        stats_text.pack(fill="x")
        stats_text.insert("end", f"Skóre: {score}\n")
        stats_text.insert("end", f"Délka A: {len(seq_a)}, Délka B: {len(seq_b)}\n")
        stats_text.insert(
            "end", f"Shody: {matches}, Neshody: {mismatches}, Mezery: {gaps}\n"
        )
        stats_text.insert("end", f"Identita: {identity:.2f}%\n")

        def protein_distribution(protein_sequence):
            AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
            total = len(protein_sequence)
            return [protein_sequence.count(aa) / total if total > 0 else 0 for aa in AMINO_ACIDS]

        dist_a = protein_distribution(seq_a)
        dist_b = protein_distribution(seq_b)
        symbols = "ACDEFGHIKLMNPQRSTVWY"
        stats_text.insert("end", "Distribuce A:\n")
        stats_text.insert(
            "end",
            ", ".join(f"{aa}:{val:.2f}" for aa, val in zip(symbols, dist_a)) + "\n",
        )
        stats_text.insert("end", "Distribuce B:\n")
        stats_text.insert(
            "end",
            ", ".join(f"{aa}:{val:.2f}" for aa, val in zip(symbols, dist_b)) + "\n",
        )
        stats_text.configure(state="disabled")

        # Frame s zarovnaným textem (barveně zvýrazněné shody/neshody)
        scroll_frame = ttk.Frame(top)
        scroll_frame.pack(fill="both", expand=True)
        scroll_frame.grid_rowconfigure(0, weight=1)
        scroll_frame.grid_columnconfigure(0, weight=1)

        text = tk.Text(scroll_frame, wrap="none", height=10)
        text.grid(row=0, column=0, sticky="nsew")
        x_scroll = ttk.Scrollbar(scroll_frame, orient="horizontal", command=text.xview)
        x_scroll.grid(row=1, column=0, sticky="ew")
        y_scroll = ttk.Scrollbar(scroll_frame, orient="vertical", command=text.yview)
        y_scroll.grid(row=0, column=1, sticky="ns")
        text.configure(xscrollcommand=x_scroll.set, yscrollcommand=y_scroll.set)

        # Definice barev pro vizuální rozlišení shod/neshod
        text.tag_configure("match", foreground="green")
        text.tag_configure("mismatch", foreground="red")
        text.tag_configure("gap", foreground="blue")
        
        if self.algo_var.get() == "needleman_wunsch":
            text.insert("end", "A: ")
            for i, (a, b) in enumerate(zip(aligned1, aligned2)):
                tag = (
                    "match"
                    if a == b
                    else ("gap" if a == "-" or b == "-" else "mismatch")
                )
                text.insert("end", a, tag)
            text.insert("end", "\nB: ")
            for i, (a, b) in enumerate(zip(aligned1, aligned2)):
                tag = (
                    "match"
                    if a == b
                    else ("gap" if a == "-" or b == "-" else "mismatch")
                )
                text.insert("end", b, tag)

        else:
            text.insert("end", "Nalezené lokální sekvence (unikátní):\n")
            text.insert("end", f"{aligned1} ")

        # Stavové proměnné pro heatmapu
        heatmap_canvas = None
        heatmap_visible = False

        # Metoda pro zobrazí/skrití heatmapy s trasováním zarovnání
        def toggle_heatmap():
            nonlocal heatmap_canvas, heatmap_visible

            if heatmap_visible:
                if heatmap_canvas:
                    heatmap_canvas.get_tk_widget().pack_forget()
                    heatmap_canvas = None
                scroll_frame.pack(fill="both", expand=True)
                heatmap_visible = False
            else:
                scroll_frame.pack_forget()
                fig = None

                if self.algo_var.get() == "needleman_wunsch":
                    score_matrix, pointer_matrix = self.sim.needleman_wunsh_matrix(
                        seq_a, seq_b
                    )
                    fig = self.sim.plot_needleman_wunsch_heatmap_with_trace(
                        score_matrix, pointer_matrix, return_fig=True
                    )

                else:
                    score_matrix, pointer_matrix, max_score, max_positions = (
                        self.sim.smith_waterman_matrix(seq_a, seq_b)
                    )
                    smith_data = self.sim.smith_waterman_alignment(
                        seq_a, seq_b, all_maxima=True
                    )

                    trace_coords = []
                    for item in smith_data:
                        print(item)
                        trace_coords.append(item[4])

                    fig = Figure(figsize=(6, 4))
                    fig = self.sim.plot_smith_waterman_heatmap_with_trace(
                        score_matrix, trace_coords, return_fig=True
                    )

                heatmap_canvas = FigureCanvasTkAgg(fig, master=top)
                heatmap_canvas.draw()
                heatmap_canvas.get_tk_widget().pack(fill="both", expand=True)
                heatmap_visible = True

        # tlačítko
        ttk.Button(info_frame, text="Zobrazit heatmapu", command=toggle_heatmap).pack(
            pady=5
        )

