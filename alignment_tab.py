import tkinter as tk
import time
from tkinter import ttk
from ProcessFiles import fetch_sequence
from processProperties import protein_distribution
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg



class AlignmentTab:
    def __init__(self, parent, db_manager, similarity_module):
        self.db = db_manager
        self.sim = similarity_module

        self.frame = ttk.Frame(parent)
        self.frame.pack(fill="both", expand=True, padx=10, pady=10)

        self.input_type_a = tk.StringVar(value="Soubor")
        self.input_type_b = tk.StringVar(value="Soubor")

        self._tree = None
        self._seqs_a = []
        self._seqs_b = []

        self._build_interface()

    def _build_interface(self):
        input_frame = ttk.LabelFrame(self.frame, text="Zarovnat: A proti B")
        input_frame.pack(fill="x", pady=5)

        row_a = ttk.Frame(input_frame)
        row_a.pack(fill="x", pady=2)
        ttk.Label(row_a, text="Vstup A:").pack(side="left", padx=5)
        self.combo_type_a = ttk.Combobox(row_a, textvariable=self.input_type_a, state="readonly", values=["Soubor", "Enzym", "Skupina", "Vlastní sekvence"])
        self.combo_type_a.pack(side="left")
        self.combo_type_a.bind("<<ComboboxSelected>>", lambda e: self._update_input_widget('A'))
        self.input_widget_a = ttk.Combobox(row_a, state="readonly")
        self.input_widget_a.pack(side="left", padx=10)

        row_b = ttk.Frame(input_frame)
        row_b.pack(fill="x", pady=2)
        ttk.Label(row_b, text="Vstup B:").pack(side="left", padx=5)
        self.combo_type_b = ttk.Combobox(row_b, textvariable=self.input_type_b, state="readonly", values=["Soubor", "Enzym", "Skupina", "Vlastní sekvence"])
        self.combo_type_b.pack(side="left")
        self.combo_type_b.bind("<<ComboboxSelected>>", lambda e: self._update_input_widget('B'))
        self.input_widget_b = ttk.Combobox(row_b, state="readonly")
        self.input_widget_b.pack(side="left", padx=10)

        algo_frame = ttk.LabelFrame(self.frame, text="Algoritmus")
        algo_frame.pack(fill="x", pady=5)
        self.algo_var = tk.StringVar(value="needleman_wunsch")
        ttk.Radiobutton(algo_frame, text="Needleman-Wunsch", variable=self.algo_var, value="needleman_wunsch").pack(side="left", padx=5)
        ttk.Radiobutton(algo_frame, text="Smith-Waterman", variable=self.algo_var, value="smith_waterman").pack(side="left", padx=5)

        button_frame = ttk.Frame(self.frame)
        button_frame.pack(pady=10)
        ttk.Button(button_frame, text="Zarovnat", command=self._on_align).pack(side="left", padx=5)
        self.cancel_requested = False
        ttk.Button(button_frame, text="Zrušit", command=self._cancel_alignment).pack(side="left", padx=5)

        self.progress_var = tk.DoubleVar(value=0)
        self.progress_bar = ttk.Progressbar(self.frame, variable=self.progress_var, maximum=100)
        self.progress_bar.pack(fill="x", padx=10, pady=5)

        self.progress_label = ttk.Label(self.frame, text="")
        self.progress_label.pack()

        self.results_frame = ttk.Frame(self.frame)
        self.results_frame.pack(fill="both", expand=True)

        self._update_input_widget('A')
        self._update_input_widget('B')

    def _update_input_widget(self, side):
        var = self.input_type_a if side == 'A' else self.input_type_b
        current_frame = self.input_widget_a if side == 'A' else self.input_widget_b
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
                new_widget["values"] = ["Skupina 1", "Skupina 2"]

        new_widget.pack(side="left", padx=10)
        if side == 'A':
            self.input_widget_a = new_widget
        else:
            self.input_widget_b = new_widget

    def _cancel_alignment(self):
        self.cancel_requested = True

    def _on_align(self):
        for widget in self.results_frame.winfo_children():
            widget.destroy()

        type_a = self.input_type_a.get()
        type_b = self.input_type_b.get()

        def fetch_sequences(widget, source_type):
            if source_type == "Vlastní sekvence":
                return [("Vlastní_A", widget.get("1.0", "end").replace("\n", "").replace(" ", "").strip())]
            elif source_type == "Soubor" or source_type == "Enzym":
                id_ = widget.get()
                return [(r[1], fetch_sequence(r[7], r[8])) for r in self.db.get_all_samples()
                        if (source_type == "Soubor" and r[7] == id_) or (source_type == "Enzym" and r[1] == id_)]
            elif source_type == "Skupina":
                return [(f"Skupina_{i}", s) for i, s in enumerate(["ATGC", "AAGC"])]
            else:
                return [("Vlastní_B", widget.get())]

        seqs_a = fetch_sequences(self.input_widget_a, type_a)
        seqs_b = fetch_sequences(self.input_widget_b, type_b)

        self._seqs_a = seqs_a
        self._seqs_b = seqs_b

        if not seqs_a or not seqs_b:
            ttk.Label(self.results_frame, text="Musíš zadat/načíst sekvence pro oba vstupy.").pack()
            return

        is_single = len(seqs_a) == 1 and len(seqs_b) == 1

        if is_single:
            name_a, a = seqs_a[0]
            name_b, b = seqs_b[0]
            if self.algo_var.get() == "needleman_wunsch":
                aligned1, aligned2, score = self.sim.needleman_wunsh_alignment(a, b)
            else:
                alignment = self.sim.smith_waterman_alignment(a, b)
                score =  alignment[0][2]
                aligned1 = []
                aligned2 = []
                for item in alignment:
                    aligned1.append(item[0])
                    aligned2.append(item[1])


            self._show_alignment_detail(name_a, name_b, a, b, aligned1, aligned2, score)
        else:
            table_frame = ttk.Frame(self.results_frame)
            table_frame.pack(fill="both", expand=True)

            x_scroll = ttk.Scrollbar(table_frame, orient="horizontal")
            y_scroll = ttk.Scrollbar(table_frame, orient="vertical")

            columns = ["enzym"] + [name for name, _ in seqs_b]
            self._tree = ttk.Treeview(
                table_frame,
                columns=columns,
                show="headings",
                xscrollcommand=x_scroll.set,
                yscrollcommand=y_scroll.set
            )

            x_scroll.config(command=self._tree.xview)
            y_scroll.config(command=self._tree.yview)
            x_scroll.pack(side="bottom", fill="x")
            y_scroll.pack(side="right", fill="y")
            self._tree.pack(side="left", fill="both", expand=True)

            for col in columns:
                self._tree.heading(col, text=col)

            total = len(seqs_a) * len(seqs_b)
            count = 0
            start_time = time.time()

            for name_a, sa in seqs_a:
                row_scores = []
                for name_b, sb in seqs_b:
                    if self.algo_var.get() == "needleman_wunsch":
                        _, _, score = self.sim.needleman_wunsh_alignment(sa, sb)
                    else:
                        alignment = self.sim.smith_waterman_alignment(sa, sb)
                        score = alignment[0][2] if alignment else 0
                    if self.cancel_requested:
                        break
                    row_scores.append(score)
                    count += 1
                    self.progress_var.set((count / total) * 100)
                    elapsed = time.time() - start_time
                    if count > 0:
                        est_total = elapsed / count * total
                        remaining = max(0, est_total - elapsed)
                        self.progress_label.config(text=f"Zpracováno: {count} / {total} | Odhadovaný čas: {int(remaining)} s")
                    self.progress_label.update()
                    self.progress_bar.update()
                self._tree.insert("", "end", values=[name_a] + row_scores)

            self.progress_var.set(0)
            total_elapsed = int(time.time() - start_time)
            self.progress_label.config(text=f"Zarovnání dokončeno. Čas výpočtu: {total_elapsed} s")
            self.cancel_requested = False

            def on_select(event):
                item = self._tree.identify_row(event.y)
                col = self._tree.identify_column(event.x)
                if not item or not col:
                    return
                row_idx = self._tree.index(item)
                col_idx = int(col.replace('#', '')) - 1  # odečíst 1 kvůli názvu sloupce, ale tabulka má navíc i první sloupec s názvy enzymů – tedy posunout o 1 ještě dále
                col_idx -= 1
                if col_idx < 0:
                    return
                name_a, seq_a = self._seqs_a[row_idx]
                if col_idx >= len(self._seqs_b): return  # ochrana proti kliknutí mimo rozsah
                name_b, seq_b = self._seqs_b[col_idx]

                if self.algo_var.get() == "needleman_wunsch":
                    aligned1, aligned2, score = self.sim.needleman_wunsh_alignment(seq_a, seq_b)
                    self._show_alignment_detail(name_a, name_b, seq_a, seq_b, aligned1, aligned2, score)
                else:
                    alignment = self.sim.smith_waterman_alignment(seq_a, seq_b)
                    score =  alignment[0][2]
                    aligned1 = []
                    aligned2 = []
                    for item in alignment:
                        aligned1.append(item[0])
                        aligned2.append(item[1])
                        
                    self._show_alignment_detail(name_a, name_b, seq_a, seq_b, aligned1, aligned2, score)


            self._tree.bind("<ButtonRelease-1>", on_select)

    def _show_alignment_detail(self, name_a, name_b, seq_a, seq_b, aligned1, aligned2, score):
        import numpy as np
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        matches = sum(1 for x, y in zip(aligned1, aligned2) if x == y)
        mismatches = sum(1 for x, y in zip(aligned1, aligned2) if x != y and x != '-' and y != '-')
        gaps = aligned1.count('-') + aligned2.count('-')
        identity = matches / max(len(aligned1), len(aligned2)) * 100 if max(len(aligned1), len(aligned2)) > 0 else 0

        top = tk.Toplevel(self.frame)
        top.title(f"Detail zarovnání: {name_a} vs {name_b}")
        top.geometry("800x500")

        info_frame = ttk.Frame(top)
        info_frame.pack(fill="x")

        stats_text = tk.Text(info_frame, height=10, wrap="word")
        stats_text.pack(fill="x")
        stats_text.insert("end", f"Skóre: {score}\n")
        stats_text.insert("end", f"Délka A: {len(seq_a)}, Délka B: {len(seq_b)}\n")
        stats_text.insert("end", f"Shody: {matches}, Neshody: {mismatches}, Mezery: {gaps}\n")
        stats_text.insert("end", f"Identita: {identity:.2f}%\n")

        dist_a = protein_distribution(seq_a)
        dist_b = protein_distribution(seq_b)
        symbols = "ACDEFGHIKLMNPQRSTVWY"
        stats_text.insert("end", "Distribuce A:\n")
        stats_text.insert("end", ", ".join(f"{aa}:{val:.2f}" for aa, val in zip(symbols, dist_a)) + "\n")
        stats_text.insert("end", "Distribuce B:\n")
        stats_text.insert("end", ", ".join(f"{aa}:{val:.2f}" for aa, val in zip(symbols, dist_b)) + "\n")
        stats_text.configure(state="disabled")

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

        text.tag_configure("match", foreground="green")
        text.tag_configure("mismatch", foreground="red")
        text.tag_configure("gap", foreground="blue")

        if self.algo_var.get() == "needleman_wunsch":
            text.insert("end", "A: ")
            for i, (a, b) in enumerate(zip(aligned1, aligned2)):
                tag = "match" if a == b else ("gap" if a == "-" or b == "-" else "mismatch")
                text.insert("end", a, tag)
            text.insert("end", "\nB: ")
            for i, (a, b) in enumerate(zip(aligned1, aligned2)):
                tag = "match" if a == b else ("gap" if a == "-" or b == "-" else "mismatch")
                text.insert("end", b, tag)
            for i, char in enumerate(aligned2):
                tag = "match" if aligned1[i] == aligned2[i] else ("gap" if char == "-" or aligned1[i] == "-" else "mismatch")
                text.insert("end", char, tag)
        else:
            text.insert("end", "Nalezené lokální sekvence (unikátní):\n")
            text.insert("end", f"{aligned1} ")

        # Stavové proměnné pro heatmapu
        heatmap_canvas = None
        heatmap_visible = False

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
                    score_matrix, pointer_matrix = self.sim.needleman_wunsh_matrix(seq_a, seq_b)
                    fig = self.sim.plot_needleman_wunsch_heatmap_with_trace(score_matrix, pointer_matrix, return_fig = True)

                else:
                    score_matrix, pointer_matrix , max_score, max_positions = self.sim.smith_waterman_matrix(seq_a, seq_b)
                    smith_data = self.sim.smith_waterman_alignment(seq_a, seq_b, all_maxima= True)

                    trace_coords = []
                    for item in smith_data:
                        print(item)
                        trace_coords.append(item[4])

                    fig = Figure(figsize=(6, 4))
                    fig = self.sim.plot_smith_waterman_heatmap_with_trace(score_matrix, trace_coords, return_fig = True)

                heatmap_canvas = FigureCanvasTkAgg(fig, master=top)
                heatmap_canvas.draw()
                heatmap_canvas.get_tk_widget().pack(fill="both", expand=True)
                heatmap_visible = True

        ttk.Button(info_frame, text="Zobrazit heatmapu", command=toggle_heatmap).pack(pady=5)
