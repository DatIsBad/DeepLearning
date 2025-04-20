import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from processProperties import ProcessProperties
from itertools import product

class ClusteringTab:
    def __init__(self, parent, db_manager, group_manager, switch_to_alignment_tab):
        self.db = db_manager
        self.groups = group_manager
        self.processor = ProcessProperties()
        self.processed_data = None
        self.switch_to_alignment_tab = switch_to_alignment_tab

        self.frame = ttk.Frame(parent)
        self.frame.pack(fill="both", expand=True)

        self._build_interface()

    def _build_interface(self):
        left_frame = ttk.Frame(self.frame)
        left_frame.pack(side="left", fill="y", padx=10, pady=10)

        ttk.Label(left_frame, text="Vyber první vstup (soubor nebo skupina):").pack(anchor="w")
        self.source_a_var = tk.StringVar()
        self.source_a_combo = ttk.Combobox(left_frame, textvariable=self.source_a_var, state="readonly")
        self.source_a_combo.pack(fill="x", pady=2)

        ttk.Label(left_frame, text="Vyber druhý vstup (soubor nebo skupina):").pack(anchor="w")
        self.source_b_var = tk.StringVar()
        self.source_b_combo = ttk.Combobox(left_frame, textvariable=self.source_b_var, state="readonly")
        self.source_b_combo.pack(fill="x", pady=2)

        ttk.Label(left_frame, text="Zarovnávací algoritmus:").pack(anchor="w", pady=(10, 0))
        self.algo_var = tk.StringVar(value="needleman_wunsch")
        algo_combo = ttk.Combobox(left_frame, textvariable=self.algo_var, state="readonly")
        algo_combo["values"] = ["needleman_wunsch", "smith_waterman"]
        algo_combo.pack(fill="x")

        ttk.Label(left_frame, text="Vyber vlastnosti:").pack(anchor="w", pady=(10, 0))
        self.feature_vars = {}
        feature_names = [
            "alignment_score", "distribution", "length", "is_fragment",
            "molecular_weight", "isoelectric_point", "gravy",
            "instability_index", "aromaticity", "aliphatic_index"
        ]
        for name in feature_names:
            var = tk.BooleanVar(value=False)
            cb = ttk.Checkbutton(left_frame, text=name, variable=var)
            cb.pack(anchor="w")
            self.feature_vars[name] = var

        ttk.Label(left_frame, text="Shlukovací algoritmus:").pack(anchor="w", pady=(10, 0))
        self.clustering_algo_var = tk.StringVar(value="kmeans")
        cluster_combo = ttk.Combobox(left_frame, textvariable=self.clustering_algo_var, state="readonly")
        cluster_combo["values"] = ["kmeans", "agglomerative"]
        cluster_combo.pack(fill="x")

        ttk.Button(left_frame, text="Analyzovat", command=self._on_analyze).pack(pady=5, fill="x")

        ttk.Label(left_frame, text="Počet shluků:").pack(anchor="w", pady=(10, 0))
        self.cluster_count_var = tk.IntVar(value=3)
        ttk.Spinbox(left_frame, from_=2, to=10, textvariable=self.cluster_count_var).pack(fill="x")

        ttk.Label(left_frame, text="Zobrazit graf:").pack(anchor="w", pady=(10, 0))
        self.plot_type_var = tk.StringVar(value="vše")
        plot_options = ["vše", "shluky", "silhouette", "délky"]
        ttk.Combobox(left_frame, textvariable=self.plot_type_var, values=plot_options, state="readonly").pack(fill="x")

        ttk.Button(left_frame, text="Shlukovat", command=self._on_cluster).pack(pady=5, fill="x")

        self.right_frame = ttk.Frame(self.frame)
        self.right_frame.pack(side="left", fill="both", expand=True, padx=10, pady=10)

        self.output_text = tk.Text(self.right_frame, height=1, state="disabled")
        self.output_text.pack(fill="x")

        self.canvas_area = ttk.Frame(self.right_frame)
        self.canvas_area.pack(fill="both", expand=True)

        self.refresh_lists()

    def refresh_lists(self):
        filenames = self.db.get_unique_filenames()
        groups = self.groups.get_all_group_names()
        combined = ["soubor: " + f for f in filenames] + ["skupina: " + g for g in groups]
        self.source_a_combo["values"] = combined
        self.source_b_combo["values"] = combined
        if combined:
            self.source_a_combo.current(0)
            self.source_b_combo.current(min(1, len(combined)-1))

    def _clear_output_text(self):
        self.output_text.configure(state="normal")
        self.output_text.delete("1.0", "end")
        self.output_text.configure(height=1)
        self.output_text.configure(state="disabled")

    def _show_text(self, content):
        self._clear_output_text()
        self.output_text.configure(state="normal")
        self.output_text.insert("end", content + "\n\n")
        self.output_text.see("end")
        new_height = int(self.output_text.index("end-1c").split(".")[0])
        self.output_text.configure(height=new_height)
        self.output_text.configure(state="disabled")

    def _clear_canvas(self):
        for widget in self.canvas_area.winfo_children():
            widget.destroy()

    def _show_figure(self, fig):
        self._clear_canvas()
        canvas = FigureCanvasTkAgg(fig, master=self.canvas_area)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True, pady=5)
        toolbar_frame = ttk.Frame(self.canvas_area)
        toolbar_frame.pack(fill="x")
        toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
        toolbar.update()

    def _on_analyze(self):
        selected = {k: v.get() for k, v in self.feature_vars.items()}
        if sum(selected.values()) < 3:
            messagebox.showwarning("Výběr vlastností", "Vyber alespoň 3 vlastnosti.")
            return

        self.processor.choose_features(selected)

        def parse_input(value):
            if value.startswith("soubor: "):
                return [row[0] for row in self.db.get_samples_by_filename(value[8:])]
            elif value.startswith("skupina: "):
                return [int(row[0]) for row in self.groups.get_group_ids(value[9:])]
            return []

        ids_a = parse_input(self.source_a_var.get())
        ids_b = parse_input(self.source_b_var.get())

        if not ids_a or not ids_b:
            messagebox.showwarning("Neplatný výběr", "Vyber platné dva zdroje dat.")
            return

        algorithm = self.algo_var.get()
        id_pairs = [(a, b) for a, b in product(ids_a, ids_b) if a != b]
        raw_data = {}
        chunk_size = 400
        for i in range(0, len(id_pairs), chunk_size):
            chunk = id_pairs[i:i + chunk_size]
            chunk_result = self.db.get_alignment_data(chunk, algorithm)
            raw_data.update(chunk_result)
            self._show_text(f"Načteno zarovnání pro {min(i + chunk_size, len(id_pairs))} z {len(id_pairs)} dvojic")

        if not raw_data:
            result = messagebox.askyesno("Chybí zarovnání", "Pro některé dvojice není dostupné zarovnání. Chceš přejít do záložky Zarovnání a vypočítat je?")
            if result:
                self.switch_to_alignment_tab()
            return

        alignment_data = []
        for (idA, idB), data in raw_data.items():
            nameA = self.db.get_samples_by_id(idA)[0][1]
            nameB = self.db.get_samples_by_id(idB)[0][1]
            score = data["score"]
            seqA = data["aligned_seq_a"].replace("-", "")
            seqB = data["aligned_seq_b"].replace("-", "")
            fragA = self.db.get_samples_by_id(idA)[0][6]
            fragB = self.db.get_samples_by_id(idB)[0][6]
            alignment_data.append((nameA, nameB, score, seqA, seqB, fragA, fragB))

        self.processed_data = self.processor.generate_vector_matrix_from_alignment(alignment_data)
        self.analysis = self.processor.process_properties(self.processed_data)

        fig = self.processor.plot_analysis_clean(self.analysis["pca"], self.analysis["names"])
        self._show_figure(fig)

        variance = self.analysis["variance"]
        summary = (
            f"PCA 1: {variance[0]*100:.1f}%\n"
            f"PCA 2: {variance[1]*100:.1f}%\n"
            f"Celkem vysvětleno: {sum(variance[:2])*100:.1f}%\n"
            f"Počet enzymových dvojic: {len(self.analysis['names'])}"
        )
        self._show_text(summary)

    def _on_cluster(self):
        self._clear_output_text()
        if not self.processed_data:
            messagebox.showwarning("Nejprve analyzuj", "Proveď nejdříve předběžnou analýzu.")
            return

        method = self.clustering_algo_var.get()
        n_clusters = self.cluster_count_var.get()

        names, labels, coords, vectors = self.processor.prepare_and_cluster(self.processed_data, n_clusters, method)
        pca_variance = self.processor.run_pca(vectors)[1]
        stats = ProcessProperties.get_cluster_statistics(labels, vectors, pca_variance)

        text = ""
        for cluster_id, count in stats["counts"].items():
            text += f"Shluk {cluster_id}: {count} enzymů\n"
        text += f"\nSilhouette skóre: {stats['silhouette_score']:.3f}"
        self._show_text(text)

        selected = self.plot_type_var.get()
        if selected in ("vše", "shluky"):
            self._show_figure(self.processor.plot_clusters_clean(coords, labels))
        if selected in ("vše", "silhouette"):
            self._show_figure(self.processor.plot_silhouette(vectors, labels))
        if selected in ("vše", "délky"):
            self._show_figure(self.processor.plot_length_histogram_by_cluster(vectors, labels))

        variance = self.processor.run_pca(vectors)[1]
        summary = (
            f"\nPCA rozptyl:\n"
            f" - První složka: {variance[0]*100:.2f}%\n"
            f" - Druhá složka: {variance[1]*100:.2f}%\n"
            f"Celkově zachyceno: {sum(variance[:2])*100:.2f}%"
        )
        self._show_text(summary)
