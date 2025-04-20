import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from processProperties import ProcessProperties

class ClusteringTab:
    def __init__(self, parent):
        self.processor = ProcessProperties()
        self.alignment_data = []  # bude nastavena později přes metodu
        self.processed_data = None

        self.frame = ttk.Frame(parent)
        self.frame.pack(fill="both", expand=True)

        self._build_interface()

    def set_alignment_data(self, alignment_data):
        self.alignment_data = alignment_data

    def _build_interface(self):
        left_frame = ttk.Frame(self.frame)
        left_frame.pack(side="left", fill="y", padx=10, pady=10)

        self.feature_vars = {}
        ttk.Label(left_frame, text="Vyber vlastnosti:").pack(anchor="w")
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
        self.algo_var = tk.StringVar(value="kmeans")
        algo_combo = ttk.Combobox(left_frame, textvariable=self.algo_var, state="readonly")
        algo_combo["values"] = ["kmeans", "agglomerative"]
        algo_combo.pack(fill="x")

        ttk.Button(left_frame, text="Analyzovat", command=self._on_analyze).pack(pady=5, fill="x")

        ttk.Label(left_frame, text="Počet shluků:").pack(anchor="w", pady=(10, 0))
        self.cluster_count_var = tk.IntVar(value=3)
        ttk.Spinbox(left_frame, from_=2, to=10, textvariable=self.cluster_count_var).pack(fill="x")

        ttk.Button(left_frame, text="Shlukovat", command=self._on_cluster).pack(pady=5, fill="x")

        self.right_frame = ttk.Frame(self.frame)
        self.right_frame.pack(side="left", fill="both", expand=True, padx=10, pady=10)

        self.output_text = tk.Text(self.right_frame, height=1, state="disabled")
        self.output_text.pack(fill="x")

        ttk.Button(self.right_frame, text="Vymazat výstup", command=self._clear_output_text).pack(pady=(0, 5), anchor="e")

        self.canvas_area = ttk.Frame(self.right_frame)
        self.canvas_area.pack(fill="both", expand=True)

    def _on_analyze(self):
        selected = {k: v.get() for k, v in self.feature_vars.items()}
        if sum(selected.values()) < 3:
            messagebox.showwarning("Výběr vlastností", "Vyber alespoň 3 vlastnosti.")
            return

        self.processor.choose_features(selected)
        self.processed_data = self.processor.generate_vector_matrix_from_alignment(self.alignment_data)
        analysis = self.processor.process_properties(self.processed_data)

        self.analysis_result = analysis  # uložit pro pozdější použití

        fig = self.processor.plot_analysis_clean(analysis["pca"], analysis["names"])
        self._clear_canvas()
        self._show_figure(fig)

        variance = analysis["variance"]
        summary = (
            f"PCA 1: {variance[0]*100:.1f}%\n"
            f"PCA 2: {variance[1]*100:.1f}%\n"
            f"Celkem vysvětleno: {sum(variance[:2])*100:.1f}%\n"
            f"Počet enzymových dvojic: {len(analysis['names'])}"
        )
        self._show_text(summary)

    def _on_cluster(self):
        if not self.processed_data:
            messagebox.showwarning("Nejprve analyzuj", "Proveď nejdříve předběžnou analýzu.")
            return

        method = self.algo_var.get()
        n_clusters = self.cluster_count_var.get()

        names, labels, coords, vectors = self.processor.prepare_and_cluster(self.processed_data, n_clusters, method)
        self.cluster_data = (labels, coords, vectors)

        stats = self.processor.get_cluster_statistics(labels, vectors, self.processor.run_pca(vectors)[1])

        cluster_text = ""
        for k, v in stats["counts"].items():
            cluster_text += f"Shluk {k}: {v} enzymů\n"
        cluster_text += f"\nSilhouette skóre: {stats['silhouette_score']:.3f}"
        self._show_text(cluster_text)

        self._clear_canvas()

        fig1 = self.processor.plot_clusters_clean(coords, labels)
        fig2 = self.processor.plot_silhouette(vectors, labels)
        fig3 = self.processor.plot_length_histogram_by_cluster(vectors, labels)

        self._show_figure(fig1)
        self._show_figure(fig2)
        self._show_figure(fig3)

        variance = self.processor.run_pca(vectors)[1]
        summary = (
            f"\nPCA rozptyl:\n"
            f" - První složka: {variance[0]*100:.2f}%\n"
            f" - Druhá složka: {variance[1]*100:.2f}%\n"
            f"Celkově zachyceno: {sum(variance[:2])*100:.2f}%"
        )
        self._show_text(summary)

    def _clear_output_text(self):
        self.output_text.configure(state="normal")
        self.output_text.delete("1.0", "end")
        self.output_text.configure(height=1)
        self.output_text.configure(state="disabled")

    def _clear_canvas(self):
        for widget in self.canvas_area.winfo_children():
            widget.destroy()

    def _show_figure(self, fig):
        canvas = FigureCanvasTkAgg(fig, master=self.canvas_area)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True, pady=5)

    def _show_text(self, content):
        self.output_text.configure(state="normal")
        self.output_text.insert("end", content + "\n\n")
        self.output_text.see("end")
        new_height = int(self.output_text.index("end-1c").split(".")[0])
        self.output_text.configure(height=new_height)
        self.output_text.configure(state="disabled")
