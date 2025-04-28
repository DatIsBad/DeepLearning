import tkinter as tk
from tkinter import ttk, messagebox
from tkinter import filedialog, simpledialog, messagebox
from tkinter.filedialog import asksaveasfilename
from ProcessFiles import read_enzyme_headers, fetch_sequence
from alignment_tab import AlignmentTab
from clustering_tab import ClusteringTab
from prediction_tab import PredictionTab
from ProcessSimilarity import *



class App:
    def __init__(self, root, db_manager, group_manager, exporter, filter_manager):
        self.root = root
        self.db = db_manager
        self.groups = group_manager
        self.exporter = exporter
        self.filterer = filter_manager
        self.filtered_data = []

        self.root.minsize(900, 600)
        self.calculate_data_ranges()

        # Hlavní horizontální dělení (sidebar + obsah)
        self.main_frame = ttk.Frame(self.root)
        self.main_frame.pack(fill="both", expand=True)

        # Sidebar vlevo
        self.sidebar = ttk.Frame(self.main_frame, width=200)
        self.sidebar.pack(side="left", fill="y")
        self._build_group_sidebar()  # Tato metoda vytvoří obsah postranního panelu

        # Content frame (obsah záložek) vpravo
        self.content_frame = ttk.Frame(self.main_frame)
        self.content_frame.pack(side="left", fill="both", expand=True)

        # Notebook obsažený v content_frame
        self.notebook = ttk.Notebook(self.content_frame)
        self.notebook.pack(fill="both", expand=True)

        # Složka "Data" Notebooku
        self.data_tab = ttk.Frame(self.notebook)
        self.build_gui(self.data_tab)
        self.notebook.add(self.data_tab, text="Data")

        # Složka "Zarovnání" Notebooku
        self.alignment_tab = AlignmentTab(self.notebook, self.db, ProcessSimilarity(), self.groups)
        self.notebook.add(self.alignment_tab.frame, text="Zarovnání")

        # Složka "Shlukování" Notebooku
        self.clustering_tab = ClusteringTab(
            self.notebook,
            self.db,
            self.groups,
            switch_to_alignment_tab=self._switch_to_alignment_tab
        )
        self.notebook.add(self.clustering_tab.frame, text="Shlukování")
        
        # Složka "Predikce" Notebooku
        self.prediction_tab  = PredictionTab(
            self.notebook,
            self.db,
            self.groups
        )
        self.notebook.add(self.prediction_tab.frame, text="Predikce sekvencí")




        self.update_filenames()
        self.populate_tree(self.db.get_all_samples())

    def _switch_to_alignment_tab(self):
        self.notebook.select(self.alignment_tab.frame)

    def build_gui(self, parent):
        self.root.title("Restrikční enzymy – databázový nástroj")

        self.control_frame = ttk.LabelFrame(parent, text="Ovládání")
        self.control_frame.pack(fill="x", padx=5, pady=5)

        #ttk.Button(self.control_frame, text="Vytvořit databázi", command=self.db.create_database).grid(row=0, column=0, padx=5)
        ttk.Button(self.control_frame, text="Načíst soubor", command=self.load_file).grid(row=0, column=1, padx=5)
        ttk.Button(self.control_frame, text="Reset databáze", command=self.reset_database).grid(row=0, column=2, padx=5)
        ttk.Button(self.control_frame, text="Export filtrované data do CSV", command=self.export_csv).grid(row=0, column=3, padx=5)
        ttk.Button(self.control_frame, text="Smazat soubory", command=self.open_delete_window).grid(row=0, column=4, padx=5)

        filter_frame = ttk.LabelFrame(parent, text="Filtry")
        filter_frame.pack(fill="x", padx=10, pady=5)

        tree_frame = ttk.Frame(parent)
        tree_frame.pack(fill="both", expand=True, padx=10, pady=5)
        tree_frame.columnconfigure(0, weight=1)
        tree_frame.rowconfigure(0, weight=1)

        self._add_group_button_to_controls()

        filter_top = ttk.Frame(filter_frame)
        filter_top.pack(fill="x")
        ttk.Label(filter_top, text="Enzym:").pack(side="left")
        self.filter_enzyme = ttk.Entry(filter_top, width=15)
        self.filter_enzyme.pack(side="left", padx=5)

        self.filter_fragment = tk.BooleanVar()
        ttk.Checkbutton(filter_top, text="Pouze fragmenty", variable=self.filter_fragment).pack(side="left", padx=5)

        ttk.Label(filter_top, text="Soubor:").pack(side="left")
        self.filter_filename = tk.StringVar()
        self.filename_combo = ttk.Combobox(filter_top, textvariable=self.filter_filename, state="readonly", width=20)
        self.filename_combo.pack(side="left", padx=5)

        filter_mid = ttk.Frame(filter_frame)
        filter_mid.pack(fill="x", pady=5)
        ttk.Label(filter_mid, text="Rozsah délky:").pack(side="left")
        self.length_range_label = ttk.Label(filter_mid, text=f"{self.min_seq_length} - {self.max_seq_length}")
        self.length_range_label.pack(side="left", padx=10)

        self.length_min = tk.IntVar(value=self.min_seq_length)
        self.length_max = tk.IntVar(value=self.max_seq_length)

        self.length_slider_min = ttk.Scale(filter_mid, from_=self.min_seq_length, to=self.max_seq_length, orient="horizontal", command=self.update_length_label, variable=self.length_min, length=200)
        self.length_slider_min.pack(side="left", padx=5)

        self.length_slider_max = ttk.Scale(filter_mid, from_=self.min_seq_length, to=self.max_seq_length, orient="horizontal", command=self.update_length_label, variable=self.length_max, length=200)
        self.length_slider_max.set(self.max_seq_length)
        self.length_slider_max.pack(side="left", padx=5)

        filter_bottom = ttk.Frame(filter_frame)
        filter_bottom.pack(fill="x")
        ttk.Label(filter_bottom, text="Min počet rozpoznávacích míst:").pack(side="left")
        self.filter_min_sites = ttk.Entry(filter_bottom, width=10)
        self.filter_min_sites.pack(side="left", padx=5)

        ttk.Button(filter_bottom, text="Použít filtr", command=self.apply_filter).pack(side="left", padx=10)

        self.tree = ttk.Treeview(
        tree_frame,
            columns=("id", "enzyme", "sample", "orf", "rec_seq", "size", "fragment", "filename", "line"),
            show="headings"
        )
        self.tree.tag_configure("evenrow", background="#f0f0ff")
        self.tree.tag_configure("oddrow", background="#ffffff")

        for col in self.tree["columns"]:
            self.tree.heading(col, text=col, command=lambda _col=col: self.sort_column(_col, False))
            self.tree.column(col, width=100)

        vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(tree_frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        tree_frame.rowconfigure(0, weight=1)
        tree_frame.columnconfigure(0, weight=1)

    def update_length_label(self, _):
        min_val = self.length_min.get()
        max_val = self.length_max.get()
        if min_val > max_val:
            self.length_min.set(max_val)
        self.length_range_label.config(text=f"{self.length_min.get()} - {self.length_max.get()}")

    def apply_filter(self):
        enzyme = self.filter_enzyme.get().strip()
        only_fragment = self.filter_fragment.get()
        filename = self.filter_filename.get()
        if filename == "(všechny)":
            filename = None

        min_size = self.length_min.get()
        max_size = self.length_max.get()
        min_sites = self.filter_min_sites.get().strip()

        data = self.filterer.filter_data(enzyme, only_fragment, filename)
        if min_size:
            data = [row for row in data if int(row[5]) >= int(min_size)]
        if max_size:
            data = [row for row in data if int(row[5]) <= int(max_size)]
        if min_sites:
            data = [row for row in data if len(row[4].split(",")) >= int(min_sites)]

        self.filtered_data = data
        self.populate_tree(self.filtered_data)

    def calculate_data_ranges(self):
        data = self.db.get_all_samples()
        sizes = [int(row[5]) for row in data if str(row[5]).isdigit()]
        if sizes:
            self.min_seq_length = min(sizes)
            self.max_seq_length = max(sizes)
        else:
            self.min_seq_length = 0
            self.max_seq_length = 2000

    def open_delete_window(self):
        top = tk.Toplevel(self.root)
        top.title("Vyber soubory k odstranění")

        filenames = self.db.get_unique_filenames()

        select_all_var = tk.BooleanVar()
        select_fragment_var = tk.BooleanVar()

        def toggle_select_all():
            if select_all_var.get():
                listbox.select_set(0, tk.END)
                select_fragment_var.set(False)
            else:
                listbox.select_clear(0, tk.END)

        def toggle_select_fragments():
            listbox.select_clear(0, tk.END)
            if select_fragment_var.get():
                for i, fname in enumerate(filenames):
                    data = self.db.get_samples_by_filename(fname)
                    if any(row[6] for row in data):  # column 6 = isFragment
                        listbox.select_set(i)
                select_all_var.set(False)

        tk.Checkbutton(top, text="Vybrat vše", variable=select_all_var, command=toggle_select_all).pack(anchor="w", padx=10, pady=(10, 0))
        
        listbox = tk.Listbox(top, selectmode=tk.MULTIPLE, width=50)
        for name in filenames:
            listbox.insert(tk.END, name)
        listbox.pack(padx=10, pady=10)

        def confirm_delete():
            selected = [listbox.get(i) for i in listbox.curselection()]
            for fname in selected:
                self.db.delete_by_filename(fname)
            self.update_filenames()
            self.apply_filter()
            top.destroy()

        ttk.Button(top, text="Smazat vybrané", command=confirm_delete).pack(pady=5)
        ttk.Button(top, text="Zavřít", command=top.destroy).pack()

    def open_warning_window(self, text):

        pass

    def load_file(self):
        path = filedialog.askopenfilename(title="Vyber soubor", filetypes=[("Textové soubory", "*.txt")])
        if not path:
            return
        try:
            data = read_enzyme_headers(path)
            self.db.insert_data(data)
            messagebox.showinfo("Hotovo", "Data byla načtena.")
            self.update_filenames()
            self.populate_tree(self.db.get_all_samples())
        except Exception as e:
            messagebox.showerror("Chyba", str(e))

    def reset_database(self):
        if messagebox.askyesno("Potvrzení", "Opravdu chcete resetovat databázi?"):
            self.db.reset_database()
            self.update_filenames()
            self.tree.delete(*self.tree.get_children())

    def export_csv(self):
        self.exporter.export_to_csv(self.filtered_data)

    def update_filenames(self):
        filenames = self.db.get_unique_filenames()
        filenames.insert(0, "(všechny)")
        self.filename_combo["values"] = filenames
        self.filename_combo.set("(všechny)")
        if hasattr(self, "clustering_tab"):
            self.clustering_tab.refresh_lists()

    def populate_tree(self, data):
        self.tree.delete(*self.tree.get_children())
        for index, row in enumerate(data):
            tag = "evenrow" if index % 2 == 0 else "oddrow"
            self.tree.insert("", "end", values=row, tags=(tag,))

    def sort_column(self, col, reverse):
        l = [(self.tree.set(k, col), k) for k in self.tree.get_children("")]
        try:
            l.sort(key=lambda t: int(t[0]), reverse=reverse)
        except ValueError:
            l.sort(key=lambda t: t[0], reverse=reverse)
        for index, (val, k) in enumerate(l):
            self.tree.move(k, "", index)
        self.tree.heading(col, command=lambda: self.sort_column(col, not reverse))


# -------------------------------------------------------------------------------------
# Methods to work with Groups
    def _add_selected_to_group(self):
        selection = self.tree.selection()
        if not selection:
            messagebox.showwarning("Žádný výběr", "Nejprve vyber enzymy v tabulce.")
            return

        if not hasattr(self, 'selected_group_name') or not self.selected_group_name:
            messagebox.showinfo("Není vybraná skupina", "Nejdříve klikni na skupinu vlevo.")
            return

        group_names = self.groups.get_all_group_names()
        if not group_names:
            messagebox.showinfo("Žádné skupiny", "Nejdříve vytvoř skupinu vlevo.")
            return

        for iid in selection:
            values = self.tree.item(iid, "values")
            enzyme, filename, line = values[1], values[7], int(values[8])
            sequence = fetch_sequence(filename, line)
            temp = list(values)
            temp.append(sequence)
            self.groups.add_to_group(self.selected_group_name, {tuple(temp)})

        self.groups.save_groups()
        self._refresh_groups()

    def _add_group_button_to_controls(self):
        btn = ttk.Button(self.control_frame, text="Přidat do skupiny", command=self._add_selected_to_group)
        btn.grid(row=0, column=5, padx=5)

    def _build_group_sidebar(self):
        sidebar_header = ttk.Frame(self.sidebar)
        sidebar_header.pack(fill="x", pady=5)

        ttk.Label(sidebar_header, text="Skupiny", font=("Arial", 10, "bold")).pack(side="left", padx=5)
        ttk.Button(sidebar_header, text="+", width=3, command=self._prompt_new_group).pack(side="right", padx=5)

        self.group_canvas = tk.Canvas(self.sidebar, borderwidth=0)
        self.group_frame = ttk.Frame(self.group_canvas)
        self.group_scroll = ttk.Scrollbar(self.sidebar, orient="vertical", command=self.group_canvas.yview)
        self.group_canvas.configure(yscrollcommand=self.group_scroll.set)

        self.group_scroll.pack(side="right", fill="y")
        self.group_canvas.pack(side="left", fill="both", expand=True)
        self.group_canvas.create_window((0, 0), window=self.group_frame, anchor="nw")

        self.group_frame.bind("<Configure>", lambda e: self.group_canvas.configure(scrollregion=self.group_canvas.bbox("all")))

        self.group_widgets = {}
        self._refresh_groups()

    def _prompt_new_group(self):
        name = simpledialog.askstring("Nová skupina", "Zadej název nové skupiny:")
        if name:
            if name in self.groups.get_all_group_names():
                messagebox.showerror("Chyba", "Skupina s tímto názvem již existuje.")
                return
            self.groups.add_to_group(name, set())
            self.groups.save_groups()
            self._refresh_groups()

    def _refresh_groups(self):
        for widget in self.group_frame.winfo_children():
            widget.destroy()

        for name in self.groups.get_all_group_names():
            self._create_collapsible_group(name)

    def _create_collapsible_group(self, group_name):
        def select_group():
            if hasattr(self, 'selected_group_label') and self.selected_group_label:
                try:
                    self.selected_group_label.config(background="SystemButtonFace", foreground="black")
                except tk.TclError:
                    # Label už neexistuje, tak ho ignorujeme
                    self.selected_group_label = None
            self.selected_group_name = group_name
            self.selected_group_label = name_label
            name_label.config(background="#00ff59", foreground="black")

        container = ttk.Frame(self.group_frame)
        container.pack(fill="x", pady=2, padx=5)

        header = ttk.Frame(container)
        header.pack(fill="x")

        toggle_btn = ttk.Button(header, text="⯆", width=2)
        toggle_btn.pack(side="left")

        name_label = tk.Label(header, text=group_name, font=("Arial", 9, "bold"), anchor="w")
        name_label.pack(side="left", fill="x", expand=True, padx=5)
        name_label.bind("<Button-1>", lambda e: select_group())

        rename_btn = ttk.Button(header, text="✎", width=2, command=lambda: self._rename_group(group_name))
        rename_btn.pack(side="right", padx=2)
        delete_btn = ttk.Button(header, text="🗑", width=2, command=lambda: self._delete_group(group_name))
        delete_btn.pack(side="right", padx=2)
        export_btn = ttk.Button(header, text="📤", width=2, command=lambda: self._export_group(group_name))
        export_btn.pack(side="right", padx=2)

        content_frame = ttk.Frame(container)
        content_frame.pack(fill="x")

        enzymes = self.groups.get_group_ids(group_name)
        headers = ["Enzym", "Sample", "ORF", "Sekvence"]
        table = ttk.Frame(content_frame)
        table.pack(fill="x", padx=10)

        for col, name in enumerate(headers):
            ttk.Label(
                table,
                text=name,
                font=("Arial", 9, "bold"),
                background="#f0f0f0",
                anchor="w",
                padding=3
            ).grid(row=0, column=col, sticky="nsew", padx=2, pady=1)

        for row_idx, enz in enumerate(enzymes, start=1):
            data = list(enz)[1:5] if isinstance(enz, tuple) and len(enz) >= 5 else [str(enz), "?", "?", "?"]
            for col_idx, value in enumerate(data):
                ttk.Label(
                    table,
                    text=value,
                    anchor="w",
                    padding=3
                ).grid(row=row_idx, column=col_idx, sticky="nsew", padx=2, pady=1)
            del_btn = ttk.Button(table, text="x", width=2, command=lambda e=enz, g=group_name: self._remove_enzyme_from_group(g, e))
            del_btn.grid(row=row_idx, column=len(headers), padx=2)

        def toggle():
            if content_frame.winfo_viewable():
                content_frame.pack_forget()
                toggle_btn.config(text="⯈")
            else:
                content_frame.pack(fill="x")
                toggle_btn.config(text="⯆")

        toggle_btn.config(command=toggle)

    def _remove_enzyme_from_group(self, group, enzyme):
        self.groups.remove_from_group(group, enzyme)
        self.groups.save_groups()
        self._refresh_groups()

    def _delete_group(self, group_name):
        if messagebox.askyesno("Smazat skupinu", f"Opravdu chceš smazat skupinu '{group_name}'?"):
            self.groups.remove_group(group_name)
            self.groups.save_groups()
            self._refresh_groups()

    def _export_group(self, group_name):
        data = self.groups.get_group_ids(group_name)
        if not data:
            messagebox.showinfo("Prázdná skupina", "Skupina neobsahuje žádné enzymy.")
            return

        file_path = asksaveasfilename(defaultextension=".txt", filetypes=[(f"DATA\\{group_name}", "*.txt")])
        if not file_path:
            return

        try:
            with open(file_path, "w", encoding="utf-8") as f:
                for item in data:
                    if isinstance(item, tuple) and len(item) >= 10:
                        enzyme_name = item[1]
                        sample = item[2]
                        orf = item[3]
                        rec_seq = item[4]
                        size = item[5]
                        fragment = item[6]
                        sequence = item[9]

                        label = f">{enzyme_name}"
                        if sample:
                            label = f"{sample}.{enzyme_name}"
                        if orf:
                            label = f"{label}ORF{orf}"

                        rec_seq = rec_seq or []
                        rec_text = ", ".join(rec_seq) if isinstance(rec_seq, (list, tuple)) else str(rec_seq)

                        label = f"{label}    {rec_text} {size} aa"
                        if fragment.lower() == "true":
                            label += " fragment"

                        formatted_seq = "\n".join([" ".join([sequence[i+j:i+j+10] for j in range(0, 50, 10)]) for i in range(0, len(sequence), 50)])
                        f.write(f"{label}\n{formatted_seq}\n\n\n")
                    else:
                        f.write(f">{str(item)}\n?\n")

            messagebox.showinfo("Hotovo", f"Skupina '{group_name}' byla exportována.")
        except Exception as e:
            messagebox.showerror("Chyba při exportu", str(e))

    def _rename_group(self, old_name):
        new_name = simpledialog.askstring("Přejmenovat skupinu", "Nový název:", initialvalue=old_name)
        if new_name and new_name != old_name:
            if new_name in self.groups.get_all_group_names():
                messagebox.showerror("Chyba", "Skupina s tímto názvem už existuje.")
                return
            members = self.groups.get_group_ids(old_name)
            self.groups.remove_group(old_name)
            self.groups.add_to_group(new_name, members)
            self.groups.save_groups()
            self._refresh_groups()

























