import tkinter as tk
from tkinter import ttk, messagebox
from tkinter import filedialog
from ProcessFiles import read_enzyme_headers
from alignment_tab import AlignmentTab
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

        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True)

        self.data_tab = ttk.Frame(self.notebook)
        self.build_gui(self.data_tab)
        self.notebook.add(self.data_tab, text="Data")

        self.alignment_tab = AlignmentTab(self.notebook, self.db, ProcessSimilarity())
        """self.alignment_tab.combo_file1['values'] = self.db.get_unique_filenames()
        self.alignment_tab.combo_file2['values'] = self.db.get_unique_filenames()
        all_data = self.db.get_all_samples()
        enzyme_names = sorted(set(row[1] for row in all_data))
        self.alignment_tab.combo_enzyme['values'] = enzyme_names
        self.alignment_tab.combo_enzyme1['values'] = enzyme_names
        self.alignment_tab.combo_enzyme2['values'] = enzyme_names"""
        self.notebook.add(self.alignment_tab.frame, text="Zarovnání")


        self.update_filenames()
        self.populate_tree(self.db.get_all_samples())

    def build_gui(self, parent):
        self.root.title("Restrikční enzymy – databázový nástroj")

        control_frame = ttk.LabelFrame(parent, text="Ovládání")
        control_frame.pack(fill="x", padx=5, pady=5)

        ttk.Button(control_frame, text="Vytvořit databázi", command=self.db.create_database).grid(row=0, column=0, padx=5)
        ttk.Button(control_frame, text="Načíst soubor", command=self.load_file).grid(row=0, column=1, padx=5)
        ttk.Button(control_frame, text="Reset databáze", command=self.reset_database).grid(row=0, column=2, padx=5)
        ttk.Button(control_frame, text="Export do CSV", command=self.export_csv).grid(row=0, column=3, padx=5)
        ttk.Button(control_frame, text="Smazat soubory", command=self.open_delete_window).grid(row=0, column=4, padx=5)

        filter_frame = ttk.LabelFrame(parent, text="Filtry")
        filter_frame.pack(fill="x", padx=10, pady=5)

        tree_frame = ttk.Frame(parent)
        tree_frame.pack(fill="both", expand=True, padx=10, pady=5)
        tree_frame.columnconfigure(0, weight=1)
        tree_frame.rowconfigure(0, weight=1)



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