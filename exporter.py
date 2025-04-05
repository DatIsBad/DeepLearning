import csv
from tkinter import filedialog, messagebox

class DataExporter:
    def export_to_csv(self, data):
        if not data:
            messagebox.showwarning("Žádná data", "Nejsou dostupná žádná data pro export.")
            return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV soubory", "*.csv")],
            title="Uložit jako..."
        )

        if not filepath:
            return

        try:
            with open(filepath, mode="w", newline="", encoding="utf-8") as file:
                writer = csv.writer(file)
                header = ["id", "enzyme", "sample", "orf", "rec_seq", "size", "fragment", "filename", "line"]
                writer.writerow(header)
                for row in data:
                    writer.writerow(row)

            messagebox.showinfo("Export dokončen", f"Data byla úspěšně exportována do {filepath}")

        except Exception as e:
            messagebox.showerror("Chyba při exportu", str(e))