from gui_app import App
from ttkthemes import ThemedTk
from database_manager import DatabaseManager
from group_manager import GroupManager
from exporter import DataExporter
from filter_manager import FilterManager

if __name__ == "__main__":
    root = ThemedTk(theme="arc")

    # Injektuj logiku tříd
    db_manager = DatabaseManager()
    group_manager = GroupManager()
    exporter = DataExporter()
    filter_manager = FilterManager(db_manager)

    # Inicializace aplikace
    app = App(root, db_manager, group_manager, exporter, filter_manager)
    root.geometry("1250x700")
    root.mainloop()