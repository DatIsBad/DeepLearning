# Nástroj pro analýzu dat z databáze restrikčních enzymů

Tato aplikace slouží k analýze a predikci restrikčních enzymů na základě dat z databáze REBASE.

## Funkce
- **Správa databáze**: import, reset, mazání a filtrování dat enzymů.
- **Zarovnání sekvencí**: Needleman-Wunsch a Smith-Waterman algoritmy.
- **Shlukování enzymů**: pomocí vlastností sekvencí a metod KMeans / Agglomerative clusteringu.
- **Predikce rozpoznávacích míst**: pomocí trénovaného Transformer-based modelu (BertLight).

## Struktura projektu
- `App/` – hlavní adresář aplikace obsahující všechny moduly:
  - `gui_app.py` – hlavní grafické rozhraní (GUI).
  - `alignment_tab.py` – záložka zarovnání sekvencí.
  - `clustering_tab.py` – záložka shlukování enzymů.
  - `prediction_tab.py` – záložka predikce restrikčních míst.
  - `processPrediction.py` – logika trénování a predikce modelu.
  - `ProcessSimilarity.py` – algoritmy zarovnání sekvencí.
  - `processProperties.py` – extrakce vlastností sekvencí a shlukování.
  - `database_manager.py` a `Database.py` – práce s SQLite databází.
  - `group_manager.py`, `filter_manager.py`, `exporter.py` – podpora pro skupiny, filtry a export.
  - `ProcessFiles.py` – nástroje pro práci se sekvencemi.
- `MODEL/` – složka s uloženými trénovanými modely (každý model má vlastní podsložku).
- `DATABASE/` – databáze enzymů (`R_Enzime.db`).
- `DATA/` – složka pro vstupní soubory se sekvencemi enzymů.

## Instalace
1. Nainstalujte požadované knihovny:
    ```
    pip install -r requirements.txt
    ```
2. Ujistěte se, že máte správně připravenou složku `DATABASE` s databází.

3. Spusťte hlavní soubor:
    ```
    python main.py
    ```

## Požadavky
- Python 3.9 nebo vyšší
- GUI rozhraní je vytvořeno pomocí `tkinter`
- Využívá knihovny `torch`, `sklearn`, `biopython`, `matplotlib` a další (viz requirements.txt)

## Poznámky
- Modely jsou uloženy v podsložkách složky `MODEL`, v každé podsložce je:
  - `model.pth`
  - `label_to_motif.json`
  - `tokenizer_vocab.json`
- Predikce nemusí být přesná kvůli omezenému množství trénovacích dat a vysoké unikátnosti sekvencí enzymů.

## Autor
Khai Dat Phan
