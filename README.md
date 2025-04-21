# Nástroj pro analýzu dat z databáze restrikčních enzymů

Tento nástroj umožňuje analýzu, vizualizaci, porovnávání a predikci restrikčních enzymů na základě dat z databáze REBASE. Je napsán v jazyce Python a využívá grafické rozhraní postavené na knihovně `tkinter`.

## Funkce aplikace

- Načítání a správa databáze enzymů (import z FASTA souborů)
- Filtrování, třídění a zobrazení dat v přehledné tabulce
- Tvorba uživatelských skupin enzymů
- Export dat do CSV
- Porovnávání enzymů pomocí Needleman-Wunsch a Smith-Waterman algoritmů (záložka **Zarovnání**)
- Shlukování enzymů dle zvolených vlastností (záložka **Shlukování**)
- Vizualizace pomocí PCA a dalších grafů
- Podpora různých vlastností: skóre zarovnání, distribuce AK, délka, fragment, fyzikální vlastnosti atd.

## Instalace

1. **Nainstaluj Python 3.10+**
2. **Nainstaluj požadované knihovny:**

```bash
pip install -r requirements.txt
```

3. **Spusť aplikaci:**

```bash
python main.py
```

## 📁 Struktura projektu

```
.
├── main.py                          # Spouštěcí bod aplikace
├── gui_app.py                      # Hlavní GUI aplikace
├── alignment_tab.py                # Záložka Zarovnání
├── clustering_tab.py               # Záložka Shlukování
├── ProcessSimilarity.py            # Zarovnávací algoritmy
├── processProperties.py            # Výpočet vlastností a shlukování
├── ProcessFiles.py                 # Načítání a zpracování enzymových sekvencí
├── database_manager.py             # Vrstva nad databází
├── Database.py                     # Práce s SQLite databází
├── exporter.py                     # Export dat do CSV
├── group_manager.py                # Správa uživatelských skupin
├── groups.json                     # Uložené skupiny enzymů
├── requirements.txt                # Seznam knihoven
├── DATA/                           # Složka pro vstupní FASTA soubory
├── DATABASE/                       # Složka s .db souborem a SQL schématem
├── README.md                       # Tento soubor
└── ThesisSpecification_PHA0051.pdf # Oficiální zadání bakalářské práce
```

## Databáze

Databáze `R_Enzime.db` obsahuje:
- `samples` – informace o jednotlivých enzymech
- `alignments` – výsledky zarovnání mezi enzymy

## Odkazy na literaturu

- Pingoud, A. *Restriction Endonucleases*, Springer, 2012  
- Gusfield, D. *Algorithms on Strings, Trees, and Sequences*, Cambridge Univ Press, 1997  
- Chollet, F. *Deep Learning with Python*, Manning, 2021  
- Tramontano, A. *Introduction to Bioinformatics*, Chapman & Hall/CRC, 2018  

##  Autor

- **Student:** Khai Dat Phan  
- **Vedoucí:** Ing. Michal Vašinek, Ph.D.  
- **Fakulta:** FEI VŠB-TUO  
- **Termín odevzdání:** 30.04.2025