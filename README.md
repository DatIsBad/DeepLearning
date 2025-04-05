# Nástroj pro analýzu restrikčních enzymů

Tato aplikace slouží k načítání, filtrování, třídění a správě dat o restrikčních enzymech (např. z databáze REBASE).

## Funkce
- Vytváření a reset databáze
- Načítání dat ze souborů
- Filtrování dat podle enzymu, fragmentu, souboru
- Správa uživatelských skupin (přidávání, uložení, načtení)
- Export vyfiltrovaných dat do CSV
- Odstraňování dat podle souboru

## Spuštění
1. Nainstaluj požadavky:
   ```
   pip install -r requirements.txt
   ```

2. Spusť aplikaci:
   ```
   python main.py
