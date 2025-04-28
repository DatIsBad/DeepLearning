import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, random_split, Subset
from sklearn.metrics import classification_report
from sklearn.model_selection import StratifiedShuffleSplit
import numpy as np
from typing import List, Tuple
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter, defaultdict
from database_manager import DatabaseManager
from ProcessFiles import fetch_sequence
import random
import json
from pathlib import Path


import os

# Třída pro zpracování predikcí restrikčních míst z proteinových sekvencí
# Obsahuje tokenizaci, dataset, model, trénování, evaluaci a predikci
class ProcessPrediction:
    # Inicializuje správce databáze, cestu k modelu a připraví kontejnerové proměnné
    def __init__(self, db: DatabaseManager, model_path: str = "MODEL\\dump\\my_model.pth"):
        self.db = db
        self.model_path = model_path
        self.tokenizer = None
        self.dataset = None
        self.label_to_motif = None
        self.model = None
      
    # Tokenizace proteinových sekvencí pomocí k-merů (např. 3-mer)
    class KmerTokenizer:
        def __init__(self, k=3):
            self.k = k
            self.vocab = {"<PAD>": 0}

                # Vytvoří slovník nejčastějších k-merů z trénovacích sekvencí
        def build_vocab(self, sequences: List[str], max_vocab_size=5000):
            kmer_counts = defaultdict(int)
            for seq in sequences:
                for i in range(len(seq) - self.k + 1):
                    kmer = seq[i:i+self.k]
                    kmer_counts[kmer] += 1

            sorted_kmers = sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)
            for idx, (kmer, _) in enumerate(sorted_kmers[:max_vocab_size-1], start=1):
                self.vocab[kmer] = idx

            print(f"[Tokenizer] Vocab size: {len(self.vocab)}")

                # Zakóduje vstupní sekvenci na číselný vektor pomocí k-mer tokenizace
        def encode(self, seq: str, max_len: int) -> List[int]:
            tokens = [seq[i:i+self.k] for i in range(len(seq) - self.k + 1)]
            encoded = [self.vocab.get(tok, 0) for tok in tokens[:max_len]]
            return encoded + [0] * (max_len - len(encoded))
        
        
    # Dataset přizpůsobený pro PyTorch - uchovává zakódované sekvence a cílové motivy
    class ProteinMotifDataset(Dataset):
        def __init__(self, sequences: List[str], motifs: List[str], tokenizer, max_len: int = 128):
            self.tokenizer = tokenizer
            self.max_len = max_len
            self.motif_to_label = {motif: idx for idx, motif in enumerate(sorted(set(motifs)))}
            self.label_to_motif = {idx: motif for motif, idx in self.motif_to_label.items()}

            self.X_seq = [self.tokenizer.encode(seq, self.max_len) for seq in sequences]
            self.y = [self.motif_to_label[m] for m in motifs]

        # Vrací délku datasetu
        def __len__(self):
            return len(self.y)

        # Vrací jeden vzorek (sekvence, štítek) pro daný index
        def __getitem__(self, idx):
            return torch.tensor(self.X_seq[idx], dtype=torch.long), torch.tensor(self.y[idx], dtype=torch.long)

    # Transformer model pro klasifikaci motivů z proteinových sekvencí
    class BertLight(nn.Module):
        def __init__(self, vocab_size=1024, embed_dim=128, num_heads=4, hidden_dim=256, num_classes=10, max_len=128):
            super().__init__()
            self.embedding = nn.Embedding(vocab_size, embed_dim, padding_idx=0)
            nn.init.xavier_uniform_(self.embedding.weight)
            self.layer_norm = nn.LayerNorm(embed_dim)
            self.pos_embedding = nn.Parameter(torch.randn(1, max_len, embed_dim))
            encoder_layer = nn.TransformerEncoderLayer(d_model=embed_dim, nhead=num_heads, dim_feedforward=hidden_dim, batch_first=True)
            self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=2)
            self.classifier = nn.Sequential(
                nn.Linear(embed_dim, 128),
                nn.ReLU(),
                nn.Linear(128, num_classes)
            )

        # Definuje průchod dat modelem (embedding → transformer → klasifikace)
        def forward(self, x):
            x = self.embedding(x) + self.pos_embedding[:, :x.size(1), :]
            x = self.layer_norm(x)
            x = self.transformer(x)
            x = x.mean(dim=1)
            return self.classifier(x)


    # Trénuje model na trénovacích datech a sleduje výkon na validačních datech
    # Argumenty:
    #   epochs – počet epoch
    #   batch_size – velikost batchů
    #   lr – learning rate
    #   patience – počet epoch bez zlepšení pro early stopping
    def train_model(self, epochs=50, batch_size=32, lr=0.0005, patience=10):
        train_loader = DataLoader(self.train_dataset, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(self.val_dataset, batch_size=batch_size)

        optimizer = torch.optim.Adam(self.model.parameters(), lr=lr)
        criterion = nn.CrossEntropyLoss()

        best_val_acc = 0.0
        patience_counter = 0
        label_counts = Counter([label for _, label in self.train_dataset])
        print("Rozložení tříd v trénovacích datech:")
        for label, count in label_counts.items():
            print(f"  Třída {label}: {count} vzorků")

            # Přepneme zpět do trénovacího režimu
            self.model.train()
        for epoch in range(epochs):
            total_loss = 0
            correct = 0
            total = 0


            # Procházíme jednotlivé batch-e trénovacích dat
            for batch_idx, (seqs, labels) in enumerate(train_loader):
                if batch_idx % 10 == 0:
                    print(f"  Batch {batch_idx+1}/{len(train_loader)}")
                # Nuluje gradienty před zpětnou propagací
                optimizer.zero_grad()
                outputs = self.model(seqs)
                # Vypočítá ztrátu mezi predikcemi a reálnými štítky
                loss = criterion(outputs, labels)
                # Zpětná propagace chyby
                loss.backward()
                # Aktualizace vah modelu
                optimizer.step()

                total_loss += loss.item()
                pred = outputs.argmax(dim=1)
                correct += (pred == labels).sum().item()
                total += len(labels)

            acc = correct / total * 100
            print(f"Epoch {epoch+1}/{epochs} | Train Loss: {total_loss:.4f} | Train Acc: {acc:.2f}%")

            # Přepneme model do evaluačního režimu
            self.model.eval()
            all_preds = []
            all_labels = []
            with torch.no_grad():
                for seqs, labels in val_loader:
                    outputs = self.model(seqs)
                    preds = outputs.argmax(dim=1)
                    all_preds.extend(preds.tolist())
                    all_labels.extend(labels.tolist())

            report = classification_report(all_labels, all_preds, zero_division=0, output_dict=True)
            val_acc = report['accuracy'] * 100
            print(f"Validation accuracy: {val_acc:.2f}%")
            if val_acc > best_val_acc:
                best_val_acc = val_acc
                patience_counter = 0
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    print("Early stopping triggered.")
                    break
            self.model.train()


    # Načte sekvence a motivy z databáze, předzpracuje je a vytvoří tokenizer a dataset
    # Argumenty:
    #   filenames – volitelný seznam názvů souborů
    #   min_len – minimální délka sekvence
    #   only_full – True = vynechá fragmenty
    #   top_n – počet nejčastějších motivů, které ponecháme
    def extract_data_from_db(self, filenames=None, min_len=50, only_full=False, top_n=30):
        data = []
        if filenames:
            for file in filenames:
                data.extend(self.db.get_samples_by_filename(file))

        else:
            data = self.db.get_all_samples()

        records = []
        for row in data:
            motif = row[4]
            if not motif or len(motif) < 3:
                continue
            seq = fetch_sequence(row[7], row[8])
            if not seq or len(seq) < min_len:
                continue
            if only_full and row[6] == 1:
                continue
            records.append((seq, motif))

        

        motifs = [m for _, m in records]
        common = set([m for m, _ in Counter(motifs).most_common(top_n)])
        records = [(s, m) for s, m in records if m in common]

        # Filtrujeme motivy s alespoň 2 výskyty
        motif_counts = Counter(m for s, m in records)
        records = [(s, m) for s, m in records if motif_counts[m] >= 2]

        
        if not records:
            raise ValueError("Nebyly nalezeny dostatečné motivy pro trénink.")


        sequences = [s for s, m in records]
        motifs = [m for s, m in records]


        self.tokenizer = ProcessPrediction.KmerTokenizer(k=3)
        self.tokenizer.build_vocab(sequences, max_vocab_size=1024)
        self.dataset = ProcessPrediction.ProteinMotifDataset(sequences, motifs, self.tokenizer, max_len=128)
        self.label_to_motif = self.dataset.label_to_motif

        return sequences, motifs


    # Načte model z uloženého .pth souboru do self.model
    # Načte model z disku a inicializuje architekturu na základě tokenizeru a počtu tříd
    def load_model(self):
        model_dir = Path(self.model_path).parent

        if not (model_dir / "model.pth").exists():
            raise FileNotFoundError(f"Model {self.model_path} neexistuje.")

        # Načti tokenizer vocab
        try:
            with open(model_dir / "tokenizer_vocab.json", "r", encoding="utf-8") as f:
                vocab = json.load(f)
            self.tokenizer = ProcessPrediction.KmerTokenizer(k=3)
            self.tokenizer.vocab = vocab
            print("Tokenizer vocab načten.")
        except Exception as e:
            raise ValueError(f"Chyba při načítání tokenizer_vocab.json: {e}")

        # Načti label_to_motif
        try:
            with open(model_dir / "label_to_motif.json", "r", encoding="utf-8") as f:
                self.label_to_motif = json.load(f)
            self.label_to_motif = {int(k): v for k, v in self.label_to_motif.items()}
            print("Label mapa načtena.")
        except Exception as e:
            raise ValueError(f"Chyba při načítání label_to_motif.json: {e}")

        # Inicializuj model
        self.model = ProcessPrediction.BertLight(
            vocab_size=len(self.tokenizer.vocab),
            num_classes=len(self.label_to_motif),
            max_len=128
        )
        self.model.load_state_dict(torch.load(model_dir / "model.pth"))
        print(f"Model načten ze souboru '{self.model_path}'.")





    # Uloží model a mapu štítků (label_to_motif) do souboru
    # Uloží model a label-to-motif mapu na disk, pro pozdější predikci nebo nasazení
    def save_model(self, model_name: str):
        save_path = Path("MODEL") / model_name
        save_path.mkdir(parents=True, exist_ok=True)
        torch.save(self.model.state_dict(), save_path / "model.pth")
        with open(save_path / "label_to_motif.json", "w", encoding="utf-8") as f:
            json.dump(self.label_to_motif, f, indent=4)
        with open(save_path / "tokenizer_vocab.json", "w", encoding="utf-8") as f:
            json.dump(self.tokenizer.vocab, f, indent=4)
        print(f"Model uložen do {save_path}")




    # Rozdělí dataset na trénovací a validační část pomocí stratifikovaného dělení
    def split_dataset(self, test_ratio: float = 0.2, random_state: int = 42):
        if self.dataset is None:
            raise ValueError("Dataset nebyl inicializován. Nejprve zavolej extract_data_from_db().")

        labels = [label for _, label in self.dataset]
        splitter = StratifiedShuffleSplit(n_splits=1, test_size=test_ratio, random_state=random_state)
        train_idx, val_idx = next(splitter.split(np.zeros(len(labels)), labels))

        self.train_dataset = Subset(self.dataset, train_idx)
        self.val_dataset = Subset(self.dataset, val_idx)
        print(f"Dataset rozdělen: {len(train_idx)} trénovacích, {len(val_idx)} validačních vzorků.")


    # Vyhodnotí model na validačních datech a vypíše metriky
    # Spustí vyhodnocení modelu na validačních datech
    # Výsledkem je klasifikační report (přesnost, F1, atd.)
    def evaluate(self):
        if self.val_dataset is None:
            raise ValueError("Validační dataset není k dispozici. Zavolej nejprve split_dataset().")

        loader = DataLoader(self.val_dataset, batch_size=32)
        self.model.eval()
        all_preds = []
        all_labels = []

        with torch.no_grad():
            for seqs, labels in loader:
                outputs = self.model(seqs)
                preds = outputs.argmax(dim=1)
                all_preds.extend(preds.tolist())
                all_labels.extend(labels.tolist())

        report = classification_report(all_labels, all_preds, zero_division=0)
        print("Vyhodnocení na validačním datasetu:")
        print(report)


    # Vrátí predikovaný motiv na základě vstupní sekvence pomocí modelu
    # Argumenty:
    #   sequence – vstupní proteinová sekvence
    #   max_len – maximální délka sekvence po zakódování
    def predict(self, sequence: str, max_len: int = 128):
        if self.model is None:
            raise ValueError("Model není načtený.")
        self.model.eval()
        encoded = self.tokenizer.encode(sequence, max_len)
        tensor = torch.tensor(encoded, dtype=torch.long).unsqueeze(0)
        with torch.no_grad():
            output = self.model(tensor)
            pred_label = output.argmax(dim=1).item()
            return self.label_to_motif.get(pred_label, "UNKNOWN")
        


    def auto_train(self, model_name = "TestModel_2" , top_n = 30):
        sequences, motifs = self.extract_data_from_db(
        min_len=30,   # minimální délka sekvence
        top_n= top_n    # ponecháme top 100 nejčastějších motivů (upraveno dle požadavku)
        )
        print(f"✅ Načteno {len(sequences)} sekvencí.")

        # 4. Rozdělení na trénovací a validační sadu
        self.split_dataset(test_ratio=0.2)

        # 5. Vytvoření modelu
        self.model = ProcessPrediction.BertLight(
            vocab_size=len(self.tokenizer.vocab),
            num_classes=len(self.label_to_motif),
            max_len=128
        )

        # 6. Trénink modelu
        print("⏳ Trénuji model...")
        self.train_model(
            epochs=100,        # počet epoch
            batch_size=32,    # velikost batchů
            lr=0.0002,        # learning rate
            patience=20        # early stopping
        )

        # 7. Vyhodnocení modelu
        print("📈 Vyhodnocení modelu:")
        self.evaluate()

        # 8. Uložení modelu
        self.save_model(model_name)
