from typing import List, Tuple, Dict  # typová kontrola
from ProcessFiles import fetch_sequence
from collections import Counter
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader



class ProcessPrediction:
    def __init__(self, db_manager):
        self.db = db_manager
        self.motif_to_label: Dict[str, int] = {}
        self.label_to_motif: Dict[int, str] = {}
        self.amino_to_idx = self._build_amino_index()

    # ------------------- DATASET -------------------

    def _build_amino_index(self) -> Dict[str, int]:
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # běžné aminokyseliny
        return {aa: idx + 1 for idx, aa in enumerate(amino_acids)}  # 0 = padding


    # Extrahuje dvojice (proteinová sekvence, DNA motiv) ze záznamů v databázi.
    def extract_protein_motif_pairs(self, filename:None) -> List[Tuple[str, str]]:
        pairs = []
        for row in self.db.get_samples_by_filename(filename):
            motif = row[4]  # rec_sequence – rozpoznávací DNA sekvence (motiv)
            if not motif or len(motif) < 3:
                continue

            file = row[7]
            line = row[8]

            try:
                protein_seq = fetch_sequence(file, line)
                if not protein_seq or len(protein_seq) < 30:
                    continue
                pairs.append((protein_seq, motif))
            except Exception as e:
                print(f"Chyba při načítání {file}:{line} - {e}")
                continue

        return pairs
    

    def prepare_classification_dataset(self, max_len: int = 500, top_n: int = 100) -> Tuple[np.ndarray, np.ndarray]:
        pairs = self.extract_protein_motif_pairs()

        # Sečti frekvence motivů
        counter = Counter(motif for _, motif in pairs)  # spočítá kolikrát se každý motiv v datech vyskytuje.
        most_common = counter.most_common(top_n)        # dostanu 100(defaunt) nejvíce objevovaných rozpoznávacích sekvencí
        selected_motifs = set(m for m, _ in most_common)

        # Mapování motiv <-> label
        self.motif_to_label = {motif: idx for idx, (motif, _) in enumerate(most_common)}
        self.label_to_motif = {idx: motif for motif, idx in self.motif_to_label.items()}

        X = []  #seznam proteinových sekvencí
        y = []

        for protein_seq, motif in pairs:
            if motif not in self.motif_to_label:
                continue

            # Sekvenci převedeme na indexy aminokyselin
            indices = [self.amino_to_idx.get(aa, 0) for aa in protein_seq[:max_len]]
            indices += [0] * (max_len - len(indices))  # doplnění nul

            X.append(indices)
            y.append(self.motif_to_label[motif])

        return (np.array(X), np.array(y))
    


    # ------------------- Pytorch -------------------
    
    def train_pytorch_model(
        self,
        X: np.ndarray,  # seznam proteinových sekvencí (zakódovaných jako čísla) z prepare_classification_dataset
        y: np.ndarray,  # seznam integer štítků odpovídajících motivům
        max_len: int,   # maximální délka vstupní sekvence (použito při paddingu)
        num_classes: int,  # počet možných výstupních tříd (motivů)
        epochs: int = 6,  # počet trénovacích epoch
        batch_size: int = 32  # velikost dávky při tréninku
    ):
        # Zjistí, zda je dostupná GPU akcelerace. Pokud ano, použije cuda, jinak cpu
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


        # -------------------------------
        # Dataset obalující vstupní data a labely (nutné pro DataLoader z Pytorch)
        class ProteinDataset(Dataset):
            def __init__(self, data, labels):
                self.data = torch.tensor(data, dtype=torch.long)        # vstupní sekvence
                self.labels = torch.tensor(labels, dtype=torch.long)    # odpovídající labely

            def __len__(self):
                return len(self.data)

            def __getitem__(self, idx):
                return self.data[idx], self.labels[idx]


        # Architektura modelu: Embedding → LSTM → Dense → Výstup
        class MotifClassifier(nn.Module):
            def __init__(self, vocab_size, embed_dim, hidden_dim, num_classes):
                super().__init__()

                # Vstupní vrstva: každé číslo (aminokyselina) se převede na vektor
                self.embedding = nn.Embedding(vocab_size, embed_dim, padding_idx=0)
                # Obousměrné LSTM – sleduje kontext sekvence z obou směrů
                self.lstm = nn.LSTM(embed_dim, hidden_dim, batch_first=True, bidirectional=True)
                # Skrytá plně propojená vrstva
                self.fc1 = nn.Linear(hidden_dim * 2, 128)
                # Výstupní vrstva: počet neuronů = počet tříd (motivů)
                self.fc2 = nn.Linear(128, num_classes)

            def forward(self, x):
                x = self.embedding(x)                   # 1) Embedding
                _, (h_n, _) = self.lstm(x)              # 2) LSTM – výstup = skrytý stav obou směrů
                h = torch.cat((h_n[0], h_n[1]), dim=1)  # 3) Spojení směrů do jednoho vektoru
                x = F.relu(self.fc1(h))                 # 4) Aktivace přes ReLU
                return self.fc2(x)                      # 5) Výstupní logity pro každou třídu
            

        # -------------------------------
        dataset = ProteinDataset(X, y)
        dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)   # Dávkovač dat pro mini-batch trénink

        model = MotifClassifier(vocab_size=len(self.amino_to_idx) + 1, embed_dim=64, hidden_dim=64, num_classes=num_classes).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)       # Optimalizátor a ztrátová funkce
        criterion = nn.CrossEntropyLoss()                               # pro vícetřídovou klasifikaci


        # Vstup do trénovacího režimu
        model.train()
        for epoch in range(epochs):
            total_loss = 0
            correct = 0

            # Mini-batch tréninková smyčka
            for batch_x, batch_y in dataloader:
                # Přesuň vstup a výstup na správné zařízení (CPU/GPU); device bylo nastaveno na začátku metody
                batch_x, batch_y = batch_x.to(device), batch_y.to(device)


                optimizer.zero_grad()               # Reset gradientů
                outputs = model(batch_x)            # Forward průchod modelem
                loss = criterion(outputs, batch_y)  # Výpočet ztráty
                loss.backward()                     # Zpětná propagace (výpočet gradientů)
                optimizer.step()                    # Aktualizace vah podle gradientů

                total_loss += loss.item()
                correct += (outputs.argmax(dim=1) == batch_y).sum().item()  # Počet správně klasifikovaných vzorků

            # Výpočet přesnosti za epochu
            acc = correct / len(dataset)
            print(f"Epoch {epoch+1}/{epochs} - Loss: {total_loss:.4f} - Accuracy: {acc:.4f}")

        return model



    # ------------------- Work -------------------
    # Predikuje motiv pro zadanou proteinovou sekvenci pomocí natrénovaného modelu
    def predict_motif(self, protein_sequence: str, model, max_len: int) -> str:
        model.eval()
        device = next(model.parameters()).device

        # Převeď sekvenci na indexy a doplň na max_len
        indices = [self.amino_to_idx.get(aa, 0) for aa in protein_sequence[:max_len]]
        indices += [0] * (max_len - len(indices))

        input_tensor = torch.tensor([indices], dtype=torch.long).to(device)
        with torch.no_grad():
            output = model(input_tensor)
            pred_label = output.argmax(dim=1).item()

        return self.label_to_motif.get(pred_label, "Neznámý")

    # Uloží PyTorch model na disk pod zadaným názvem souboru
    def save_model(self, model, path: str):
        torch.save(model.state_dict(), path)

    # Načte model ze souboru a vrátí ho připravený k použití
    def load_model(self, path: str, max_len: int, num_classes: int):
        class MotifClassifier(nn.Module):
            def __init__(self, vocab_size, embed_dim, hidden_dim, num_classes):
                super().__init__()
                self.embedding = nn.Embedding(vocab_size, embed_dim, padding_idx=0)
                self.lstm = nn.LSTM(embed_dim, hidden_dim, batch_first=True, bidirectional=True)
                self.fc1 = nn.Linear(hidden_dim * 2, 128)
                self.fc2 = nn.Linear(128, num_classes)

            def forward(self, x):
                x = self.embedding(x)
                _, (h_n, _) = self.lstm(x)
                h = torch.cat((h_n[0], h_n[1]), dim=1)
                x = F.relu(self.fc1(h))
                return self.fc2(x)

        model = MotifClassifier(vocab_size=len(self.amino_to_idx) + 1, embed_dim=64, hidden_dim=64, num_classes=num_classes)
        model.load_state_dict(torch.load(path, map_location=torch.device('cpu')))
        model.eval()
        return model




















