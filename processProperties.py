import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_distances
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import silhouette_samples, silhouette_score
import matplotlib.pyplot as plt
from collections import Counter
from Bio.SeqUtils.ProtParam import ProteinAnalysis



class ProcessProperties:
    def __init__(self):
        self.features = {
            "alignment_score": False,       # data se získají ze záložky Zarovnání
            "distribution": False,          # protein_distribution
            "length": False,                # délka sekvence enzymu
            "is_fragment": False,           # pokud je daný enzym fragment

            # extract_physical_properties
            "molecular_weight": False,   
            "isoelectric_point": False,   
            "gravy": False,  
            "instability_index": False,   
            "aromaticity": False,
            "aliphatic_index": False,    
        }

    # Metoda pro nastavení které vlastnosti se využijí pro shlukování
    def choose_features(self, features: dict):
        for key in self.features:
            if key in features:
                self.features[key] = features[key]

    # Výpošet distribuce aminokyselin
    def protein_distribution(self, protein_sequence):
        AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
        total = len(protein_sequence)
        return [protein_sequence.count(aa) / total if total > 0 else 0 for aa in AMINO_ACIDS]


    # PCA (Principal Component Analysis):
    # Redukuje dimenzionalitu dat na 2D prostor při zachování co největší možné variability.
    # Výsledné souřadnice se používají pro vizualizaci shluků.
    def run_pca(self, vectors, n_components=2):
        pca = PCA(n_components=n_components)
        transformed = pca.fit_transform(vectors)
        return transformed, pca.explained_variance_ratio_


    # Kosinová vzdálenost:
    # Měří úhel mezi dvěma vektory – čím menší úhel, tím jsou si vektory (a tedy enzymy) podobnější.
    # Používá se pro vyhodnocení podobnosti ve vektorovém prostoru.
    def compute_cosine_distance_matrix(self, vectors):
        return cosine_distances(vectors)


    def aliphatic_index(self, sequence):
        sequence = sequence.upper()
        length = len(sequence)
        if length == 0:
            return 0.0

        ala = sequence.count("A") / length * 100
        val = sequence.count("V") / length * 100
        ile = sequence.count("I") / length * 100
        leu = sequence.count("L") / length * 100

        #a = 2.9, b = 3.9 jsou empirické koeficienty viz https://www.cymobase.org/cymobase/help_popup?for_window=true&topic=aliphatic_index/#:~:text=where%20X(Ala)%2C%20X,the%20side%20chain%20of%20alanine.
        return ala + 2.9 * val + 3.9 * (ile + leu) 


    # Fizikální vlastnosti získané z knihovny biopython
    def extract_physical_properties(self, sequence):
        analysis = ProteinAnalysis(sequence)
        props = []
        # Součet hmotností všech aminokyselin v proteinové sekvenci; Udává velikost (hmotnost) celého enzymu
        if self.features["molecular_weight"]:
            props.append(analysis.molecular_weight())

        # pH, při kterém má protein nulový celkový náboj; Důležitý pro chování enzymu v různých pH prostředích
        if self.features["isoelectric_point"]:
            props.append(analysis.isoelectric_point())

        # Průměrná hydropatie; Záporné hodnoty znamenají, že protein je spíš hydrofilní, kladné že je hydrofobní
        if self.features["gravy"]:
            props.append(analysis.gravy())

        # Empirický index předpovídající, zda je protein in vitro stabilní (>40 = nestabilní, <40 = stabilní); Ukazuje, jak dlouho enzym vydrží
        if self.features["instability_index"]:
            props.append(analysis.instability_index())

        # Podíl aromatických aminokyselin (Phe, Trp, Tyr)
        if self.features["aromaticity"]:
            props.append(analysis.aromaticity())

        # Poměr alifatických AK (Ala, Val, Ile, Leu) – tedy necyklických, hydrofobních zbytků
        if self.features["aliphatic_index"]:
            props.append(self.aliphatic_index(sequence))
        return props
    
        
    # Metoda pro generaci vektoru z tabulky v záložce Zarovnání
    def generate_vector_matrix_from_alignment(self, alignment_results):
        result = []
        for nameA, nameB, score, seqA, seqB, isFragA, isFragB in alignment_results:
            alignment_vec = [score] if self.features["alignment_score"] else []
            vecA = self.enrich_vector(seqA, alignment_vec, isFragA)
            vecB = self.enrich_vector(seqB, alignment_vec, isFragB)
            combined = vecA + vecB
            result.append((nameA, nameB, combined))
        return result


    # metoda pro získání statistik o vypracovaném výslledku shlukovacího algiritmu 
    def get_cluster_statistics(labels, vectors, pca_variance):
        counts = Counter(labels)
        silhouette = silhouette_score(vectors, labels)
        total_explained = sum(pca_variance[:2])
        return {
            "counts": dict(counts),                 # počet enzymů v každém shluku
            "silhouette_score": silhouette,         # průměrné silhouette skóre
            "pca_variance": pca_variance,           # Variance jednotlivých PCA složek
            "pca_explained_total": total_explained  # Součet prvních dvou složek – kolik procent variability zachycuje 2D projekce
        }




    #---------------------------------------------------------------------------------------------------------------
    # Vector pro jeden enzym
    def enrich_vector(self, protein_seq, alignment_vec, is_fragment):
        vector = []
        if self.features["alignment_score"]:
            vector += alignment_vec
        if self.features["length"]:
            vector.append(len(protein_seq))
        if self.features["distribution"]:
            vector += self.protein_distribution(protein_seq)
        vector += self.extract_physical_properties(protein_seq)
        if self.features["is_fragment"]:
            vector.append(1 if is_fragment else 0)
        return vector


    # Tato metoda provádí základní zpracování vlastností enzymů pro účely analýzy.
    # 1. Vstupem je seznam dat, kde každý prvek obsahuje jména dvou enzymů a jejich spojený vektor vlastností.
    # 2. Metoda provede PCA redukci vektorů na 2D (pro pozdější vizualizaci).
    # 3. Vypočítá kosinovou vzdálenost mezi všemi vektory (matice podobnosti).
    # Výstupem je slovník se jmény, PCA výsledky a vzdálenostmi – může se použít např. pro porovnání před shlukováním.
    def process_properties(self, data):
        # Očekává data ve formátu: [(nameA, nameB, [v1, v2, v3, v4, v5]), ...]
        vectors = [entry[2] for entry in data]
        names = [(entry[0], entry[1]) for entry in data]

        pca_result, explained_variance = self.run_pca(vectors)
        distance_matrix = self.compute_cosine_distance_matrix(vectors)

        return {
            "names": names,
            "vectors": vectors,
            "pca": pca_result,
            "variance": explained_variance,
            "cosine_distances": distance_matrix,
        }


    # Tato metoda kombinuje výpočet shluků a PCA projekci:
    # 1. Vstupem je seznam enzymových dvojic a jejich vektorů vlastností.
    # 2. Metoda podle volby 'method' provede buď KMeans nebo hierarchické shlukování.
    # 3. Výsledkem jsou:
    #    - jména porovnávaných enzymů
    #    - přiřazení každého vzorku do shluku
    #    - PCA souřadnice (2D) pro vizualizaci
    #    - původní vektory vlastností (např. pro další analýzy)
    def prepare_and_cluster(self, data_rows, n_clusters=5, method="kmeans"):
        vectors = [entry[2] for entry in data_rows]
        names = [(entry[0], entry[1]) for entry in data_rows]

        if method == "kmeans":
            model = KMeans(n_clusters=n_clusters, random_state=42)
        elif method == "agglomerative":
            model = AgglomerativeClustering(n_clusters=n_clusters)
        else:
            raise ValueError("Nepodporovaný algoritmus shlukování.")

        labels = model.fit_predict(vectors)
        pca = PCA(n_components=2)
        coords = pca.fit_transform(vectors)

        return names, labels, coords, vectors




    #---------------------------------------------------------------------------------------------------------------
    def plot_clusters_clean(coords, labels, title="Shlukování enzymů (čisté)"):
        fig, ax = plt.subplots(figsize=(10, 8))
        for cluster_id in np.unique(labels):
            cluster_coords = coords[labels == cluster_id]
            ax.scatter(
                cluster_coords[:, 0],
                cluster_coords[:, 1],
                label=f"Shluk {cluster_id}",
                alpha=0.7
            )
        ax.set_title(title)
        ax.set_xlabel("PCA 1")
        ax.set_ylabel("PCA 2")
        ax.legend()
        ax.grid(True)
        fig.tight_layout()
        return fig


    # -------------------------------------
    def plot_silhouette(vectors, labels):
        silhouette_vals = silhouette_samples(vectors, labels)
        avg_score = silhouette_score(vectors, labels)

        fig, ax = plt.subplots(figsize=(10, 6))
        y_lower = 10
        for i in np.unique(labels):
            cluster_vals = silhouette_vals[labels == i]
            cluster_vals.sort()
            y_upper = y_lower + len(cluster_vals)
            ax.fill_betweenx(np.arange(y_lower, y_upper), 0, cluster_vals, label=f'Shluk {i}')
            y_lower = y_upper + 10

        ax.axvline(avg_score, color="red", linestyle="--", label=f"Průměrné skóre: {avg_score:.2f}")
        ax.set_xlabel("Silhouette hodnota")
        ax.set_ylabel("Vzorky")
        ax.set_title("Silhouette analýza shluků")
        ax.legend()
        fig.tight_layout()
        return fig



    # -------------------------------------
    def plot_length_histogram_by_cluster(vectors, labels):
        lengths = [v[0] for v in vectors]  # Délka je na pozici 0

        fig, ax = plt.subplots(figsize=(10, 6))
        for cluster_id in np.unique(labels):
            cluster_lengths = [lengths[i] for i in range(len(labels)) if labels[i] == cluster_id]
            ax.hist(cluster_lengths, bins=15, alpha=0.6, label=f'Shluk {cluster_id}')
        ax.set_xlabel("Délka sekvence")
        ax.set_ylabel("Počet enzymů")
        ax.set_title("Rozložení délek sekvencí podle shluku")
        ax.legend()
        fig.tight_layout()
        return fig



