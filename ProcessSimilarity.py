import numpy as np
import matplotlib.pyplot as plt
from Bio.Align import substitution_matrices


class ProcessSimilarity:
    def __init__(self):
        # načtení BLOSUM62;
        # toto načtení je velice náročné
        # lepší načíst jen jednou a umožnit algoritmům použít BLOSUM62 bez čekání
        self.blosum62 = substitution_matrices.load("BLOSUM62")
        self.gap_penalty = -4

    def get_blosum_score(self, a, b, default=-4):
        matrix = self.blosum62
        try:
            return matrix[(a, b)]
        except KeyError:
            try:
                return matrix[(b, a)]
            except KeyError:
                return default

    # Needleman-Wunsh (global alignment)
    # return: score_matrix : [], pointer_matrix : []
    def needleman_wunsh_matrix(self, enzyme_1, enzyme_2):
        len_1 = len(enzyme_1) + 1
        len_2 = len(enzyme_2) + 1
        score_matrix = [[0] * len_2 for _ in range(len_1)]
        pointer_matrix = [[None] * len_2 for _ in range(len_1)]

        for i in range(len_1):
            score_matrix[i][0] = i * self.gap_penalty
            pointer_matrix[i][0] = (i - 1, 0) if i > 0 else None
        for j in range(len_2):
            score_matrix[0][j] = j * self.gap_penalty
            pointer_matrix[0][j] = (0, j - 1) if j > 0 else None

        for i in range(1, len_1):
            for j in range(1, len_2):
                match_score = self.get_blosum_score(enzyme_1[i - 1], enzyme_2[j - 1])
                diag = score_matrix[i - 1][j - 1] + match_score
                up = score_matrix[i - 1][j] + self.gap_penalty
                left = score_matrix[i][j - 1] + self.gap_penalty
                max_score = max(diag, up, left)

                score_matrix[i][j] = max_score
                if max_score == diag:
                    pointer_matrix[i][j] = (i - 1, j - 1)
                elif max_score == up:
                    pointer_matrix[i][j] = (i - 1, j)
                else:
                    pointer_matrix[i][j] = (i, j - 1)

        return (score_matrix, pointer_matrix)

    # Smith-Waterman (local alignment)
    # Return: score_matrix : [], pointer_matrix : [], max_score : int, max_positions : []
    def smith_waterman_matrix(self, enzyme_1, enzyme_2):
        len_1 = len(enzyme_1) + 1
        len_2 = len(enzyme_2) + 1
        score_matrix = [[0] * len_2 for _ in range(len_1)]
        pointer_matrix = [[None] * len_2 for _ in range(len_1)]

        max_score = 0
        max_positions = []  # ← změna: list všech maxim

        for i in range(1, len_1):
            for j in range(1, len_2):
                match_score = self.get_blosum_score(enzyme_1[i - 1], enzyme_2[j - 1])
                diag = score_matrix[i - 1][j - 1] + match_score
                up = score_matrix[i - 1][j] + self.gap_penalty
                left = score_matrix[i][j - 1] + self.gap_penalty
                max_cell = max(0, diag, up, left)

                score_matrix[i][j] = max_cell
                if max_cell == 0:
                    pointer_matrix[i][j] = None
                elif max_cell == diag:
                    pointer_matrix[i][j] = (i - 1, j - 1)
                elif max_cell == up:
                    pointer_matrix[i][j] = (i - 1, j)
                else:
                    pointer_matrix[i][j] = (i, j - 1)

                if max_cell > max_score:
                    max_score = max_cell
                    max_positions = [(i, j)]  # ← nový nejlepší bod, resetujeme seznam
                elif max_cell == max_score:
                    max_positions.append((i, j))  # ← další bod se stejným skóre

        return (score_matrix, pointer_matrix, max_score, max_positions)


    # ----------------------------------------------------------------------------------------------------
    # compute scores of enzyme_1 and enzyme_2 with needleman method
    # return: (aligned_seq1, aligned_seq2, score_matrix[M][N])
    def needleman_wunsh_alignment(self, enzyme_1, enzyme_2):
        enzyme_1 = enzyme_1.replace(" ", "").replace("\n", "")
        enzyme_2 = enzyme_2.replace(" ", "").replace("\n", "")

        M, N = len(enzyme_1), len(enzyme_2)
        score_matrix, pointer_matrix = self.needleman_wunsh_matrix(enzyme_1, enzyme_2)
        aligned_seq1, aligned_seq2 = "", ""
        i, j = M, N

        # Calculate alignment
        while pointer_matrix[i][j] is not None:
            prev_i, prev_j = pointer_matrix[i][j]
            if prev_i == i - 1 and prev_j == j - 1:
                aligned_seq1 = enzyme_1[i - 1] + aligned_seq1
                aligned_seq2 = enzyme_2[j - 1] + aligned_seq2
            elif prev_i == i - 1:
                aligned_seq1 = enzyme_1[i - 1] + aligned_seq1
                aligned_seq2 = "-" + aligned_seq2
            elif prev_j == j - 1:
                aligned_seq1 = "-" + aligned_seq1
                aligned_seq2 = enzyme_2[j - 1] + aligned_seq2
            i, j = prev_i, prev_j

        # scores were saved from the end to beginning, thus list need to be reversed

        return (aligned_seq1, aligned_seq2, score_matrix[M][N])

    # compute scores of enzyme_1 and enzyme_2 with smith method
    # return: [(aligned_seq1, aligned_seq2, max_score, max_score, path, coords)] vypočer kosinovy vzdalenosti; PCA - vizualizace ; PAMBOSOM - podobnost; lstm (6); keras
    def smith_waterman_alignment(self, enzyme_1, enzyme_2, all_maxima=True):
        enzyme_1 = enzyme_1.replace(" ", "").replace("\n", "")
        enzyme_2 = enzyme_2.replace(" ", "").replace("\n", "")

        score_matrix, pointer_matrix, max_score, max_positions = (
            self.smith_waterman_matrix(enzyme_1, enzyme_2)
        )

        def traceback(start_i, start_j):
            aligned_seq1, aligned_seq2 = "", ""
            i, j = start_i, start_j
            traceback_path = []
            traceback_coords = [(i, j)]
            while pointer_matrix[i][j] is not None:
                prev_i, prev_j = pointer_matrix[i][j]
                traceback_coords.append((prev_i, prev_j))
                if prev_i == i - 1 and prev_j == j - 1:
                    aligned_seq1 = enzyme_1[i - 1] + aligned_seq1
                    aligned_seq2 = enzyme_2[j - 1] + aligned_seq2
                    traceback_path.append("diag")
                elif prev_i == i - 1:
                    aligned_seq1 = enzyme_1[i - 1] + aligned_seq1
                    aligned_seq2 = "-" + aligned_seq2
                    traceback_path.append("up")
                elif prev_j == j - 1:
                    aligned_seq1 = "-" + aligned_seq1
                    aligned_seq2 = enzyme_2[j - 1] + aligned_seq2
                    traceback_path.append("left")
                i, j = prev_i, prev_j
            return (
                aligned_seq1,
                aligned_seq2,
                list(reversed(traceback_path)),
                list(reversed(traceback_coords)),
            )

        if all_maxima:
            result = []
            for position in max_positions:
                aligned_seq1, aligned_seq2, path, coords = traceback(*position)
                result.append((aligned_seq1, aligned_seq2, max_score, path, coords))
            return result

        else:
            i, j = max_positions[0]
            aligned_seq1, aligned_seq2, path, coords = traceback(i, j)
            return [(aligned_seq1, aligned_seq2, max_score, path, coords)]

    # Levenshteinova vzdálenost (edit distance)
    # Return: edit_distance : int
    def levenshtein_distance_alignment(self, enzyme_1: str, enzyme_2: str) -> int:
        enzyme_1 = enzyme_1.replace(" ", "").replace("\n", "")
        enzyme_2 = enzyme_2.replace(" ", "").replace("\n", "")

        M, N = len(enzyme_1), len(enzyme_2)
        distance_matrix = [[0] * (N + 1) for _ in range(M + 1)]

        for i in range(M + 1):
            distance_matrix[i][0] = i
        for j in range(N + 1):
            distance_matrix[0][j] = j

        for i in range(1, M + 1):
            for j in range(1, N + 1):
                cost = 0 if enzyme_1[i - 1] == enzyme_2[j - 1] else 1
                distance_matrix[i][j] = min(
                    distance_matrix[i - 1][j] + 1,  # odstranění
                    distance_matrix[i][j - 1] + 1,  # vložení
                    distance_matrix[i - 1][j - 1] + cost,  # substituce
                )

        return (distance_matrix, distance_matrix[M][N])

    # ----------------------------------------------------------------------------------------------------
    # Calculates alignment statistics based on aligned sequences
    # dict with keys: matches, mismatches, gaps, identity_percent, alignment_length, score
    def needleman_wunsch_stats(
        self, aligned_seq1: str, aligned_seq2: str, score: int
    ) -> dict:
        matches = mismatches = gaps = 0
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a == "-" or b == "-":
                gaps += 1
            elif a == b:
                matches += 1
            else:
                mismatches += 1

        length = len(aligned_seq1)
        identity_percent = matches / length * 100 if length > 0 else 0.0

        return {
            "matches": matches,
            "mismatches": mismatches,
            "gaps": gaps,
            "identity_percent": identity_percent,
            "alignment_length": length,
            "score": score,
        }

    # Calculates alignment statistics based on aligned sequences
    # dict with keys: matches, mismatches, gaps, identity_percent, alignment_length, score
    def smith_waterman_stats(
        self, aligned_seq1: str, aligned_seq2: str, score: int
    ) -> dict:
        matches = mismatches = gaps = 0
        for a, b in zip(aligned_seq1, aligned_seq2):
            if a == "-" or b == "-":
                gaps += 1
            elif a == b:
                matches += 1
            else:
                mismatches += 1
        length = len(aligned_seq1)
        identity_percent = matches / length * 100 if length > 0 else 0.0
        return {
            "matches": matches,
            "mismatches": mismatches,
            "gaps": gaps,
            "identity_percent": identity_percent,
            "alignment_length": length,
            "score": score,
        }

    # Vrací statistiku podobnosti na základě Levenshteinovy vzdálenosti
    # Return: dict with keys: edit_distance, identity_percent
    def levenshtein_stats(self, enzyme_1: str, enzyme_2: str) -> dict:
        edit_distance = self.levenshtein_distance(enzyme_1, enzyme_2)
        max_len = max(len(enzyme_1), len(enzyme_2))
        identity_percent = (1 - edit_distance / max_len) * 100 if max_len > 0 else 0.0

        return {
            "edit_distance": edit_distance,
            "identity_percent": identity_percent,
        }

    # ----------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------

    # Displays a heatmap of the score matrix with the traceback path overlaid (for Needleman-Wunsch).
    def plot_needleman_wunsch_heatmap_with_trace(
        self, score_matrix, pointer_matrix, return_fig=False
    ):
        matrix = np.array(score_matrix)
        fig, ax = plt.subplots(figsize=(6, 6))
        cax = ax.imshow(matrix, cmap="viridis", origin="upper")
        fig.colorbar(cax)

        # Reconstruct the traceback path from bottom right to top left
        i, j = len(score_matrix) - 1, len(score_matrix[0]) - 1
        trace_coords = [(i, j)]
        while pointer_matrix[i][j] is not None:
            i, j = pointer_matrix[i][j]
            trace_coords.append((i, j))

        if trace_coords:
            y_coords, x_coords = zip(*trace_coords)
            ax.plot(
                x_coords,
                y_coords,
                color="red",
                marker="o",
                markersize=4,
                linewidth=2,
                label="Traceback path",
            )
            ax.legend()

        ax.set_title("Needleman-Wunsch Score Matrix with Traceback Path")
        ax.set_xlabel("Sequence 2")
        ax.set_ylabel("Sequence 1")

        if return_fig:
            return fig
        else:
            plt.show()

    # Displays a heatmap of the score matrix with the traceback path overlaid (for Smith-Waterman).
    def plot_smith_waterman_heatmap_with_trace(
        self, score_matrix, trace_coords_list, return_fig=False
    ):
        matrix = np.array(score_matrix)
        fig, ax = plt.subplots(figsize=(6, 6))
        cax = ax.imshow(matrix, cmap="viridis", origin="upper")
        fig.colorbar(cax)

        if trace_coords_list:
            for i, trace_coords in enumerate(trace_coords_list):
                y_coords, x_coords = zip(*trace_coords)
                ax.plot(
                    [x for x in x_coords],
                    [y for y in y_coords],
                    marker="o",
                    markersize=3,
                    linewidth=2,
                    label=f"Traceback {i+1}",
                    alpha=0.8,
                )
            ax.legend()

        ax.set_title("Smith-Waterman Score Matrix with Traceback Path")
        ax.set_xlabel("Sequence 2")
        ax.set_ylabel("Sequence 1")

        if return_fig:
            return fig
        else:
            plt.show()

    def plot_levenshtein_heatmap(self, distance_matrix, return_fig=False):
        import numpy as np
        import matplotlib.pyplot as plt

        matrix = np.array(distance_matrix)
        fig, ax = plt.subplots(figsize=(6, 6))
        cax = ax.imshow(matrix, cmap="plasma", origin="upper")
        fig.colorbar(cax)

        ax.set_title("Levenshtein Distance Matrix")
        ax.set_xlabel("Sequence 2")
        ax.set_ylabel("Sequence 1")

        if return_fig:
            return fig
        else:
            plt.show()

    def format_alignment_blast_style(self, seq1: str, seq2: str) -> str:
        """
        Formats the alignment of two sequences in BLAST-style alignment:
        Query:   ACTG-ACTG
                || | ||||
        Sbjct:   ACTGCACTG
        """
        match_line = ""
        for a, b in zip(seq1, seq2):
            if a == b:
                match_line += "|"
            elif a == "-" or b == "-":
                match_line += " "
            else:
                match_line += "."

        return f"Query:  {seq1}\n        {match_line}\nSbjct:  {seq2}"



