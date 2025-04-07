import numpy as np
import matplotlib.pyplot as plt


class ProcessSimilarity:
    def __init__(self):
        pass

    # Needleman-Wunsh (global alignment)
    # return: score_matrix : [], pointer_matrix : []
    def needleman_wunsh_matrix(self, enzyme_1, enzyme_2, match=1, mismatch=-1, gap=-2):
        enzyme_1 = enzyme_1.replace(" ", "").replace("\n", "")
        enzyme_2 = enzyme_2.replace(" ", "").replace("\n", "")

        M, N = len(enzyme_1), len(enzyme_2)
        score_matrix = [[0] * (N + 1) for _ in range(M + 1)]
        pointer_matrix = [[None] * (N + 1) for _ in range(M + 1)]

        for i in range(M + 1):
            score_matrix[i][0] = i * gap
            pointer_matrix[i][0] = (i - 1, 0) if i > 0 else None
        for j in range(N + 1):
            score_matrix[0][j] = j * gap
            pointer_matrix[0][j] = (0, j - 1) if j > 0 else None

        for i in range(1, M + 1):
            for j in range(1, N + 1):
                diag = score_matrix[i - 1][j - 1] + (
                    match if enzyme_1[i - 1] == enzyme_2[j - 1] else mismatch
                )
                delete = score_matrix[i - 1][j] + gap
                insert = score_matrix[i][j - 1] + gap
                score = max(diag, delete, insert)
                score_matrix[i][j] = score

                if score == diag:
                    pointer_matrix[i][j] = (i - 1, j - 1)
                elif score == delete:
                    pointer_matrix[i][j] = (i - 1, j)
                else:
                    pointer_matrix[i][j] = (i, j - 1)

        return (score_matrix, pointer_matrix)

    # Smith-Waterman (local alignment)
    # Return: score_matrix : [], pointer_matrix : [], max_score : int, max_positions : []
    def smith_waterman_matrix(self, enzyme_1, enzyme_2, match=1, mismatch=-1, gap=-2):
        enzyme_1 = enzyme_1.replace(" ", "").replace("\n", "")
        enzyme_2 = enzyme_2.replace(" ", "").replace("\n", "")

        M, N = len(enzyme_1), len(enzyme_2)
        score_matrix = np.zeros((M + 1, N + 1))
        pointer_matrix = [[None] * (N + 1) for _ in range(M + 1)]
        max_score = 0
        max_positions = []

        for i in range(1, M + 1):
            for j in range(1, N + 1):
                diag = score_matrix[i - 1][j - 1] + (
                    match if enzyme_1[i - 1] == enzyme_2[j - 1] else mismatch
                )
                delete = score_matrix[i - 1][j] + gap
                insert = score_matrix[i][j - 1] + gap
                score = max(0, diag, delete, insert)
                score_matrix[i][j] = score

                if score == 0:
                    pointer_matrix[i][j] = None
                elif score == diag:
                    pointer_matrix[i][j] = (i - 1, j - 1)
                elif score == delete:
                    pointer_matrix[i][j] = (i - 1, j)
                else:
                    pointer_matrix[i][j] = (i, j - 1)

                if score > max_score:
                    max_score = score
                    max_positions = [(i, j)]
                elif score == max_score:
                    max_positions.append((i, j))

        return (score_matrix, pointer_matrix, max_score, max_positions)

    # compute scores of enzyme_1 and enzyme_2 with needleman method
    # return: (aligned_seq1, aligned_seq2, score_matrix[M][N])
    def needleman_wunsh_alignment(
        self, enzyme_1, enzyme_2, match=1, mismatch=-1, gap=-2
    ):
        enzyme_1 = enzyme_1.replace(" ", "").replace("\n", "")
        enzyme_2 = enzyme_2.replace(" ", "").replace("\n", "")

        M, N = len(enzyme_1), len(enzyme_2)
        score_matrix, pointer_matrix = self.needleman_wunsh_matrix(
            enzyme_1, enzyme_2, match, mismatch, gap
        )
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
    # return: [(aligned_seq1, aligned_seq2, max_score, max_score, path, coords)]
    def smith_waterman_alignment(
        self, enzyme_1, enzyme_2, match=1, mismatch=-1, gap=-2, all_maxima=True
    ):
        enzyme_1 = enzyme_1.replace(" ", "").replace("\n", "")
        enzyme_2 = enzyme_2.replace(" ", "").replace("\n", "")

        score_matrix, pointer_matrix, max_score, max_positions = (
            self.smith_waterman_matrix(enzyme_1, enzyme_2, match, mismatch, gap)
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

    # ----------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------

    # Displays a heatmap of the score matrix with the traceback path overlaid (for Needleman-Wunsch).
    def plot_needleman_wunsch_heatmap_with_trace(self, score_matrix, pointer_matrix, return_fig=False):
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
    def plot_smith_waterman_heatmap_with_trace(self, score_matrix, trace_coords_list, return_fig=False):
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

    # visualise score matrix as heatmap
    def visualise_hm_score_matrix(self, score_matrix):
        # Vizualizace matice skóre
        plt.figure(figsize=(8, 6))
        plt.imshow(score_matrix, cmap="Blues", interpolation="nearest")
        plt.colorbar(label="Score")
        plt.title("Score Matrix")
        plt.xlabel("Sequence 2")
        plt.ylabel("Sequence 1")
        plt.show()
