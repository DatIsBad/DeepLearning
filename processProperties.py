import re

# map to convert degenerated DNA symbols to regex syntax
degenerate_dna_map = {
    "A": "A",
    "T": "T",
    "G": "G",
    "C": "C",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ATGC]",
}


# convert degenerated DNA (recognition sequence) symbols to regex (regular expression)
def convert_dna_to_regex(seq: str) -> str:
    return "".join(degenerate_dna_map.get(base.upper(), base) for base in seq)


# count all recognition sites in dna sequence using recognition sequence, that can be degenerated
def count_recognition_sites_in_dna(dna_sequence: str, recognition_seq: str) -> int:
    dna_sequence = dna_sequence.replace("\n", "").replace(" ", "")
    print(dna_sequence)
    pattern = convert_dna_to_regex(recognition_seq)
    matches = re.findall(f"(?={pattern})", dna_sequence.upper())
    print(matches)
    return len(matches)




# return protein distribution count
def protein_distribution(protein_sequence: str):
    protein_sequence = protein_sequence.replace("\n", "").replace(" ", "")
    sum = len(protein_sequence)
    return (
        protein_sequence.count("A") / sum,
        protein_sequence.count("G") / sum,
        protein_sequence.count("M") / sum,
        protein_sequence.count("S") / sum,
        protein_sequence.count("C") / sum,
        protein_sequence.count("H") / sum,
        protein_sequence.count("N") / sum,
        protein_sequence.count("T") / sum,
        protein_sequence.count("D") / sum,
        protein_sequence.count("I") / sum,
        protein_sequence.count("P") / sum,
        protein_sequence.count("V") / sum,
        protein_sequence.count("E") / sum,
        protein_sequence.count("K") / sum,
        protein_sequence.count("Q") / sum,
        protein_sequence.count("W") / sum,
        protein_sequence.count("F") / sum,
        protein_sequence.count("L") / sum,
        protein_sequence.count("R") / sum,
        protein_sequence.count("Y") / sum,
    )


def process_properties(data):
    def openFile(filePath):
        #with open(filePath, 'r') as file:



        pass
    properties = []
    # item = tuple (id, enzyme, sample, orf, rec_sequence, size, fragment, filename, line)
    for item in data:
        id, enzyme_name, sample, orf, rec_sequence, size, fragment, filename, line = (
            item
        )

        distribution = ()
        idk = 0
        properties.append(
            (
                id,
                size,
                distribution,
                idk,
                fragment,
            )
        )
    return properties


