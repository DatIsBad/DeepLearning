import numpy as np
from pathlib import Path

genetic_code = {
    "UUU": "Phe",
    "UCU": "Ser",
    "UAU": "Tyr",
    "UGU": "Cys",
    "UUC": "Phe",
    "UCC": "Ser",
    "UAC": "Tyr",
    "UGC": "Cys",
    "UUA": "Leu",
    "UCA": "Ser",
    "UAA": "STOP",
    "UGA": "STOP",
    "UUG": "Leu",
    "UCG": "Ser",
    "UAG": "STOP",
    "UGG": "Trp",
    "CUU": "Leu",
    "CCU": "Pro",
    "CAU": "His",
    "CGU": "Arg",
    "CUC": "Leu",
    "CCC": "Pro",
    "CAC": "His",
    "CGC": "Arg",
    "CUA": "Leu",
    "CCA": "Pro",
    "CAA": "Gln",
    "CGA": "Arg",
    "CUG": "Leu",
    "CCG": "Pro",
    "CAG": "Gln",
    "CGG": "Arg",
    "AUU": "Ile",
    "ACU": "Thr",
    "AAU": "Asn",
    "AGU": "Ser",
    "AUC": "Ile",
    "ACC": "Thr",
    "AAC": "Asn",
    "AGC": "Ser",
    "AUA": "Ile",
    "ACA": "Thr",
    "AAA": "Lys",
    "AGA": "Arg",
    "AUG": "Met",
    "ACG": "Thr",
    "AAG": "Lys",
    "AGG": "Arg",
    "GUU": "Val",
    "GCU": "Ala",
    "GAU": "Asp",
    "GGU": "Gly",
    "GUC": "Val",
    "GCC": "Ala",
    "GAC": "Asp",
    "GGC": "Gly",
    "GUA": "Val",
    "GCA": "Ala",
    "GAA": "Glu",
    "GGA": "Gly",
    "GUG": "Val",
    "GCG": "Ala",
    "GAG": "Glu",
    "GGG": "Gly",
}

translation = {
    "Ala": "A",
    "Gly": "G",
    "Met": "M",
    "Ser": "S",
    "Cys": "C",
    "His": "H",
    "Asn": "N",
    "Thr": "T",
    "Asp": "D",
    "Ile": "I",
    "Pro": "P",
    "Val": "V",
    "Glu": "E",
    "Lys": "K",
    "Gln": "Q",
    "Trp": "W",
    "Phe": "F",
    "Leu": "L",
    "Arg": "R",
    "Tyr": "Y",
}


# using data from file (header and dna sequence) we can convert them to RNA and finally to protein; output would be protein data that can be saved into a file
def convert_DNA_to_protein(header, dna_sequence):
    # transcription from DNA to RNA
    # replaces thymin ('T') for uracil ('U')
    def transcribe(dna):
        return dna.replace("T", "U")

    # ------------------------------------------------

    # translation from RNA to Protein
    def translate(rna):
        protein = []
        output = ""

        for i in range(0, len(rna) - 2, 3):  # Process every 3 bases (codon)
            codon = rna[i : i + 3]
            amino_acid = genetic_code.get(codon, "")  # if it's not there, return empty

            if amino_acid == "STOP":
                break
            else:
                protein.append(amino_acid)

        for i in range(0, len(protein)):
            output = output + translation.get(protein[i], "")

        return output

    # ------------------------------------------------

    parts = header.split()
    newHeader = ""
    isFragment = False

    for number in range(0, len(parts)):
        if number == 0:
            newHeader = parts[0] + "  "
        elif parts[number].isdigit() == False:
            newHeader += " " + parts[number]
        else:
            if parts[len(parts) - 1] == "fragment":
                isFragment = True
            break

    newData = translate(transcribe(dna_sequence))

    newHeader = newHeader + " " + str(len(newData)) + " aa"
    if isFragment:
        newHeader += " fragment\n"
    else:
        newHeader += "\n"

    result = ""
    for i in range(0, len(newData), 10):
        # After every 50 characters (i.e., every 5 chunks of 10), add a newline
        if (i // 10) % 5 == 0 and result != "":
            result += "\n"
        # Add a chunk of 10 characters followed by a space
        result += newData[i : i + 10] + " "

    result = newHeader + result

    return result


# open a DNA file and convers it into Protein file
def convert_file_DNA_to_protein(filePath):
    first10 = 0

    with open(filePath, "r") as fileDNA:
        with open("output.txt", "w") as fileProtein:
            print(f"file {filePath} was opened")
            lines = fileDNA.readlines()

            header = ""
            gene_data = ""

            # go through each line
            for line in lines:
                if len(line) > 1:
                    # if header is detected, put it into header
                    if line[0] == ">":
                        header = line.rstrip("\n")
                    elif len(header) > 1:
                        gene_data = gene_data + line.replace(" ", "").replace("\n", "")

                elif header != "":
                    first10 += 1
                    result = convert_DNA_to_protein(header, gene_data)
                    fileProtein.write(result)
                    fileProtein.write("\n\n\n")
                    header = ""
                    gene_data = ""


# enzyme_data_file will return list of enzymes [name, sample number, orf, [recognition sequence], size, isFragment, file name, line in file]
def read_enzyme_headers(filePath: str):
    def processHeader(value: str):
        data = value.replace(",", " ").split()
        sample = ""
        orf = ""
        recognition_seq = []
        id = data[0].replace(">", "")
        size = 0
        fragment = False
        temp = ""

        # check if it's sample of enzyme (R1.BbrUI | R2.BbrUI = both are same enzyme, different "sample")
        if "." in id:
            temp = id.split(".")
            sample = temp[0]
            id = temp[1]

        # check if enzyme as orf
        if "ORF" in id:
            temp = id.split("ORF")
            orf = temp[1]
            id = temp[0]

        # check for recognition sequences and append them into recognition_seq
        i = 1
        while data[i].isdigit() == False:
            recognition_seq.append(data[i])
            i = i + 1

        # after recognition sequence is always sequence size
        size = int(data[i])

        # check if this sample is fragment
        for item in data:
            if item == "fragment":
                fragment = True
                break

        result = [id, sample, orf, recognition_seq, size, fragment]

        return result

    # --------------------------------------
    result = []

    with open(filePath, "r") as file:
        lines = file.readlines()
        header = ""
        num_line = 1
        num_header = -1
        path_parts = Path(filePath).parts

        for line in lines:
            # line with text
            if len(line) > 1:
                if line[0] == ">":
                    header = line.rstrip("\n")
                    num_header = num_line

            # line without text
            else:
                if len(header) > 0:
                    temp1 = processHeader(header)
                    temp1.append(path_parts[len(path_parts) - 1])
                    temp1.append(num_header)
                    result.append(temp1)
                    header = ""
                    num_header = -1

            num_line = num_line + 1

        # if data was not inputed yet
        if len(header) > 0:
            temp1 = processHeader(header)
            temp1.append(path_parts[len(path_parts) - 1])
            temp1.append(num_header)
            result.append(temp1)
            header = ""
            num_header = -1

    return result


# compare 2 files if they match
def compare_files(filePath_1, filePath_2):
    data_1 = read_enzyme_headers(filePath_1)
    data_2 = read_enzyme_headers(filePath_2)

    for i in range(0, len(data_1)):
        if (
            data_1[i][0] == data_2[i][0]
            and data_1[i][1] == data_2[i][1]
            and data_1[i][2] == data_2[i][2]
            and data_1[i][3] == data_2[i][3]
            and data_1[i][4] == data_2[i][4]
            and data_1[i][5] == data_2[i][5]
        ):
            pass
        else:
            print(
                f"Datas of {data_1[i][0]} and {data_2[i][0]} are not matching ({data_1[i][7]} : {data_2[i][7]})"
            )

#fetch sequence of enzyme
def fetch_sequence(fileName, line_start):
    fileName = "DATA\\" + fileName
    with open(fileName, "r") as file:
        lines = file.readlines()
        i = line_start
        sequence = ""

        while(i < len(lines)):
            if len(lines[i]) == 1:
                break

            sequence = sequence + lines[i].strip().replace(" ", "")
            i = i + 1

        return sequence





path1 = "DATA\\Type_I_restriction_enzymes_Gold_Standards_Protein.txt"
path2 = "DATA\\Type_II_restriction_enzymes_Gold_Standards_Protein.txt"
path3 = "DATA\\Type_III_restriction_enzymes_Gold_Standards_Protein.txt"
path4 = "DATA\\Type_VI_restriction_enzymes_Gold_Standards_Protein.txt"

