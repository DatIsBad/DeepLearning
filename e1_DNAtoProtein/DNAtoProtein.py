genetic_code = {
    'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
    'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys', 
    'UUA': 'Leu', 'UCA': 'Ser', 'UAA': 'STOP', 'UGA': 'STOP', 
    'UUG': 'Leu', 'UCG': 'Ser', 'UAG': 'STOP', 'UGG': 'Trp', 
    'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
    'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
    'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
    'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
    'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser', 
    'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser', 
    'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
    'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
    'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly', 
    'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
    'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
    'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
}

translation = {
    'Ala' :'A',  'Gly' :'G',  'Met' :'M',  'Ser' :'S',
    'Cys' :'C',  'His' :'H',  'Asn' :'N',  'Thr' :'T', 
    'Asp' :'D',  'Ile' :'I',  'Pro' :'P',  'Val' :'V', 
    'Glu' :'E',  'Lys' :'K',  'Gln' :'Q',  'Trp' :'W', 
    'Phe' :'F',  'Leu' :'L',  'Arg' :'R',  'Tyr' :'Y'
}


# transcription from DNA to RNA
# replaces thymin ('T') for uracil ('U')
def transcribe(dna):
    return dna.replace('T', 'U')



# translation from RNA to Protein
def translate(rna):
    protein = []
    output = ""

    for i in range(0, len(rna) - 2, 3): # Process every 3 bases (codon)
        codon = rna[i:i+3]
        amino_acid = genetic_code.get(codon, '') #if it's not there, return empty

        if amino_acid == 'STOP':
            break
        else:
            protein.append(amino_acid)

    for i in range(0, len(protein)):
        output = output + translation.get(protein[i], '')

    return output



def dnaToProtein(header,data):
    parts = header.split()
    newHeader = ""
    if(len(parts) == 4):
        newHeader = parts[0] + '   ' + parts[1]
    else:
        newHeader = parts[0]

    newData = translate(transcribe(data))
    newHeader = newHeader + ' ' + str(len(newData)) + ' aa\n'

    result = ''
    for i in range(0, len(newData), 10):
        # Add a chunk of 10 characters followed by a space
        result += newData[i:i+10] + ' '
        # After every 50 characters (i.e., every 5 chunks of 10), add a newline
        if (i // 10 + 1) % 5 == 0:
            result += '\n'

    result = newHeader + result

    return result



def openFile(filePath):
    first10 = 0

    with open(filePath, 'r') as fileDNA:
        with open('output.txt', 'w') as fileProtein:
            print(f"file {filePath} was opened")
            lines = fileDNA.readlines()

            header = ""
            gene_data = ""

            #go through each line
            for line in lines:
                if len(line) > 1:
                    if line[0] == '>':
                        header = line.rstrip('\n')
                    else:
                        gene_data = gene_data + line.replace(' ','').replace('\n', '')
                
                elif header != "":
                    first10 += 1

                    result = dnaToProtein(header, gene_data)
                    fileProtein.write(result)
                    fileProtein.write('\n\n')
                    header = ""
                    gene_data = ""