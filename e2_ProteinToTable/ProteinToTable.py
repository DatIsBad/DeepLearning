import numpy as np
from pathlib import Path

def processHeader(value:str):
    data = value.replace(",", " ").split()
    variation = ""
    orf = ""
    modifier = []
    id = data[0].replace(">", "")
    size = 0
    fragment = False
    temp = ""

    if '.' in id:
        temp = id.split('.')
        variation = temp[0]
        id = temp[1]

    if 'ORF' in id:
        temp = id.split('ORF')
        orf = temp[1]
        id = temp[0]

    i = 1
    while data[i].isdigit() == False:
        modifier.append(data[i])
        i = i + 1

    size = int(data[i])



    if data[ len(data) - 1 ] == 'fragment':
        fragment = True

    result = [id, variation, orf, modifier, size, fragment]
    
    return result

def openFile(filePath:str):
    result = []
    
    with open(filePath, 'r') as file:
        lines = file.readlines()
        header = ""
        num_line = 0
        num_header = -1
        path_parts = Path(filePath).parts
        #gene_data = ""

        for line in lines:
            if len(line) > 1:
                if line[0] == '>':
                    header = line.rstrip('\n')
                    num_header = num_line

                else:
                    #gene_data += line.rstrip('\n').replace(' ', '')

                    pass
            
            else:
                if len(header) > 0:
                    temp1 = processHeader(header)
                    #temp1.append(gene_data)
                    temp1.append( path_parts[len(path_parts) - 1] )
                    temp1.append(num_header)
                    result.append(temp1)
                    header = ""
                    num_header = -1
                    #gene_data = ""

            num_line = num_line + 1
    return result

path1 = "DeepLearning\\DATA\\Type_I_restriction_enzymes_Gold_Standards_Protein.txt"
path2 = "DeepLearning\\DATA\\Type_II_restriction_enzymes_Gold_Standards_Protein.txt"
path3 = "DeepLearning\\DATA\\Type_III_restriction_enzymes_Gold_Standards_Protein.txt"
path4 = "DeepLearning\\DATA\\Type_VI_restriction_enzymes_Gold_Standards_Protein.txt"


'''data = openFile(path1)

i = 0
for item in data:
    print(item)
    i = i + 1

    if i > 10:
        break
'''
