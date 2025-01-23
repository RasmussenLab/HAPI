#!/usr/bin/env python3
import sys, re, csv
import os

def readInt(str):
    if str == "nan":
        return -1
    else:
        return int(str)
        
DelSNPname = "rs66552573"

outFile_path = rf"C:\Users\k_rav\OneDrive\Skrivebord\Speciale\1NOV22120\Cell\Cell_reviewers_comments\Deletions_background\Results_background\PyScriptDelres\Py{DelSNPname}\res.top4_{DelSNPname}.tsv"
highlightOutFile_path = rf"C:\Users\k_rav\OneDrive\Skrivebord\Speciale\1NOV22120\Cell\Cell_reviewers_comments\Deletions_background\Results_background\PyScriptDelres\Py{DelSNPname}\resStar.top4_{DelSNPname}.tsv"
artefact_path = rf"C:\Users\k_rav\OneDrive\Skrivebord\Speciale\1NOV22120\Cell\Cell_reviewers_comments\Deletions_background\Results_background\PyScriptDelres\Py{DelSNPname}\resArtefact.{DelSNPname}.tsv"
Hapi_res_path = rf"C:\Users\k_rav\OneDrive\Skrivebord\Speciale\1NOV22120\Cell\Cell_reviewers_comments\Deletions_background\Results_background\Top4\top4_{DelSNPname}.tsv"

# Open the output files using the defined variables
outFile = open(outFile_path, "w")
highlightOutFile = open(highlightOutFile_path, "w")
artefactOutFile = open(artefact_path, "w")

# Process the Hapi_res file
firstLine = True
for line in open(Hapi_res_path, "r"):
    if firstLine:
        headerColumns = line.strip().split("\t")
        outFile.write(line)
        highlightOutFile.write(line)
        artefactOutFile.write(line)
        firstLine = False
        continue

    columns = line.strip().split("\t")

#    if columns[0] != 'NEO516':
#        continue

#rs13077753_ref	rs13077753_alt	rs12053824_ref	rs12053824_alt	rs13093624_ref	rs13093624_alt	rs67007922_ref	rs67007922_alt

    N_reads_ref = readInt(columns[headerColumns.index("N_reads_ref")])
    print(f"N_reads_ref: {N_reads_ref}")
    ref1 = readInt(columns[headerColumns.index("rs4073957_ref")])
    print(f"ref1: {ref1}")
    ref2 = readInt(columns[headerColumns.index("rs34901767_ref")])
    print(f"ref2: {ref2}")
    ref3 = readInt(columns[headerColumns.index("rs7642436_ref")]) 
    print(f"ref3: {ref3}")
    ref4 = readInt(columns[headerColumns.index("rs4508714_ref")])# #Not ancient damage 
    print(f"ref4: {ref4}")

    N_reads_del = readInt(columns[headerColumns.index("N_reads_del")])
    print(f"N_reads_del: {N_reads_del}")
    alt1 = readInt(columns[headerColumns.index("rs4073957_alt")])
    print(f"alt1: {alt1}")
    alt2 = readInt(columns[headerColumns.index("rs34901767_alt")])
    print(f"alt2: {alt2}")
    alt3 = readInt(columns[headerColumns.index("rs7642436_alt")])
    print(f"alt3: {alt3}")
    alt4 = readInt(columns[headerColumns.index("rs4508714_alt")])  #Not ancient damage 
    print(f"alt4: {alt4}")

    sum_x = readInt(columns[headerColumns.index("sum_x")])
    print(f"sum_x: {sum_x}")

    referencecount = readInt(columns[headerColumns.index("referencecount")])
    print(f"referencecount: {referencecount}")
    haplocount = readInt(columns[headerColumns.index("haplocount")])
    print(f"haplocount: {haplocount}")

    if N_reads_del > 0:
        if N_reads_ref == 0 and N_reads_del >=2 and sum_x >= 1 and ref1 < 1 and ref2 < 1 and ref3 < 1 and ref4 < 1:
            permissive = "DD"
        elif N_reads_ref == 0 and N_reads_del >=1 and sum_x >= 1 and ref1 < 1 and ref2 < 1 and ref3 < 1 and ref4 < 1 and referencecount < 10:
            permissive = "DD"
    
        elif alt4 > 0 and ref1 < 1 and ref2 < 1 and ref3 < 1 and ref4 < 1:
            permissive = "RD"
        elif sum_x >= 1 or alt4 >= 1: 
            permissive = "RD"
        elif N_reads_del >= 2:
            permissive = "RD"
        elif haplocount >= 30: 
           permissive = "RD"
        else:
            permissive = "RR"

    elif N_reads_del == 0  and (sum_x >= 2 or alt4 >= 1) and haplocount >= 18:
        permissive = "RD"
    else:
        permissive = "RR"

    if permissive == "DD" and haplocount >= 18:
        strict = "DD"

    elif N_reads_del > 0 and (sum_x >= 2 or alt4 >= 1) and haplocount >= 18: 
        strict = "RD"
    
    elif N_reads_del > 0 and sum_x >= 1 and haplocount >= 30: 
        strict = "RD"

    elif N_reads_del > 1 and (sum_x >= 1 or haplocount >= 18):
        strict = "RD"
  
    elif N_reads_del == 0 and sum_x >= 2 or alt4 >= 1 and haplocount >= 30:#or and 
        strict = "RD"
#elif N_reads_del == 0 and sum_x >= 2 and alt4 >= 1 and haplocount >= 18:#or and 
        #strict = "RD"
    else:
        strict = "RR"



    threshold = 10  
    if N_reads_del != 0 and N_reads_ref / N_reads_del > threshold:
        permissive = "RR"
        strict = "RR"
        artefactOutFile.write(line)



    print(f"strict: {strict}")
    print(f"permissive: {permissive}")

    outFile.write(columns[0])
    outFile.write("\t" + columns[1])
    outFile.write("\t" + permissive)
    outFile.write("\t" + strict)
    for i in range(4, len(columns)):
        outFile.write("\t" + columns[i])
    outFile.write("\n")

    highlightOutFile.write(columns[0])
    highlightOutFile.write("\t" + columns[1])
    highlightOutFile.write("\t" + permissive)
    if columns[1] != permissive:
        highlightOutFile.write("*")
    highlightOutFile.write("\t" + strict)
    if columns[1] != strict:
        highlightOutFile.write("*")
    for i in range(4, len(columns)):
        highlightOutFile.write("\t" + columns[i])
    highlightOutFile.write("\n")

# Close the files
outFile.close()
highlightOutFile.close()
import pandas as pd
# read outFile
res = pd.read_csv(outFile_path, sep='\t')
mapping = pd.read_csv(f"{basepath}/PyScriptDelres_no_duplicated_RISE/mapping_Sample_Assigned.tsv", sep='\t', header=None, names=['Assigned', 'Sample'])
# left join mapping to res and put Assigned as first column
res = pd.merge(res, mapping, on='Sample', how='left')

res = res[['Assigned'] + [col for col in res.columns if col != 'Assigned']]

res.to_csv(f"{basepath}/PyScriptDelres_no_duplicated_RISE/Py{DelSNPname}/AssignedDataset.txt", sep='\t', index=False)

