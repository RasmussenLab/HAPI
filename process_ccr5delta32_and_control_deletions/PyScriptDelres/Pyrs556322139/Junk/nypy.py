#!/usr/bin/env python3
import sys, re, csv
import os

def readInt(str):
    if str == "nan":
        return -1
    else:
        return int(str)
        
DelSNPname = "rs556322139"

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

    #if columns[0] != 'NEO942':
        #continue
    N_reads_del = readInt(columns[headerColumns.index("N_reads_del")])
    print(f"N_reads_del: {N_reads_del}")
    N_reads_ref = readInt(columns[headerColumns.index("N_reads_ref")])
    print(f"N_reads_ref: {N_reads_ref}")

    ref1 = readInt(columns[headerColumns.index("rs3806694_ref")])
    print(f"ref1: {ref1}")
    ref2 = readInt(columns[headerColumns.index("rs34368826_ref")])
    print(f"ref2: {ref2}")
    ref3 = readInt(columns[headerColumns.index("rs2276850_ref")])
    print(f"ref3: {ref3}")
    ref4 = readInt(columns[headerColumns.index("rs13324142_ref")])
    print(f"ref4: {ref4}")
    referencecount = readInt(columns[headerColumns.index("referencecount")])
    print(f"referencecount: {referencecount}")

    alt1 = readInt(columns[headerColumns.index("rs3806694_alt")])  #None acient
    print(f"alt1: {alt1}")
    alt2 = readInt(columns[headerColumns.index("rs34368826_alt")]) #None acient
    print(f"alt2: {alt2}")
    alt3 = readInt(columns[headerColumns.index("rs2276850_alt")])
    print(f"alt3: {alt3}")
    alt4 = readInt(columns[headerColumns.index("rs13324142_alt")])
    print(f"alt4: {alt4}")
    haplocount = readInt(columns[headerColumns.index("haplocount")])
    print(f"haplocount: {haplocount}")
    sum_x = readInt(columns[headerColumns.index("sum_x")])
    print(f"sum_x: {sum_x}")


    if N_reads_del > 0:
        if N_reads_del > 0 and N_reads_ref == 0 and referencecount == 0 and haplocount >= 5:
            permissive = "DD"
        elif N_reads_del >= 3:
            permissive = "RD"
        elif sum_x >=1 or haplocount >= 5:
            permissive = "RD"
        else:
            permissive = "RR"
    else:
        permissive = "RR"

    if permissive == "DD":
        strict = "DD"
           
    elif permissive == "RD":
        strict = "RD"

    elif permissive == "RR":
        strict = "RR"    
        
    threshold = 3  
    if N_reads_del != 0 and N_reads_ref / N_reads_del > threshold and haplocount <=5:
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


threshold = 3  
print(N_reads_del,0,N_reads_del != 0)
if N_reads_del != 0 and N_reads_ref / N_reads_del > threshold and haplocount <=5:
        print("IM HERE")
        print(columns[0])
        permissive = "RR"
        strict = "RR"
        artefactOutFile.write(line)
