#!/usr/bin/env python3
# Author: Kirstine Ravn
import sys, re, csv
import os

def readInt(str):
    if str == "NA":
        return -1
    else:
        return int(str)
        

import argparse
parser = argparse.ArgumentParser(description="Process results files from HAPI, after manual curation")
parser.add_argument("--deletion_rsid", type=str, help="The name of the deletion SNP")
parser.add_argument("--path", type=str, help="The path to the directory containing the results files", default=".")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

# add help
args = parser.parse_args()

DelSNPname = args.deletion_rsid
basepath = args.path


basepath = "/Users/lmz306/Library/CloudStorage/OneDrive-UniversityofCopenhagen/delta_ccr5/2024_07_15_latest_results_8del/2024_07_18_Newfiles_Kirstine"
outFile_path = f"{basepath}/PyScriptDelres_no_duplicated_RISE/Py{DelSNPname}/res.top4_{DelSNPname}.tsv"
highlightOutFile_path = f"{basepath}/PyScriptDelres_no_duplicated_RISE/Py{DelSNPname}/resStar.top4_{DelSNPname}.tsv"
artefact_path = f"{basepath}/PyScriptDelres_no_duplicated_RISE/Py{DelSNPname}/resArtefact.{DelSNPname}.tsv"
Hapi_res_path = f"{basepath}/Top4/top4_{DelSNPname}.tsv"

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

    #if columns[0] != 'RISE493':
        #continue

# rs17787940_ref	rs17787940_alt	rs59909742_ref	rs59909742_alt	rs75045913_ref	rs75045913_alt	rs17195637_ref	rs17195637_alt
    N_reads_ref = readInt(columns[headerColumns.index("N_reads_ref")])
    print(f"N_reads_ref: {N_reads_ref}")
    ref1 = readInt(columns[headerColumns.index("rs17787940_ref")]) 
    print(f"ref1: {ref1}")
    ref2 = readInt(columns[headerColumns.index("rs59909742_ref")])
    print(f"ref2: {ref2}")
    ref3 = readInt(columns[headerColumns.index("rs75045913_ref")])
    print(f"ref3: {ref3}")
    ref4 = readInt(columns[headerColumns.index("rs17195637_ref")])
    print(f"ref4: {ref4}")

    N_reads_del = readInt(columns[headerColumns.index("N_reads_del")])
    print(f"N_reads_del: {N_reads_del}")
    alt1 = readInt(columns[headerColumns.index("rs17787940_alt")])  #NOT affected by ancient damage
    print(f"alt1: {alt1}")
    alt2 = readInt(columns[headerColumns.index("rs59909742_alt")])  #NOT affected by ancient damage
    print(f"alt2: {alt2}")
    alt3 = readInt(columns[headerColumns.index("rs75045913_alt")]) #NOT affected by ancient damage
    print(f"alt3: {alt3}")
    alt4 = readInt(columns[headerColumns.index("rs17195637_alt")]) #ancient damaged
    print(f"alt4: {alt4}")


    sum_x = readInt(columns[headerColumns.index("sum_x")])
    print(f"sum_x: {sum_x}")

    referencecount = readInt(columns[headerColumns.index("referencecount")])
    print(f"referencecount: {referencecount}")
    haplocount = readInt(columns[headerColumns.index("haplocount")])
    print(f"haplocount: {haplocount}")

    if N_reads_del > 0:
        if N_reads_del >= 2 and N_reads_ref == 0 and referencecount == 0:
            permissive = "DD" #sum_x >= 1 dob allels
        elif N_reads_del >= 1 and N_reads_ref == 0 and referencecount == 0 and haplocount >= 2: #New line added >= 2
            permissive = "DD"
        elif N_reads_del >= 2:
            permissive = "RD"
        elif N_reads_del == 1 and haplocount >=1:
            permissive = "RD"
        else:
            permissive = "RR"  

    elif N_reads_del == 0 and haplocount >= 2 and (referencecount >=1 or N_reads_ref >=1):
         #(N_reads_ref >=1 or referencecount >= 2): 
        permissive = "RD"

    else:
        permissive = "RR"


    if permissive == "DD":
        if haplocount >= 2:
            strict = "DD"
        else:
            strict = "RR"

    elif permissive == "RD":
        if N_reads_del == 0 and (alt1 >=1 or alt2 >=1 or alt3 >=1) and haplocount >= 3: 
            strict = "RD"
        
        elif N_reads_del >= 3:
            strict = "RD"

        elif N_reads_del >= 2 and haplocount >= 1: #(alt1 >=1 or alt2 >=1 or alt3 >=1)
            strict = "RD"

        elif N_reads_del >= 1 and haplocount >= 2:
            strict = "RD"
        else:
            strict = "RR"     

    else:
        strict = "RR"


    threshold = 5
    if N_reads_del != 0 and N_reads_ref / N_reads_del > threshold:
        permissive = "RR"
        strict = "RR"
    elif N_reads_del == 0 and N_reads_ref >= 7:
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

