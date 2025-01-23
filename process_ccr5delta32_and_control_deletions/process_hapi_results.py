#!/usr/bin/env python3
# Author: Kirstine Ravn
import sys, re, csv
import os
import argparse

def readInt(str):
    if str == "NA":
        return -1
    else:
        return int(str)


#Script after runned on data from Martin list 
#Filepath: outFile_path = rf"C:\Users\k_rav\OneDrive\Skrivebord\Speciale\1NOV22120\Cell\Cell_reviewers_comments\Deletions_background\Results_background\New_files150724\
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

# DelSNPname = "rs333"
# basepath = "/Users/lmz306/Library/CloudStorage/OneDrive-UniversityofCopenhagen/delta_ccr5/2024_07_15_latest_results_8del/2024_07_18_Newfiles_Kirstine"

folderpath = f"{basepath}/Py{DelSNPname}/"


#-----------------------------
# Datasetassigned.py
#-----------------------------
#!/usr/bin/env python3
import sys, re, csv
import os

# Initialize dictionaries to store counts for each dataset and filter combination
data_group_filter_counts = {
    "VK": {"No_filter": {"RR": 0, "RD": 0, "DD": 0}, "Permissive_filter": {"RR": 0, "RD": 0, "DD": 0}, "Strict_filter": {"RR": 0, "RD": 0, "DD": 0}},
    "RISE": {"No_filter": {"RR": 0, "RD": 0, "DD": 0}, "Permissive_filter": {"RR": 0, "RD": 0, "DD": 0}, "Strict_filter": {"RR": 0, "RD": 0, "DD": 0}},
    "NEO": {"No_filter": {"RR": 0, "RD": 0, "DD": 0}, "Permissive_filter": {"RR": 0, "RD": 0, "DD": 0}, "Strict_filter": {"RR": 0, "RD": 0, "DD": 0}},
    "BOTAI": {"No_filter": {"RR": 0, "RD": 0, "DD": 0}, "Permissive_filter": {"RR": 0, "RD": 0, "DD": 0}, "Strict_filter": {"RR": 0, "RD": 0, "DD": 0}}
}

# Initialize dictionaries to store total counts for each dataset
total_counts = {"VK": 0, "RISE": 0, "NEO": 0, "BOTAI": 0}

# Open the dataset file for reading


with open(f"{basepath}/Py{DelSNPname}/AssignedDataset.txt", "r") as datasetAssigned:
    # Read the dataset line by line
    for line in datasetAssigned:
        # Skip the header line
        if line.startswith("Assigned"):
            continue

        # Split the line into columns
        columns = line.strip().split("\t")  # Use strip() to remove any trailing newline characters or spaces

        # Ensure the line has the expected number of columns
        if len(columns) < 5:
            print(f"Skipping line with unexpected number of columns: {line.strip()}")
            continue

        # Extract values from specific columns
        data_group = columns[0]
        No_filter = columns[2].strip()  # Remove any leading/trailing whitespace
        Permissive_filter = columns[3].strip()  # Remove any leading/trailing whitespace
        Strict_filter = columns[4].strip()  # Remove any leading/trailing whitespace

        # Update total count for the data group
        total_counts[data_group] += 1

        # Count occurrences based on conditions for No_filter
        if No_filter in data_group_filter_counts[data_group]["No_filter"]:
            data_group_filter_counts[data_group]["No_filter"][No_filter] += 1
        else:
            print(f"Unexpected value '{No_filter}' in No_filter column for data group '{data_group}'")

        # Count occurrences based on conditions for Permissive_filter
        if Permissive_filter in data_group_filter_counts[data_group]["Permissive_filter"]:
            data_group_filter_counts[data_group]["Permissive_filter"][Permissive_filter] += 1
        else:
            print(f"Unexpected value '{Permissive_filter}' in Permissive_filter column for data group '{data_group}'")

        # Count occurrences based on conditions for Strict_filter
        if Strict_filter in data_group_filter_counts[data_group]["Strict_filter"]:
            data_group_filter_counts[data_group]["Strict_filter"][Strict_filter] += 1
        else:
            print(f"Unexpected value '{Strict_filter}' in Strict_filter column for data group '{data_group}'")

# Write the counts to the output file
output_file_path = f"{folderpath}/Count.top4_{DelSNPname}.tsv"
with open(output_file_path, "w") as CountOutFile:
    for data_group in ["VK", "RISE", "NEO", "BOTAI"]:
        CountOutFile.write(f"Data Group: {data_group}\n")
        CountOutFile.write(f"Total samples: {total_counts[data_group]}\n")
        CountOutFile.write(f"No_filter - RR count: {data_group_filter_counts[data_group]['No_filter']['RR']}\n")
        CountOutFile.write(f"No_filter - RD count: {data_group_filter_counts[data_group]['No_filter']['RD']}\n")
        CountOutFile.write(f"No_filter - DD count: {data_group_filter_counts[data_group]['No_filter']['DD']}\n")
        CountOutFile.write(f"Permissive_filter - RR count: {data_group_filter_counts[data_group]['Permissive_filter']['RR']}\n")
        CountOutFile.write(f"Permissive_filter - RD count: {data_group_filter_counts[data_group]['Permissive_filter']['RD']}\n")
        CountOutFile.write(f"Permissive_filter - DD count: {data_group_filter_counts[data_group]['Permissive_filter']['DD']}\n")
        CountOutFile.write(f"Strict_filter - RR count: {data_group_filter_counts[data_group]['Strict_filter']['RR']}\n")
        CountOutFile.write(f"Strict_filter - RD count: {data_group_filter_counts[data_group]['Strict_filter']['RD']}\n")
        CountOutFile.write(f"Strict_filter - DD count: {data_group_filter_counts[data_group]['Strict_filter']['DD']}\n")

print(os.getcwd())
print (f"File {output_file_path} written.")
import csv

input_file_path = f'{folderpath}/Count.top4_{DelSNPname}.tsv'
output_file_path = f"{folderpath}/MAF/resMAF.tsv"

# Initialize the counters for each data group and filter type
counts = {}

# Open and read the input file
with open(input_file_path, 'r') as file:
    reader = file.readlines()
    
    current_group = None
    for line in reader:
        line = line.strip()
        if line.startswith("Data Group:"):
            current_group = line.split(":")[1].strip()
            counts[current_group] = {'No_filter': {}, 'Permissive_filter': {}, 'Strict_filter': {}, 'Total_samples': 0}
        elif line.startswith("Total samples:"):
            counts[current_group]['Total_samples'] = int(line.split(":")[1].strip())
        elif 'count' in line:
            parts = line.split(" - ")
            filter_type = parts[0]
            count_type, count_value = parts[1].split(": ")
            count_value = int(count_value.strip())
            if filter_type not in counts[current_group]:
                counts[current_group][filter_type] = {'RR': 0, 'RD': 0, 'DD': 0}
            count_key = count_type.split()[0]
            counts[current_group][filter_type][count_key] = count_value

# Open the output file and write the results
with open(output_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')

    # Write the header
    writer.writerow(['Data Group', 'Filter Type', 'RR', 'RD', 'DD', 'Total Samples', 'Total Alleles', 'MAF'])

    # Write the results for each data group and filter type
    for data_group, filter_counts in counts.items():
        total_samples = filter_counts.pop('Total_samples', 0)
        for filter_type, count_dict in filter_counts.items():
            RD_number = count_dict['RD']
            DD_number = count_dict['DD']
            total_alleles = 2 * total_samples
            MAF = (RD_number + DD_number) / total_alleles if total_alleles != 0 else 0
            # Formatting MAF to display only four decimal places
            MAF_formatted = "{:.4f}".format(MAF)
            writer.writerow([
                data_group, filter_type, count_dict['RR'], RD_number, DD_number, total_samples, total_alleles, MAF_formatted
            ])

print("done")

import pandas as pd
import openpyxl
import os

# os.getcwd()
# basepath = "/Users/lmz306/Library/CloudStorage/OneDrive-UniversityofCopenhagen/delta_ccr5/2024_07_15_latest_results_8del/2024_07_18_Newfiles_Kirstine/PyScriptDelres_no_duplicated_RISE"
# DelSNPname = "rs333"
# folderpath = f"{basepath}/Py{DelSNPname}/"


resHAP = pd.read_csv(f"{folderpath}/res.top4_{DelSNPname}.tsv", sep="\t")
Assigned= pd.read_csv(f"{folderpath}/AssignedDataset.txt", sep="\t")
Counts = pd.read_csv(f"{folderpath}/Count.top4_{DelSNPname}.tsv", sep="\t", names=["Column"])
MAF = pd.read_csv(f"{folderpath}/MAF/resMAF.tsv", sep="\t")


with pd.ExcelWriter(f"{folderpath}/results_{DelSNPname}.xlsx", engine='openpyxl') as writer:
    resHAP.to_excel(writer, sheet_name='resHAP', index=False)
    Assigned.to_excel(writer, sheet_name='Assigned', index=False)
    Counts.to_excel(writer, sheet_name='Counts', index=False)
    MAF.to_excel(writer, sheet_name='MAF', index=False)

# select columns Assigned, Sample, No_filter, Permissive_filter, Strict_filter, pRR_Data_n, pRD_Data_n, pDD_Data_n
Assigned_forevan = Assigned.filter(['Assigned', 'Sample', 'No_filter', 'Permissive_filter', 'Strict_filter', 'pRR_Data_n', 'pRD_Data_n', 'pDD_Data_n'])

# order by Assigned
Assigned_forevan = Assigned_forevan.sort_values(by=['Assigned'])
with pd.ExcelWriter(f"{basepath}/{DelSNPname}_E.xlsx", engine='openpyxl') as writer:
    Assigned_forevan.to_excel(writer, sheet_name=f"{DelSNPname}", index=False)
    

# with pd.ExcelWriter(f"{basepath}/summary_8del.xlsx", engine='openpyxl', mode='a') as writer:
#     MAF.to_excel(writer, sheet_name=f"{DelSNPname}_MAF", index=False)

print("Excel file created successfully.")

