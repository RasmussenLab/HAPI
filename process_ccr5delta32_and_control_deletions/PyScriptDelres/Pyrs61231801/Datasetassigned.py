#!/usr/bin/env python3
import sys, re, csv
import os
# Author: Kirstine Ravn
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

DelSNPname = "rs61231801"

with open(rf"C:\Users\k_rav\OneDrive\Skrivebord\Speciale\1NOV22120\Cell\Cell_reviewers_comments\Deletions_background\Results_background\New_files150724\PyScriptDelres\Py{DelSNPname}\AssignedDataset.txt", "r") as datasetAssigned:
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
output_file_path = f"Count.top4_{DelSNPname}.tsv"
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
