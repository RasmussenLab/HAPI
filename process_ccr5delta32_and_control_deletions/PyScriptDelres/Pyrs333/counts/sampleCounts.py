import csv

# Initialize the counters for each filter type
counts = {
    'No_filter': {'RR': 0, 'RD': 0, 'DD': 0},
    'Permissive_filter': {'RR': 0, 'RD': 0, 'DD': 0},
    'Strict_filter': {'RR': 0, 'RD': 0, 'DD': 0}
}

# File paths
input_file_path = 'counts/table.txt'
output_file_path = 'counts/XresCounts.tsv'

# Open and read the input file
with open(input_file_path, 'r') as file:
    reader = csv.reader(file, delimiter='\t')  # Assuming tab-separated values
    header = next(reader)  # Skip the header row

    # Update counts for each filter
    for row in reader:
        counts['No_filter'][row[1]] += 1
        counts['Permissive_filter'][row[2]] += 1
        counts['Strict_filter'][row[3]] += 1

# Open the output file and write the results
with open(output_file_path, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')

    # Write the header
    writer.writerow(['Filter Type', 'RR', 'RD', 'DD', 'Total Samples', 'Total Alleles', 'MAF'])

    # Write the results for each filter type
    for filter_type, count_dict in counts.items():
        total_samples = sum(count_dict.values())
        RD_number = count_dict['RD']
        DD_number = count_dict['DD']
        total_alleles = 2 * total_samples
        MAF = (RD_number + DD_number) / total_alleles
        # Formatting MAF to display only four decimal places
        MAF_formatted = "{:.4f}".format(MAF)
        writer.writerow([filter_type, count_dict['RR'], RD_number, DD_number, total_samples, total_alleles, MAF_formatted])
        print("done")