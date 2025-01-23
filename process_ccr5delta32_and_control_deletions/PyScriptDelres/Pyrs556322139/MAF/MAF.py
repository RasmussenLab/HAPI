import csv

DelSNPname = "rs556322139"

# File paths using the variable
input_file_path = rf'C:\Users\k_rav\OneDrive\Skrivebord\Speciale\1NOV22120\Cell\Cell_reviewers_comments\Deletions_background\Results_background\New_files150724\PyScriptDelres\Py{DelSNPname}\Count.top4_{DelSNPname}.tsv'
output_file_path = r'MAF/resMAF.tsv'

# Initialize the counters for each data group and filter type
counts = {}

# Open and read the input file
try:
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
except FileNotFoundError:
    print(f"Error: The file {input_file_path} does not exist.")
    exit()
except Exception as e:
    print(f"An error occurred while reading the file: {e}")
    exit()

# Open the output file and write the results
try:
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

    print(f"Output saved to {output_file_path}")
except Exception as e:
    print(f"An error occurred while writing the file: {e}")
