import os
import csv
import argparse

def generate_tsv(root_dir, output_tsv='sample_sheet.tsv', lab='ucl'):
    rows = []

    for subdir in sorted(os.listdir(root_dir)):
        print(subdir)
        subdir_path = os.path.join(root_dir, subdir)
        if not os.path.isdir(subdir_path):
            continue

        fq1 = fq2 = None
        for file in os.listdir(subdir_path):
            file_path = os.path.join(subdir_path, file)
            if 'R1' in file and os.path.isfile(file_path):
                fq1 = file
            elif 'R2' in file and os.path.isfile(file_path):
                fq2 = file

        if fq1 and fq2:
            rows.append([subdir, fq1, fq2, lab])
        else:
            print(f"[WARNING] Skipping {subdir}: R1 or R2 file not found.")

    with open(output_tsv, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['SAMPLE_NAME', 'FQ1', 'FQ2', 'LAB'])
        writer.writerows(rows)

    print(f"[INFO] TSV file saved as: {output_tsv}")

def main():
    parser = argparse.ArgumentParser(description="Generate TSV sample sheet from subdirectories.")
    parser.add_argument('root_dir', help='Root directory containing sample subdirectories')
    parser.add_argument('--output', default='sample_sheet.tsv', help='Output TSV filename')
    parser.add_argument('--lab', default='cordlife', help='Lab name (default: cordlife)')

    args = parser.parse_args()
    generate_tsv(args.root_dir, args.output, args.lab)

if __name__ == '__main__':
    main()
