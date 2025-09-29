# Genome-Wide Protein Physicochemical Properties Analysis Tool

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python tool for the high-throughput computation of protein physicochemical properties directly from genome-scale FASTA files. This utility is designed for proteome analysis, functional genome annotation, and to support structural bioinformatics pipelines.

-----

## Key Features

- **FASTA Input**: Efficiently reads and processes multi-protein FASTA files.
- **Comprehensive Analysis**: Calculates key physicochemical properties for each protein sequence using Biopython's `ProtParam` module:
    - Molecular Weight
    - Theoretical Isoelectric Point (pI)
    - Instability Index
    - Aliphatic Index
    - Grand Average of Hydropathicity (GRAVY)
    - Net Charge at neutral pH
- **Flexible Output**: Saves the results in a clean, tabular format, supporting `.csv`, `.tsv`, and Excel (`.xlsx`).
- **Handles empty sequences gracefully with warnings**
- **Validates input files and output formats**

-----

## Requirements

- Python 3.x
- Biopython
- Pandas
- Openpyxl (optional, for `.xlsx` export functionality)

You can install all the necessary libraries with a single command:

```bash
pip install biopython pandas openpyxl
```

-----

## Usage

Run the script from your terminal, providing the input FASTA file and the desired output filename.

```bash
python protein_properties.py <input_fasta_file> <output_file.csv>
```

**Example:**

```bash
python protein_properties.py proteome.fasta results.csv
```

### Example Output

The script generates a file with the following structure:

| Gene ID  | Length | MW (kDa) | pI   | Instability | Aliphatic Index | GRAVY  | Charge at pH7 |
| :------- | :----- | :------- | :--- | :---------- | :-------------- | :----- | :------------ |
| Protein1 | 350    | 38.92    | 6.77 | 42.13       | 88.29           | -0.12  | 1.00          |
| Protein2 | 280    | 29.72    | 5.49 | 33.67       | 76.25           | 0.05   | -1.00         |

-----

## Applications

This tool is useful for a variety of bioinformatics tasks, including:

- Genome-wide protein annotation and characterization
- Comparative proteomics and genomics
- Proteome-wide clustering and classification based on physical properties
- Data preparation for structural bioinformatics and machine learning pipelines

-----

## Project Status

This project is currently a work in progress.

- ✅ **Core Script**: Main analysis and file I/O functionality is complete.
- 🔄 **Upcoming Improvements**:
    - Basic data visualization plots.
    - Visualization of protein property distributions.

-----

## License

This project is licensed under the **MIT License**. You are free to use, modify, and distribute it. See the `LICENSE` file for more details.
