# Translocation Detection Pipeline

This repository contains a Python-based pipeline for detecting chromosomal translocations from Hi-C data. The pipeline identifies potential translocations by analyzing contact matrices and comparing sample data against an average reference.

## Overview

Chromosomal translocations are genetic abnormalities where segments of chromosomes are transferred between non-homologous chromosomes. This pipeline uses Hi-C data to identify such translocations by:

1. Processing Hi-C contact matrices to detect triangular patterns
2. Comparing sample data against an average reference
3. Identifying statistically significant deviations that indicate potential translocations

## Features

- Process multiple Hi-C files in batch
- Calculate trans contacts for normalization
- Detect triangular patterns in contact matrices
- Create average reference from all samples
- Identify statistically significant translocations
- Output results in tab-separated format

## Requirements

- Python 3.7+
- hicstraw
- pandas
- numpy
- scikit-learn
- tqdm

Install dependencies using pip:

```bash
pip install hicstraw pandas numpy scikit-learn tqdm
```

## Usage

Run the pipeline with the following command:

```bash
python trans_detection_with_dot_cords.py --hic_dir /path/to/hic_files --resolution 2500000 --output_dir /path/to/results
```

### Arguments

- `--hic_dir`: Directory containing .hic files (required)
- `--resolution`: Resolution in base pairs (default: 2500000)
- `--output_dir`: Output directory for results (required)

### Output Structure

The pipeline creates the following output structure:

```
output_dir/
├── sample_contacts.csv          # Total trans contacts for each sample
├── average_triangles.txt        # Average triangle values across samples
├── triangles/                   # Individual triangle detection results per sample
│   └── sample_name_triangles.txt
└── translocations/              # Individual translocation results per sample
    └── sample_name_translocations.txt
```

## Pipeline Steps

1. **Trans Contact Calculation**: Calculate total trans contacts for each Hi-C sample
2. **Triangle Detection**: Identify triangular patterns in contact matrices for each sample
3. **Average Reference Creation**: Create an average reference from all samples
4. **Translocation Calling**: Compare each sample against the average to detect significant translocations

## Output Format

The final translocation results file (`detected_translocations.txt`) contains columns:
- `sample`: Sample name
- `chrom1`: First chromosome in the translocation
- `chrom2`: Second chromosome in the translocation
- `point_count`: Number of significant data points supporting the translocation

## License

[Specify your license here if applicable]

## Citation

[Include citation information if this is based on published work]