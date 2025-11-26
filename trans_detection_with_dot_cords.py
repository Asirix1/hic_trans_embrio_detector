#!/usr/bin/env python3
"""
Translocation Detection Pipeline
Usage: python translocation_pipeline.py --hic_dir /path/to/hic_files --resolution 2500000 --output_dir /path/to/results
"""

import os
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm
import hicstraw
import itertools
import math
from sklearn.linear_model import LinearRegression
import glob

def calculate_trans_contacts(hic_file, resolution):
    """Calculate total trans contacts for a sample"""
    hic = hicstraw.HiCFile(hic_file)
    valid_chroms = [
        chrom.name for chrom in hic.getChromosomes()
        if chrom.name not in ["All", "chrX", "chrY", "chrM"]
    ]
    
    total_contacts = 0
    for chr1, chr2 in itertools.combinations(valid_chroms, 2):
        try:
            matrix_data = hic.getMatrixZoomData(
                chr1, chr2, "observed", "NONE", "BP", resolution
            )
            chr1_len = hic.getChromosome(chr1).length
            chr2_len = hic.getChromosome(chr2).length
            records = matrix_data.getRecordsAsMatrix(0, chr1_len, 0, chr2_len)
            total_contacts += np.nansum(records)
        except Exception as e:
            print(f"Warning: Error processing {chr1}-{chr2} in {os.path.basename(hic_file)}: {str(e)}")
    return total_contacts

def zero_values_above_percentile(data, percentile=95):
    """Replace values above percentile with zero"""
    if not data:
        return data
    threshold = np.percentile(data, percentile)
    return [0 if x >= threshold else x for x in data]

def calculate_triangle_sum(matrix, start_i, start_j, height, width, triangle_type):
    """Calculate sum of values in a triangular region of the matrix"""
    rows, cols = matrix.shape
    triangle_values = []
    
    if triangle_type == 1:  # Left-top triangle
        for row_offset in range(height):
            i = start_i - height + 1 + row_offset
            num_cols = min(row_offset + 1, width)
            for col_offset in range(num_cols):
                j = start_j - col_offset
                if 0 <= i < rows and 0 <= j < cols and not np.isnan(matrix[i, j]):
                    triangle_values.append(matrix[i, j])
    elif triangle_type == 2:  # Right-bottom triangle
        for row_offset in range(height):
            i = start_i + row_offset
            num_cols = max(0, width - row_offset)
            for col_offset in range(num_cols):
                j = start_j + col_offset
                if 0 <= i < rows and 0 <= j < cols and not np.isnan(matrix[i, j]):
                    triangle_values.append(matrix[i, j])
    elif triangle_type == 3:  # Left-bottom triangle
        for row_offset in range(height):
            i = start_i + row_offset
            cols_in_row = max(1, width - row_offset)
            start_col = start_j - cols_in_row + 1
            for col_offset in range(cols_in_row):
                j = start_col + col_offset
                if 0 <= i < rows and 0 <= j < cols and not np.isnan(matrix[i, j]):
                    triangle_values.append(matrix[i, j])
    elif triangle_type == 4:  # Right-top triangle
        for row_offset in range(height):
            i = start_i - height + 1 + row_offset
            num_cols = min(row_offset + 1, width)
            for col_offset in range(num_cols):
                j = start_j + col_offset
                if 0 <= i < rows and 0 <= j < cols and not np.isnan(matrix[i, j]):
                    triangle_values.append(matrix[i, j])
    
    triangle_values = zero_values_above_percentile(triangle_values)
    return np.sum(triangle_values)

def detect_triangles(hic_file, resolution, total_contacts, output_dir):
    """Detect triangles in Hi-C data and save results"""
    sample_name = os.path.basename(hic_file).replace('.hic', '')
    output_file = os.path.join(output_dir, f"{sample_name}_triangles.txt")
    
    # Skip if already processed
    if os.path.exists(output_file):
        return output_file
    
    hic = hicstraw.HiCFile(hic_file)
    valid_chroms = [
        chrom.name for chrom in hic.getChromosomes()
        if chrom.name not in ["All", "chrX", "chrY", "chrM"]
    ]
    chrom_pairs = list(itertools.combinations(valid_chroms, 2))
    
    results = []
    triangle_height = 10
    triangle_base = 10
    
    for chri, chrj in tqdm(chrom_pairs, desc=f"Processing {sample_name}"):
        try:
            matrix_data = hic.getMatrixZoomData(chri, chrj, "observed", "NONE", "BP", resolution)
            chri_len = hic.getChromosome(chri).length
            chrj_len = hic.getChromosome(chrj).length
            matrix = matrix_data.getRecordsAsMatrix(0, chri_len, 0, chrj_len)
            
            if matrix is None or matrix.size == 0:
                continue
                
            max_val = 0
            best_coords = None
            best_type = None
            
            # Find highest value triangle
            for i in range(matrix.shape[0]):
                for j in range(matrix.shape[1]):
                    for t_type in range(1, 5):
                        tri_sum = calculate_triangle_sum(
                            matrix, i, j, triangle_height, triangle_base, t_type
                        )
                        if tri_sum > max_val:
                            max_val = tri_sum
                            best_coords = (i, j)
                            best_type = t_type
            
            if best_coords and max_val > 0:
                i, j = best_coords
                start1 = max(0, i - 10) * resolution
                end1 = min(matrix.shape[0], i + 10) * resolution
                start2 = max(0, j - 10) * resolution
                end2 = min(matrix.shape[1], j + 10) * resolution
                
                # Normalize by total contacts
                norm_value = (max_val * 500000) / total_contacts if total_contacts > 0 else 0
                
                results.append({
                    'chrom1': chri,
                    'start1': start1,
                    'end1': end1,
                    'chrom2': chrj,
                    'start2': start2,
                    'end2': end2,
                    'value': norm_value,
                    'triangle_number': best_type
                })
        except Exception as e:
            print(f"Warning: Error processing {chri}-{chrj} in {sample_name}: {str(e)}")
    
    # Save results
    if results:
        df = pd.DataFrame(results)
        df = df.sort_values('value', ascending=False)
        df.to_csv(output_file, sep='\t', index=False)
    
    return output_file

def create_average_file(triangle_files, output_dir):
    """Create average triangle values file from sample files"""
    avg_file = os.path.join(output_dir, "average_triangles.txt")
    
    # Skip if already exists
    if os.path.exists(avg_file):
        return avg_file
    
    values_dict = {}
    
    for file in tqdm(triangle_files, desc="Calculating averages"):
        sample_name = os.path.basename(file).replace('_triangles.txt', '')
        df = pd.read_csv(file, sep='\t')
        
        for _, row in df.iterrows():
            key = (
                row['chrom1'], row['start1'], row['end1'],
                row['chrom2'], row['start2'], row['end2'],
                row['triangle_number']
            )
            if key not in values_dict:
                values_dict[key] = []
            values_dict[key].append(row['value'])
    
    # Calculate averages
    avg_results = []
    for key, values in tqdm(values_dict.items(), desc="Computing averages"):
        valid_values = [v for v in values if not math.isnan(v)]
        if valid_values:
            avg_value = sum(valid_values) / len(valid_values)
            avg_results.append(list(key) + [avg_value])
    
    # Save average file
    columns = [
        'chrom1', 'start1', 'end1',
        'chrom2', 'start2', 'end2',
        'triangle_number', 'value'
    ]
    avg_df = pd.DataFrame(avg_results, columns=columns)
    avg_df.to_csv(avg_file, sep='\t', index=False)
    
    return avg_file

def detect_translocations(sample_file, avg_file, output_dir):
    """Detect translocations using sample and average data"""
    sample_name = os.path.basename(sample_file).replace('_triangles.txt', '')
    output_file = os.path.join(output_dir, f"{sample_name}_translocations.txt")
    
    # Parameters optimized from original code
    params = {
        'offset': 0.1,
        'sig': 6.6,
        'slop_coef': 0.1,
        'chrom_ratio_threshold': 0.2,
        'min_points_threshold': 230
    }
    
    # Load data
    sample_df = pd.read_csv(sample_file, sep='\t')
    avg_df = pd.read_csv(avg_file, sep='\t')
    merged_df = pd.merge(
        sample_df, avg_df,
        on=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'triangle_number'],
        suffixes=('_sample', '_avg')
    )
    
    # Filter low values
    merged_df = merged_df[merged_df['value_sample'] >= 2].copy()
    if merged_df.empty:
        return []
    
    # Linear regression setup
    X = merged_df['value_sample'].values.reshape(-1, 1)
    y = merged_df['value_avg'].values
    reg = LinearRegression().fit(X, y)
    
    # Calculate thresholds
    x_99 = np.percentile(X, 99.5)
    clipped_X = np.clip(X, None, x_99)
    mean_x = np.mean(clipped_X)
    std_clipped = np.std(clipped_X)
    x_quantile = mean_x + params['sig'] * std_clipped
    y_at_x_quantile = reg.predict([[x_quantile]])[0] + params['offset']
    
    # Find potential translocations
    slope_line = reg.coef_[0] * params['slop_coef']
    y_threshold = slope_line * (merged_df['value_sample'] - x_quantile) + y_at_x_quantile
    potential_trans = merged_df[
        (merged_df['value_avg'] < y_threshold) & 
        (merged_df['value_sample'] > x_quantile)
    ]
    
    # Filter by chromosome pairs
    if not potential_trans.empty:
        chrom_counts = potential_trans.groupby(['chrom1', 'chrom2']).size()
        chrom_counts_filtered = chrom_counts[chrom_counts > params['min_points_threshold']]
        
        if not chrom_counts_filtered.empty:
            chrom_ratios = chrom_counts_filtered / len(potential_trans)
            significant_pairs = chrom_ratios[chrom_ratios > params['chrom_ratio_threshold']].index
            
            filtered_trans = potential_trans[
                potential_trans[['chrom1', 'chrom2']].apply(tuple, axis=1).isin(significant_pairs)
            ]
            
            # Get unique chromosome pairs
            detected_pairs = filtered_trans[['chrom1', 'chrom2']].drop_duplicates()
            detected_pairs['point_count'] = detected_pairs.apply(
                lambda row: chrom_counts_filtered.get((row['chrom1'], row['chrom2']), 0), 
                axis=1
            )
            
            # Save results
            detected_pairs.to_csv(output_file, sep='\t', index=False)
            return detected_pairs.to_dict('records')
    
    return []

def main():
    parser = argparse.ArgumentParser(description="Translocation detection pipeline")
    parser.add_argument("--hic_dir", required=True, help="Directory containing .hic files")
    parser.add_argument("--resolution", type=int, default=2500000, 
                        help="Resolution in base pairs (default: 2500000)")
    parser.add_argument("--output_dir", required=True, help="Output directory for results")
    args = parser.parse_args()
    
    # Create output directories
    os.makedirs(args.output_dir, exist_ok=True)
    triangle_dir = os.path.join(args.output_dir, "triangles")
    translocation_dir = os.path.join(args.output_dir, "translocations")
    os.makedirs(triangle_dir, exist_ok=True)
    os.makedirs(translocation_dir, exist_ok=True)
    
    # Get Hi-C files
    hic_files = glob.glob(os.path.join(args.hic_dir, "*.hic"))
    if not hic_files:
        raise ValueError(f"No .hic files found in {args.hic_dir}")
    
    print(f"Found {len(hic_files)} Hi-C files")
    
    # Step 1: Calculate trans contacts for each sample
    print("\nStep 1: Calculating trans contacts...")
    contacts_data = []
    for hic_file in tqdm(hic_files, desc="Calculating contacts"):
        sample_name = os.path.basename(hic_file).replace('.hic', '')
        total_contacts = calculate_trans_contacts(hic_file, args.resolution)
        contacts_data.append({'sample': sample_name, 'total_contacts': total_contacts})
    
    contacts_df = pd.DataFrame(contacts_data)
    contacts_file = os.path.join(args.output_dir, "sample_contacts.csv")
    contacts_df.to_csv(contacts_file, index=False)
    print(f"Saved contact data to {contacts_file}")
    
    # Step 2: Detect triangles for each sample
    print("\nStep 2: Detecting triangles...")
    triangle_files = []
    for hic_file in tqdm(hic_files, desc="Detecting triangles"):
        sample_name = os.path.basename(hic_file).replace('.hic', '')
        total_contacts = contacts_df[contacts_df['sample'] == sample_name]['total_contacts'].values[0]
        triangle_file = detect_triangles(hic_file, args.resolution, total_contacts, triangle_dir)
        triangle_files.append(triangle_file)
    
    # Step 3: Create average triangle file
    print("\nStep 3: Creating average triangle file...")
    avg_file = create_average_file(triangle_files, args.output_dir)
    print(f"Created average file: {avg_file}")
    
    # Step 4: Detect translocations
    print("\nStep 4: Detecting translocations...")
    all_translocations = []
    
    for triangle_file in tqdm(triangle_files, desc="Calling translocations"):
        sample_name = os.path.basename(triangle_file).replace('_triangles.txt', '')
        translocations = detect_translocations(triangle_file, avg_file, translocation_dir)
        
        for trans in translocations:
            all_translocations.append({
                'sample': sample_name,
                'chrom1': trans['chrom1'],
                'chrom2': trans['chrom2'],
                'point_count': trans['point_count']
            })
    
    # Save final results
    if all_translocations:
        results_df = pd.DataFrame(all_translocations)
        results_file = os.path.join(args.output_dir, "detected_translocations.txt")
        results_df.to_csv(results_file, sep='\t', index=False)
        print(f"\nTranslocation detection complete! Results saved to {results_file}")
        print(f"Found {len(all_translocations)} translocations across {len(hic_files)} samples")
    else:
        print("\nNo translocations detected in any samples")

if __name__ == "__main__":
    main()