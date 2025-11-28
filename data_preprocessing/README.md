# Preparation of Input Hi-C Files

For correct operation, this pipeline requires Hi-C files that have been filtered to remove blacklist regions (Ogata, J. D. et al., 2023). We recommend generating these files using the Juicer package, following the same procedure used in our study. Specifically, we used a modified version of Juicer available at:
https://github.com/genomech/juicer1.6_compact

# Requirements
We recommend reading the instructions in the official Juicer repository (https://github.com/aidenlab/juicer?ysclid=mijdkz1ume687884273).

The environment with the necessary dependencies for data preprocessing can be installed by running:
```bash
conda create -n juicer -c bioconda bwa samtools openjdk=17.0.8 python=3.10
```

# Running

## Step 1: Running Juicer on FASTQ Files
We processed the FASTQ files for each sample using the following command:

```bash
for SAMPLE in sample_1 sample_2 sample_n; do 
  bash juicer1.6_compact-main/scripts/juicer.sh \
    -g T2T \
    -d folder_for_${SAMPLE} \
    -s DpnII \
    -y T2T-CHM13v2.0_DpnII.txt \
    -p hgT2T.chromsizes \
    -z hgT2T.fa \
    -D juicer1.6_compact-main \
    -t 10
done
```

This step produces the merged_nodups.txt file for each sample.

## Step 2: Filtering Blacklist Regions
Next, the `merged_nodups.txt` files were filtered to remove blacklist regions using the `merged_nodups_filter.py` script:
```bash
for SAMPLE in sample_1 sample_2 sample_n ; do 
  python3 merged_nodups_filter.py \
    folder_for_${SAMPLE}/aligned/merged_nodups.txt \
    T2T.excluderanges_high_signal.bed \
    folder_for_${SAMPLE}/aligned/filtered_merged_nodups.txt
done
```

This step generates filtered_merged_nodups.txt for each sample.

## Step 3: Generating Final .hic Files
Finally, the filtered files were converted into `.hic` format using `juicer_tools.jar` from the Juicer package:
```bash
for SAMPLE in sample_1 sample_2 sample_n ; do 
  java -jar /data/amivanov/Tools/juicer1.6_compact-main/scripts/common/juicer_tools.jar pre \
    -q 30 \
    folder_for_${SAMPLE}/aligned/filtered_merged_nodups.txt \
    folder_for_${SAMPLE}/aligned/${SAMPLE}_filtered_inter_30.hic \
    hgT2T.chromsizes
done
```

The resulting `.hic` files are the final inputs used in the pipeline.
