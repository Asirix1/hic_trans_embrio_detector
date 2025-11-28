import sys

def load_bed_regions(bed_file):
    from collections import defaultdict
    regions = defaultdict(list)
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end = line.strip().split()[:3]
            regions[chrom].append((int(start), int(end)))
    return regions

def is_in_regions(chrom, pos, regions):
    if chrom not in regions:
        return False
    for start, end in regions[chrom]:
        if start <= pos <= end:
            return True
    return False

def filter_nodups(nodups_file, bed_regions, output_file):
    with open(nodups_file) as f_in, open(output_file, "w") as f_out:
        for line in f_in:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split()
            chr1, pos1 = fields[1], int(fields[2])
            chr2, pos2 = fields[5], int(fields[6])
            if not is_in_regions(chr1, pos1, bed_regions) and not is_in_regions(chr2, pos2, bed_regions):
                f_out.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_nodups_by_bed.py merged_nodups regions.bed output.txt")
        sys.exit(1)

    nodups_file = sys.argv[1]
    bed_file = sys.argv[2]
    output_file = sys.argv[3]

    bed_regions = load_bed_regions(bed_file)
    filter_nodups(nodups_file, bed_regions, output_file)
