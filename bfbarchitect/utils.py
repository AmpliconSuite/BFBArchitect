import logging
import pandas as pd
import pysam
import re

def create_logger(name, log_file):
    """Create a logger"""
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    
    logger.handlers.clear()
    
    # Create file handler
    handler = logging.FileHandler(log_file, mode='w')
    handler.setLevel(logging.DEBUG)
    
    # Create formatter
    formatter = logging.Formatter('[%(name)s:%(levelname)s]\t%(message)s')
    handler.setFormatter(formatter)
    
    logger.addHandler(handler)
    return logger

def get_normal_coverage(cns_fn, bam_fn):
    # Get normal genome regions
    cns = pd.read_csv(cns_fn, sep="\t")
    sex_chr = re.compile(r"^(chr)?[XY]$", re.IGNORECASE)
    cns = cns[~cns.chromosome.apply(lambda x: bool(sex_chr.match(x)))]
    segments = cns.sort_values(by='log2').reset_index(drop=True) # sort all segments by log2
    l = int(len(segments) / 2.4)
    r = l + 1
    total_length = 0
    log2_cn = []
    regions = []
    while total_length < 10_000_000:
        l_chr, l_start, l_end = segments.loc[l].chromosome, int(segments.loc[l].start), int(segments.loc[l].end)
        regions.append((l_chr, l_start, l_end))
        total_length += l_end - l_start
        log2_cn.append(segments.loc[l].log2)
        r_chr, r_start, r_end = segments.loc[r].chromosome, int(segments.loc[r].start), int(segments.loc[r].end)
        regions.append((r_chr, r_start, r_end))
        total_length += r_end - r_start
        log2_cn.append(segments.loc[r].log2)
        l -= 1
        r += 1
    # Get normal coverage
    total_length, total_bases = 0, 0
    for region in regions:
        (coverage, _) = get_coverage_and_rc(bam_fn, region)
        region_length = region[2] - region[1] + 1
        total_bases += coverage * region_length
        total_length += region_length
    normal_cov = total_bases / total_length if total_length > 0 else 0
    return normal_cov

def get_chrom_length(bam_fn, chrom):
    """Return the length of chrom from the BAM header, or None if not present."""
    with pysam.AlignmentFile(bam_fn, 'rb') as bam:
        lengths = dict(zip(bam.references, bam.lengths))
    return lengths.get(chrom)


def get_coverage_and_rc(bam_fn, interval, qc_threshold=0):
    total_length, total_bases = 0, 0
    bam = pysam.AlignmentFile(bam_fn, "rb")
    read_count = 0
    chrom, start, end = interval
    total_length += end - start + 1
    for read in bam.fetch(chrom, start, end):
        if read.mapping_quality < qc_threshold or read.seq == None:
            continue
        read_count += 1
        for block_start, block_end in read.get_blocks():
            if block_end < start or block_start > end:
                continue
            total_bases += min(block_end, end) - max(block_start, start)
    coverage = total_bases / total_length
    return (coverage, read_count)