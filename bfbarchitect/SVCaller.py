import pysam
from collections import defaultdict
import re
import argparse

try:
    from bfbarchitect.datatypes import CigarAlignment, SV, REVERSE_STRAND, CHR_TO_IDX
    from bfbarchitect.utils import get_normal_coverage, get_coverage_and_rc
except:
    from datatypes import CigarAlignment, SV, REVERSE_STRAND, CHR_TO_IDX
    from utils import get_normal_coverage, get_coverage_and_rc

def query_ends_from_cigar(cigar_str: str, strand: str) -> tuple[int, int, int]:
    """
    Retrieve alignment ends and alignment length from cigar string

    Args:
        cigar_str: CIGAR string
        strand: Strand of the read ("+" or "-")

    Returns:
        The start and end position on the read on positive strand, and the alignment length on
        the reference genome.
    """
    # Consumable ops on the reference genome (M, D, N, =, X)
    ref_consumable_ops = {"M", "D", "N", "=", "X"}
    query_consumable_ops = {"M", "I", "=", "X"}  # Consumable ops on the query

    query_start = 0
    query_consumed = 0
    ref_consumed = 0  # This will track the length on the reference genome
    
    # Parse the CIGAR string using regex to extract (length, operation) tuples
    cigar = re.findall(r"(\d+)([A-Z=])", cigar_str)

    # If the strand is negative, reverse the CIGAR operations
    if strand == "-":
        cigar = cigar[::-1]

    # Process CIGAR to find query and reference alignment lengths
    for length, cigar_op in cigar:
        length = int(length)  # Convert length to an integer

        # Handle soft/hard clipping which affects query_start but not reference
        if cigar_op == "S" or cigar_op == "H":  # Soft/Hard clip
            if query_consumed == 0:  # Before any alignment operation
                query_start += length

        # Handle query alignment (including matches and insertions)
        if cigar_op in query_consumable_ops:
            query_consumed += length

        # Handle reference alignment (e.g., matches, deletions, skipped regions)
        if cigar_op in ref_consumable_ops:
            ref_consumed += length

    query_end = query_start + query_consumed
    ref_alignment_length = ref_consumed  # Reference length consumed

    return (query_start, query_end, ref_alignment_length)

def cluster_SVs(SVs, max_bp_distance = 100):
    # sort all foldbacks by # support reads
    sv_list = sorted(SVs.items(), key=lambda item: len(item[1]), reverse=True)
    clustered_sv = {}
    for i in range(len(sv_list)):
        sv1, read_list1 = sv_list[i]
        if len(read_list1) == 0:
            continue
        clustered_sv[sv1] = read_list1
        for j in range(i+1, len(sv_list)):
            sv2, read_list2 = sv_list[j]
            if sv1.is_equal(sv2, max_bp_distance = max_bp_distance) == True:
                clustered_sv[sv1] += read_list2
                sv_list[j][1].clear()
    clustered_sv = dict(sorted(clustered_sv.items(), key=lambda item: (CHR_TO_IDX[item[0].chrom1], item[0].bp1)))
    return clustered_sv

def call_SVs(bam_file, region, min_mapq=20, min_ref_length=100, output_fn=None, normal_cov=10, min_cn=0.75):
    """Find SVs in a specific region"""
    bam = pysam.AlignmentFile(bam_file, 'rb')
    reads = bam.fetch() if region == None else bam.fetch(region[0], region[1], region[2])
    alignments: dict[str, list[CigarAlignment]] = defaultdict(list)

    # Parse primary and supplementary alignments
    for read in reads:
        # Skip if no supplementary alignments or low quality
        if not read.has_tag('SA'):
            continue
        
        # Parse primary alignment
        primary_chrom = read.reference_name
        primary_pos = read.reference_start
        primary_strand = '-' if read.is_reverse else '+'
        (read_start, read_end, ref_length) = query_ends_from_cigar(read.cigarstring, primary_strand)
            
        if primary_strand == "+":
            primary_start = primary_pos
            primary_end = primary_pos + ref_length - 1
        else:
            primary_start = primary_pos + ref_length - 1
            primary_end = primary_pos
        mapping_quality, edit_dist = float(read.mapping_quality), float(read.get_tag("NM"))/(read_end-read_start)
        if primary_chrom in CHR_TO_IDX and mapping_quality >= min_mapq and ref_length >= min_ref_length:
            alignment = CigarAlignment(primary_chrom, primary_start, primary_end, primary_strand, \
                                   ref_length, read.query_name, read_start, read_end, mapping_quality, edit_dist)
            if alignment not in alignments[read.query_name]:
                alignments[read.query_name].append(alignment)
        # Parse supplementary alignments
        for sa_tag in read.get_tag('SA').split(';'):
            if not sa_tag:
                continue
                
            sa_info = sa_tag.split(',')
            sa_chrom = sa_info[0]
            sa_pos = int(sa_info[1]) - 1 # 1-based to 0-based
            sa_strand = sa_info[2]
            sa_cigar = sa_info[3]
            (read_start, read_end, ref_length) = query_ends_from_cigar(sa_cigar, sa_strand)
            mapping_quality, edit_dist = float(sa_info[4]), float(sa_info[-1])/(read_end-read_start)
            if sa_chrom not in CHR_TO_IDX or mapping_quality < min_mapq or ref_length < min_ref_length:
                continue
            if sa_strand == "+":
                sa_start = sa_pos
                sa_end = sa_pos + ref_length - 1
            else:
                sa_start = sa_pos + ref_length - 1
                sa_end = sa_pos
            alignment = CigarAlignment(sa_chrom, sa_start, sa_end, sa_strand, \
                                       ref_length, read.query_name, read_start, read_end, mapping_quality, edit_dist)
            if alignment not in alignments[read.query_name]:
                alignments[read.query_name].append(alignment)
    # Get SVs from read alignments
    SVs: dict[SV, list] = defaultdict(list)
    TST_SVs = []
    for read_name, alignment_list in alignments.items():
        alignment_list.sort(key=lambda a: a.read_start)
        supported_SVs = []
        for i in range(len(alignment_list)-1):
            j = i + 1
            alignment1, alignment2 = alignment_list[i], alignment_list[j]
            query_gap = alignment2.read_start - alignment1.read_end
            mapping_quality = min(alignment1.mapping_quality, alignment2.mapping_quality)

            chrom1, chrom2 = alignment1.chrom, alignment2.chrom
            bp1, bp2 = alignment1.end, alignment2.start
            strand1, strand2 = alignment1.strand, REVERSE_STRAND[alignment2.strand]
            sv = SV(chrom1, bp1, strand1, chrom2, bp2, strand2)
            supported_SVs.append(sv)
            SVs[sv].append((read_name, query_gap, mapping_quality))
        if len(supported_SVs) > 1:
            alignment1, alignment2 = alignment_list[0], alignment_list[-1]
            query_gap = alignment2.read_start - alignment1.read_end
            mapping_quality = min(alignment1.mapping_quality, alignment2.mapping_quality)
            
            # if -1000 < query_gap and query_gap < 1000:
            chrom1, chrom2 = alignment1.chrom, alignment2.chrom
            bp1, bp2 = alignment1.end, alignment2.start
            strand1, strand2 = alignment1.strand, REVERSE_STRAND[alignment2.strand]
            TST_sv = SV(chrom1, bp1, strand1, chrom2, bp2, strand2)
            TST_sv.TST = True
            SVs[TST_sv].append((read_name, query_gap, mapping_quality))
            TST_SVs.append(TST_sv)
            
    # Cluster foldback inversions
    SVs = cluster_SVs(SVs)
    for sv, read_list in SVs.items():
        if sv in TST_SVs:
            sv.TST = True
        # remove duplicate read_names
        # unique_reads = defaultdict(list)
        # for read, gap in read_list:
        #     unique_reads[read].append(gap)
        # SVs[sv] = [(read, round(sum(gaps)/len(gaps))) for read, gaps in unique_reads.items()]
    # Print results
    if output_fn != None:
        output_file = open(output_fn, 'a')
        output_file.write('SV\tSV_type\tTST\t#Support_reads\tMapping_quality\tQuery_gaps\tSupport_reads\n')
        for sv, info in SVs.items():
            if 2*len(info)/normal_cov < min_cn:
                continue
            reads, query_gaps, mapping_qualities = [read for (read, _, _) in info], [str(gap) for (_, gap, _) in info], [str(mq) for (_, _, mq) in info]
            output_file.write(f"{sv}\t{sv.type}\t{sv.TST}\t{len(info)}\t{','.join(mapping_qualities)}\t{','.join(query_gaps)}\t{','.join(reads)}\n")
        output_file.close()
    output_inversions = {}
    for sv, info in SVs.items():
        if 2*len(info)/normal_cov < min_cn:
            continue
        output_inversions[sv] = len(info)
    return output_inversions

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Call foldback inversions from bam.")
    parser.add_argument("--bam", help = "Path to a sorted bam file", required = True)
    parser.add_argument("--cns", help = "Path to a cns file for normal coverage estimation", required = True)
    parser.add_argument("--region", help = "Call SVs from a specific region (e.g. chr1:1000000-2000000)", required = True)
    parser.add_argument("--output_prefix", help = "Output file name", required = True)
    parser.add_argument("--min_sv_cn", type=float, default=0.75, help="Minimum copy number for SV calling (default: 0.75)")
    parser.add_argument("--min_mapq", type=int, default=20, help="Minimum mapping quality for SV calling (default: 20)")
    args = parser.parse_args()

    normal_cov = get_normal_coverage(args.cns, args.bam)

    region = None
    if args.region != None:
        chrom = args.region.split(':')[0]
        start = int(args.region.split(':')[1].split('-')[0])
        end = int(args.region.split('-')[1])
        region = (chrom, start, end)

    output_read_fn = None if args.output_prefix == None else f'{args.output_prefix}_reads.txt'
    if output_read_fn != None:
        output_file = open(output_read_fn, 'w')
        output_file.close()
    SVs = call_SVs(args.bam, region, normal_cov=normal_cov, output_fn=output_read_fn, min_cn=args.min_sv_cn, min_mapq=args.min_mapq)
    print(f'Saved structural variants to {args.output_prefix}_reads.txt.')