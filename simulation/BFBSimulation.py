import argparse
import random
import pysam

def simulate_BFB_sequence(n_segments, n_foldbacks, is_reverse, no_rep):
    BFB_sequence = []
    if is_reverse:
        for seg in range(n_segments, 0, -1):
            BFB_sequence.append(-seg)
    else:
        for seg in range(1, n_segments + 1):
            BFB_sequence.append(seg)
    deletions = [0 for _ in range(n_segments)]
    deletions[random.randint(0, n_segments - 1)] = 1
    foldback_count = 0
    breakpoints = [BFB_sequence[-1]]
    while foldback_count < n_foldbacks and len(breakpoints) < 2*(n_segments - 1):
        new_BFB_sequence, new_deletions = BFBOneCycle(BFB_sequence, deletions)
        breakpoint = new_BFB_sequence[-1]
        if no_rep and breakpoint in breakpoints: # If foldback is duplicate, rerun one BFB cycle
            continue

        breakpoints.append(breakpoint)
        if breakpoint > 0:
            breakpoints.append(-(breakpoint+1))
        else:
            breakpoints.append(-breakpoint-1)
        BFB_sequence = new_BFB_sequence
        deletions = new_deletions
        _, L, R = count_vectors(n_segments, BFB_sequence)
        foldback_count = sum(L) + sum(R)
    return BFB_sequence, deletions

def BFBOneCycle(input_sequence, deletions):
    breakpoint = random.randint(1, len(input_sequence)-1)
    duplicate = input_sequence[breakpoint:len(input_sequence)]
    new_sequence = input_sequence + ReverseComplementSeq(duplicate)

    duplicate = deletions[breakpoint:len(input_sequence)][::-1]
    for i in range(len(duplicate)):
        if duplicate[i] > 0:
            duplicate[i] += 1
    new_deletions = deletions + duplicate
    return new_sequence, new_deletions

def ReverseComplementSeq(sequence):
    complement = [-seg for seg in sequence]
    reverse_complement = complement[::-1]
    return reverse_complement

def count_vectors(n_segments, amplicon_sequence):
    C, L, R = [0] * n_segments, [0] * n_segments, [0] * n_segments
    for seg1, seg2 in zip(amplicon_sequence[:-1], amplicon_sequence[1:]):
        idx = abs(seg1) - 1
        C[idx] += 1
        if seg1 * seg2 < 0 and abs(seg1) == abs(seg2):
            if seg1 > 0:
                R[idx] += 1
            else:
                L[idx] += 1
    C[abs(amplicon_sequence[-1]) - 1] += 1
    return C, L, R

def get_segment_children(BFB_sequence, i, children):
    parent = BFB_sequence[i]
    for j in range(i + 1, len(BFB_sequence)):
        child = BFB_sequence[j]
        mid = (i + j) // 2 + 1
        if parent == -child and BFB_sequence[i:mid] == ReverseComplementSeq(BFB_sequence[mid:j+1]):
            children[i].append(j)

def generate_segment_coordinates(chr, start, end, n, d):  
    if (end - start) < (n - 1) * (d + 1):  
        raise ValueError("Impossible to generate n integers with given constraints")  
    
    result = []  
    # Calculate minimum total span needed  
    min_span = (n - 1) * (d + 1)  
    
    # Start with random first position with enough room at the end  
    start_max = end - min_span  
    if start_max < start:  
        raise ValueError("Range too small for constraints")  
    
    current = random.randint(start, start_max)  
    result.append(current)  
    
    for _ in range(n - 1):  
        # Next number must be at least current + d + 1  
        # And leave enough room for remaining numbers  
        remaining = n - len(result)  
        min_next = current + d + 1  
        max_next = end - (remaining - 1) * (d + 1)  
        
        if min_next > max_next:  
            # Backtrack or raise error  
            raise ValueError("Cannot satisfy constraints with current random choice")  
        
        current = random.randint(min_next, max_next)  
        result.append(current)
    
    result.insert(0, start)
    result[-1] = end
    seg_coords = [(chr, result[i], result[i + 1]) if i == 0 \
                  else (chr, result[i] + 1, result[i + 1]) for i in range(len(result) - 1)]
    return seg_coords

def GetSequencesFromGenome(seg_coords, fasta_file_path):
    ref = pysam.FastaFile(fasta_file_path)
    res = {}
    for i, seg in enumerate(seg_coords):
        print(seg)
        seg_id = i + 1
        res[seg_id] = ref.fetch(seg[0], seg[1], seg[2])
        res[-seg_id] = ReverseComplementStr(res[seg_id])
    ref.close()
    return res

def ReverseComplementStr(dna_str):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 
                  'N': 'N', 'n': 'n'}
    # Complement the sequence
    complemented_sequence = ''.join([complement[base] for base in dna_str])
    # Reverse the complementary sequence
    reverse_complemented_sequence = complemented_sequence[::-1]
    return reverse_complemented_sequence

def generate_graph_file(output_fn, segments, CNs, SVs):
    # Generate graph.txt
    out_file = open(output_fn, 'w')
    out_file.write('SequenceEdge: StartPosition, EndPosition, PredictedCN, AverageCoverage, Size, NumberOfLongReads\n')
    for i, seg in enumerate(segments):
        size = seg[2] - seg[1] + 1
        entry = f'sequence	{seg[0]}:{seg[1]}-	{seg[0]}:{seg[2]}+	{CNs[i]}	0	{size}	0\n'
        out_file.write(entry)
    out_file.write('BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfLongReads\n')
    for i in range(len(segments) - 1):
        seg1, seg2 = segments[i], segments[i + 1]
        cn = min(CNs[i], CNs[i + 1])
        entry = f'concordant	{seg1[0]}:{seg1[2]}+->{seg2[0]}:{seg2[1]}-	{cn}	0\n'
        out_file.write(entry)
    for sv, cn in SVs.items():
        entry = f'discordant	{sv}	{cn}	0\n'
        out_file.write(entry)
    out_file.close()

def generate_cycle_file(output_fn, segments, BFB_sequence, length):
    # Generate cycle.txt
    out_file = open(output_fn, 'w')
    out_file.write(f'Interval	1	{segments[0][0]}	{segments[0][1]}	{segments[-1][2]}\n')
    out_file.write('List of cycle segments\n')
    for i, segment in enumerate(segments):
        out_file.write(f'Segment	{i+1}	{segment[0]}	{segment[1]}	{segment[2]}\n')
    sequence = ','.join([str(seg) for seg in BFB_sequence])
    out_file.write(f'Path=1;Copy_count=1.0;Length={length};Segments={sequence}\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Simulate a BFB sequence.")
    parser.add_argument("--reference", help = "Reference genome fasta file.", required = True)
    parser.add_argument("--regions", help = "BFB region (e.g. chr1:1000000-2000000", required = True)
    parser.add_argument("--output_prefix", help = "Prefix of output files.", required = True)
    parser.add_argument("--sequence", help = "Optional: specify BFB sequence.", required = False)
    args = parser.parse_args()

    # simulate a BFB sequence with deletions
    n_segments = random.randint(5, 12)
    n_foldbacks = random.randint(5, 15)
    is_reverse = False
    no_rep = True

    BFB_sequence, deletions = simulate_BFB_sequence(n_segments, n_foldbacks, is_reverse, no_rep)

    n_deletions = 8
    level = 1
    deletion_indices = [i for i in range(len(deletions)) if deletions[i] >= level]
    while len(deletion_indices) > n_deletions:
        level += 1
        deletion_indices = [i for i in range(len(deletions)) if deletions[i] >= level]
    
    if args.sequence:
        BFB_sequence = [int(seg[:-1]) if seg[-1] == '+' else -int(seg[:-1]) for seg in args.sequence.split(',')]
        n_segments = max([abs(seg) for seg in BFB_sequence])
        deletion_indices = []
    print('BFB sequence', BFB_sequence)


    # simulate segment coordinates
    chrom = args.regions.split(':')[0]
    start, end = int(args.regions.split(':')[1].split('-')[0]), int(args.regions.split('-')[1])
    seg_coords = generate_segment_coordinates(chrom, start, end, n_segments, 100000)
    print('Segment coordinates:', seg_coords)
    print(n_segments, len(seg_coords))

    # generate SVs
    foldback_length = [random.randint(10, 2000) for _ in range(n_segments)]
    SVs = {}
    for i in range(len(BFB_sequence) - 1):
        seg1 = BFB_sequence[i]
        seg2 = BFB_sequence[i + 1]
        if seg1 == -seg2: # foldback inversion
            i = abs(seg1) - 1
            chrom, start, end = seg_coords[i]
            if seg1 > 0:
                sv = f'{chrom}:{end-foldback_length[i]}+->{chrom}:{end}+'
            else:
                sv = f'{chrom}:{start}-->{chrom}:{start+foldback_length[i]}-'
            if sv not in SVs:
                SVs[sv] = 0
            SVs[sv] += 1
    
    if args.sequence is None:
        del_idx = abs(BFB_sequence[deletion_indices[0]]) - 1
        chrom, start, end = seg_coords[del_idx]
        deletion_length = random.randint((end - start) // 4, (end - start) // 2)
        start = random.randint(start + 3000, end - deletion_length)
        end = start + deletion_length
        deletion_sv = f'{chrom}:{start}+->{chrom}:{end}-'
        SVs[deletion_sv] = len(deletion_indices)

    # generate fasta
    segment_strs = GetSequencesFromGenome(seg_coords, args.reference)
    if args.sequence is None:
        segment_strs[f'{del_idx+1}_del1'] = segment_strs[del_idx+1][:start - seg_coords[del_idx][1]+1] + \
                                    segment_strs[del_idx+1][end - seg_coords[del_idx][1]:]
        segment_strs[f'{del_idx+1}_del2'] = ReverseComplementStr(segment_strs[f'{del_idx+1}_del1'])
    with open(f'{args.output_prefix}.fa', 'w') as fasta_out:
        final_str = ""
        for i, seg in enumerate(BFB_sequence):
            segment_str = ""
            if i in deletion_indices:
                if seg > 0:
                    segment_str = segment_strs[f'{abs(seg)}_del1']
                else:
                    segment_str = segment_strs[f'{abs(seg)}_del2']
            else:
                segment_str = segment_strs[seg]
            if i > 0 and seg == -BFB_sequence[i - 1]:
                print('Simulate deletion in foldback:', seg, BFB_sequence[i - 1], foldback_length[abs(seg) - 1])
                segment_str = segment_str[foldback_length[abs(seg) - 1]:]
            final_str += segment_str
        fasta_out.write(f">{args.output_prefix}" + "\t" + "Length:"+str(len(final_str))+"\n")
        fasta_out.write(final_str)
        fasta_out.write("\n")
        fasta_out.close()
    print(BFB_sequence)
    # generate graph and cycle files
    C, L, R = count_vectors(n_segments, BFB_sequence)
    print('Count vectors', C, L, R)
    generate_graph_file(f'{args.output_prefix}_graph.txt', seg_coords, C, SVs) 
    generate_cycle_file(f'{args.output_prefix}_cycles.txt', seg_coords, BFB_sequence, -1)
