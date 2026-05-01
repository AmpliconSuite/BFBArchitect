# Map a strand to its opposite strand
REVERSE_STRAND = {"+": "-", "-": "+"}


def chrom_sort_key(chrom: str) -> tuple:
    """Return a sort key for a chromosome name, prefix-agnostic.

    Works for any organism: numeric autosomes sort by number, sex/mito
    chromosomes sort after, unrecognised contigs sort last.
    Handles both 'chr7' and '7' identically.
    """
    name = chrom[3:] if chrom.startswith('chr') else chrom
    _tail = {'X': 9000, 'Y': 9001, 'M': 9002, 'MT': 9002}
    if name in _tail:
        return (_tail[name],)
    try:
        return (int(name),)
    except ValueError:
        return (99999,)


def chrom_in_dict(d: dict, chrom: str) -> bool:
    """Return True if chrom (or its chr-prefix-toggled form) is a key in d.

    Use this to check chromosome validity against reference data rather than
    maintaining a separate hard-coded list of valid names.
    """
    if chrom in d:
        return True
    alt = chrom[3:] if chrom.startswith('chr') else 'chr' + chrom
    return alt in d

CHR_SIZES = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
}

CHR_CENTRO = {
    "chr1": 122026459,
    "chr2": 92188145,
    "chr3": 90772458,
    "chr4": 49712061,
    "chr5": 46485900,
    "chr6": 58553888,
    "chr7": 58169653,
    "chr8": 44033744,
    "chr9": 43389635,
    "chr10": 39686682,
    "chr11": 51078348,
    "chr12": 34769407,
    "chr13": 16000000,
    "chr14": 16000000,
    "chr15": 17083673,
    "chr16": 36311158,
    "chr17": 22813679,
    "chr18": 15460899,
    "chr19": 24498980,
    "chr20": 26436232,
    "chr21": 10864560,
    "chr22": 12954788,
    "chrX": 58605579,
    "chrY": 10316944,
}


def load_centromere_bed(bed_path):
    """Parse a centromere BED file; return {chrom: midpoint}."""
    from collections import defaultdict
    spans = defaultdict(lambda: [float('inf'), 0])
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            spans[chrom][0] = min(spans[chrom][0], start)
            spans[chrom][1] = max(spans[chrom][1], end)
    return {chrom: (lo + hi) // 2 for chrom, (lo, hi) in spans.items()}


def build_centromere_dict(bed_path=None):
    """Return centromere dict: CHR_CENTRO overlaid with values from bed_path if provided."""
    if bed_path is None:
        return CHR_CENTRO
    merged = dict(CHR_CENTRO)
    merged.update(load_centromere_bed(bed_path))
    return merged


class CigarAlignment:
    def __init__(self, chrom:str, start:int, end:int, strand:str, ref_length:int, read_name:str, read_start:int, read_end:int, \
                 mapping_quality:float, edit_dist:float):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.ref_length = ref_length  # Length on reference genome (!= read query length above)
        self.read_name = read_name
        self.read_start = read_start
        self.read_end = read_end
        self.mapping_quality = mapping_quality
        self.edit_dist = edit_dist
    
    def __hash__(self):
        id = ','.join([self.chrom, self.start, self.end, self.strand, self.read_name, self.read_start, self.read_end])
        return hash(id)  # Must match __eq__
    
    def __eq__(self, other):
        if not isinstance(other, CigarAlignment):
            return False
        id1 = ','.join([self.chrom, str(self.start), str(self.end), self.strand, self.read_name, str(self.read_start), str(self.read_end)])
        id2 = ','.join([other.chrom, str(other.start), str(other.end), other.strand, other.read_name, str(other.read_start), str(other.read_end)])
        return id1 == id2

class SV:
    def __init__(self, chrom1:str, bp1:int, strand1:str, chrom2:str, bp2:int, strand2:str):
        self.chrom1 = chrom1
        self.bp1 = bp1
        self.strand1 = strand1
        self.chrom2 = chrom2
        self.bp2 = bp2
        self.strand2 = strand2
        self.type = self.get_SV_type()
        self.TST = False
        self.sort_breakpoints()
    
    def sort_breakpoints(self) -> None:
        if chrom_sort_key(self.chrom1) > chrom_sort_key(self.chrom2) or \
            (chrom_sort_key(self.chrom1) == chrom_sort_key(self.chrom2) and self.bp1 > self.bp2):
            self.chrom1, self.chrom2 = self.chrom2, self.chrom1
            self.bp1, self.bp2 = self.bp2, self.bp1
            self.strand1, self.strand2 = self.strand2, self.strand1
    
    def is_equal(self, other, max_bp_distance = 100) -> bool:
        if isinstance(other, SV) == False:
            return False
        if self.chrom1 != other.chrom1 or self.chrom2 != other.chrom2:
            return False
        if self.strand1 != other.strand1 or self.strand2 != other.strand2:
            return False
        if abs(self.bp1-other.bp1) > max_bp_distance or abs(self.bp2-other.bp2) > max_bp_distance:
            return False
        return True
    
    def get_SV_type(self) -> str:
        if self.is_foldback():
            return 'FBI'
        if self.chrom1 != self.chrom2:
            return 'TRA'
        if self.strand1 == self.strand2:
            return 'INV'
        if self.bp1 < self.bp2:
            return 'DEL' if self.strand1 == '+' else 'DUP'
        else:
            return 'DUP' if self.strand1 == '+' else 'DEL'

    def is_foldback(self, max_distance=50000) -> bool:
        if self.chrom1 != self.chrom2:
            return False
        if self.strand1 != self.strand2:
            return False
        if abs(self.bp1 - self.bp2) > max_distance:
            return False
        return True
    
    def is_in_region(self, region: tuple[str, int, int], flanking_length=1000000) -> bool:
        chrom, start, end = region
        flag1, flag2 = False, False
        if self.chrom1 == chrom and start-flanking_length <= self.bp1 <= end+flanking_length:
            flag1 = True
        if self.chrom2 == chrom and start-flanking_length <= self.bp2 <= end+flanking_length:
            flag2 = True
        return flag1, flag2
    
    def __str__(self):
        return f'{self.chrom1}:{self.bp1}{self.strand1}->{self.chrom2}:{self.bp2}{self.strand2}'
    def __repr__(self):
        return f'({self.chrom1}:{self.bp1}:{self.strand1})->({self.chrom2}:{self.bp2}:{self.strand2})'
    def __hash__(self):
        return hash(str(self))  # Hash based on the property
    def __eq__(self, other):
        return isinstance(other, SV) and str(self) == str(other)