# BFBArchitect

BFBArchitect is a tool for detecting and reconstructing breakage-fusion-bridge (BFB) cycles from
long-read sequencing (currently support Oxford Nanopore). 

## Prerequisites
- CNVkit>=0.9.10 (https://cnvkit.readthedocs.io/en/master/quickstart.html)
- pandas>=2.3.3 (https://pandas.pydata.org/docs/whatsnew/index.html)
- PuLP>=3.3.0 (https://coin-or.github.io/pulp/main/includeme.html)
- pysam>=0.23.3 (https://pysam.readthedocs.io/en/latest/release.html)
- Python>=3.12.8 (https://www.python.org/downloads/release/python-3128/)

## Installation
BFBArchitect can be installed and run on most modern Unix-like operating systems (e.g. Ubuntu 18.04+, CentOS 7+, macOS). It requires python>=3.8 and the above dependencies.

First, pull the source code:
```
git clone git@github.com:AmpliconSuite/BFBArchitect.git
cd /path/to/BFBArchitect
```

### Option A: conda (recommended)
```bash
conda env create -f environment.yml
conda activate bfbarchitect
pip install -e . --no-deps    # -e (editable) makes changes take effect without reinstalling
```

### Option B: pip + virtual environment
```bash
python3 -m venv BFBArchitect_venv
source BFBArchitect_venv/bin/activate
pip install -e .
```

After installation, `BFBArchitect.py` is available as a system-wide command. You can call it from any directory:
```bash
BFBArchitect.py --help
```

### Gurobi license (recommended for efficient ILP solving)
Download a Gurobi optimizer license ([free for academic use](https://support.gurobi.com/hc/en-us/articles/360040541251-How-do-I-obtain-a-free-academic-license)) and place the ```gurobi.lic``` file at ```$HOME/gurobi.lic```. BFBArchitect autodetects the available solver: Gurobi is used when the license file exists and `gurobipy` is installed; otherwise it falls back to the open-source CBC solver (slower, no solution pool). The solver can also be set explicitly via `--solver gurobi|cbc`.
   

## Running
Before running BFBArchitect, genome-wide copy number (CN) calls must be generated from the aligned long-read data by running the follow script:
```
python /path/to/BFBArchitect/scripts/call_CNV.py <input.bam> /path/to/BFBArchitect/scripts/hg38full_ref_5k.cnn <output_dir> <threads>
```
This will create a file called ```[input].cns```, which is a required argument in BFBArchitect. 

Then run BFBArchitect to reconstruct potential BFB sequences for any genomic region ```chrom:start-end``` with copy number amplification. (The amplicon region can be detected by standard pipelines like [CoRAL](https://github.com/AmpliconSuite/CoRAL).)
### Usage
```
python /path/to/BFBArchitect/bfbarchitect/BFBArchitect.py --bam <input.bam> --cns <input.cns> --region <chrom:start-end> --output_prefix <dir/output_prefix> [--segmentation] [--no-deletion] [--coverage <sequencing coverage>]
```
BFBArchitect supports reconstructing BFB sequences at the whole-genome level, given CoRAL results at ```CoRAL_output_directory```: 
```
python /path/to/BFBArchitect/scripts/batch_run.py --directory <CoRAL_output_directory> --bam <input.bam> --cns <input.cns> --output_prefix <dir/output_prefix> [--segmentation] [--no-deletion] [--coverage <sequencing coverage>]
```
BFBArchitect also supports reconstructing BFB sequences directly from an [AmpliconArchitect](https://github.com/AmpliconSuite/AmpliconArchitect) (or BFBArchitect) `_graph.txt` file, with no BAM or CNS file required. Three modes are supported:
```
# Auto-detect BFB candidate regions (default) — one output set per region
BFBArchitect.py --graph <AA_graph.txt> --output_prefix <dir/output_prefix>

# Process a specific region only — single output at <output_prefix>_BFB_*
BFBArchitect.py --graph <AA_graph.txt> --region chr7:120000000-125000000 --output_prefix <dir/output_prefix>

# Treat all segments as one region — single output at <output_prefix>_BFB_*
BFBArchitect.py --graph <AA_graph.txt> --whole_graph --output_prefix <dir/output_prefix>

# Disable deletion handling, if needed for a control run
BFBArchitect.py --graph <AA_graph.txt> --no-deletion --output_prefix <dir/output_prefix>
```
`--region` and `--whole_graph` are mutually exclusive.

### Required arguments (BAM mode)
- --bam <.bam file>: Aligned long reads
- --cns <.cns file>: The .cns file from genome-wide copy number calling
- --region <string>: A string that represents the amplified genomic region (e.g. chr1:1-1000000)
- --output_prefix <string>: The directory and prefix for all output files

### Required arguments (graph mode)
- --graph <_graph.txt file>: An AA-format breakpoint graph file
- --output_prefix <string>: The directory and prefix for all output files

### Optional arguments (BAM mode)
- --segmentation: Consider copy number variation when segmenting the amplicon region.
- --coverage <integer>: Sequencing coverage (if provided, estimating coverage from cns will be skipped)

### Optional arguments (deletion handling)
- --no-deletion: Disable deletion handling. Deletion handling is enabled by default in both BAM and graph modes. In BAM mode, deletion-support evidence is added back into affected segment copy-number estimates. In graph mode, same-chromosome deletion-edge CN is added back to sequence segments skipped by those deletion edges before constructing the BFB CN vector.
- --deletion: Enable deletion handling explicitly. This is the default and is retained for compatibility with older commands.

### Optional arguments (both modes)
- --multiple: Reconstruct multiple optimal BFB candidate sequences (requires Gurobi)
- --solver gurobi|cbc: ILP solver to use (default: autodetect)
- -t / --threads <int>: Number of threads for the ILP solver (default: 8)
- --region <string>: (graph mode) process a specific region only, bypassing auto-detection (e.g. chr7:120000000-125000000). Mutually exclusive with --whole_graph.
- --whole_graph: (graph mode only) treat all segments as a single region instead of auto-detecting BFB regions
- --max-graph-segments <int>: (graph mode only) maximum number of graph segments allowed per graph-mode region (default: 100; use 0 to disable)
- -g / --gene <gtf_file>: Gene annotation for visualization (graph mode only)
- --centromere <.bed file>: Path to a BED file of centromere regions (≥3 tab-separated columns: chrom, start, end). Multiple rows per chromosome are merged to a single midpoint. Falls back to built-in hg38 defaults if omitted. An example GRCh38 file is provided at `resources/GRCh38_centromere.bed`.

### Output
- graph.txt: A text file describing the segment and structural variant information of a breakpoint graph constructed from the amplicon region. 
- cycles.txt: A text file including the reconstructed BFB sequences. 
- reads.txt: A text file storing information of supporting reads for structural variants. 

## Sample run
Please download the sample input from this [link](https://drive.google.com/file/d/1OVAKD8kiH3vK9e2hE6YecMIoAulS_oId/view?usp=sharing), 
which includes sample.sorted.bam, sample.sorted.cns, and sample.sorted.cnr (for visualization). Run the following command:
```
python /path/to/BFBArchitect/bfbarchitect/BFBArchitect.py --bam BFBArchitect_input/sample.sorted.bam --cns BFBArchitect_input/sample.sorted.cns --region chr7:120000000-125000000 --output_prefix sample --coverage 15.0
```
The following files will be output:
1. [sample_graph.txt](https://github.com/AmpliconSuite/BFBArchitect/blob/main/sample/sample_graph.txt)
   ```
    SequenceEdge: StartPosition, EndPosition, PredictedCN, AverageCoverage, Size, NumberOfLongReads
    sequence	chr7:120000001-	chr7:123649137+	1	15.953404599498457	3649137	11764
    sequence	chr7:123649138-	chr7:124186212+	10	84.07632453567938	537075	8843
    sequence	chr7:124186213-	chr7:124347510+	5	48.4874021996553	161298	1607
    sequence	chr7:124347511-	chr7:124740114+	3	31.59089311367179	392604	2494
    sequence	chr7:124740115-	chr7:124999999+	7	61.75353714142794	259885	3388
    BreakpointEdge: StartPosition->EndPosition, PredictedCN, NumberOfLongReads
    concordant	chr7:123649137+->chr7:123649138-	1	50
    concordant	chr7:124186212+->chr7:124186213-	5	66
    concordant	chr7:124347510+->chr7:124347511-	3	40
    concordant	chr7:124740114+->chr7:124740115-	3	46
    discordant	chr7:123649138-->chr7:123650665-	4	33
    discordant	chr7:124184685+->chr7:124186212+	2	14
    discordant	chr7:124346947+->chr7:124347510+	1	7
    discordant	chr7:124740115-->chr7:124740275-	2	17
    discordant	chr7:124807844+->chr7:124840986-	7	51
    discordant	chr7:124998462+->chr7:124999999+	4	29
   ```
2. [sample_cycles.txt](https://github.com/AmpliconSuite/BFBArchitect/blob/main/sample/sample_cycles.txt)
   ```
    Interval	1	chr7	120000001	124999999
    List of cycle segments
    Segment	1	chr7	120000001	123649137
    Segment	2	chr7	123649138	124186212
    Segment	3	chr7	124186213	124347510
    Segment	4	chr7	124347511	124740114
    Segment	5	chr7	124740115	124999999
    List of longest subpath constraints
    Path=1;Copy_count=1;Segments=1+,2+,3+,4+,5+,5-,5+,5-,5+,5-,4-,3-,2-,2+,2-,2+,3+,3-,2-,2+,2-,2+,3+,4+,5+,5-;Path_constraints_satisfied=;Score=0.30515965787039095;Multiplicity=1
   ```
3. [sample_reads.txt](https://github.com/AmpliconSuite/BFBArchitect/blob/main/sample/sample_reads.txt)
   ```
    SV	SV_type	TST	#Support_reads	Query_gaps	Support_reads
    ...
   ```

After reconstructing BFB sequences, visualization can be generated by running the following command:
```
python ~/BFBArchitect/bfbarchitect/BFBVisualizer.py --graph sample_graph.txt --cycle sample_cycles.txt --cnr BFBArchitect_input/sample.sorted.cnr --output_prefix sample
```
![Visualization generated by BFBArchitect](https://github.com/AmpliconSuite/BFBArchitect/blob/main/sample/sample_1.png)

## Library usage

BFBArchitect can be imported directly to reconstruct BFB sequences from Python without invoking the CLI.

### From a graph file

```python
from bfbarchitect import reconstruct_bfb_from_graph, write_bfb_graph, write_bfb_cycles, visualize_BFB
from bfbarchitect import build_centromere_dict

# hg38 defaults; supply a BED file for another assembly:
#   centromere_dict = build_centromere_dict('/path/to/centromere.bed')
centromere_dict = build_centromere_dict()

results = reconstruct_bfb_from_graph(
    'path/to/sample_graph.txt',
    centromere_dict=centromere_dict,
    deletion=True,  # Optional; default is True. Set False for a no-deletion control.
    threads=8,    # Optional: number of ILP solver threads (default: 8)
    silent=True   # Optional: suppress terminal output/logs
)
# results is a list of dicts, one per detected BFB region

for r in results:
    chrom, start, end = r['region']
    prefix = f'{chrom}_{start}'
    print(f"Region {chrom}:{start}-{end}")
    
    # Optionally write output files and visualize:
    write_bfb_graph(f'{prefix}_graph.txt', r['new_segments'], r['svs'], r['sv_info'])
    write_bfb_cycles(f'{prefix}_cycles.txt', r['new_segments'],
                     r['bfb_strings'], r['scores'], r['multiplicity'])
    
    # Generate the visualization PDF/PNG
    visualize_BFB(
        cycle_file=f"{prefix}_cycles.txt",
        graph_file=f"{prefix}_graph.txt",
        cnr_file=None,  # or path to .cnr for scatter plot
        output_prefix=f"{prefix}_BFB"
    )
```

Each result dict contains:
| Key | Type | Description |
|---|---|---|
| `region` | `(chrom, start, end)` | Detected BFB candidate region coordinates |
| `new_segments` | `list[(chrom, start, end, cn_float, coverage, read_count)]` | Segmented amplicon |
| `bfb_strings` | `list[list[int]]` | Reconstructed BFB paths (segment indices) |
| `scores` | `list[float]` | Per-candidate score (lower is better) |
| `multiplicity` | `int` | ILP multiplicity factor |
| `svs` | `list[SV]` | Discordant SVs in the region |
| `sv_info` | `dict[SV, (float, int)]` | Per-SV predicted CN and read count |

### From pre-segmented data

```python
from bfbarchitect import reconstruct_bfb, write_bfb_graph, write_bfb_cycles, visualize_BFB
from bfbarchitect import build_centromere_dict

centromere_dict = build_centromere_dict()

# new_segments: list of (chrom, start, end, cn_float, coverage, read_count)
# cn, lf, rf:   per-segment integer vectors aligned to new_segments
chrom = new_segments[0][0]
BFB_strings, scores, multiplicity = reconstruct_bfb(
    new_segments, cn, lf, rf,
    centromere_dict.get(chrom, 0),
    threads=8,  # Optional: number of ILP solver threads (default: 8)
)

# Optionally write output files and visualize:
prefix = "my_sample"
write_bfb_graph(prefix + '_graph.txt', new_segments, SVs, sv_info)
write_bfb_cycles(prefix + '_cycles.txt', new_segments, BFB_strings, scores, multiplicity)
visualize_BFB(prefix + '_cycles.txt', prefix + '_graph.txt', None, prefix)
```

## Contact
BFBArchitect is developed and maintained by Bafna Lab at UC San Diego. Please raise an issue or reach out to chl221@ucsd.edu if you have any questions. 

## License
This project is is licensed under the BSD 3-Clause License - see the [LICENSE.txt](https://github.com/AmpliconSuite/BFBArchitect/blob/main/LICENSE) file for details
