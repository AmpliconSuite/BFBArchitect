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
BFBArchitect can be installed and run on most modern Unix-like operating systems (e.g. Ubuntu 18.04+, CentOS 7+, macOS). It requires python>=3.8 and the above dependencies. Please follow the instructions to install (more installation options will be provided soon):
1. Pull the source code
    ```
    git clone git@github.com:AmpliconSuite/BFBArchitect.git
    cd /path/to/BFBArchitect
    ```
2. Create a virtual environment (optional)
    ```
    python3 -m venv BFBArchitect_venv
    source BFBArchitect_venv/bin/activate
    ```
3. Install dependencies locally
    ```
    pip install .
    ```
4. Recommended for efficient ILP solving: Download a Gurobi optimizer license ([free for academic use](https://support.gurobi.com/hc/en-us/articles/360040541251-How-do-I-obtain-a-free-academic-license))
   * Place the ```gurobi.lic``` file in ```$HOME/gurobi.lic```.
   

## Running
Before running BFBArchitect, genome-wide copy number (CN) calls must be generated from the aligned long-read data by running the follow script:
```
python /path/to/BFBArchitect/scripts/call_CNV.py <input.bam> /path/to/BFBArchitect/scripts/hg38full_ref_5k.cnn <output_dir> <threads>
```
This will create a file called ```[input].cns```, which is a required argument in BFBArchitect. 

Then run BFBArchitect to reconstruct potential BFB sequences for any genomic region ```chrom:start-end``` with copy number amplification. (The amplicon region can be detected by standard pipelines like [CoRAL](https://github.com/AmpliconSuite/CoRAL).)
### Usage
```
python /path/to/BFBArchitect/src/BFBArchitect.py --bam <input.bam> --cns <input.cns> --region <chrom:start-end> --output_prefix <dir/output_prefix> [--segmentation] [--deletion] [--coverage <sequencing coverage>]
```
BFBArchitect also supports reconstructing BFB sequences at the whole-genome level, given CoRAL results at ```CoRAL_output_directory```: 
```
python /path/to/BFBArchitect/scripts/batch_run.py --directory <CoRAL_output_directory> --bam <input.bam> --cns <input.cns> --output_prefix <dir/output_prefix> [--segmentation] [--deletion] [--coverage <sequencing coverage>]
```
### Required arguments
- --bam <.bam file>: Aligned long reads
- --cns <.cns file>: The .cns file from genome-wide copy number calling
- --region <string>: A string that represents the amplified genomic region (e.g. chr1:1-1000000)
- --output_prefix <string>: The directory and prefix for all output files
### Optional arguments
- --segmentation: Consider copy number variation when segmenting the amplicon region. 
- --deletion: Handle deletion when reconstructing BFB sequences. 
- --coverage <integer>: Sequencing coverage (if provided, estimating coverage from cns will be skipped)

### Output
- graph.txt: A text file describing the segment and structural variant information of a breakpoint graph constructed from the amplicon region. 
- cycles.txt: A text file including the reconstructed BFB sequences. 
- reads.txt: A text file storing information of supporting reads for structural variants. 

## Sample run
Please download the sample input from this [link](https://drive.google.com/file/d/1OVAKD8kiH3vK9e2hE6YecMIoAulS_oId/view?usp=sharing), 
which includes sample.sorted.bam, sample.sorted.cns, and sample.sorted.cnr (for visualization). Run the following command:
```
python /path/to/BFBArchitect/src/BFBArchitect.py --bam BFBArchitect_input/sample.sorted.bam --cns BFBArchitect_input/sample.sorted.cns --region chr7:120000000-125000000 --output_prefix sample --coverage 15.0
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
python ~/BFBArchitect/src/BFBVisualizer.py --graph sample_graph.txt --cycle sample_cycles.txt --cnr BFBArchitect_input/sample.sorted.cnr --output_prefix sample
```
![Visualization generated by BFBArchitect](https://github.com/AmpliconSuite/BFBArchitect/blob/main/sample/sample_1.png)

## Contact
BFBArchitect is developed and maintained by Bafna Lab at UC San Diego. Please raise an issue or reach out to chl221@ucsd.edu if you have any questions. 

## License
This project is is licensed under the BSD 3-Clause License - see the [LICENSE.txt](https://github.com/AmpliconSuite/BFBArchitect/blob/main/LICENSE) file for details