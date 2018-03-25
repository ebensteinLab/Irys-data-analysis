# Irys-data-analysis
Scripts and information regarding Irys data analysis

**If you use this software, please cite: [Genome-wide epigenetic profiling of 5-hydroxymethylcytosine by long-read optical mapping (https://doi.org/10.1101/260166)](https://doi.org/10.1101/260166)**

## Prerequisites
1. Python 3.4 or higher
2. Output files from BioNano MQR

## Scripts
### Filter xmap & q_cmap files by confidence and percent alignment
Use this script to filter an xmap file and its corresponding q_cmap file according to minimum alignment confidence and minimum percent of the molecule's length aligned to the reference (avoid ambiguous alignments).

_Input files:_ xmap file, q_cmap file

_Output files:_ **filtered** xmap file, **filtered** q_cmap file

**Usage:**
```
filter_xmap_by_confidence_and_percent_alignment.py [-h] [-o OUTPUT DIR]
                                                        [-c MIN CONFIDENCE]
                                                        [-a MIN PERCENT ALIGNMENT]
                                                        [-p PREFIX]
                                                        xmap q_cmap

Filter BioNano xmap and q_cmap by alignment confidence and percent alignment out of molecule length

positional arguments:
  xmap                  input xmap file
  q_cmap                input q_cmap file

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT DIR, --output OUTPUT DIR
                        output files directory. default current working directory
  -c MIN CONFIDENCE, --conf MIN CONFIDENCE
                        min alignment confidence. default 12
  -a MIN PERCENT ALIGNMENT, --alignment MIN PERCENT ALIGNMENT
                        min percent of molecule length aligned. default 60
  -p PREFIX, --prefix PREFIX
                        output files prefix
```

**Example:**

Filter on minimum confidence 20 and minimum aligned length percentage 75:
```
filter_xmap_by_confidence_and_percent_alignment.py -o output_dir -c 20 -a 75 -p min.conf.20.min.alignment.75. xmap_file.xmap cmap_file_q.cmap
```

### Convert xmap & cmap files to stretched BED file
This script converts single molecule optical mapping alignment data stored in xmap and cmap files to bioinformatics compatible BED files.

_Input files:_ xmap file, q_cmap file, r_cmap file, chromosome key file (produced when converting the reference fasta file to the refernce cmap)

_Output options:_
1. BED file of labels in QUERY channel

   The reported labels take up 3 bp in the BED file. You can then extend them separately to account for optical resolution.

2. BED file of alignment regions of molecules that had at least 1 label in QUERY channel

   **Note** that in order to get a BED file of **all** aligned molecules, you need to run the script twice, once for each channel, then merge the files using the _Merge two molecule regions BED files_ script.

**Usage:**
```
xmap_cmap_to_BED_stretched.py [-h] [-o OUTPUT DIR]
                                   [-c MIN_CONFIDENCE]
                                   [-a ALIGNMENT LABEL CHANNEL]
                                   [-q QUERY LABEL CHANNEL] [-lb] [-mb]
                                   [-p PREFIX]
                                   xmap q_cmap r_cmap key

Convert BioNano xmap and cmap to labels and molecules BED files

positional arguments:
  xmap                  input xmap file
  q_cmap                input q_cmap file
  r_cmap                input r_cmap file
  key                   ref cmap key file

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT DIR, --output OUTPUT DIR
                        output files directory. default current working directory
  -c MIN_CONFIDENCE, --conf MIN_CONFIDENCE
                        min alignment confidence. default 12
  -a ALIGNMENT LABEL CHANNEL, --alignment ALIGNMENT LABEL CHANNEL
                        alignment label channel. default 1
  -q QUERY LABEL CHANNEL, --query QUERY LABEL CHANNEL
                        query label channel. default 2
  -lb, --label_bed      save query labels BED file. default False
  -mb, --molecule_bed   save molecules containing query labels. default False
  -p PREFIX, --prefix PREFIX
                        output files prefix
```

**Exmaple:**

Running:
```
xmap_cmap_to_BED_stretched.py -o output_dir -lb -mb -p output_file_prefix xmap_file.xmap mol_q.cmap ref_r.cmap ref_key.txt
```
will output two two BED files, one of all labels in channel 2, and another of all molecules that had at least one label in channel 2. Only molecules with alignment confidence >= 12 will be considered.

### Merge two BED files containing molecule regions
Use this script to merge two molecule regions BED files from two different channels produced by the _xmap_cmap_to_BED_stretched_ script. Running this script ensures that you get the most complete range possible for each molecule.

_Input files:_ Two molecule regions BED files, one for each channel.

_Output files:_ One merged molecule regions BED file.

**Usage:**
```
merge_2_mol_region_files.py [-h] [-o OUTPUT DIR] [-p PREFIX] regions_ch1 regions_ch2

Merge channel 1 and channel 2 molecule region BED files

positional arguments:
  regions_ch1           molecule regions channel 1 input path
  regions_ch2           molecule regions channel 2 input path

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT DIR, --output OUTPUT DIR
                        output file directory. default current working directory
  -p PREFIX, --prefix PREFIX
                        output file prefix
```

### Convert xmap & cmap files to GFF3 format
A collection of scripts to convert single molecule optical mapping data to GFF3 format, in order to view single molecule alignments in JBrowse.
