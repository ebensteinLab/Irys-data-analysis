# Irys-data-analysis
Scripts and information regarding Irys data analysis

**If you use this software, please cite: [Genome-wide epigenetic profiling of 5-hydroxymethylcytosine by long-read optical mapping](https://doi.org/10.1101/260166)**

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

### Convert xmap & cmap files to stretched BED file

### Merge two BED files containing molecule regions

### Convert xmap & cmap files to GFF3 format
A collection of scripts to convert single molecule optical mapping data to GFF3 format, in order to view single molecule alignments in JBrowse.
