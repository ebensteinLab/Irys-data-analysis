# Scripts for wavelet decomposition analysis
Scripts for automated correlation analysis between optical mapping and sequencing data using ENCODE's wavelet decomposition pipeline.

## Scripts
### automated_wavelet_pipe.py
This is the main script used in the analysis, which converts coverage data in bedGraph format to a continuous vector of values to be used
for wavelet decomposition analysis.

**input:** optical mapping bedgraph, sequencing coverage bedgraph

**output:** bedGraph containing Pearson correlation values between wavelet decomposition coefficients.

As the analysis is computationally intensive, it may be necessary to break the genome into multiple non-overlapping regions, and analyse
each region seperately.

This main script calls two R scripts for each region: _MODWT_calc.R_ is called first to reduce the resolution of the sequencing data to
fit the inherent lower resolution of optical mapping, to facilitate the comparison between the two data sets. _CWT_DOG_calc.R_ is then
called to perform wavelet decomposition and to save the decomposition coefficients.

The ```pearsonr``` function from the ```scipy.stats.stats``` module is used to calculate the correlation between the decomposition
coefficients. For each genomic position, the correlation is calculated based on the values in a window centered at that position. The size
of the window can be determined by the user.

### min_cor_regions.py
Script to go over correlation bedGraph produced for each genomic region and produce one genome-wide correlation track with a specific
minimum correlation value.

As the minimal Pearson correlation value is -1, a value of -2 (default) will result in a genome-wide correlation track.

This script can be used to identify high correlation regions along the genome.

### merge_same_val_in_bedgraph.py
Script used to reduce the size of resulting bedGraph files on disk.
**input:** bedGraph file to be reduced
**output:** reduced size bedGraph

This script merges together consecutive regions in a bedGraph file which have the same score.
