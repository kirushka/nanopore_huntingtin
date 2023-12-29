# Verification of CAG repeat length in Htt plasmid using ONT reads

## ONT reads analysis

Bash-script `01_ONT-reads-analysis.sh` containes all steps of ONT reads preprocessing (QC, alignment, and repeat length estimation).

## Plots

The results of ONT reads analysis were visualized using R. R-scripts `02_read_length.R`, `03_read_number.R`, and `04_repeat_length.R` contain code for plotting reads length distribution (based on NanoPlot output), read number per group distribution (based on NanoPlot output), and CAG repeat length distribution in Q138 group (or Q138 + unclassified reads), respectively.

## Supplementary files

- `htt_repeat_linear.txt` - BED files containing CAG repeat coordinates in the reference plasmid.