input: text file of proteomics data on linear scale (calculated in MS Excel)

hierarchical clustering: performs and plots hierarchical clustering analysis (Figure 2c)

input: Excel file reanalysed_mouse_vs_invitro, annotation file protein_annotation.tsv

dataset_overview: calculates detected protein numbers (Figure 2b)

data_cleanup: reads and cleans data (filtering out proteins measured in only one replicate), calculates normalised data

followed by: multiple t-tests in Graphpad Prism

volcano plots all species: reads in t-tests from Prism file and plots volcano plots with all species (Figure 2d and Supplementary Figure 7)

enrichments all species: performs enrichment analysis across all species for all comparisons (Figure 2e and Supplementary Figure 8)

enrichment APF vs. mouse per species: performs enrichment analysis per species for the comparison APF vs. mouse (Extended Data Figure 5)

