## "R_merge_proteomics_volcano"
To use for merging proteomics results data with subcellular annotations and making volcano plots.

The input file "PTS3886_AHA_Surf4_protein_filtered.csv" contains post-processed SILAC proteomics data comparing WT
and SURF4 KO cell media. P-values using a one-sample t-test and other statistics were calculated
in Perseus software. The other input file "PTS3886_Subcell_annot.csv" contains subcellular annotations Uniprot.

1) merges results with protein Uniprot annotations for sub-cellular localisation
2) adjusts p-value for multiple testing (FDR)
3) colours protein hits by sub-cellular localisation
4) produces a volcano plot with statistically significantly changed proteins in the KO labelled
5) saves the merged csv file which now also contains adjusted p-value and significance tresholds

## "PNAS_2022_HEK_SURF4KO"
To use for making a volcano plot for further annotation in Illustrator.

The input file "Perseus-one-sample-t-testR.csv" contains post-processed SILAC proteomics data comparing WT
and SURF4 KO cell media. P-values using a one-sample t-test and other statistics were calculated
in Perseus software. Subcellular localisation data for each hit were added from Uniprot.

1) adjusts p-value for multiple testing (FDR)
2) colours protein hits by sub-cellular localisation
3) produces volcano plots, saves the last one with statistically significantly changed proteins in the KO labelled
4) saves the updated csv file which now also contains adjusted p-value and significance tresholds
