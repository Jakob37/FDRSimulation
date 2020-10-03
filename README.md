# Simulation of the impact of using different p-value and FDR cutoffs

Access the running Shiny application at:
https://jakobshinyapps.shinyapps.io/fdr_exploration

In -omics (the study of biological molecules on a large scale), the analysis is often carried out for hundreds or thousands of molecules at once. 
When performing statistics for all these molecules, care needs to be taken to manage the false positives produced. If using a simple p-value threshold, 
this will, by definition, lead to a large part of false positives, 
as for instance, for p < 0.05, 5% of random features will still pass this threshold. 

A common tool to tackle this is the false discovery rate correction (FDR). These attempts to correct for the many tests to ensure that only a 
fraction of the molecules found as different are false positives. If applied correctly, this reduces the issue with false positives but instead 
induces a large number of false negatives - missed features regulated.

Here, you can explore how p-value thresholds and FDR-thresholds impact what you see in the data analysis. I illustrate this using p-value 
histograms, an efficient way to overview many performed tests. I demonstrate two FDR methods - the in omics commonly performed Benjamini Hochberg, 
and the even stricter Bonferroni correction.

