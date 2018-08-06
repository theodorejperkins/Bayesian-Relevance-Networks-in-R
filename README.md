# Bayesian-Relevance-Networks-in-R

This R source code is in support of the publication "Uncovering Robust Patterns of MicroRNA Co-Expression across Cancers using Bayesian Relevance Networks" by Ramachandran, Sanchez-Taltavull, and Perkins, PLoS ONE, Vol. 12, No. 8, Art. e0183103 (also appeared at GLBIO 2017, where it won "Outstanding Presentation" prize).

Sourcing the R code will provide for you five functions: BayesianCorrelation_Grouped, BayesianPermutation_Grouped, PearsonCorrelation_Grouped and PearsonPermutation_Grouped, and FDRAnalysis.

The first and third functions are different ways of computing correlations between different entities in a count matrix from a set of sequencing experiments. Typically, these entities would be genes or microRNAs whose expression is assessed by RNA-seq or single-cell RNA-seq. But they could be other things as well. The second and fourth functions are ways of estimating a null distribution of correlations under the hypothesis of no true correlation. The final function calcultes empirical false discovery rates, based on the outputs of the first and second functions, or the second and third functions. 

The basic steps of a complete analysis are to 1) Compute Bayesian (or Pearson) correlations between entities (genes, microRNAs, etc.), 2) Compute a null distribution of those correlations by permuting the data, and 3) threshold the correlations based on the permuted distribution and a desired false discovery rate threshold.
