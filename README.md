**tfbsincbeta:** Software for a Bayesian method for estimating the number of 
binding sites for a transcription factor (TF), based on a known position-frequency
matrix (PFM) for the TF, within a set of promoter sequences. The method
incorporates a prior probability distribution on the number of transcription 
factor binding sites (TFBS) that is similar to the beta-binomial distribution
but with modifications to account for the double-stranded nature of DNA. Samples
from the posterior probability distribution of the number of TFBS are generated
using a Metropolis-Hastings algorithm with a proposal generator that is weighted
based on the Shannon entropies of the probabilities for presence/absence of
a binding site each possible TFBS position. The software accompanies the 
manuscript "An empirical prior improves accuracy for Bayesian estimation of 
transcription factor binding site frequencies within gene promoters" by
Stephen Ramsey, which has been submitted to the journal *Bioinformatics and
Biology Insights.*

**Author:**  Stephen Ramsey, Oregon State University (lab.saramsey.org)

**Date:**  Sept. 10, 2016

**License:** This software is distributed under the Apache Software License 2.0.
Please see the file LICENSE for details on the software licensing
agreement.

**Usage notes:** The R script, "tfbsincbeta.R", reads a data file "Matrices.txt" of
TFBS PFMs in tab-delimited format (see header comment for "tfbsincbeta.R" and an
example file in "data"). In order to generate the empirical performance results
in the above-referenced article, the "Matrices.txt" file contains all TF PFMs
from the TRANSFAC Professional database version 2015.1. That database PFMs can
be obtained from QIAGEN but the author is not permitted to redistribute the
database. The R script makes use of the R package "parallel".


