Codon usage indices 

This document describes the indices calculated by CodonW, by default only 
the G+C content of the sequence is reported. The others being dependent on 
the genetic code selected. More than one index may be calculated at the same 
time.  

Codon Adaptation Index (CAI) (Sharp and Li 1987). 
CAI is a measurement of the relative adaptiveness of the codon usage of a 
gene towards the codon usage of highly expressed genes. The relative 
adaptiveness (w) of each codon is the ratio of the usage of each codon, to 
that of the most abundant codon for the same amino acid. The relative 
adaptiveness of codons for albeit a limited choice of species, can be 
selected from Menu 3. The user can also input a personal choice of values. 
The CAI index is defined as the geometric mean of these relative 
adaptiveness values. Non-synonymous codons and termination codons (dependent 
on genetic code) are excluded. 
 
To prevent a codon absent from the reference set but present in other genes 
from having a relative adaptiveness value of zero, which would cause CAI to 
evaluate to zero for any genes which used that codon; it was suggested that 
absent codons should be assigned a frequency of 0.5 when estimating ? (Sharp 
and Li 1987). An alternative suggestion was that ? should be adjusted to 
0.01 where otherwise it would be less than this value (Bulmer 1988). CodonW 
does not adjust the ? value if a non-zero-input value is found; zero values 
are assigned a value of 0.01. 

Frequency of Optimal codons (Fop) (Ikemura 1981). 
This index, is the ratio of optimal codons to synonymous codons (genetic 
code dependent). Optimal codons for several species are in-built and can be 
selected using Menu 3. By default, the optimal codons of E. coli are 
assumed. The user may also enter a personal choice of optimal codons. If 
rare synonymous codons have been identified, there is a choice of 
calculating the original Fop index or a modified Fop index. Fop values for 
the original index are always between 0 (where no optimal codons are used) 
and 1 (where only optimal codons are used). When calculating the modified 
Fop index, negative values are adjusted to zero. 

Codon Bias Index (CBI) (Bennetzen and Hall 1982). 
Codon bias index is another measure of directional codon bias, it measures 
the extent to which a gene uses a subset of optimal codons. CBI is similar 
to Fop as used by Ikemura, with expected usage used as a scaling factor. In a 
gene with extreme codon bias, CBI will equal 1.0, in a gene with random 
codon usage CBI will equal 0.0. Note that it is possible for the number of 
optimal codons to be less than expected by random change. This results in a 
negative value for CBI.

The effective number of codons (NC) (Wright 1990).
This index is a simple measure of overall codon bias and is analogous to the 
effective number of alleles measure used in population genetics. Knowledge 
of the optimal codons or a reference set of highly expressed genes is 
unnecessary. Initially the homozygosity for each amino acid is estimated 
from the squared codon frequencies (see Equation 5).

	
If amino acids are rare or missing, adjustments must be made. When 
there are no amino acids in a synonymous family, Nc is not calculated 
as the gene is either too short or has extremely skewed amino acid 
usage (Wright 1990). An exception to this is made for genetic codes 
where isoleucine is the only 3-fold synonymous amino acid, and is not 
used in the protein gene. The reported value of Nc is always between 20 
(when only one codon is effectively used for each amino acid) and 61 
(when codons are used randomly). If the calculated Nc is greater than 
61 (because codon usage is more evenly distributed than expected), it 
is adjusted to 61.

G+C content of the gene. 
The frequency of nucleotides that are guanine or cytosine.

G+C content 3rd position of synonymous codons (GC3s).
This the fraction of codons, that are synonymous at the third codon 
position, which have either a guanine of cytosine at that third codon 
position. 

Silent base compositions. 
Selection of this option calculates four separate indices, i.e. G3s, C3s, 
A3s & T3s. Although correlated with GC3s, this index is not directly 
comparable. It quantifies the usage of each base at synonymous third codon 
positions. When calculating GC3s each synonymous amino acid has at least one 
synonym with G or C in the third position. Two or three fold synonymous 
amino acids do not have an equal choice between bases in the synonymous 
third position. The index A3s is the frequency that codons have an A at their 
synonymous third position, relative to the amino acids that could have a 
synonym with A in the synonymous third codon position. The codon usage 
analysis of Caenorhabditis elegans identified a trend correlated with the 
frequency of G3s. Though it was not clear whether it reflected variation in 
base composition (or mutational biases) among regions of the C. elegans 
genome, or another factor (Stenico et al. 1994).

Length silent sites (Lsil). 
Frequency of synonymous codons.

Length  amino acids (Laa). 
Equivalent to the number of translatable codons.

Hydropathicity of protein. 
The general average hydropathicity or (GRAVY) score, for the hypothetical 
translated gene product. It is calculated as the arithmetic mean of the sum 
of the hydropathic indices of each amino acid (Kyte and Doolittle 1982). 
This index has been used to quantify the major COA trends in the amino acid 
usage of E. coli genes (Lobry and Gautier 1994). 

Aromaticity score
The frequency of aromatic amino acids (Phe, Tyr, Trp) in the hypothetical 
translated gene product. The hydropathicity and aromaticity protein scores 
are indices of amino acid usage. The strongest trend in the variation in the 
amino acid composition of E. coli genes is correlated with protein 
hydropathicity, the second trend is correlated with gene expression, while 
the third is correlated with aromaticity (Lobry and Gautier 1994). The 
variation in amino acid composition can have applications for the analysis 
of codon usage. If total codon usage is analysed, a component of the 
variation will be due to differences in the amino acid composition of genes. 



Bennetzen, J. L., and B. D. Hall, (1982). Codon selection in yeast. Journal 
of Biological Chemistry 257: 3026-3031.
Bulmer, M., (1988). Are codon usage patterns in unicellular organisms 
determined by selection-mutation balance. Journal of Evolutionary 
Biology 1: 15-26.
Ikemura, T., (1981). Correlation between the abundance of Escherichia coli 
transfer RNAs and the occurrence of the respective codons in its 
protein genes: a proposal for a synonymous codon choice that is 
optimal for the E. coli system. Journal of Molecular Biology 151: 389-
409.
Kyte, J., and R. Doolittle, (1982). A simple method for displaying the 
hydropathic character of a protein. Journal of Molecular Biology 157: 
105-132.
Lobry, J. R., and C. Gautier, (1994). Hydrophobicity, expressivity and 
aromaticity are the major trends of amino acid usage in 999 
Escherichia coli chromosome encoded genes. Nucleic Acids Research 22: 
3174-3180.
Sharp, P. M., and W. H. Li, (1987). The codon adaptation index a measure of 
directional synonymous codon usage bias, and its potential 
applications. Nucleic Acids Research 15: 1281-1295.
Stenico, M., A. T. Lloyd and P. M. Sharp, (1994). Codon usage in 
Caenorhabditis elegans delineation of translational selection and 
mutational biases. Nucleic Acids Research 22: 2437-2446.
Wright, F., (1990). The effective number of codons used in a gene. Gene  87 
: 23-29.

