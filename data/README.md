[![DOI](https://zenodo.org/badge/7281/RodrigoGM/Mmu2QTL.png)](http://dx.doi.org/10.5281/zenodo.12793)

data/
---
This directory contains mouse body fat mass data and its covariates used to perform a QTL mapping. You may use the data for self learning and teaching purposes. However, we would appreciate that you notify us at <jfmedrano at ucdavis dot edu> or <rodrigo dot gularte at ulg dot ac dot be> if you are distributing the data to a class or colleagues, or using it for your own research.

Plese don't hesitate to send us an email if you have any questions.

* *SBC_hg2d_pheno.csv* : At the moment contains only body fat data and the covariates used to perform an in-depth QTL analysis on distal mouse choromosome 2.  This is the raw, curated data. The curation process consisted of double checking all the phenotypes for every single individual with our paper records at the moment of the dissection and after finishing with the data input process.

* *SBC_hg2d_geno.csv* : This file contains the genotypes of the mice used for the QTL analysis. There are two types of markers microsatellites and SNP.  For microsatellites the midpoint between the start and end positions were used.  All locations are based on Genome build mm9/NCBI_b37.
 
* *abv.dat* : This file contains the abreviations, the phenotype's full name, and any information pertaining the phenotypes

* *Fatq2b - Atp5e_Ctsz.sds*, *Fatq2b - Gnas_Rab22a.sds*, *Fatq2b - Gus_GusRanin.sds*, *Fatq2b - Stx16_SDHA.sds* : These files contain the raw qPCR expression data for each individual mouse.  Files are in *.sds format and can be read/used with ABI Data File Converter or ABI Prism(R). Please vist ABI <appliedbiosystems.com> for more information on how to access these files.
