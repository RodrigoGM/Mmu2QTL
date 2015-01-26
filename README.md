[![DOI](https://zenodo.org/badge/7281/RodrigoGM/Mmu2QTL.png)](http://dx.doi.org/10.5281/zenodo.12793)


Mmu2 QTL
==========
This repository contains the raw data and the analysis scripts for mapping QTL on mouse chromosome 2 using CAST derived subcongenic strains.

**Please not that all the files contained in the data/ directory are excluded from the BSD license**. To use them for any purpose e.g. research or teaching, we would appreciate that you notify us by email.

Thank you in advance

### data/

This directory contains mouse body fat mass data and its covariates used to perform a QTL mapping.\
 You may use the data for self learning and teaching purposes. However, we would appreciate that you\
 notify us at <jfmedrano at ucdavis dot edu> or <rodrigo dot gularte at ulg dot ac dot be>\
 if you are distributing the data to a class or colleagues.

Plese don't hesitate to send us an email if you have any questions.

### figures/

Contains an SVG file with the published figures.

### analysis/

Contains shell scripts for deplyoing the analysis scripts in R.

*System requirements* : for the analysis you will require ```R:Language and Environment```, and the packages ```r/qtl```, ```boa```, ```coda```, ```hdrcde```, and ```RColorBrewer```.  You can install these as :
   ```install.packages(c("qtl", "boa", "coda", "hdrcde", "RColorBrewer"))```

All the scripts run on a single core, and most should work with ± 4GB of RAM.  However, scripts for the repicated CIM analyses will require more RAM ± 16 GB or more.

R/qtl allows for mutlithreading with the ```n.cluster``` argument (requires packages ```snow``` and ```rlecuyer```).  In addition ```lapply``` and ```sapply``` loops can be run in paralell with ```snow``` using ```parLapply``` or ```parSapply```. , For compatibility with most systems these options were not implemented.

### src/

Scripts and other source code to perform the data analysis.


