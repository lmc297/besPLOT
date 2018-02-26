# besPLOT

Interactive PCA and NMDS plots

## Overview

besPLOT takes any numerical matrix as input, runs various statistical analyses, and outputs interactive and downloadable plots of the results.
Users can upload metadata to incorporate into their plots.

### Input

#### Required

* **Matrix file** *(n x m matrix)*: a whitespace-, tab-, space-, comma-, semicolon-, or colon-separated file
with n rows and m columns
    + *Rows* represent individuals or samples
    + *Columns* represent measured attributes
    + *Row 1* should contain column names (i.e. the file should have a header)
    + *Column 1* should contain sample names
    + *Data* in row 2 to row n and column 2 to column m should be numerical
    
#### Optional

* **Metadata file** *(n x 2 matrix)*: a whitespace-, tab-, space-, comma-, semicolon-, or colon-separated file
with n rows and 2 columns
    + *Rows* represent individuals or samples
    + *Column 1* contains sample names and should contain sample names identical to those in column 1 of the matrix file (order doesn't matter)
    + *Column 2* contains metadata associated with each sample (e.g. host, country)
    + *No header* should be included (i.e. sample names and metadata start in the first row)
    
The current version of besPLOT supports the following analyses:

* Principal component analysis (PCA) using the ```prcomp``` function in R's ```stats``` package

* Non-metric multidimensional scaling (NMDS) using the ```metaMDS``` function in R's ```vegan``` package

The current version of besPLOT supports the following types of metadata:

* Categorical metadata

besPLOT is a shiny application that can be run locally on your computer through R Studio. 
The application can be downloaded from https://github.com/lmc297/besPLOT.

Post issues at https://github.com/lmc297/besPLOT/issues

As always, some statistical analyses may not be appropriate for certain data sets. If you're unsure of what you're doing, consult a statistician.

### Citation

#### If you found the besPLOT app and/or its source code to be useful, please cite:

Carroll, Laura M. 2018. besPLOT: Interactive PCA and NMDS plots. Version 1.0.0. https://github.com/lmc297/besPLOT.

#### If you are interested in how the statistical methods implemented in besPLOT can be applied to genomic data (or if you just want to see some cool plots), check out the following papers:

Carroll, Laura M., Jasna Kovac, Rachel A. Miller, Martin Wiedmann. 2017. Rapid, high-throughput identification of anthrax-causing and emetic *Bacillus cereus* group genome assemblies using BTyper, a computational tool for virulence-based classification of *Bacillus cereus* group isolates using nucleotide sequencing data. *Applied and Environmental Microbiology* 2017 Jun 16. pii: AEM.01096-17. doi: 10.1128/AEM.01096-17.

Carroll, Laura M., Martin Wiedmann, Henk den Bakker, Julie Siler, Steven Warchocki, David Kent, Svetlana Lyalina, Margaret Davis, William Sischo, Thomas Besser, Lorin D. Warnick, Richard V. Pereira. Whole-Genome Sequencing of Drug-Resistant *Salmonella enterica* Isolated from Dairy Cattle and Humans in New York and Washington States Reveals Source and Geographic Associations. *Applied and Environmental Microbiology* 2017 May 31;83(12).


------------------------------------------------------------------------

## Launching besPLOT (Latest Version)

1. Download R, if necessary: https://www.r-project.org/

2. Dowload R Studio, if necessary: https://www.rstudio.com/products/rstudio/download/

3. Open R Studio, and install the following packages, if necessary, by typing the following commands into R Studio's console:

```
install.packages("shiny")
install.packages("ggplot2")
install.packages("vegan")
install.packages("plyr")
install.packages("dplyr")
install.packages("cluster")
install.packages("ggrepel")
```

4. Load shiny package

```
library(shiny)
```

5. Launch the app by typing the following command into R Studio's console:
```
runGitHub("besPLOT","lmc297")
```

You're ready to go!

------------------------------------------------------------------------

## Launching besPLOT (Older Versions)

**Currently, there is only 1 version of besPLOT; however, as more versions are added, this is how to run them:**

1. Download R, if necessary: https://www.r-project.org/

2. Dowload R Studio, if necessary: https://www.rstudio.com/products/rstudio/download/

3. Open R Studio, and install the following packages, if necessary, by typing the following commands into R Studio's console:

```
nstall.packages("shiny")
install.packages("ggplot2")
install.packages("vegan")
install.packages("plyr")
install.packages("dplyr")
install.packages("cluster")
install.packages("ggrepel")
```

4. Load shiny package

```
library(shiny)
```

5. Launch version X.Y.Z of the app by typing the following command into R Studio's console:
```
runGitHub(repo = "BMiner", username = "lmc297", subdir = "archive/besPLOT-X.Y.Z")
```

You're ready to go!

------------------------------------------------------------------------

# Uploading Data

## Required Steps

**For an example matrix file, see section "Matrix File Examples" below**

0a. Make sure that the text file containing your matrix has columns separated by one of the following: 
whitespace, tab, space, comma, colon, or semicolon

0b. Make sure your samples/individuals are in rows, and the attributes you are measuring for each sample are in columns

0c. Make sure the first column of your matrix contains the name of each sample

0d. Make sure the first row of your matrix contains the name of each attribute (i.e. your matrix has a header: ```read.table(header = TRUE)```

0e. Make sure all other rows/columns contain numerical data


1. Upload your file  by clicking the "Browse" button under "Matrix file (n x m matrix)" in the left panel.

2. Select the delimiter for your matrix in the drop-down menu under "Matrix delimiter", i.e.
tell besPLOT which character is used to separate columns of your matrix (whitespace, tab, space, comma, colon, or semicolon)

## Additional Steps (Categorical Metadata)

**For an example metadata file, see section "Metadata File Examples" below**

0a. Make sure that the text file containing your metadata has columns separated by one of the following: 
whitespace, tab, space, comma, colon, or semicolon

0b. Make sure your metadata file has 2 columns

0c. Make sure the first column of your metadata file contains the name of each sample, and make sure they match the sample
names in column 1 of your matrix file (the order doesn't matter)

0d. Make sure the second column of your metadata file contains categorical metadata (besPLOT currently does not support
continuous metadata)

0e. Make sure your matrix does not contain a header (i.e. the first row contains the sample name and metadata for your first sample)


1. If you have a categorical metadata file, you can upload it by selecting the "Browse" 
button under "Metadata file (n x 2 matrix)" in the left panel

2. Select the delimiter for your metadata file in the drop-down menu under "Metadata delimiter", i.e.
tell besPLOT which character is used to separate columns of your metadata (whitespace, tab, space, comma, colon, or semicolon)

------------------------------------------------------------------------

# Analyses Using besPLOT

## Principal Component Analysis (PCA)

This performs PCA using the ```prcomp``` function in R's ```stats``` package. An  interactive plot is produced using ```ggplot2```.

* By default, besPLOT plots the first, second, and third principal components (PC1, PC2, and PC3, respectively) on
the x-axis, y-axis, and as the point size, respectively; if you want to change which PCs are viewed, you can select them
from the lists in the left panel.

* You can determine which samples are associated with which points by clicking on them.

* You can zoom in by dragging your mouse and double-clicking on an area of the plot (double-click to zoom back out). 

* You can color your points using your categorical metadata by checking the "Overlay Metadata" box.

* To download the plot, click "Download this plot" in the left panel.

Note: I have hard-coded `prcomp` to scale and center the data.

## Non-Metric Multidimensional Scaling (NMDS)

This performs NMDS using the ```metaMDS``` function in R's ```vegan``` package and produces an interactive, 2-dimensional plot 
using ```ggplot2```.

NMDS is a method of ordination in which input data are fit to *k* dimensions (in besPLOT, k is set to 2 so that the results
can be plotted in 2 dimensions). 

NMDS is valuable because it makes few assumptions about the nature of your data, and it allows any distance/dissimilarity measure to be used. 

besPLOT has a variety of dissimilarity metrics available, so consult the documentation for ```vegdist```
to learn which is appropriate for your data (type ```?vegdist``` in your R Studio console)

Potential drawbacks of NMDS are that it is more compuationally intensive and, thus, slower than other ordination methods. 
For very large data sets, I would recommend using PCA (described above), as it is much faster. 
Also, depending on the nature of your data, `metaMDS` may never be able to find an optimal solution. 
Note: I have hard-coded `metaMDS` as run through besPLOT to try a maximum of 10000 random starts (`trymax=10000`). 
`metaMDS` will stop once this limit is reached.

* In besPLOT, NMDS plots display NMDS axis 1 on the x-axis (NMDS1) and NMDS axis 2 on the y (NMDS2). 

* Points represent your samples/individuals , and clicking a point will give you a list of the isolates associated with it 
under the plot. 

* Dragging your mouse and double-clicking allows you to zoom in on the plot (double-click to zoom back out). 

* If you have input categorical data, you can draw convex hulls around corresponding assemblies by checking the 
"Overlay Metadata" box. 

* To export your plot, click "Download this plot" in the left panel.


------------------------------------------------------------------------

# Matrix File Examples

Example of a comma-separated matrix file (species counts from ```data(dune)```, with "site" column added as column 1):
```
site,Achimill,Agrostol,Bromhord,Juncarti,Planlanc,Salirepe
1,1,0,0,0,0,0
2,3,0,4,0,0,0
3,0,4,0,0,0,0
4,0,8,3,0,0,0
5,2,0,2,0,5,0
```

Example of a tab-separated/whitespace-separated matrix file (presence/absence data):
```
subject  cuz_im_fly  i_dont_gotta_nap  i_could_sell_a_mill_saying_nothing_on_the_track i_represent_new_york
Mims  1 1 1 1 
You 0 0 0 0 
Snoop 1 0 0 0
```
------------------------------------------------------------------------

# Metadata File Examples

Example of a semicolon-separated matrix file ("Management" variable associated with species counts from ```data(dune)```):
```
1;SF
2;BF
3;SF
4;SF
5;HF
```

Example of a whitespace-separated matrix file:
```
Mims  Hot
You   Not
Snoop Unknown
```

------------------------------------------------------------------------


# References

### R Packages

H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.

Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.

Hadley Wickham and Romain Francois (2016). dplyr: A Grammar of Data Manipulation. R package version 0.5.0. https://CRAN.R-project.org/package=dplyr

Jari Oksanen, F. Guillaume Blanchet, Michael Friendly, Roeland Kindt, Pierre Legendre, Dan McGlinn, Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens, Eduard Szoecs and Helene Wagner (2017). vegan: Community Ecology Package. R package version 2.4-2. https://CRAN.R-project.org/package=vegan

Kamil Slowikowski (2016). ggrepel: Repulsive Text and Label Geoms for 'ggplot2'. R package version 0.6.5. https://CRAN.R-project.org/package=ggrepel

Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2017).  cluster: Cluster Analysis Basics and Extensions. R package version 2.0.6.

Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2017). shiny: Web Application Framework for R. R package version 1.0.1. https://CRAN.R-project.org/package=shiny

### Journals

J.B. Kruskal. Multidimensional scaling by optimizing goodness
of fit to a nonmetric hypothesis. 1964. *Psychometrika* vol. 29, no. 1.

J.B. Kruskal. Nonmetric multidimensional scaling: A numerical method. 1964. *Psychometrika* vol. 29, issue 2, 115-129.

Jonathan M. Chase, Nathan J. B. Kraft, Kevin G. Smith, Mark Vellend, Brian D Inouye. Using null models to disentangle variation in community dissimilarity from variation in α-diversity. 2011. *Ecosphere* vol. 2, issue 2, 1-11.

### Other Materials

Mims, Shawn. 2007. <a href="https://www.youtube.com/watch?v=TwyE3WJ4AWo">This is Why I'm Hot.</a> *Music Is My Savior* Prod. by Blackout Movement & Mike Moosh for American King Music, Capitol Records.

Steven M. Holland. Non-Metric Multidimensional Scaling (MDS). May 2008. http://strata.uga.edu/software/pdf/mdsTutorial.pdf

RHG Jongman, CJF ter Braak, & OFR van Tongeren. (1987). Data Analysis in Community and Landscape Ecology. Pudoc, Wageningen.

------------------------------------------------------------------------

Disclaimer: besPLOT is pretty neat! However, no tool is perfect, and you should always interpret your results with caution. I am not responsible for misinterpretations (statistical or otherwise) of besPLOT results.
