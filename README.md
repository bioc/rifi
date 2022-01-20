# Rifi framework

## About
Rifi framework is an open source R package, attempted to estimate decay by probe or by bin while using microarrays or RNA-seq. The estimation of decay is a result of fit of intensities upon time serie points. Rifi automizes the processes of fitting and utilizes a dynamic programming attempt for the clustering of the genome by the coefficients extracted from the fit. 

# Installation 

### Dependencies

Rifi framework has the following dependencies, whereas the brackets indicate the version Rifi has been build and tested on. Make sure the requirements are satisfied by your system. 

  [devtools](https://www.rdocumentation.org/packages/devtools/versions/1.13.6) (>= 1.13.6)
  
  [roxygen2](https://cran.r-project.org/web/packages/roxygen2/index.html) (>= 7.1.12)

  [car](https://www.rdocumentation.org/packages/car/versions/3.0-11) (>= 3.0-11)
  
  [cowplot](https://www.rdocumentation.org/packages/cowplot/versions/1.1.1) (>= 1.1.1)
  
  [doMC](https://cran.r-project.org/web/packages/doMC/index.html) (>= 1.3.7)
  
  [parallel](https://rdocumentation.org/packages/parallel/versions/3.6.2) (>= 3.6.2)
  
  [dplyr](https://www.rdocumentation.org/packages/dplyr/versions/0.7.8) (>= 1.0.7)
  
  [ggplot2](https://ggplot2.tidyverse.org/) (>= 3.3.5)
  
  [foreach](https://rdocumentation.org/packages/foreach/versions/1.5.1) (>= 1.5.1)
  
  [nls2](https://www.rdocumentation.org/packages/nls2/versions/0.1-2/topics/nls2) (>= 0.1-2)
  
  [nnet](https://cran.r-project.org/web/packages/nnet/index.html) (>= 7.3-16)
  
  [rlang](https://cran.r-project.org/web/packages/rlang/index.html) (>= 0.4.12)
  
  [scales](https://www.rdocumentation.org/packages/scales/versions/0.4.1) (>= 0.4.1)
    
  [stringr](https://www.rdocumentation.org/packages/stringr/versions/1.4.0)(>= 1.4.0)
    
  [egg](https://www.rdocumentation.org/packages/egg/versions/0.4.5/topics/ggarrange) (>= 0.4.5) 
  
  [graphics](https://rdocumentation.org/packages/graphics/versions/3.6.2) (>= 3.6.2) 
  
  [stats](https://rdocumentation.org/packages/stats/versions/3.6.2) (>= 3.6.2) 
  
  [methods](https://rdocumentation.org/packages/methods/versions/3.6.2) (>= 3.6.2) 
  
  [grid](https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/grid) (>= 3.6.2) 

  [stringr](https://rdocumentation.org/packages/stringr/versions/1.4.0) (>= 1.4.0) 

  [tibble](https://rdocumentation.org/packages/tibble/versions/3.1.6) (>= 3.1.6) 

  [utils](https://rdocumentation.org/packages/utils/versions/3.6.2) (>= 3.6.2) 

  [grDevices](https://rdocumentation.org/packages/grDevices/versions/3.6.2) (>= 3.6.2) 

### To install from Bioconductor, use the following code:
```html
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
    
BiocManager::install("rifi")
```

### To install directly from github:
```html
install_github('rifi')
```

# Usage

## Overview

Rifi consists of five top level steps for data processing, at the core of which, the fitting of the input data and the fragmentation of the resulting data are conducted. The output consists of the per bin data frame, a data frame containing all segment information and an event-based data frame. Additionally, the final data can be visualized in respect to a given annotation. 

In a first step, general data preprocessing e.g., filtering, segmentation based on filtered regions, and the assignment to one of two fitting models, with which the time series data is fitted in the subsequent step, is performed. A second step evaluates suiting penalties for the fragmentation by dynamic programming in the following step, at which the bins are fragmented by the values calculated at the fitting step e.g., delay, half-life, intensity, and termination factor, and grouped into transcriptional units. Finally, potential events such as pausing sites and iTSS are found and evaluated by statistical testing. An overview of the steps conducted in the full rifi data analysis can be seen in figure 1.

<p align="center">
  <img src="https://github.com/CyanolabFreiburg/rifi/blob/main/vignettes/principle.png"/>
</p>

<p align="center">
Figure 1: The general workflow of rifi. Blue boxes are data structures, green boxes are visual outputs, orange boxes are function. Incoming arrows represent the required input, outgoing arrows represent the output produced by the function. The orange frames state the main tasks of the function they refer to with dotted arrows
</p>


## Positional Arguments

Rifi provides different functional arguments (subcalls) for individual procedures. These include `rifi_preprocess`, `rifi_fit`, `rifi_penalties`, `rifi_fragmentation`, `rifi_stats`, `rifi_summary` and `rifi_visualization`


## Input 

The input data consists of a data frame containing the RNA-seq or microarray time series data from
rifampicin treated microorganisms, including all potential replicates, column wise as relative
intensities, e.g, Table 1 shows an example input.

t<sub>0</sub> | ... | t<sub>n</sub> | ID   | position   | strand |
  :---:   | :--:| :---:    | :--: | :---:      | :---:  |
  8.673   | ... | 3.460    | 1    | 50         | +      |
  8.974   | ... | 3.060    | 1    | 50         | +      |
  18.612  | ... | 7.650    | 2    | 50         | +      |
  15.476  | ... | 7.022    | 2    | 50         | +      |
  100.674 | ... | 100.460  | 300  | 4500       | -      |
  250.645 | ... | 50.460   | 300  | 4500       | -      |
  200.541 | ... | 56.460   | 300  | 4500       | -      |

Table 1. Rifi input data table. The first n columns contain the relative intensity measurements, with the first column referring to timepoint zero, at or before the addition of rifampicin. The third to last column contains the unique ID that is identical for each replicate. The second to last column contains the position information and the last column hold the information for the strand (“+”,”-“).

## Parameters

Rifi runs in R and the command line for the workflow are: 

`prepro <- rifi_preprocess(...)`

`probe_df <- rifi_fit(...)`

`probe_pen <- rifi_penalties(...)`

`probe_frag <- rifi_fragmentation(...)`

`probe_sta <- rifi_stats(...)`

`probe_summary <- rifi_summary(...)`

`rifi_visualization(...)`

## Results

The results of the analysis are stored in two convenient data frames, one for segments (Table .2) and one for the events (Table. 3).

ID    | gene  | ...  | strand| TU     | half-life |...    |p_value_TI|
:---: | :--:  | :---:| :--:  | :---:  | :---:     | :---: | :---:    |   
1     | gene1 | ...  | +     | TU_1   | 0.63      |...    |0.002     |
2     | gene1 | ...  | +     | TU_1   | 0.63      | ...   |0.002     |  
3     | gene1 | ...  | +     | TU_2   | 1.22      |...    |0.002     |
5     | NA    | ...  | -     | TU_501 | 2.03      | ...   | NA       |  
8     | NA    | ...  | -     | TU_501 | 0.93      |...    | NA       |
12    | gene2 | ...  | -     | TU_502 | 1.55      | ...   | NA       |  


Table 2: The segment data frame output. All data given is the value averaged over the given fragment. Additional columns represented by the dots are additional information.

ID    | gene  | ...  | strand| event  | p_value  |
:---: | :--:  | :---:| :--:  | :---:  | :---:    |   
1     | gene1 | ...  | +     | ps     | 0.0003   |
2     | gene1 | ...  | +     | ps     | 0.5      | 
3     | gene1 | ...  | +     | iTSS   | 0.01     |
5     | NA    | ...  | -     | ps     | 0.006    | 
8     | NA    | ...  | -     | iTSS   | 0.007    |
12    | gene2 | ...  | -     | iTSS   | 0.007    |

Table 3: The event data frame output. All data given is the value averaged over the given fragment. Additional columns represented by the dots are additional information.

<p align="center">
  <img src="https://github.com/CyanolabFreiburg/rifi/blob/main/vignettes/genome_fragments_plot.png"/>
</p>

<p align="center">
Figure 2: Figure 2: fragments and events plot. The plot shows 4 sections, annotation, delay, half-life (HL) and intensity. The annotation contains transcription units from rifi (TUs) and features of genome annotation. The fragments with different color are correspondingly plotted upon delay, HL and intensity. The events are indicated by different color (see help more detail). Significant events are assigned with an '*'. 
</p>

## Testing

Two data frames are provided to test Rifi: `example_input_e_coli` and `example_input_synechocystis_6803`

# Troubleshooting

contact jens.georg@biologie.uni-freiburg.de or create an issue