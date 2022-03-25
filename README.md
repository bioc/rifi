# Rifi framework

## About
Rifi is an open source R package, attempted to estimate decay by probe or by bin, based on high resolution microarray or RNA-seq data. The decay is fitted based on Rifampicin time serie intensity data with consideration of the typical delay of decay onset when using Rifampicin. Subsequently, probes/bins of equal properties, i.e halflife, intensity or polymerase velocity, are combined into segments by dynamic programming, independent of an existing genome annotation. This allows to detect transcript segments of different stability or transcriptional events within one annotated gene. In addition to the classic decay constant/half-life analysis, 'rifi' detects processing sites, transcription pausing sites, internal transcription start sites in operons, sites of partial transcription termination in operons, identifies areas of likely transcriptional interference by the collision mechanism and gives an estimate of the transcription velocity. All data are integrated to give an estimate of continous transcriptional units, i.e. operons. Comprehensive output tables and visualizations of the full genome result and the individual fits for all probes/bins are produced.
  
<br/>

<p align="center">
  <img src="https://github.com/CyanolabFreiburg/rifi/blob/main/vignettes/genome_fragments_plot.png"/>
</p>


<sub>  
<b>Figure 1:</b> Example visualization of an Rifampicin microarray experiment from <i>Synechocystis</i> PCC6803. A segment of the forward strand with its GenBank annotation is shown. 
The first track shows the delay of the onset of the decay for the individual probes. The delay should be linearly increasing for continuous transcripts, which are clustered by dynamic programming and indicated by matching colors and a trendline. A sudden delay increase between two segments indicates a transcription polymerase pausing site (PS), while a sudden decrease indicates a new (internal) transcriptional start site (iTSS). The slope of the delay segment allows to estimate the speed of the RNA Polymerase. Changes in the velocity are indicated by a "V". 
The second track shows the fitted half-life of the probes and the clustered half-life segments. If two segments within the same transcriptional unit have different half-life (HL) a processing/stabilization site for one or the other segment can be assumed.
The third track shows the intensity and the intensity segments at timepoint 0 (before Rifampicin addition). If two segments within the same transcriptional unit have different intensities this (FC) this could be either due to a partial termination (Ter) or a new transcriptional start site (NS). 
Significant events are assigned with an '*'. 
</sub>


# Installation 

### Dependencies

Rifi framework has the following dependencies, whereas the brackets indicate the version Rifi has been build and tested on. Make sure the requirements are satisfied by your system. 

  [devtools](https://www.rdocumentation.org/packages/devtools/versions/1.13.6) (>= 1.13.6)
  
  [roxygen2](https://cran.r-project.org/web/packages/roxygen2/index.html) (>= 7.1.12)
  
  [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)(>= 1.24.0)

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

<sub>
<b>Figure 2:</b> The general workflow of rifi. Blue boxes are data structures, green boxes are visual outputs, orange boxes are function. Incoming arrows represent the required input, outgoing arrows represent the output produced by the function. The orange frames state the main tasks of the function they refer to with dotted arrows
</sub>

## Positional Arguments

Rifi provides different functional arguments (subcalls) for individual procedures. These include `rifi_preprocess`, `rifi_fit`, `rifi_penalties`, `rifi_fragmentation`, `rifi_stats`, `rifi_summary` and `rifi_visualization`


## Input 

The input is a summarizedExperiment data format containing the RNA-seq or microarray time series data from rifampicin treated microorganisms, including all potential replicates, column wise as relative intensities (assay), e.g, Table 1 shows an example input.


  0       | ... | 20       | 0      | 1   |    |
  :---:   |:--: | :---:    | :--:   |:---:|:--:|
  8.673   | ... | 3.460    |53.855  | NA  | ...|
  8.974   | ... | 3.060    |188.61  | NA  | ...|
  18.612  | ... | 7.650    |163.83  | NA  | ...|
  15.476  | ... | 7.022    |243.12  | NA  | ...|
  100.674 | ... | 100.460  |198.52  | NA  | ...|
  250.645 | ... | 50.460   |165.122 | NA  | ...|
  200.541 | ... | 56.460   |87.235  | NA  | ...|


<sub>Table 1. Rifi input data table (assay). The columns contain the relative intensity measurements of all replicates, with the first column referring to timepoint zero, at or before the addition of rifampicin.
</sub>

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

The results of the analysis are stored on rowRanges and metaData. rowRanges contain fragments, events and statistical tests. MetaData combine several output from different functions and genome annotation as data frame format e.g.(Table .2 and Table. 3).

ID    | gene  | ...  | strand| TU     | half-life |...    |p_value_TI|
:---: | :--:  | :---:| :--:  | :---:  | :---:     | :---: | :---:    |   
1     | gene1 | ...  | +     | TU_1   | 0.63      |...    |0.002     |
2     | gene1 | ...  | +     | TU_1   | 0.63      | ...   |0.002     |  
3     | gene1 | ...  | +     | TU_2   | 1.22      |...    |0.002     |
5     | NA    | ...  | -     | TU_501 | 2.03      | ...   | NA       |  
8     | NA    | ...  | -     | TU_501 | 0.93      |...    | NA       |
12    | gene2 | ...  | -     | TU_502 | 1.55      | ...   | NA       |  


<sub>Table 2: The segment data frame output. All data given is the value averaged over the given fragment. Additional columns represented by the dots are additional information.</sub>

ID    | gene  | ...  | strand| event  | p_value  |
:---: | :--:  | :---:| :--:  | :---:  | :---:    |   
1     | gene1 | ...  | +     | ps     | 0.0003   |
2     | gene1 | ...  | +     | ps     | 0.5      | 
3     | gene1 | ...  | +     | iTSS   | 0.01     |
5     | NA    | ...  | -     | ps     | 0.006    | 
8     | NA    | ...  | -     | iTSS   | 0.007    |
12    | gene2 | ...  | -     | iTSS   | 0.007    |


<sub>Table 3: The event data frame output. All data given is the value averaged over the given fragment. Columns represented by the dots are additional information.</sub>

## Quick start

Two examples with their corresponding annotation gff file are available. One from RNA-seq (E.coli) and the other from microarrays (Synechocystis 6803). To run E.coli example, you would need to download the data cited below from rifi/data and the annotation file from rifi/inst/extdata/ on the same directory of your work. 

#### 1. Annotation file
gff_e_coli.gff3.gz
gff_synechocystis_6803.gff.gz

#### 2. SummarizedExperiment input data
example_input_e_coli.RData
example_input_synechocystis_6803.RData

```
wrapper_e.coli <- rifi_wrapper(inp = example_input_e_coli, cores = 2, path = "gff_e_coli.gff3.gz", bg = 0, restr = 0.01)
```

```
wrapper_synechocystis <- rifi_wrapper(inp = example_input_synechocystis_6803, 
cores = 20, path = "./gff_synechocystis_6803.gff.gz", bg = 4000, restr = 0.01)
```
#### check the result
rifi_output and a plot are generated in your directory.
rifi_output is a summarizedExperiment format containing assay, rowRanges, colData and metaData. More details could be found on vignette.

## Testing

Two data frames are provided to test Rifi: `example_input_e_coli` and `example_input_synechocystis_6803`

# Troubleshooting

contact jens.georg@biologie.uni-freiburg.de or create an issue
