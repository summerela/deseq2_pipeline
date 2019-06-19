# DESeq for Differential Chromatin Accessibility

This pipeline can be used to conduct a preliminary analysis of differences in chromatin accessibility between two user-defined groupings. The pipeline begins with a list of aggregation ID's and user-specified groupings to compare. The script
takes care of all pre-processing steps, runs a DESeq2 comparison and then outputs bed files of all differences, significant differences, and visualizations.

[More information on DESeq2](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf)

![workflow](/dhs_workflow_diagram.pdf?raw=true "Workflow")

## Run-time assumptions
- User has LIMS API token [setup](https://lims.stamlab.org/apitoken) and stored in ~/.bash_profile as LIMS_API_TOKEN
- Aggregation level cutcounts and called peak files are available in LIMS directory

## Installation/Requirements

    module purge
    module load deseq_pipeline

Note: If you are having issues loading modules, check $LD_LIBRARY_PATH

## To run:

#### 1. Create a text file of aggregate id's to analyze and their grouping metadata

    agg_id   condition   
    AG5536     WT          
    AG5537     DNMT3+      
    
If your sample does not have a LIMS-generated aggregation ID, you can also just provide file paths to the relevant cutcounts and peaks files. You can also provide a combination of aggregate ID's and file paths: 

    agg_id    condition    cutcounts                                          peaks
    sample1    WT          /home/cutcount_dir/sample1_cutcounts.bed.starch    /home/peaks_dir/sample1_peaks.bed.starch
    AG5536     HET 
    sample2    HET         /home/cutcount_dir/sample2_cutcounts.bed.starch    /home/peaks_dir/sample2_peaks.bed.starch

#### 2. Create a config file with following categories/variables:

    [paths]  
    
    # path to aggregation id's, one per line  
    agg_list=./metadata.txt  

    [transforms]  

    # choose custom for John Lazar's normalization script  
    # choose none to use DESeq filtering  
    norm_method=custom  

    # if using John Lazar's normalization method,
    # select true to use geomean, false for loess
    geomean=true  

    [deseq_options]

    # experimental design for creation of DESeq exp object
    experimental_design="condition + time"  

    # supply required contrast argument to DESeq
    contrasts="condition,WT,DNMT3+"  


#### 3. Run as:

    deseq.py -c config.txt

## Required Input
- metadata.txt: tsv format, one aggregate ID per line, w/header, ex: agg_id | condition | any_other_grouping
    - [Example metadata.txt](/deseq/examples/agg_id.txt)
- config.txt: Config file containing the following options.
    - [Simple example:](/deseq/examples/simple_config.txt) The minimum information required using default options.
    - [Complex example:](/deseq/examples/advanced_config.txt) Example with all configurable options.

## Options

### Pipeline options:
- Input Master (default= None):
    Do you already have a master DHS reference, or should it be created from your input samples?
    - "none" = create from samples
    - path/to/master_dhs.bed = use this file for a reference
- Threshold (default = 0):
    Remove any rows from master DHS reference that have loci in less than X samples
- Input Options: choose between using peaks or hotspots to count overlaps between your master list    
- Normalization Method (default = Custom):
    DESeq2 normalization is tailored for RNAseq data and does not fit well to differential DHS loci data because many regions will have no DHS sites, thus making sites with DHS seem over-important. To deal with this issue, John Lazar has created a normalization method that uses the Geometric Mean or Loess smoothing methods to provide a much better fit for DHS data. This method is recommended.
    - "custom" = Use John Lazar's custom normalization method
    - "none" = Use DESeq normalization
- Geomean (default = True):
    If using custom normalization method:
    - "true" = use Geometric Mean
    - "false" = use Loess smoothing

### DESeq2 Options
- Experimental design: Model used in DESeq2
    Design indicates the effects we want to measure. For example the following would examine the effect of condition and time on the data.
    - ex: "condition + time"
    The following would examine the effect of condition on the data, controlling for batch effect:
    - ex: "~batch + condition"
    The experimental design should only contain column names and values found in your agg_id.txt file
    - [More information](/deseq/examples/experimental_design.R)
- Contrasts: Allow you to examine the differences in the experimental design between to groupings. These should be specified as a comma separated list containing (in order): "column_name,group1,group2".
    - ex: "condition,treatment,time"
-Save R objects: Save intermediate R objects for reference, and to assist in running contrasts with the same experiemental design more quickly. TRUE to save, FALSE to skip.
- FDR: Set your FDR threshold. By default this is set to 0.01. 

## Example
Click [here](https://github.com/Altius/deseq_pipeline/blob/dev/deseq/examples/Walkthrough.md) to view an example walk-thru.
