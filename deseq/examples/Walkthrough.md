# DHS Analysis Pipeline
The DHS Analysis pipeline uses a combination of python, bedops and R’s DESeq2 package to search for differences in DHS between input samples. To make this easy to use, the pipeline has been loaded on the server as a module.

Because Bill has turned this into a module for you, you will never need to worry about loading any of its requirements/dependencies individually.

The following example uses a test dataset built for the Informatics 101 training courses and will help you learn how to use the pipeline. 

## 1. Copy the RCC data set to your directory
To walk you through how to use the pipeline, we will use the RCC dataset that you worked with in previous lessons.

From the command line, let’s get this information into your home directory and walk through the options.

    cp -R /home/selasady/public_html/RCC/* ./
    
Change into the your newly copied directory: 
   
    cd ./RCC/
    
## 2. Ensure your LIMS token is properly setup
The first time you run the pipeline, you'll need to make sure your LIMS api url and token are properly setup. Follow [this tutorial](https://lims.altiusinstitute.org/apitoken).

## 3. Setup your environment variables
First, we’ll clean out all the modules you might already have loaded to ensure we don’t run into any conflicts:
    
    module purge # unload all currently loaded modules
    module load deseq_pipeline # load DHS script dependencies
    
## 4. Setup your Experiment Metadata
The pipeline requires two input files:   
- A text file containing the aggregation ID and grouping information for your experiment 
- A config file that let’s you specify experimental paramaters  

### Sample Metadata
In your RCC folder, the metadata.txt file contains a tab-separated list of aggregate ID’s to analyze, and their associated grouping information. 

Let’s take a look at that file:
    
    cat ./metadata.txt  
    agg_id condition AG3819 RCC
        AG3898    TUB
        AG3899    TUB
        AG6993    RCC
        
<b> Note: The variables in your condition column are case-sensitive </b>  

You can add more complexity to your experimental design by adding additional columns, and you can even add a separate column for labels: 
- If you create a label column, you must specify it in your config file  
- Any label column you create will appear on the graphs   
- The metadata file must be tab separated   
- Your metadata file must contain “agg_id” and “condition” columns  


    agg_id  condition    time   label
    AG3819  RCC          day1   rccDay1 
    AG3898  TUB          day2   tubDay2 
    AG3899  TUB          day1   tubDay1 
    AG6993  RCC          day2   rccDay2  
    
### Config File
In the RCC directory, config.txt contains the minimum amount of options required to run the pipeline.

You can always find an example on the [Altius GitHub repo](https://github.com/Altius/DHS_pipeline/blob/dev/deseq/examples/simple_config.txt).
      
    [paths]  
      
    # path to metadata file with agg_id's
    agg_list=./metadata.txt
      
    [transforms]
    # if a DHS is not found in this many samples, remove it
      
    # enter 0 for no filter
    sample_threshold=1
      
    [deseq_options]
      
    # contrast argument for DESeq (column,var2,var1) # example "condition,dnmt3, wt"    
    # note that these variables must match the case and spelling in your metadata
    # the base variable (wt, day0, etc) must be listed last
    contrasts=condition,RCC,TUB

The DHS pipeline is designed to take all input observations and pull out a pair-wise analysis of two variables specfied with the contrasts argument. [Click here to read more about experimental design and contrasts](https://github.com/Altius/DHS_pipeline/blob/dev/deseq/examples/experimental_design.R). 

<b>Additional options </b> The config file above leaves a lot of options set at default. All advanced options can be found [here](https://github.com/Altius/DHS_pipeline/blob/dev/deseq/examples/advanced_config.txt).
    
# Run the pipeline
Now that we have our experiment setup, we can run the pipeline by pointing the module to our configuration file:

    deseq.py -c ./config.txt
    
# View the Results
Results are output into the following directories:  

### input_files
• symlinks to peak and count files  
• metadata.txt with input agg id’s, condition, sample ID, library number, taxonomy, assay,
view URL and file locations  

### Counts
• master_dhs.bed  
• individual overlap count files  
• annotated, combined counts files   
• raw counts and normalized counts   
• normalization factors

### Results
- deseq2_results.bed: all results found for the comparison between your contrasts  
- deseq2_significant_results.bed: significant results (with your selected FDR threshold) for your contrasts  
- Saved R objects  
- Pics directory containing PDF and PNG versions of heatmaps, ma plot, counts of signifi- cantly different DHS sites, and top 8 loci with greatest fold change  

Log files are created and time-stamped for python and R scripts; check these if anything went wrong, or to view input parameters. 

# Getting information and help
The Altius GitHub repo for the project contains:   
- [Runtime instructions](https://github.com/Altius/DHS_pipeline)   
- Example [metadata](https://github.com/Altius/DHS_pipeline/blob/dev/deseq/examples/agg_id.txt) and [config](https://github.com/Altius/DHS_pipeline/blob/dev/deseq/examples/advanced_config.txt) files 
- Additional information on [experimental design](https://github.com/Altius/DHS_pipeline/blob/dev/deseq/examples/experimental_design.R)  

The [DESeq2 manual](https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) provides additional resources for setting up your experimental design and interpreting the results.

If you need help, please contact me (Summer) or ask compbio!