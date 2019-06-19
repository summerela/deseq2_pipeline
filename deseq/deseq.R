#!/usr/bin/env R

# Dependencies:
#    - module load R/3.2.5
#    - module laod cairo/1.14.2

# DESeq2 Options:

# - Experimental design: Model used in DESeq2
#    Design indicates the effects we want to measure. For example the following would examine the effect
#    of condition and time on the data.
#    - ex: "condition + time"
#    The following would examine the effect of condition on the data, controlling for batch effect:
#    - ex: "~batch + condition"
#    The experimental design should only contain column names and values found in your agg_id.txt file
#    - More information [here](/examples/experimental_design.R)
#
# - Contrasts: Allow you to examine the differences in the experimental design between to groupings.
# These should be specified as a comma separated list containing (in order): "column_name,group1,group2".
#    - ex: "condition,treatment,time"
#
# -Save R objects: Save intermediate R objects for reference, and to assist in running contrasts with the
# same experiemental design more quickly. TRUE to save, FALSE to skip.

# To run standalone:
# load modules, then run:
#Rscript /home/selasady/deseq_pipeline/deseq/deseq.R
# -d "condition"
# -m master_dhs.bed
# -c counts directory
# -i agg_ids.txt
# -f combined_counts.bed.size_factors
# -o ./results/
# -s true
# -a contrasts
# -p .05


#####################
## Parse User Args ##
#####################

# check user args first
library(optparse)

### Fetch command line arguments ###
option_list = list(
    make_option(c("-d", "--design"), type="character", default=NULL,
        help="enter experimental design as string 'condition + taxonomy'.
        Put the variable of interest at the end of the formula, control level
        should be first level.", metavar="exp_design"),
    make_option(c("-m", "--master"), type="character", default=NULL,
        help="file containing master dhs hotspots", metavar="master_dhs"),
    make_option(c("-c", "--counts"), type="character", default=NULL,
        help="Path to file containing hotspot counts", metavar="counts"),
    make_option(c("-i", "--info"), type="character", default=NULL,
        help="metadata file for input samples in tsv file_path | sample_name | condition/group", metavar="input_metadata"),
    make_option(c("-f", "--factors"), type="character", default=NULL,
        help="enter path to custom norm .size_factors file, or 'none' to skip normalization",
        metavar="norm_factors"),
    make_option(c("-o", "--outdir"), type="character", default=NULL,
        help="path to output directory", metavar="output_dir"),
    make_option(c("-a", "--contrasts"), type="character", default=NULL,
        help="Enter contrasts in the format 'colname,contrast1,contrast1'.",
        metavar="contrasts"),
    make_option(c("-s", "--saverds"), type="character", default='yes',
        help="Enter 'no' to prevent script from saving deseq and deseq results objects.",
        metavar="save_objects"),
    make_option(c("-l", "--label_col"), type="character", default="condition",
        help="Column name to use as graphing labels. Default='condition'", metavar="label_col"),
    make_option(c("-p", "--pval_cutoff"), type="double", default=0.01,
        help="Enter ajdusted Pvalue cutoff (FDR). Default='0.01'", metavar="pval_cutoff"),
    make_option(c("-t", "--threshold"), type="double", default=0,
        help="If a hotspot is not found in at least x samples, that loci will be filtered out.'", metavar="threshold"),
    make_option(c("-n", "--input_method"), type="character", default="hotspots",
        help="Select baseline comparison of either hotspots or peaks", metavar="input_method")
	 )

# capture user arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# setup global variables
exp_design = unlist(strsplit(gsub(" |\\~", "", opt$design), "\\+"))
contrast_col = unlist(strsplit(as.character(gsub(" ", "", opt$contrast)), ","))[1]
contrast1 = unlist(strsplit(as.character(gsub(" ", "", opt$contrast)), ","))[2]
contrast2 = unlist(strsplit(as.character(gsub(" ", "", opt$contrast)), ","))[3]
label_col = opt$label_col
combined_counts = opt$counts
pval_cutoff = as.numeric(opt$pval_cutoff)

#####################
## Load Libraries  ##
#####################

# check if package is is_installed
is_installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

# load libraries, install if needed
load_or_install<-function(package_names) {
     for(package_name in package_names){
         if(!is_installed(package_name)){
             install.packages(package_name,repos="http://cran.fhcrc.org/")
             }
     library(package_name,character.only=TRUE,quietly=TRUE,verbose=FALSE)
     }
}

# load libraries
load_or_install(c("optparse", "DESeq2", "ggplot2","readr", "getopt",
    "BiocParallel", "parallel", "gplots", "futile.logger", "optparse",
    "RColorBrewer", "gtools", "genefilter", "gridExtra", "grid"))

###################
## Setup Logging ##
###################

# capture today's date
today = format(Sys.Date(), "%Y%m%d")

# create log file path
log_file = paste0(as.character(getwd()), "/", today, "_deseq2_R.log")

# script log created as today_deseq2_R.log in working directory

# if there is an existing log file, erase it
if (file.exists(log_file)){
    file.remove(log_file)
}

# setup a logger to handle TRACE level and greater stdout (ERRORS, WARNINGS, etc.)
# tee outputs to both file and stdout
flog.logger(name="deseq_log", TRACE, appender=appender.tee(log_file))
# format output to include -l log level, -n namespace, -f function, -m message
layout <- layout.format('[~l] [~n.~f] ~m')
# apply layout
flog.layout(layout)

flog.info("\n \n ############   Loaded Packages  ################## \n", name="deseq_log")

# output attached packages and versions
packinfo <- installed.packages(fields = c("Package", "Version"))
flog.info(capture.output(sessionInfo()), name="deseq_log")

# backup: catch all exceptions and pass to log file
error_handle <- file(log_file, open = "a")
sink(error_handle, type =  "message", append=TRUE)

# function to stop execution on error and log error message
log_stop <- function(error_message){
    stop(call=TRUE, geterrmessage())
    flog.error(error_message, name="deseq_log")
    flog.error(geterrmessage(), name="deseq_log")
}

###########################
## Setup mutliprocessing ##
###########################

# use 35% of the cores
use_cores = round(detectCores() * .35)
register(MulticoreParam(use_cores), default=TRUE)
flog.info("Using %i cores.", use_cores, name="deseq_log")

###################################
## Check input file paths exist  ##
###################################

# check if out_dir exists, if not create
dir.create(file.path(opt$outdir), showWarnings = FALSE)
flog.info(paste0("Output directory: ", opt$outdir), name="deseq_log")

# check that user specified files exist
check_files <- function(input_file){
    if (file.exists(input_file)){
        flog.info(paste0(input_file, " path exists."), name="deseq_log")
    } else {
        log_stop(paste0(input_file, " path does not exist."))
    }
}
# create list of input files and check each
infile_list = c(opt$master, combined_counts, opt$info)
lapply(infile_list, check_files)

############################
#### Log input parameters ##
############################
flog.info("\n \n ############   Experimental Parameters   ################## \n", name="deseq_log")
flog.info(paste0("Master DHS file: ", opt$master), name="deseq_log")
flog.info(paste0("Metadata: ", opt$info), name="deseq_log")
flog.info(paste0("Counts file: ", combined_counts), name="deseq_log")
flog.info(paste0("Experimental Design: ", opt$design), name="deseq_log")
flog.info(paste0("Contrasts: ", opt$contrasts), name="deseq_log")
flog.info(paste0("Label column: ", opt$label_col), name="deseq_log")
flog.info(paste0("FDR cutoff: ", pval_cutoff), name="deseq_log")
flog.info(paste0("Analysis Background: ", opt$input_method), name="deseq_log")

if (opt$factors != "none"){
    check_files(opt$factors)
    flog.info(paste0("Custom normalization factors: ", opt$factors), name="deseq_log")
    } else {
        flog.info(paste0("Using DESeq normalization"), name="deseq_log")
}

#########################
## load required files ##
#########################

flog.info("\n \n ############   Loading Files   ################## \n", name="deseq_log")

read_files <- function(infile_path, delim_val, col_val){
    cat(paste("Reading input file:", infile_path, "\n", sep=" "))
    if (file.exists(infile_path)){
        out_df = as.data.frame(read_delim(infile_path, trim_ws=TRUE, delim=delim_val, col_names=col_val),
        stringsasfactors = FALSE)

        flog.info("Loaded: %s.", infile_path, name="deseq_log")
        return(out_df)
    } else {
        flog.error("File %s does not exist.", infile_path, name="deseq_log")
        log_stop("File", infile_path, "does not exist.", sep=' ')
        }
    }
# read in master dhs list
master = read_files(opt$master, delim_val='\t', col_val=F)
colnames(master) <- c("chrom", "start", "stop")

# read in metadata
coldata_df = read_files(opt$info, delim_val='\t', col_val=T)

## read in counts table
counts_tbl = read_files(combined_counts, delim_val='\t', col_val=F)
rownames(counts_tbl) <- paste0(master$chrom, ":", master$start, "_", master$stop)

###############################
## setup experimental design ##
###############################
flog.info("\n ############   Formatting Input   ################## \n", name="deseq_log")


# check that experimental design column(s) exists in input metadata
check_design <- function(input_var){
    x = gsub("~", "", input_var)

    flog.info(x, name="deseq_log")

    if (x %in% colnames(coldata_df)){
        flog.info(paste0("Column ", x, " found in ", opt$info), name="deseq_log")
    }else {
        flog.error(paste0("Column ", x, " not found in ", opt$info), name="deseq_log")
        log_stop(paste0("Column ", x, " not found in ", opt$info))
    }
}

flog.info("\n Checking experimental design: \n", name="deseq_log")
lapply(exp_design, check_design)


# format input metadata
set_exp_design <- function(){

    flog.info("Setting up experimental design parameters...", name="deseq_log")

    # if label_col != condition, add label col to exp_df for labelling graphs
    if (tolower(label_col) != "condition"){

        #subset columns scpecified in experimental design, coerce to factors and add label
        design_cols = as.vector(paste("coldata_df$", exp_design, sep=""), mode="list")
        exp_df_cmd <- paste0("data.frame(", paste("coldata_df$agg_id,", design_cols, collapse=","), ")" )

        exp_df = eval(parse(text=exp_df_cmd))

        df_label = paste0("coldata_df$", label_col)
        exp_df$label = eval(parse(text=df_label))

        # create row ids based on agg_id and label_label_col
        names(exp_df) <- gsub("coldata_df.", "", names(exp_df))
        rownames(exp_df) <- paste(exp_df$agg_id, exp_df$label, sep="_")

        # set as factors
        exp_df[] = lapply(exp_df, as.factor)

        return(exp_df)

    } else {

        # create row ids based on chrom/start/stop positions
        design_lbl = as.vector(paste("coldata_df$", exp_design, sep=""), mode="list")
        row_ids = paste("paste(", paste("coldata_df$agg_id", paste(design_lbl, collapse=","), sep=","), ", sep='_')", sep="")

        #subset columns scpecified in experimental design and coerce to factors
        exp_df_cmd <- paste0("data.frame(", paste(design_lbl, collapse=","), ")" )

        exp_df = eval(parse(text=exp_df_cmd))

        # set rownames and colnames
        rownames(exp_df) <- eval(parse(text=row_ids))
        names(exp_df) <- gsub("coldata_df.", "", names(exp_df))

        # set as factors
        exp_df[] = lapply(exp_df, as.factor)

        return(exp_df)
    }

}

exp_df = ftry(set_exp_design(), error=log_stop("Error creating experimental design dataframe"))

#######################
#### check contrasts ##
#######################

# check that contrast column exists in metadata file
flog.info("\n Checking contrast columns... \n", name="deseq_log")
check_design(contrast_col)

# check that contrast values exist in metadata file
con_col = as.character(paste0("coldata_df$", contrast_col))
val_list = as.list(eval(parse(text=con_col)))

if (contrast1 %in% val_list){
    flog.info(paste0("Contrasts: ", contrast1, " variable found in ", con_col), name="deseq_log")
    }else {
        log_stop(paste0("Contrasts: ", contrast1, " variable not found in ", con_col))
    }

if (contrast2 %in% val_list){
    flog.info(paste0("Contrasts: ", contrast2, " variable found in ", con_col), name="deseq_log")
    }else {
        log_stop(paste0("Contrasts: ", contrast2, " variable not found in ", con_col))
    }


####################################
###### label, format counts table ##
####################################

format_counts <- function() {

    # convert any na to 0
    counts_tbl[is.na(counts_tbl)] <- 0

    # set colnames to match metadata table
    names(counts_tbl) <- rownames(exp_df)

    return(counts_tbl)
}

counts_df = ftry(format_counts(), error=log_stop("Error formatting counts table."))

##########################
### Create DESeq object ##
##########################

create_deseq <- function(){

    flog.info("\n ############   Creating DESeq Object   ################## \n", name="deseq_log")

    print ("Creating DESeq object.")

    # create design formula from user arg opt$design
    design_string = paste("~", as.character(opt$design), sep="")

    # create deseq object
    DESeq_tbl <- DESeqDataSetFromMatrix(countData=counts_df,
        colData=exp_df, design=as.formula(design_string))

    flog.info("DEseq object created as DESeq_tbl", name="deseq_log")

    colnames(DESeq_tbl) <- colnames(counts_df)

    return(DESeq_tbl)
}

# create DESeq experiment object
DESeq_tbl <- ftry(create_deseq(),error=log_stop("Error running DESeqDataSetFromMatrix."))


#################
#### Normalize ##
#################
# normalize results
load_norms <- function(){

    norm_factors <- as.matrix(read_delim(opt$factors, delim='\t', col_names=FALSE))
    normalizationFactors(DESeq_tbl) <- norm_factors

    flog.info("'Applied normalization factors from %s", opt$factors, name="deseq_log")

   return(DESeq_tbl)

}

# if custom norm factors provided, add to DESeq_tbl
if (tolower(opt$factors) != 'none'){
    DESeq_tbl = ftry(load_norms(), error=log_stop("Error applying custom normalization factors."))
}

##########################
######### Run DESeq  #####
##########################

## run DESeq
run_deseq <- function(){
    flog.info("\n ############   Running DESeq   ################## \n", name="deseq_log")
    print("Running DESeq")

    # if filtering selected, filter before running DESeq
    if (opt$threshold > 0){
        log_msg = paste0("\n Filtered out loci that don't appear in at least ", opt$threshold, " sample(s). \n")
        flog.info(log_msg, name="deseq_log")
        dds <- estimateSizeFactors(DESeq_tbl)
        idx <- rowSums( counts(dds, normalized=FALSE) > 0 ) >= opt$threshold
        dds <- dds[idx,]
        DEout <- DESeq(dds, parallel=TRUE)

    } else {
        DEout <- DESeq(DESeq_tbl, parallel=TRUE)
    }

    flog.info("Created DESeq table object as DEout.", name="deseq_log")
    return(DEout)
}

# create DESeq object
DEout <- ftry(run_deseq(), error=log_stop("Error running DESeq."))

# check that DEout contains results
if (length(DEout) < 1) {
    print("No results found.")
    flog.info("No results found.", name="deseq_log")
    quit()
} else{
    # save normalized counts to a file
    norm_counts_file = paste0(opt$outdir, "/normalized_counts.tsv")
    norm_counts = data.frame(counts(DEout, normalized=TRUE))
    write.table(norm_counts,file=norm_counts_file, quote = F, sep = '\t', row.names = F, col.names = T)
}

# save DEseq object to results directory
if (tolower(opt$saverds) != "no") {

    deseq_out = paste0(opt$outdir, "/deseq_object.rds")
    print(paste0("Saving deseq object as ", deseq_out))

    saveRDS(DEout, deseq_out)

    flog.info(paste0("Saving deseq object as ", deseq_out), name="deseq_log")
}



# uncomment for testing
# DEout <- readRDS("./results/deseq_object.rds")
#################################
######### Running Contrasts #####
#################################

flog.info("\n ############   Calculating Results   ################## \n", name="deseq_log")

# Run DESeq results for input contrasts, independant filtering enabled for multiple test correction
contrast_cmd = paste0("DEres <- results(DEout, independentFiltering = FALSE, cooksCutoff = Inf,
                      format='DataFrame', parallel=TRUE, contrast=c('",
                      paste(sapply(strsplit(opt$contrasts, ','), paste0), collapse="','"), "'))")

# pull out results for the contrasts
DEres = ftry(eval(parse(text=contrast_cmd)), error=log_stop("Error running results(DEout)."))

####################################
######### Save Results Table(s) #####
####################################

# saving all contrast results to bed file format
results_to_table <- function(){

    flog.info("\n ############   Saving Results   ################## \n", name="deseq_log")

    results_name = paste0(as.character(opt$outdir), "/deseq_results.bed9")

    # create annotated results table from contrast results
    withStats <- data.frame((cbind(chrom=master$chrom, start=master$start, stop=master$stop,
           log2FC=as.numeric(as.character(DEres$log2FoldChange)),
           abs_fc=abs(as.numeric(as.character(DEres$log2FoldChange))),
           pval=as.numeric(as.character(DEres$pvalue)),
           adjP=as.numeric(as.character(DEres$padj)),
           mean=as.numeric(as.character(DEres$baseMean)),
           stdError=as.numeric(as.character(DEres$lfcSE)))), stringsAsFactors=FALSE)

    # add rownames and comment to column names
    rownames(withStats) <- paste0(withStats$chrom, ":", withStats$start, "_", withStats$stop)
    names(withStats)[1] <- "#chrom"

    # reorder table by lowest adjusted pvalue and greatest absolute fc
    all_results_tbl <-withStats[order(withStats$adjP, -as.numeric(withStats$abs_fc), decreasing=FALSE),]
    # save results to deseq_results.bed9 file in output directory
    write.table(all_results_tbl,file=results_name, quote = F, sep = '\t', row.names = F, col.names = T)
    flog.info("All results saved to %s", results_name, name="deseq_log")

    return(all_results_tbl)
}

# create dataframe of all contrast results ordered by smallest adjP and greatest fc
all_results_tbl = results_to_table()

# save significant results to bed file
subset_results <- function(){

    sig_file_name = paste0(opt$outdir, "/deseq_sig_results.bed9")

    # subset for significant results adjP < FDR cutoff, order, and save to a bed file
    sig_tbl = all_results_tbl[as.numeric(as.character(all_results_tbl$adjP)) <= pval_cutoff,]

    # order table by lowest FDR
    sig_results = sig_tbl[order(as.numeric(sig_tbl$adjP), decreasing=FALSE),]
    sig_results = na.omit(sig_results)

    if (nrow(sig_results) > 0){
        write.table(sig_results, file=sig_file_name, quote = F, sep = '\t', row.names = F, col.names = T)
        flog.info("Significant results saved to %s", sig_file_name, name="deseq_log")

        return(sig_results)
    # stop execution if there are no results
    }else {
    flog.info("No significant results found!", name="deseq_log")
    stop("\n No significant results found! \n")
}

}

# save all results to a file and return significant results as dataframe
sig_results = ftry(subset_results(), error=log_stop("Error saving results to bed file."))

############################
##### Log Transformation  ##
############################
flog.info("\n ############   Log Transforming Results   ################## \n", name="deseq_log")

# remove hash from column names for graphing
names(sig_results) = c("chrom", "start", "stop", "log2FC", "abs_fc", "pval", "adjP", "mean", "stdError")

# subset the counts of DEout for just the top 50 loci found in the significant results table based on contrats
DEres_ordered <- rownames(head(sig_results[order(as.numeric(sig_results$abs_fc), decreasing=TRUE),], 50))

# normalized count data
ncounts <- 1.0 + counts(DEout, normalized=TRUE)[DEres_ordered,]
# rlog transformed data
rld = rlog(DEout)[DEres_ordered,]

# save log transformed object
if (tolower(opt$saverds) != 'no') {
    outfile = paste(opt$outdir, "log_transform.rds", sep="/")
    saveRDS(rld, outfile)
    flog.info(paste("Log transform saved as ", outfile, sep=""), name="deseq_log")
}

##########################
##Create pics directory ##
##########################

# # check if pics dir exists, if not create
pic_dir = as.character(paste0(opt$outdir, "/pics"))
dir.create(file.path(pic_dir))
flog.info(paste0("Created pic directory at: ", pic_dir), name="deseq_log")

################
### Heatmaps ###
################

# create rlog matrix for plotting heatmaps
rld_sig <- assay(rld)

flog.info("\n ############   Creating Visualizations   ################## \n", name="deseq_log")

# setup color palettes
color_range = as.numeric(max(ncounts) - min(ncounts))
hmcol <- colorRampPalette(brewer.pal(9, "Blues"))(color_range)
sdcol <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(DEout$condition)))[DEout$condition]
sdcol_sidebar <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(DEout$condition)))

# heatmap function
plot_heatmap <- function(in_data, hm_title){

    heatmap.2(in_data,
        col=hmcol,
        ColSideColors=sdcol,
        trace="none",
        srtCol=315,
        adjCol=c(0,1),
        key=F,
        labRow=rownames(in_data),
        cexCol = .85,
        margin=c(10, 11),
        dendrogram="row",
        # symbreaks=TRUE,
        Colv=FALSE,
        na.rm=TRUE,
        main=hm_title)

    par(xpd=TRUE)      # plot outside the main plotting area , mar = par()$mar + c(0,0,0,7)
    legend("topright",      # location of the legend on the heatmap plot
        legend = levels(DEout$condition), # category labels
        col = sdcol_sidebar,  # color key
        lty = 1,            # line style
        lwd = 3,
        cex = 0.75,
        ncol = round(length(levels(DEout$condition))/2, digits=0)
)
}

# plot normalized count heatmap
plot_count_heatmap <- function(){

    flog.info("Plotting heatmap of normalized counts for top 50 DHS by fold change", name="deseq_log")
    plot_title = paste("Counts for top 50 Loci\n between", contrast1, "and",
        contrast2, "with the greatest Fold-Change", sep=" ")
    plot_heatmap(ncounts, plot_title)
}

# plot log count heatmap
plot_log_count_heatmap <- function(){

    flog.info("Plotting heatmap of Log Counts for for top 50 DHS by fold change", name="deseq_log")
    plot_title = paste("Top 50 Most Differential Loci\n between", contrast1, "and",
        contrast2, "by Rlog Fold-Change", sep=" ")
    plot_heatmap(rld_sig, plot_title)
}

################
### MA Plots ###
################

# MA plots
plot_ma <- function(input_table) {
    plotMA(input_table, alpha=pval_cutoff, xlab="mean DNase1 signal level", main=paste0(nrow(sig_results),
        " significantly differential loci \n between ", contrast1, " and ", contrast2, " (FDR < ", pval_cutoff, ")"))
}

# ma plot results table
plot_ma_results <- function(){
    plot_ma(DEres)
}


# plot p-value and foldchange distribution
plot_sig <- function(input_table){
    pval_ymax = nrow(input_table[!is.na(input_table$pvalue) & (input_table$pvalue < .05),]) * 1.1
    adjp_ymax = nrow(input_table[!is.na(input_table$padj) & (input_table$padj < .05),]) * 1.1
    fc_pos_ymax = nrow(input_table[!is.na(input_table$log2FoldChange) & (input_table$log2FoldChange > 1),]) * 1.1
    fc_neg_ymax = nrow(input_table[!is.na(input_table$log2FoldChange) & (input_table$log2FoldChange < -1),]) * 1.1

    par(mar=c(10,4,4,3))

    # create dataframe of pvalue counts
    pvals = data.frame(c(nrow(input_table[!is.na(input_table$pvalue) & (input_table$pvalue < .01),]),
        nrow(input_table[!is.na(input_table$pvalue) & (input_table$pvalue < .05),])),
        c("p < 01", "p < 05"), stringsAsFactors=FALSE
        )
    colnames(pvals) = c("counts", "label")

    # plot pvals
    pval_plot <- ggplot(pvals, aes(x=label, y=counts)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = sprintf("%.2f%%", counts/nrow(input_table) * 100)),
            vjust = -.5) +
        xlab("") +
        ylim(0, pval_ymax) +
        ggtitle("Nominal p-value \n")

    # create dataframe of adjusted pvalue counts
    adj_pvals = data.frame(c(nrow(input_table[!is.na(input_table$padj) & (input_table$padj < .01),]),
        nrow(input_table[!is.na(input_table$padj) & (input_table$padj < .05),])),
        c("adjP < 01", "adjP < 05"), stringsAsFactors=FALSE)
    colnames(adj_pvals) <- c("counts", "label")

    # plot adjusted pvals
    adjpval_plot <- ggplot(adj_pvals, aes(x=label, y=counts)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = sprintf("%.2f%%", counts/nrow(input_table) * 100)),
            vjust = -.5) +
        xlab("") +
        ylim(0, adjp_ymax) +
        ggtitle("False Discovery Rate (FDR) \n")

    # create dataframe of positive log2FoldChange counts
    logfc_pos = data.frame(c(nrow(input_table[!is.na(input_table$log2FoldChange) & (input_table$log2FoldChange > 2),]),
        nrow(input_table[!is.na(input_table$log2FoldChange) & (input_table$log2FoldChange > 1),])),
        c("log2fc > 2", "log2fc > 1"), stringsAsFactors=FALSE)
    colnames(logfc_pos) <- c("counts", "label")


     # create dataframe of log2FoldChange counts
     logfc_neg = data.frame(c(nrow(input_table[!is.na(input_table$log2FoldChange) & (input_table$log2FoldChange < -2),]),
         nrow(input_table[!is.na(input_table$log2FoldChange) & (input_table$log2FoldChange < -1),])),
         c("log2fc < -2", "log2fc < -1"), stringsAsFactors=FALSE)
     colnames(logfc_neg) <- c("counts", "label")

    # plot positive log foldchange
    fc_pos_plot <- ggplot(logfc_pos, aes(x=label, y=counts)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = sprintf("%.2f%%", counts/nrow(input_table) * 100)),
            vjust = -.5) +
        xlab("") +
        ylim(0, fc_pos_ymax) +
        ggtitle("Positive log2(fold-change) \n")

    # plot negative log foldchange
    fc_neg_plot <- ggplot(logfc_neg, aes(x=label, y=counts)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = sprintf("%.2f%%", counts/nrow(input_table) * 100)),
            vjust = -.5) +
        xlab("") +
        ylim(0, fc_neg_ymax) +
        ggtitle("Negative log2(fold-change) \n")

    # create overall plot title
    overall_title = paste("Number of differential loci between", contrast1,
        "and", contrast2, sep=" ")

    # plot significance
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(3,2, heights = unit(c(1, 4, 4), "null"))))
    grid.text(overall_title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
    print(pval_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(adjpval_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
    print(fc_pos_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
    print(fc_neg_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
}

# significance plot for results table
plot_sig_results <- function(input_table){

    flog.info("Plotting Significance Distribution for Results Table", name="deseq_log")
    plot_sig(DEres)
}

# plot top 30 peaks with greatest fold change
plot_counts <- function() {

    # pull out peaks with lowest adjusted pvals
    peaks <- data.frame(peak_names = rownames(head(sig_results, 30)), stringsAsFactors=FALSE)
    data = counts_df[which(rownames(counts_df) %in% peaks$peak_names),]
    pval_df = sig_results[which(rownames(sig_results) %in% peaks$peak_names),]

    # create matrix of labels for each loci
    loci_labs = as.matrix(rownames(data))
    # create matrix of labels for each input sample
    group_labs = as.matrix(colnames(data))


    for (i in 1:length(loci_labs)) {
        plot(t(data[i,]),
        pch=16,
        type="h",
        lwd = 10/(length(group_labs)/5), # decrease size of bar with increase in input samples
        lend=1,
        col = sdcol,
        xaxt="n",
        yaxt="n",
        xlab= "",
        ylab="",
        xaxs="r",
        yaxs="r",
        xlim = c(1, length(group_labs)),
        ylim= range(t(data[i,]), na.rm=T),
        axes=FALSE)

        legend("topright",
        legend = round(-log10(as.numeric(pval_df$adjP))[i], digits=2), # category labels
        lty = 0,
        bty = "n",
        lwd = 3,
        cex = .8)

        # plot labels on outer margins of either column
        if (i %% 2 == 0){mtext(loci_labs[i], side=4, line=.6, outer=FALSE, cex=.6, las=2)
            } else {
                mtext(loci_labs[i], side=2, line=.6, outer=FALSE, cex=.6, las=2)
            }
        # draw a grey box around each plot
        box(col="grey")
        # plot x axis labels on bottom plots in each column
        if (i == 29 | i == 30){text(1:length(group_labs), y=0,
            labels=group_labs,
            srt=45,
            las=2,
            pos=1,
            offset=2,
            cex=1,
            xpd = NA)
        }
        title_text = paste0("Loci with lowest adjusted P-values between ", contrast1, " and ", contrast2)
        mtext(title_text, side=3, outer=TRUE, line=.5, cex=1)
        mtext("Normalized Counts", side=2, outer=TRUE, line=12, cex=1)

    }
}

#######################
####### Plots to SVG ##
#######################

# output each plot to a svg file
run_svg_plot <- function(input_plot){

    # create output file name
    outfile = paste(pic_dir, paste0(input_plot, ".svg"), sep="/")

    svg(file=outfile, width=10,height=8)

    # run the plot
    plot_name = paste(paste("plot", input_plot, sep="_"), "()", sep="")
    eval(parse(text=plot_name))
    dev.off()
}

plot_svgs <- function() {

    flog.info("\n ############   Generating SVG Plots   ################## \n", name="deseq_log")
    flog.info("SVG plots output to %s", pic_dir, name="deseq_log")

    ftry(run_svg_plot("ma_results"), error=log_stop("Error saving MA plot to svg."))
    ftry(run_svg_plot("sig_results"), error=log_stop("Error saving significance plot to svg."))

    # the following plots should be created only if significant results are found
    ftry(run_svg_plot("count_heatmap"), error=log_stop("Error saving heatmap to svg."))
    ftry(run_svg_plot("log_count_heatmap"), error=log_stop("Error saving log heatmap to svg."))

    outfile= paste0(pic_dir, "/dhs_counts.svg")
    svg(file=outfile, width=12,height=16)
    par(mfrow=c(15,2), mar=c(0,0,1.5,.5), oma=c(6,14,2,10) + 0.1)
    ftry(plot_counts(), error=log_stop("Error saving Count Distribution plot to svg."))


}

plot_svgs()

#######################
####### Plots to PDF ##
#######################

# output each plot to a pdf file
run_pdf_plot <- function(input_plot){

    # create output file name
    outfile = paste(pic_dir, paste0(input_plot, ".pdf"), sep="/")

    pdf(file=outfile, width=10,height=8)

    # run the plot
    plot_name = paste(paste("plot", input_plot, sep="_"), "()", sep="")
    eval(parse(text=plot_name))
    dev.off()
}

plot_pdfs <- function() {

    flog.info("\n ############   Generating PDF Plots   ################## \n", name="deseq_log")
    flog.info("PDF plots output to %s", pic_dir, name="deseq_log")

    ftry(run_pdf_plot("sig_results"), error=log_stop("Error saving significance plot to PDF."))
    ftry(run_pdf_plot("ma_results"), error=log_stop("Error saving MA plot to PDF."))

    # the following plots should be created only if significant results are found
    ftry(run_pdf_plot("count_heatmap"), error=log_stop("Error saving heatmap to PDF."))
    ftry(run_pdf_plot("log_count_heatmap"), error=log_stop("Error saving log heatmap to PDF."))

    outfile= paste0(pic_dir, "/dhs_counts.pdf")
    pdf(file=outfile, width=12,height=16)
    par(mfrow=c(15,2), mar=c(0,0,1.5,.5), oma=c(6,14,2,10) + 0.1)
    ftry(plot_counts(), error=log_stop("Error saving Count Distribution plot to PDF."))

}

plot_pdfs()

####################
### Close handles ##
####################

last <- function(){
    print("DESeq2 analysis script complete. Check log to view any possible errors.")
    for(i in seq_len(sink.number())){
        sink(NULL)
    }
    q()
}
last()