#!/usr/bin/env Rscript
#' @title Command-Line Interface for gimmecpg
#' @description Main script for gimmecpg
#' @import optparse
#' @import polars
#' @import h2o


# library(itertools)
# library(future.apply)


##########################
# Command line arguments #
##########################

gimmecpg <- function() {
  option_list <- list(
    optparse::make_option(c("-i", "--input"),
      action = "store",
      default = NA,
      type = "character",
      help = "Path to bed files"
    ),
    optparse::make_option(c("-n", "--name"),
      action = "store",
      default = NA,
      type = "character",
      help = "Name of specific files. Multiple filenames can be separated with a comma"
    ),
    optparse::make_option(c("-e", "--exclude"),
      action = "store",
      default = NA,
      type = "character",
      help = "Path to a list of CpG sites to exclude"
    ),
    optparse::make_option(c("-o", "--output"),
      action = "store",
      default = NA,
      type = "character",
      help = "Path to output directory"
    ),
    optparse::make_option(c("-r", "--ref"),
      action = "store",
      default = NA,
      type = "character",
      help = "Path to reference methylation file"
    ),
    optparse::make_option(c("-c", "--minCov"),
      action = "store",
      default = 10,
      type = "integer",
      help = "Minimum coverage to consider methylation site as present. Default = 10"
    ),
    optparse::make_option(c("-d", "--maxDistance"),
      action = "store",
      default = 1000,
      type = "integer",
      help = "Maximum distance between missing site and each neighbour for the site to be imputed. Default = 1000 bp"
    ),
    optparse::make_option(c("-k", "--collapse"),
      action = "store_false",
      type = "logical",
      default = TRUE,
      help = "Choose whether to merge methylation sites on opposite strands together. Default = True"
    ),
    optparse::make_option(c("-x", "--machineLearning"),
      action = "store_true",
      type = "logical",
      default = FALSE,
      help = "Choose whether to use machine learning for imputation. Default = False (no machine learning)"
    ),
    optparse::make_option(c("-t", "--runTime"),
      action = "store",
      default = 3600,
      type = "integer",
      help = "Time (seconds) to train model. Default = 3600s (2h)"
    ),
    optparse::make_option(c("-m", "--maxModels"),
      action = "store",
      default = 20,
      type = "integer",
      help = "Maximum number of models to train within the time specified under --runTime. Excludes Stacked Ensemble models"
    ),
    optparse::make_option(c("-s", "--streaming"),
      action = "store_true",
      type = "logical",
      default = FALSE,
      help = "Choose if streaming is required (for files that exceed memory). Default = False"
    )
  )

  opt <<- optparse::parse_args(optparse::OptionParser(option_list = option_list))

#   # check mandatory options
  if (is.na(opt$input)) {
    print("ERROR: No input directory given. GIMMEcpg terminating.")
    quit()
  }

  if (is.na(opt$ref)) {
    print("ERROR: No reference file given. GIMMEcpg terminating.")
    quit()
  }

  if (is.na(opt$output)) {
    print("ERROR: No output directory given. GIMMEcpg terminating.")
    quit()
  }


  ##########################################################
  # Read file in as LazyFrame, collapse strands if needed. #
  ##########################################################


  bed_paths <- list.files(path = opt$input, pattern = "*.bed", full.names = TRUE) # list all the paths to a bed files

  # select files
  if (!is.na(opt$name)) {
    names <- gsub(",", "|", opt$name)
    bed_paths <- grep(names, bed_paths, value = TRUE)
  }

  # check files exist
  if (identical(bed_paths, character(0))) {
    print("ERROR: No matching Bed file(s) found. GIMMEcpg terminating.")
    quit()
  }


  print(paste("Merge methylation sites on opposite strands =", opt$collapse))
  print(paste("Coverage cutoff at", opt$minCov))

  lf_list <- lapply(bed_paths, read_files, collapse = opt$collapse)


  ##########################
  # Identify missing sites #
  ##########################

  if (!is.na(opt$exclude)) {
    print("Blacklisted regions will be excluded")
  } else {
    print("No blacklisted regions provided; all autosomal CG sites considered")
  }

  missing <- lapply(lf_list, missing_sites, ref = opt$ref, blacklist = opt$exclude, mincov = opt$minCov, ml = opt$machineLearning)

  print("Identified missing sites")

  ################################
  # Imputation (default is fast) #
  ################################

  if (opt$maxDistance > 0) {
    print(paste0("Imputing methylation for missing sites less than ", opt$maxDistance, " bases from each neighbour"))
  }

  results <- list()

  if (opt$machineLearning == FALSE) {
    print("Fast Imputation mode")
    imputed_lfs <- lapply(missing, fast_impute, dist = opt$maxDistance)
    results <- imputed_lfs
  } else {
    print("machineLearning mode: prepare for H2O AutoML training")
    lead_prediction <- lapply(missing, h2oTraining, maxTime = opt$runTime, maxModels = opt$maxModels, dist = opt$maxDistance, streaming = opt$streaming)
    results <- lead_prediction
    h2o::h2o.shutdown(prompt = FALSE)
  }

  batch_limit <- 10

  if (length(results) <= batch_limit) {
    print("Batch mode OFF")
    if (opt$streaming == TRUE) {
      print("Collecting results in streaming mode")
      dfs <- future_lapply(results, function(result) result$collect(streaming = TRUE)) ## NEED TO TEST future_lapply() against lapply()
      future_lapply(dfs, save_files, outpath = opt$output)
      print("All files saved")
    } else {
      print("Collecting results")
      dfs <- lapply(results, function(result) result$collect()) # future_lapply
      lapply(dfs, save_files, outpath = opt$output) # future_lapply
      print("All files saved")
    }
  } else {
    print(paste0("Batches of ", batch_limit)) ### Needs re-working ###

    # batch <- ihasNext(isplitVector(results, chunkSize = batch_limit))
    # while (ihasNext(batch)) {
    #     if (opt$streaming == TRUE) {
    #         print("Collecting batches of results in streaming mode")
    #         dfs <- future_lapply(batch, function(result) result$collect(streaming = TRUE))
    #         future_lapply(dfs, save_files, outpath = opt$output)
    #         print("All files saved")
    #     } else {
    #         print("Collecting batches of results")
    #         dfs <- future_lapply(batch, function(result) result$collect())
    #         future_lapply(dfs, save_files, outpath = opt$output)
    #         print("All files saved")
    #     }
    # }
  }
}

gimmecpg()