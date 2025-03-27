
#' @param lf a list of lazyframes
#' @param dist numeric value for minimum distance to next neighbouring CpG site
#' @returns a list of lazyframes
#' @examples fast_impute(lazyframes, 1000)
#' @keywords internal
#' @noRd 
#' @import polars
fast_impute <- function(lf, dist) {
  if (dist > 0) {
    lf <- lf$filter((pl$col("f_dist") <= dist) & (pl$col("b_dist") <= dist))
  }

  imputed <- (lf$with_columns(
    sample = pl$when(pl$col("avg")$is_null())$then(pl$lit("imputed"))$otherwise("sample")
  )
  $with_columns(
    pl$col("avg")$fill_null(
      (pl$col("b_meth") * pl$col("f_dist") + pl$col("f_meth") * pl$col("b_dist"))
      / pl$sum_horizontal("f_dist", "b_dist")
    )
  )
  $with_columns(
    methylation = pl$col("avg")
  )
  )

  imputed <- imputed$select(list("chr", "start", "end", "strand", "sample", "methylation", "total_coverage"))

  return(imputed)
}


#' @param lf a list of lazyframes
#' @param dist numeric value for minimum distance to next neighbouring CpG site
#' @param streaming boolean; use streaming or not
#' @returns a list of lazyframes
#' @examples h2oPrep(lazyframes, 1000, FALSE)
#' @keywords internal
#' @noRd 
#' @import polars
h2oPrep <- function(lf, dist, streaming) {
  features_lf <- (
    lf$filter((pl$col("avg")$is_not_null())) # any CpG site that is present in the original data, regardless of coverage, is considered "known". "Methylation" column will only be null if CpG site is in reference but not in original data
    $with_columns(
      pl$col("start")$shift(-1)$over("chr")$alias("f_start"),
      pl$col("start")$shift()$over("chr")$alias("b_start"),
      pl$col("methylation")$shift(-1)$over("chr")$alias("f_meth"),
      pl$col("methylation")$shift()$over("chr")$alias("b_meth")
    )
    $with_columns(
      (pl$col("start") - pl$col("b_start"))$alias("b_dist"),
      (pl$col("f_start") - pl$col("start"))$alias("f_dist")
    )

  )

  to_predict_lf <- lf$filter(pl$col("avg")$is_null()) # "avg" is null when either site from reference not present in sample data, or if coverage < 10

  if (dist > 0) {
    to_predict_lf <- to_predict_lf$filter((pl$col("f_dist") <= dist) & (pl$col("b_dist") <= dist))
    features_lf <- features_lf$filter((pl$col("f_dist") <= dist) & (pl$col("b_dist") <= dist))
  }


  if (opt$streaming == TRUE) {
    features <- features_lf$collect(streaming = TRUE)
    to_predict <- to_predict_lf$select(list("methylation", "b_dist", "f_dist", "b_meth", "f_meth"))$collect(streaming = TRUE)
  } else {
    features <- features_lf$select(list("methylation", "b_dist", "f_dist", "b_meth", "f_meth"))$collect()
    to_predict <- to_predict_lf$select(list("methylation", "b_dist", "f_dist", "b_meth", "f_meth"))$collect()
  }

  prepped <- list(features$to_data_frame(), to_predict$to_data_frame(), to_predict_lf)
  return(prepped)
}


#' @param lf a list of lazyframes
#' @param maxTime numeric value (seconds) for maximum model training time
#' @param maxModels numeric value for maximum models built 
#' @param dist numeric value for minimum distance to next neighbouring CpG site
#' @param streaming boolean; use streaming or not
#' @returns a list of lazyframes
#' @examples h2oTraining(lazyframes, 18000, 20, 1000, FALSE)
#' @keywords internal
#' @noRd 
#' @import polars
#' @import h2o
h2oTraining <- function(lf, maxTime, maxModels, dist, streaming) {
  print("Starting H2O AutoML training")

  training <- h2oPrep(lf, dist, streaming)[[1]] # features
  test <- h2oPrep(lf, dist, streaming)[[2]] # to_predict
  to_predict_lf <- h2oPrep(lf, dist, streaming)[[3]] # to_predict_lf

  h2o::h2o.init(port = 54321, nthreads = 24, max_mem_size = "100G")

  trainingFrame <- h2o::as.h2o(training)

  testingFrame <- h2o::as.h2o(test)

  y <- "methylation"
  x <- c("b_dist", "f_dist", "b_meth", "f_meth")
  aml <- h2o::h2o.automl(x = x, y = y, max_runtime_secs = maxTime, max_models = maxModels, seed = 1234, training_frame = trainingFrame, sort_metric = "MAE")

  lb <- h2o::h2o.get_leaderboard(object = aml, extra_columns = "ALL")
  print(lb, n = nrow(lb))

  prediction <- h2o::h2o.predict(aml@leader, testingFrame)
  prediction_df <- as.data.frame(prediction) # results

  to_predict_df <- as.data.frame(to_predict_lf$select(list("chr", "start"))$collect()) # the coords and chr information of predicted sites

  imputed_lf <- polars::as_polars_lf(cbind(to_predict_df, prediction_df))$cast(list(start = pl$UInt64))$with_columns(imputed = pl$lit("yes"))

  lf <- lf$drop("methylation") # remove the old methylation column

  res <- (
    lf$join(imputed_lf, on = c("chr", "start"), how = "full", coalesce = TRUE)
    $with_columns(
      pl$col("avg")$fill_null(pl$col("predict")),
      sample = pl$when(pl$col("imputed") == "yes")$then(pl$lit("imputed"))$otherwise("sample")
    )
    $with_columns(
      methylation = pl$col("avg")$clip(0, 100)
    )
  )

  res <- res$select(list("chr", "start", "end", "strand", "sample", "methylation", "total_coverage"))

  h2o::h2o.removeAll()

  return(res)
}
