library(h2o)
library(polars)


fast_impute <- function(lf, dist) {

    if (dist > 0) {
        lf <- lf$filter((pl$col("f_dist") <= dist) & (pl$col("b_dist") <= dist))
    }

    imputed <- lf$with_columns(
        pl$col("avg")$fill_null(
            (pl$col("b_meth") * pl$col("f_dist") + pl$col("f_meth") * pl$col("b_dist"))
            / pl$sum_horizontal("f_dist", "b_dist")
        ),
        pl$col("sample")$fill_null(pl$lit("imputed"))
    )

    imputed <- imputed$select(list("chr", "start", "end", "strand", "sample", "avg"))

    return(imputed)

}


# distBins <- function(lf, col) {

#     corrMat <- (
#         lf$with_columns(
#             pl$when(pl$col(col) <= 200)
#                 $then(pl$col(col))
#             $alias("dist_bins")
#         )
#         $with_columns(
#             pl$when(pl$col(col) > 200)
#                 $then((((pl$col(col)-1) %/% 10) + 1) * 10)
#                 $otherwise(pl$col("dist_bins"))
#             $alias("dist_bins")
#         )
#         $with_columns(
#             pl$when(pl$col(col) > 500)
#                 $then((((pl$col(col)-1) %/% 100) + 1) * 100)
#                 $otherwise(pl$col("dist_bins"))
#             $alias("dist_bins")
#         )
#         $with_columns(
#             pl$when(pl$col(col) > 2000)
#                 $then(-1)
#                 $otherwise(pl$col("dist_bins"))
#             $alias("dist_bins")
#         )
#     )

#     return(corrMat)
# }


h2oPrep <- function(lf, dist, streaming) {

    known_sites <- (
        lf$filter(pl$col("avg")$is_not_null())
        $select(list("chr", "start", "avg")) #"end", 
        $with_columns(
            pl$col("start")$shift(-1)$over("chr")$alias("f_start"),
            pl$col("start")$shift()$over("chr")$alias("b_start"),
            pl$col("avg")$shift(-1)$over("chr")$alias("f_meth"),
            pl$col("avg")$shift()$over("chr")$alias("b_meth")
        )
        $with_columns(
            (pl$col("start") - pl$col("b_start"))$alias("b_dist"), 
            (pl$col("f_start") - pl$col("start"))$alias("f_dist")
        )
        # $with_columns(
        #     ((pl$col("b_meth") * pl$col("f_dist") + pl$col("f_meth") * pl$col("b_dist"))
        #     / pl$sum_horizontal("f_dist", "b_dist"))$alias("fast_res")
        # )
        # $with_columns(
        #     (pl$col("fast_res") - pl$col("avg"))$abs()$alias("error")
        # )
        # $with_columns(
        #     (pl$col("error")$quantile(0.99, "nearest"))$alias(limit)
        # )
        # $filter(pl$col("error") < pl$col("limit"))
        
    )
   
    if (dist > 0) {
        known_sites <- known_sites$filter((pl$col("f_dist") <= dist) & (pl$col("b_dist") <= dist))
    }

    # known_sites <- distBins(known_sites, "f_dist")$rename(dist_bins = "f_dist_bins")
    # known_sites <- distBins(known_sites, "b_dist")$rename(dist_bins = "b_dist_bins")

    # corr <- (known_sites$select(list("avg", "f_meth", "f_dist", "f_dist_bins"))
    #         $group_by("f_dist_bins")$agg(pl$corr("f_meth","avg", method = "pearson")$alias("corr"))
    #         $with_columns(
    #             pl$when(pl$col("f_dist_bins") == -1)$then(0)$otherwise(pl$col("corr"))$alias("corr")
    #         )
    # )$drop_nulls()


    features_lf <- (
        known_sites
        # $join(corr, right_on = "f_dist_bins", left_on = "b_dist_bins", how = "left")$rename(corr = "b_corr")
        # $join(corr, on = "f_dist_bins", how = "left")$rename(corr = "f_corr")
        # $select(list("avg", "b_meth", "f_meth", "b_dist", "f_dist", "b_corr", "f_corr"))
        # $with_columns(
        #     ((pl$col("b_meth") + 1)$log())$alias("lg_b_meth"),
        #     ((pl$col("f_meth") + 1)$log())$alias("lg_f_meth"),
        #     ((pl$col("b_dist") + 1)$log())$alias("lg_b_dist"),
        #     ((pl$col("f_dist") + 1)$log())$alias("lg_f_dist")
        # )
    )

    missing_sites <- lf$filter(pl$col("avg")$is_null())
    
    if (dist > 0) {
        missing_sites <- missing_sites$filter((pl$col("f_dist") <= dist) & (pl$col("b_dist") <= dist))
    }
    # missing_sites <- distBins(missing_sites, "f_dist")$rename(dist_bins = "f_dist_bins")
    # missing_sites <- distBins(missing_sites, "b_dist")$rename(dist_bins = "b_dist_bins")

    to_predict_lf <- (
        missing_sites
        # $join(corr, right_on = "f_dist_bins", left_on = "b_dist_bins", how = "left")$rename(corr = "b_corr")
        # $join(corr, on = "f_dist_bins", how = "left")$rename(corr = "f_corr")
        # $with_columns(
        #     ((pl$col("b_meth") + 1)$log())$alias("lg_b_meth"),
        #     ((pl$col("f_meth") + 1)$log())$alias("lg_f_meth"),
        #     ((pl$col("b_dist") + 1)$log())$alias("lg_b_dist"),
        #     ((pl$col("f_dist") + 1)$log())$alias("lg_f_dist")
        # )
    )
    # write.csv(as.data.frame(known_sites$fetch(10000)), "/home/nchai/projects/gimmecpg-r/train2.csv") 
    # quit()

    if (opt$streaming == TRUE) {
        features <- features_lf$collect(streaming=TRUE)
        to_predict <- to_predict_lf$select(list("avg", "b_dist", "f_dist", "b_meth", "f_meth"))$collect(streaming=TRUE) #, "b_corr", "f_corr"
    } else {
        features <- features_lf$select(list("avg", "b_dist", "f_dist", "b_meth", "f_meth"))$collect() #$collect()
        to_predict <- to_predict_lf$select(list("avg", "b_dist", "f_dist", "b_meth", "f_meth"))$collect() #, "b_corr", "f_corr"
    }

   prepped <- list(features$to_data_frame(), to_predict$to_data_frame(), to_predict_lf)
   return(prepped)

}


h2oTraining <- function(lf, maxTime, maxModels, dist, streaming) {

    print("Starting H2O AutoML training")
    
    training <- h2oPrep(lf, dist, streaming)[[1]] # features
    test <- h2oPrep(lf, dist, streaming)[[2]] # to_predict
    to_predict_lf <- h2oPrep(lf, dist, streaming)[[3]] # to_predict_lf

    # write.csv(test, "/home/nchai/NiuzhengChai/test.csv")
    # quit()

    h2o.init(port = 54321, nthreads = -1, max_mem_size = "100G")

    cols <- c("avg", "b_dist", "f_dist", "b_meth", "f_meth")

    trainingFrame <- as.h2o(training)

    # splits <- h2o.splitFrame(data = fullFrame, ratios = 0.80)

    # # Create a training set from the 1st dataset in the split
    # trainingFrame <- splits[[1]]

    # # Create a leaderboard set from the 2nd dataset in the split
    # lbFrame <- splits[[2]]

    # trainingFrame[c("avg", "lg_b_dist", "lg_f_dist", "lg_b_meth", "lg_f_meth", "b_corr", "f_corr")] <- sapply(trainingFrame[
    #     c("avg", "lg_b_dist", "lg_f_dist", "lg_b_meth", "lg_f_meth", "b_corr", "f_corr")
    # ], as.numeric)

    testingFrame <- as.h2o(test)

    # for (col in 1:length(cols)) {
    #     testingFrame[cols[col]] <- h2o::as.numeric(testingFrame[cols[col]])
    #     trainingFrame[cols[col]] <- h2o::as.numeric(trainingFrame[cols[col]])
    # }

    # print(h2o.isnumeric(testingFrame["avg"]))
    # testingFrame <- apply(testingFrame, 2, h2o::as.numeric)
    # testingFrame[c("avg", "b_dist", "f_dist", "b_meth", "f_meth")] <- h2o::as.numeric(testingFrame[c("avg", "b_dist", "f_dist", "b_meth", "f_meth")])
    # print(head(testingFrame))
    print(str(trainingFrame))
    print(str(testingFrame))
    # quit()
    # h2o.shutdown(prompt = FALSE)
    # quit()
    y <- "avg"
    x <- c("b_dist", "f_dist", "b_meth", "f_meth") #, "b_corr", "f_corr"

    aml <- h2o.automl(x = x, y = y, max_runtime_secs=maxTime, max_models=maxModels, seed=1, training_frame = trainingFrame, nfolds = 5, sort_metric = "MAE", verbosity = "debug", distribution = "gamma") # nfolds = 5, ,  , leaderboard_frame = lbFrame, stopping_metric = "MAE", monotone_constraints = list(f_meth = 1) , stopping_tolerance = 0.01

    lb <- h2o.get_leaderboard(object = aml, extra_columns = 'ALL')
    print(lb, n = nrow(lb)) 

    prediction <- h2o.predict(aml@leader, testingFrame)
    prediction_df <- as.data.frame(prediction) # results
    to_predict_df <- as.data.frame(to_predict_lf$select(list("chr", "start"))$collect()) # the coords and other site information
    imputed_lf <- as_polars_lf(cbind(to_predict_df, prediction_df))$cast(list(start = pl$UInt64))
    # print(head(prediction_df))
    # print(head(to_predict_df))
    # print(imputed_lf$fetch(10))

    # print(class(prediction))
    # print(class(prediction_df))
    # print(head(prediction_df))
    # h2o.shutdown(prompt = FALSE)
    # quit()
    # prediction_lf <- as_polars_df(prediction)$lazy()
    # prediction_lf <- as.vector(prediction)
    # print(to_predict_lf)
    # imputed_lf <- to_predict_lf$join(prediction_lf, how="full", left_on = "avg", right_on = "predict", join_nulls = TRUE)
    # imputed_lf <- to_predict_lf$with_columns(predictions = prediction_lf$select("predict"))
    # print(head(a))
    # h2o.shutdown(prompt = FALSE)
    # quit()
    res <- (
        lf$join(imputed_lf, on = c("chr", "start"), how="full", coalesce=TRUE)
        $with_columns(pl$col("avg")$fill_null(pl$col("predict")), pl$col("sample")$fill_null(pl$lit("imputed")))
        $with_columns(
            avg=pl$col("avg")$clip(0, 100)
        )
    )
    
    if (dist > 0) {
        res <- res$filter((pl$col("f_dist") <= dist) & (pl$col("b_dist") <= dist))
    }

    res <- res$select(list("chr", "start", "end", "strand", "sample", "avg"))

    # print(res$fetch(10000))
    # quit()

    h2o.removeAll()

    return(res)

}
