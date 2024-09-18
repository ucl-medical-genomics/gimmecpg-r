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


distBins <- function(lf, col) {

    corrMat <- (
        lf$with_columns(
            pl$when(pl$col(col) <= 200)
                $then(pl$col(col))
            $alias("dist_bins")
        )
        $with_columns(
            pl$when(pl$col(col) > 200)
                $then((((pl$col(col)-1) %/% 10) + 1) * 10)
                $otherwise(pl$col("dist_bins"))
            $alias("dist_bins")
        )
        $with_columns(
            pl$when(pl$col(col) > 500)
                $then((((pl$col(col)-1) %/% 100) + 1) * 100)
                $otherwise(pl$col("dist_bins"))
            $alias("dist_bins")
        )
        $with_columns(
            pl$when(pl$col(col) > 2000)
                $then(-1)
                $otherwise(pl$col("dist_bins"))
            $alias("dist_bins")
        )
    )

    return(corrMat)
}


h2oPrep <- function(lf, dist, streaming) {

    known_sites <- (
        lf$filter(pl$col("avg")$is_not_null())
        $select(list("chr", "start", "end", "avg"))
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
        $with_columns(
            ((pl$col("b_meth") * pl$col("f_dist") + pl$col("f_meth") * pl$col("b_dist"))
            / pl$sum_horizontal("f_dist", "b_dist"))$alias("fast_res")
        )
        $with_columns(
            (pl$col("fast_res") - pl$col("avg"))$abs()$alias("error")
        )
        $with_columns(
            (pl$col("error")$quantile(0.99, "nearest"))$alias(limit)
        )
        $filter(pl$col("error") < pl$col("limit"))
        
    )
   
    if (dist > 0) {
        known_sites <- known_sites$filter((pl$col("f_dist") <= dist) & (pl$col("b_dist") <= dist))
    }

    known_sites <- distBins(known_sites, "f_dist")$rename(dist_bins = "f_dist_bins")
    known_sites <- distBins(known_sites, "b_dist")$rename(dist_bins = "b_dist_bins")

    corr <- (known_sites$select(list("avg", "f_meth", "f_dist", "f_dist_bins"))
            $group_by("f_dist_bins")$agg(pl$corr("f_meth","avg", method = "pearson")$alias("corr"))
            $with_columns(
                pl$when(pl$col("f_dist_bins") == -1)$then(0)$otherwise(pl$col("corr"))$alias("corr")
            )
    )$drop_nulls()


    features_lf <- (
        known_sites
        $join(corr, right_on = "f_dist_bins", left_on = "b_dist_bins", how = "left")$rename(corr = "b_corr")
        $join(corr, on = "f_dist_bins", how = "left")$rename(corr = "f_corr")
        $select(list("avg", "b_meth", "f_meth", "b_dist", "f_dist", "b_corr", "f_corr"))
        $with_columns(
            ((pl$col("b_meth") + 1)$log())$alias("lg_b_meth"),
            ((pl$col("f_meth") + 1)$log())$alias("lg_f_meth"),
            ((pl$col("b_dist") + 1)$log())$alias("lg_b_dist"),
            ((pl$col("f_dist") + 1)$log())$alias("lg_f_dist")
        )
    )

    missing_sites <- lf$filter(pl$col("avg")$is_null())

    missing_sites <- distBins(missing_sites, "f_dist")$rename(dist_bins = "f_dist_bins")
    missing_sites <- distBins(missing_sites, "b_dist")$rename(dist_bins = "b_dist_bins")

    to_predict_lf <- (
        missing_sites
        $join(corr, right_on = "f_dist_bins", left_on = "b_dist_bins", how = "left")$rename(corr = "b_corr")
        $join(corr, on = "f_dist_bins", how = "left")$rename(corr = "f_corr")
        $with_columns(
            ((pl$col("b_meth") + 1)$log())$alias("lg_b_meth"),
            ((pl$col("f_meth") + 1)$log())$alias("lg_f_meth"),
            ((pl$col("b_dist") + 1)$log())$alias("lg_b_dist"),
            ((pl$col("f_dist") + 1)$log())$alias("lg_f_dist")
        )
    )

    if (opt$streaming == TRUE) {
        features <- features_lf$collect(streaming=TRUE)
        to_predict <- to_predict_lf$select(list("avg", "lg_b_dist", "lg_f_dist", "lg_b_meth", "lg_f_meth", "b_corr", "f_corr"))$collect(streaming=TRUE)
    } else {
        features <- features_lf$collect()
        to_predict <- to_predict_lf$select(list("avg", "lg_b_dist", "lg_f_dist", "lg_b_meth", "lg_f_meth", "b_corr", "f_corr"))$collect()
    }

   prepped <- list(features, to_predict, to_predict_lf)
   return(prepped)

}


h2oTraining <- function(lf, maxTime, maxModels, dist, streaming) {

    print("Starting H2O AutoML training")

    training <- h2oPrep(lf, dist, streaming)$known
    test <- h2oPrep(lf, dist, streaming)$to_predict
    to_predict_lf <- h2oPrep(lf, dist, streaming)$to_predict_lf

    h2o.init()

    trainingFrame <- as.h2o(as.data.frame(training))
    trainingFrame[c("avg", "lg_b_dist", "lg_f_dist", "lg_b_meth", "lg_f_meth", "b_corr", "f_corr")] <- as.numeric(trainingFrame[
        c("avg", "lg_b_dist", "lg_f_dist", "lg_b_meth", "lg_f_meth", "b_corr", "f_corr")
    ])

    testingFrame <- as.h2o(as.data.frame(test))
    testingFrame[c("avg", "lg_b_dist", "lg_f_dist", "lg_b_meth", "lg_f_meth", "b_corr", "f_corr")] <- as.numeric(testingFrame[
        c("avg", "lg_b_dist", "lg_f_dist", "lg_b_meth", "lg_f_meth", "b_corr", "f_corr")
    ])

    y <- "avg"
    x <- c("avg", "lg_b_dist", "lg_f_dist", "lg_b_meth", "lg_f_meth", "b_corr", "f_corr")

    aml <- h2o.automl(x = x, y = y, max_runtime_secs=maxTime, max_models=maxModels, nfolds = 5, stopping_rounds = 3, sort_metric = "deviance", seed=1, training_frame = trainingFrame)

    lb <- aml@leaderboard

    prediction <- h2o.predict(aml@leader, testingFrame)
    prediction_lf <- pl$as.data.frame(prediction)$lazy()

    imputed_lf <- pl$concat(to_predict_lf, prediction_lf, how="horizontal")

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

    print(lb, n = nrow(lb)) 

    h2o.removeAll()

    return(res)

}
