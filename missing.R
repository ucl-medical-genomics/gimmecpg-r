library(polars)

align_to_reference <- function(bed, ref, blacklist) {

    if (!is.na(blacklist)) {

        blacklist <- (
            pl$scan_parquet(blacklist, parallel="row_groups")
            $with_columns(
                pl$col("chr")$cast(pl$dtypes$String),
                pl$col("start")$cast(pl$dtypes$UInt64)
            )
            $select(list("chr", "start"))
        )

        ref <- (
            pl$scan_parquet(ref, parallel="row_groups")
            $with_columns(
                pl$col("chr")$cast(pl$dtypes$String),
                pl$col("start")$cast(pl$dtypes$UInt64),
                pl$col("end")$cast(pl$dtypes$UInt64)
            )
            $select(list("chr", "start", "end"))
        )$filter(!pl$col("chr")$is_in(c("Y", "X")))$join(blacklist, how = "anti", on = c("chr", "start"))
        
    } else {

        ref <- (
            pl$scan_parquet(ref, parallel="row_groups")
            $with_columns(
                pl$col("chr")$cast(pl$dtypes$String),
                pl$col("start")$cast(pl$dtypes$UInt64),
                pl$col("end")$cast(pl$dtypes$UInt64)
            )
            $select(list("chr", "start", "end"))
        )$filter(!pl$col("chr")$is_in(c("Y", "X")))
        
    }

    aligned <- ref$join(bed, on=c("chr", "start"), how="left")

    return(aligned)

}


missing_sites <- function(bed, ref, blacklist, mincov, ml) {

    aligned <- align_to_reference(bed, ref, blacklist)
    
    if (opt$machineLearning == FALSE) {
        filtered <- aligned$with_columns(
            avg = pl$when(pl$col("total_coverage") < mincov)$then(NA)$when(pl$col("over") >= 0)$then(NA)$otherwise("avg") # 
            # "avg" is the column where we check if imputation needed or not; "methylation" column is final methylation value
        )
        
        missing_sites_defined <- (
            filtered$with_columns(
                pl$when(pl$col("avg")$is_not_null())$then(pl$col("start"))$alias("b_start"),
                pl$when(pl$col("avg")$is_not_null())$then(pl$col("start"))$alias("f_start"),
                pl$when(pl$col("avg")$is_not_null())$then(pl$col("avg"))$alias("b_meth"),
                pl$when(pl$col("avg")$is_not_null())$then(pl$col("avg"))$alias("f_meth")
            )
            $with_columns(
                pl$col(c("f_start", "f_meth"))$backward_fill()$over("chr"),
                pl$col(c("b_start", "b_meth"))$forward_fill()$over("chr")
            )
            $with_columns(
                (pl$col("start") - pl$col("b_start"))$alias("b_dist"),
                (pl$col("f_start") - pl$col("start"))$alias("f_dist")
            )
        )

    } else { # machine learning mode can handle low-quality sites better
        neighbours_added <- (
            aligned$with_columns(
                pl$when(pl$col("avg")$is_not_null())$then(pl$col("start"))$alias("b_start"),
                pl$when(pl$col("avg")$is_not_null())$then(pl$col("start"))$alias("f_start"),
                pl$when(pl$col("avg")$is_not_null())$then(pl$col("avg"))$alias("b_meth"),
                pl$when(pl$col("avg")$is_not_null())$then(pl$col("avg"))$alias("f_meth")
            )
            $with_columns(
                pl$col(c("f_start", "f_meth"))$backward_fill()$over("chr"),
                pl$col(c("b_start", "b_meth"))$forward_fill()$over("chr")
            )
            $with_columns(
                (pl$col("start") - pl$col("b_start"))$alias("b_dist"),
                (pl$col("f_start") - pl$col("start"))$alias("f_dist")
            )
        )

        missing_sites_defined <- neighbours_added$with_columns(
            methylation = pl$col("avg"),
            avg = pl$when(pl$col("total_coverage") < mincov)$then(NA)$when(pl$col("over") >= 0)$then(NA)$otherwise("avg") # 
            # "avg" is the column where we check if imputation needed or not; "methylation" column is final methylation value
        )

    }

    
    return(missing_sites_defined)

}