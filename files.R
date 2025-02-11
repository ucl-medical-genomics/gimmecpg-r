library(polars)


collapse_strands <- function(bed) {

    pos <- bed$filter(pl$col("strand") == "+")
    neg <- bed$filter(pl$col("strand") == "-")$with_columns(
            (pl$col("start"))$alias("cStart")
    ) # add column for start site on complementary strand

    joint <- pos$join(
            neg, left_on = c("chr", "end"), right_on = c("chr", "cStart"), how = "full", coalesce = TRUE
    )$with_columns(pl$concat_str(pl$col("strand"), pl$col("strand_right"), separator = "/", ignore_nulls = TRUE))

    merged <- (
        joint$with_columns(
            pl$min_horizontal("start", "start_right")$alias("start"),
            pl$col(c("coverage_right", "coverage")) 
            $fill_null(0)
            $cast(pl$UInt64)
        )$with_columns(
            (pl$col("coverage") + pl$col("coverage_right"))$alias("total_coverage"),
            pl$when(pl$col("strand") == "-")$then(pl$col("start") - 1)$otherwise(pl$col("start"))$alias("start")
        )$filter(pl$col("total_coverage") > 0) # remove sites with no coverage
        $with_columns(
            (
                (
                    pl$col("coverage") * pl$col("percent_methylated")
                    + pl$col("coverage_right") * pl$col("percent_methylated_right")
                )
                / pl$col("total_coverage")
            )$alias("avg")
        )
    )
    return(merged)
}


read_files <- function(file, mincov, collapse) {

    name <- basename(file)
    name <- gsub(pattern = ".bed", "", name)
    print(paste0("Scanning ", name))
    
    bed <- (
        pl$scan_csv(
            file, 
            separator = "\t", 
            skip_rows = 1, 
            has_header = FALSE, 
            dtypes = list(
                'column_1' = pl$String, 
                'column_2' = pl$UInt64, 
                'column_3'= pl$UInt64, 
                'column_6'= pl$String, 
                'column_10'= pl$UInt64, 
                'column_11'= pl$UInt64
            )
        )$select(
            list("column_1", "column_2", "column_3", "column_6", "column_10", "column_11") # only select relevant columns
        )$rename(
            column_1 = "chr", 
            column_2 = "start", 
            column_3 = "end", 
            column_6 = "strand", 
            column_10 = "coverage", 
            column_11 = "percent_methylated" # rename to something that makes more sense
        )$with_columns(
            pl$col("chr")$str$replace("(?i)Chr", "") # remove "chr" from Chr column to match reference
        )
    )

    if (collapse == TRUE) {
        data <- collapse_strands(bed)
    } else {
        data <- bed$with_columns(pl$col("percent_methylated")$alias("avg"), 
                pl$col("coverage")$alias("total_coverage")
        ) 
    }

    data_cov_filt <- (
        data$with_columns(pl$lit(name)$alias("sample"))
        $select(list("chr", "start", "strand", "avg", "sample", "total_coverage"))
        $cast(list(chr = pl$String, start = pl$UInt64, strand = pl$String, avg = pl$Float64, sample = pl$String, total_coverage = pl$UInt64))
    )
    
    return(data_cov_filt)
}


save_files <- function(df, outpath) {
    
    filename <- (
        df$unique(subset="sample", keep="any")
        $select(pl$col("sample")$filter(pl$col("sample") != "imputed")$first())
        $item()
    )
    outfile <- file.path(outpath, paste0("imputed_", filename, ".bed"))
    print(paste0("Saving ", filename))
    df$write_csv(outfile, separator = "\t")

}
