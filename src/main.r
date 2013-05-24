
"%notin%" <- function(x,y) !(x %in% y)

bring_to_front <- function(df, names) {
    names_df <- rownames(df)
    indeces_df <- which(names_df %in% names)
    df <- rbind(df[-indeces_df,], df[indeces_df,])
    df
}

get_labels <- function(x) {
    if (is.matrix(x) || is.data.frame(x)){
        if (is.null(rownames(x))){
            labels <- letters[1:nrow(x)] # BUG or FEATURE? if no labels provided, and no rownames, label a maximum of 26 rows with letters
        } else {
            labels <- rownames(x)
        }
    } else {
        if (is.null(names(x))){
            labels <- letters[1:length(x)]
        } else {
            labels <- names(x)
        }
    }

    labels
}

transparent <- function(color, alpha=.5) {
    color <- col2rgb(color)
    lc <- as.list(color)
    names(lc) <- rownames(color)
    rgb(green = lc$green, red = lc$red, blue = lc$blue, alpha = alpha * 255, maxColorValue = 255)
}

annotate_points <- function(x, y=NULL,
                            labels=get_labels(x),
                            pos="above",
                            font="bold",
                            cex=1.2,
                            col="red",...){
    pos <- switch(pos, "below"=1, "leftside"=2, "above"=3, "rightside"=4)
    font <- switch(font, "plain"=1, "bold"=2, "italic"=3, "bold_italic"=4, "symbol"=5)

    if (!is.null(labels)){
        text(x,y,labels=labels, pos=pos, font=font, xpd=NA, cex=cex, col=col, ...)
    }
}

add_ma <- function(data_set){
    for(timepoint in colnames(data_set$normalized_counts)[-1]){
        data_set[[timepoint]]$ma <- get_ma(data_set$normalized_counts, c("hpi_0", timepoint))
    }
    data_set
}

# Calculate the log2 fold-change (M) and the log10 average magnitude of that fold-change (A) between two columns of a data set.
get_ma <- function(data_set, timepoints){
    x <- data_set[, timepoints[1]]
    y <- data_set[, timepoints[2]]

    ma <- data.frame(m = log2(y) - log2(x), a = (log10(x) + log10(y))/2)
    rownames(ma) <- rownames(data_set)
    ma
}

plot_ma <- function(ma, highlight_names = NULL, threshold = NULL, names = FALSE, highlight_color = "black", ...){
    plot(ma$a, ma$m, pch = 16, col = "grey70", xlab="log10 normalized counts", ylab="log2 fold changes",...)
    if (!is.null(highlight_names)){
        points(ma[highlight_names, "a"], ma[highlight_names, "m"], col = "black", pch = 16)
        if (names){
            annotate_points(ma[highlight_names, "a"], ma[highlight_names, "m"], highlight_names, col= "black", ...)
        }
    }
}

# plot_ma <- function(probe = NULL, ...){
#     plot(lassa$d3$visual_pvalue, lassa$d3$logFC, col = ifelse(lassa$d3_condition, "black", "grey70"), main = paste0("Lassa\n", probe), ylab = "3 dpe / uninfected (log2)" , xlab = "Significance (p-value)", xaxt = "n", ...)
#     axis(1, at = 0:14, labels = c("1", NA, NA,".001", rep(NA,11)), family = "Rockwell", cex.axis=1.8)
#     segments(x0=3, y0=c(-20,1.5), y1=c(-1.5,20), col = "black")
#     segments(x0=3, y0=c(-1.5,1.5), x1=16, col = "black")
#     if (!is.null(probe)){
#         points(lassa$d3[probe, "visual_pvalue"], lassa$d3[probe, "logFC"], bg = "white", col = palette$red, pch = 24, cex = 2, lwd = 4)
#     }
# }

get_dge <- function(counts){
    dge <- DGEList(counts=counts, group=factor(colnames(counts), levels=colnames(counts)))
    dge <- calcNormFactors(dge)
    dge$normalized_counts <- equalizeLibSizes(dge)$pseudo.counts

    dge
}

to_gene_symbol <- function(ensembl_ids, ensembl_info){
    gene_symbols <- ensembl_info[match(ensembl_ids, ensembl_info$gene_id), "gene_symbol"]
    gene_symbols[is.na(gene_symbols)] <- ""
    gene_symbols[which(gene_symbols == "")] <- ensembl_ids[which(gene_symbols == "")]
    gene_symbols
}

gene_symbols_as_rownames <- function(data_set, ensembl_info_path = "data/hsa_ensembl_69.rds"){
    ensembl_info <- readRDS(ensembl_info_path)
    candidate_rownames <- to_gene_symbol(rownames(data_set), ensembl_info)
    dups <- duplicated(candidate_rownames)
    candidate_rownames[dups] <- paste0(candidate_rownames[dups], "_dup_", 1:sum(dups))
    rownames(data_set) <- candidate_rownames
    data_set
}

filter_low_counts <- function(counts, min = 10) {
    filtered_counts <- counts[rowSums(counts) >= min, ]
    filtered_counts
}

get_counts <- function(path){
    counts <- read.table(path, header=TRUE, check.names=FALSE)
    counts <- counts[rownames(counts) %notin% c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique"),]
    counts
}

get_data_set <- function(path){
    data_set <- list()
    data_set$counts <- get_counts(path)
    data_set$filtered_counts <- filter_low_counts(data_set$counts)

    # use gene names instead of Ensembl gene ids
    data_set$gene_ids <- rownames(data_set$filtered_counts)
    data_set$filtered_counts <- gene_symbols_as_rownames(data_set$filtered_counts)

    # normalize counts using edgeR's Trimmed Mean of M-values (TMM)
    data_set$dge <- get_dge(data_set$filtered_counts)
    data_set$normalized_counts <- data.frame(data_set$dge$normalized_counts)
    data_set$normalized_counts[data_set$normalized_counts < .1] <- 1 # small fix, so the log10(data_set$normalized_counts) doesn't choke with zeros

    data_set$log2_fold_changes <- log2(data_set$normalized_counts/data_set$normalized_counts[,1])

    data_set
}