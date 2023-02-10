library(BiocParallel)
# library(FRASER)
library(data.table)
library(magrittr)
library(stringr)
register(MulticoreParam(snakemake@threads))

resultsByGenes <- function (res, geneColumn = "hgncSymbol", method = "BY") {
    res <- res[order(res$pValue)]
    if (is(res, "GRanges")) {
        ans <- as.data.table(mcols(res)[, c(geneColumn, "pValue", 
            "sampleID")])
        colnames(ans) <- c("features", "pval", "sampleID")
    }
    else {
        ans <- featureNames <- res[, .(features = get(geneColumn), 
            pval = pValue, sampleID = sampleID)]
    }
    naIdx <- ans[, is.na(features)]
    ansNoNA <- ans[!is.na(features)]
    ansNoNA[, `:=`(pByFeature, min(p.adjust(pval, method = "holm"))), 
        by = "sampleID,features"]
    dupIdx <- duplicated(ansNoNA[, .(features, sampleID)])
    ansGenes <- ansNoNA[!dupIdx]
    ansGenes[, `:=`(fdrByFeature, p.adjust(pByFeature, method = method)), 
        by = "sampleID"]
    finalAns <- res[!naIdx][!dupIdx]
    finalAns$pValueGene <- ansGenes$pByFeature
    finalAns$padjustGene <- ansGenes$fdrByFeature
    finalAns
}


delta_psi_cutoff <- as.double(snakemake@params[['delta_psi_cutoff']])
outlier_type <- unlist(strsplit(snakemake@params[['outlier_type']], "__"))

df <- fread(snakemake@input[['results']])

rj_dt <- df[
  df$type %in% outlier_type
  & abs(df$deltaPsi) >= delta_psi_cutoff
]

fwrite(rj_dt, snakemake@output[['junction_level']], sep = '\t', quote = F)

res_genes_dt <- resultsByGenes(rj_dt) %>% as.data.table()
fwrite(res_genes_dt, snakemake@output[['gene_level']], sep = '\t', quote = F)
