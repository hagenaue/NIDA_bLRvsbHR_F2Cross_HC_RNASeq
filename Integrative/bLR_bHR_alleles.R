# Assigning direction of effect of cis-eQTLs based on the bLR vs. bHR distribution of alleles
# Daniel Munro
# August 2023

# Load F0 genotypes, cluster and label 2 groups, infer group names from correlation to DE FC.
# Then save bLR and bHR allele frequencies to use for orienting eQTL effects.

library(tidyverse)

load_geno <- function(filename) {
    gt <- VariantAnnotation::readGT(filename)
    geno <- apply(gt, 2, function(x) c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2,
                                       "0/0" = 0, "0/1" = 1, "1/0" = 1, "1/1" = 2)[x])
    rownames(geno) <- rownames(gt)
    geno
}

geno <- do.call(
    rbind,
    lapply(c(1:20, "X"),
           function(chrom) load_geno(str_glue("data/genotypes/F0_genotypes/chr{chrom}_flt.vcf.gz")))
)
geno <- geno[apply(geno, 1, sd) > 0, ]
pca <- prcomp(t(geno), scale. = TRUE)$x
bLR <- rownames(pca)[pca[, 'PC1'] < 0]
bHR <- rownames(pca)[pca[, 'PC1'] > 0]
stopifnot("5739-JL-0021" %in% bLR) # Labeling inferred from correlation with DEG FC

df <- bind_rows(tibble(sample_id = bLR, group = "bLR"),
                tibble(sample_id = bHR, group = "bHR"))
write_tsv(df, "data/F0_samples.tsv")

snps <- tibble(variant_id = rownames(geno)) |>
    separate(variant_id, into = c("variant_id", "alleles"), sep = "_") |>
    separate(alleles, into = c("ref", "alt"), sep = "/") |>
    mutate(
        bLR_mean = rowMeans(geno[, bLR]),
        bHR_mean = rowMeans(geno[, bHR]),
        # bLR_allele = if_else(bLR_mean < 1, ref, alt),
        # bHR_allele = if_else(bHR_mean < 1, ref, alt),
        alt_group = case_when(
            bLR_mean > bHR_mean ~ "bLR",
            bLR_mean < bHR_mean ~ "bHR",
            .default = "equal"
        )
    )
write_tsv(snps, "data/genotypes/bLR_bHR_alleles.tsv.gz")
