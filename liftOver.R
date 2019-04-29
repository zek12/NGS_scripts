# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("vcfR", version = "3.8")
library(rtracklayer)
library(dplyr)
library(VariantAnnotation)
library(vcfR)
# detach("package:vcfR", unload=TRUE)



# example
# biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# 
# chain <- import.chain("hg19ToHg18.over.chain")
# tx_hg19 <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
# tx_hg18 <- liftOver(tx_hg19, chain)
# end of example

input_vcf <- "2.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.vcf"
output_vcf <- "2.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.vcf.gz"


vcf <- readVcf(input_vcf, genome = "hg38")
granges <- rowRanges(vcf)
chain <- import.chain("various_scripts/liftOver/hg38ToHg19.over.chain")


lifted <- liftOver(granges, chain)
lifted_df <- lifted %>% as.data.frame()
duplicated(lifted_df$group_name) %>% table() # all different and matches the number of variants in VCF (5795)

lifted_df$ALT <- lapply(lifted_df$ALT, function(x) {
    return(x[[1]] %>% as.character())
}) %>% unlist()

# if it's a deletion or an insertion, POS = start + 1
# else, POS = start
# lifted_df %>% rowwise() %>%  mutate(POS = ifelse(nchar(REF) == 1 & nchar(ALT) == 1, start, start + 1 )) %>% ungroup() %>% as.data.frame() -> lifted_df
lifted_df %>% rowwise() %>%  mutate(POS = start ) %>% ungroup() %>% as.data.frame() -> lifted_df

# geno(vcf)
# geno(header(vcf))
# samples(header(vcf))
# update VCF coordinates with GRCh37
# rowRanges(vcf)[1] %>% ranges() %>% start() <- 8858646 # change POS, original was 8858645
# ref(vcf) <- lifted_df$REF
# qual(vcf)[1:5]


# using vcfR
vcf2 <- read.vcfR(input_vcf, verbose = F)
vcf2@fix[,'POS'] <- as.character(lifted_df$POS)
write.vcf(vcf2, output_vcf)
# then run: gzip -d 2.biallelic.no_star.GRCh38_pos.rsIDannotated.contigs_corrected.repeat_flagged.clinvar.lifted_rtracklayer.vcf.gz

# table(vcf2@fix[,'CHROM'] == lifted_df$seqnames)
# table(vcf2@fix[,'REF'] == lifted_df$REF)
# table(vcf2@fix[,'ALT'] == lifted_df$ALT)

lifted_df %>% filter(group_name == "chr2_73448097_TCTC_T") # this variant couldn't be lifted over using UCSC online tool
lifted_df %>% filter(group >= 1358 & group <= 1360)



# The following lines are to compare rtracklayer vs Ensembl and rtracklayer vs UCSC online tools

## Compare vs Ensembl
# temp <- read.table("temp.txt", h = F)
# colnames(temp) <- c("CHROM", "POS", "ID", "REF", "ALT")
# write.table(temp, "temp_with_header.txt", quote = F)

temp <- read.table("temp_with_header.txt", h = T)
lifted_df$ID <- info(vcf)$GRCh38_pos
i <- intersect(temp$ID, lifted_df$ID)
comparison <- data.frame(ID = i, r_lifted_pos = lifted_df %>% filter(ID %in% i) %>% .$POS, ensembl_lifted_pos = temp %>% filter(ID %in% i) %>% .$POS)
all_equal(comparison$r_lifted_pos, comparison$ensembl_lifted_pos) # all matches with Ensembl

temp %>% filter(!ID %in% i)
lifted_df %>% filter(!ID %in% i) # Ensembl couldn't lift variant chr2_85554195_G_T

## Compare vs UCSC
temp <- read.table("2.biallelic.no_star.with_ID.lifted_UCSC.bed")
colnames(temp) <- c("CHROM", "START", "END", "ID", "QUAL")
i <- intersect(temp$ID, lifted_df$ID)
comparison <- data.frame(ID = i, r_lifted_pos = lifted_df %>% filter(ID %in% i) %>% .$POS, ucsc_lifted_pos = temp %>% filter(ID %in% i) %>% .$START)
all_equal(comparison$r_lifted_pos, comparison$ucsc_lifted_pos) # all matches with UCSC

temp %>% filter(!ID %in% i)
lifted_df %>% filter(!ID %in% i) # UCSC couldnt lift variant chr2_73448097_TCTC_T








