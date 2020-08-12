dat <- readRDS("data/omm_claudia_new.rds")
dat <- dat[which(dat$chr == "Akkermansia_muciniphila"),]
# long to wide

dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
df <- data.frame(id = dat$snp_id, AF = dat$AF, sample = dat$sample)

library(tidyr)
require(reshape2)

# prepare data
data_wide <- spread(df, sample, AF)
data_wide[is.na(data_wide)] <- 0
rownames(data_wide) <- data_wide$id
data_wide$id <- NULL


# how many haplotypes?
gof <- assessNumberHaplotyes(data_wide,2:30, nReplicates = 4)
plotNumberHaplotyes(gof)


# run
res <- findSignatures(data_wide, 9)
plotHaplotypeMap(res)

