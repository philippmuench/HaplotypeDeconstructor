# Haplotype reconstruction using non-negative matrix factorization 


## usage

```r
# claudia
dat <- readRDS("data/omm_claudia_new.rds")
dat <- dat[which(dat$chr == "Enterococcus_faecalis"),]
dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
dat$sample <-  paste0("claudia ", dat$sample)
df_claudia <- data.frame(id = dat$snp_id, AF = dat$AF, sample = dat$sample)


dat <- readRDS("data/omm_ab.rds")
dat <- dat[which(dat$chr == "Enterococcus_faecalis"),]
dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
dat$sample <- paste0("AB study ",dat$mouse.id, " d",dat$day, " ", dat$mouse.group)
df_ab <- data.frame(id = dat$snp_id, AF = dat$AF, sample = dat$sample)

# reseq
dat <- readRDS("data/reseq.rds")
dat <- dat[which(dat$chr == "Enterococcus_faecalis"),]
dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
dat$sample <- paste0("reseq study ",dat$mouse.id, " d",dat$day, "", dat$mouse.group)
df_reseq <- data.frame(id = dat$snp_id, AF = dat$AF, sample = dat$sample)

df <- rbind(df_claudia, df_ab, df_reseq )

library(tidyr)
require(reshape2)

# prepare data
data_wide <- spread(df, sample, AF)
data_wide[is.na(data_wide)] <- 0
rownames(data_wide) <- data_wide$id
data_wide$id <- NULL


# how many haplotypes?
gof <- assessNumberHaplotyes(data_wide,2:30, nReplicates = 1)
plotNumberHaplotyes(gof)

# run
res <- findSignatures(data_wide, 10)
plotHaplotypeMap(res)
```
