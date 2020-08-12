
# claudia
dat <- readRDS("data-raw/omm_claudia_new.rds")
dat <- dat[which(dat$chr == "Enterococcus_faecalis"),]
dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
dat$sample <-  paste0("claudia ", dat$sample)
df_claudia <- data.frame(id = dat$snp_id, AF = dat$AF, sample = dat$sample)

dat <- readRDS("data-raw/omm_ab.rds")
dat <- dat[which(dat$chr == "Enterococcus_faecalis"),]
dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
dat$sample <- paste0("AB study ",dat$mouse.id, " d",dat$day, " ", dat$mouse.group)
df_ab <- data.frame(id = dat$snp_id, AF = dat$AF, sample = dat$sample)

# reseq
dat <- readRDS("data-raw/reseq.rds")
dat <- dat[which(dat$chr == "Enterococcus_faecalis"),]
dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
dat$sample <- paste0("reseq study ",dat$mouse.id, " d",dat$day, "", dat$mouse.group)
df_reseq <- data.frame(id = dat$snp_id, AF = dat$AF, sample = dat$sample)

df <- rbind(df_claudia, df_ab, df_reseq)

# prepare data
omm <- tidyr::spread(df, sample, AF)
omm[is.na(omm)] <- 0
rownames(omm) <- omm$id
omm$id <- NULL

usethis::use_data(omm, overwrite = TRUE)