bug <- "Akkermansia_muciniphila"
# claudia
dat <- readRDS("data-raw/omm_claudia_new.rds")
dat <- dat[which(dat$chr == bug),]
dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
dat$sample <-  paste0("claudia ", dat$sample)
df_claudia <- data.frame(id = dat$snp_id, AF = dat$AF, sample = dat$sample)

dat <- readRDS("data-raw/omm_ab.rds")
dat <- dat[which(dat$chr == bug),]
dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
dat$sample <- paste0("AB study ",dat$mouse.id, " d",dat$day, " ",
                     group = dat$group, sample = dat$mouse.group)
df_ab <- data.frame(id = dat$snp_id, AF = dat$AF, sample = dat$sample)

# reseq
#dat <- readRDS("data-raw/reseq.rds")
#dat <- dat[which(dat$chr == bug),]
#dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
#dat$sample <- paste0("reseq study ",dat$mouse.id, " d",dat$day, "", dat$mouse.group)
#df_reseq <- data.frame(id = dat$snp_id, AF = dat$AF, sample = dat$sample)

df <- rbind(df_claudia, df_ab)

# prepare data
omm2 <- tidyr::spread(df, sample, AF)
omm2[is.na(omm2)] <- 0
rownames(omm2) <- omm2$id
omm2$id <- NULL

usethis::use_data(omm2, overwrite = TRUE)