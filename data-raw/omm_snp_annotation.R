# claudia
dat <- readRDS("data-raw/omm_claudia_new.rds")
dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
df_claudia <- data.frame(id = dat$snp_id, description = dat$feature)


# ab
dat <- readRDS("data-raw/omm_ab.rds")
dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
df_ab <- data.frame(id = dat$snp_id, description = dat$feature)

# reseq
dat <- readRDS("data-raw/reseq.rds")
dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
df_reseq <- data.frame(id = dat$snp_id, description = dat$feature)

df <- rbind(df_claudia, df_ab, df_reseq)

omm_snp_annotation <- df[!duplicated(df), ]

# prepare data
usethis::use_data(omm_snp_annotation, overwrite = TRUE)