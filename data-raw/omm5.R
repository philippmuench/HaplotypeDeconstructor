bug <- "Muribaculum_intestinale"
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


df <- rbind(df_claudia, df_ab)

# prepare data
omm5 <- tidyr::spread(df, sample, AF)
omm5[is.na(omm5)] <- 0
rownames(omm5) <- omm5$id
omm5$id <- NULL

usethis::use_data(omm5, overwrite = TRUE)