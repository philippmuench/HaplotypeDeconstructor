bug <- "Akkermansia_muciniphila"
# claudia
dat <- readRDS("data-raw/omm_claudia_new.rds")
dat <- dat[which(dat$chr == bug),]
dat$sample <-  paste0("claudia ", dat$sample)
dat$group <- dat$mouse.group
df_claudia <- data.frame(sample = dat$sample,
                         group = dat$group, 
                         day = NA)

dat <- readRDS("data-raw/omm_ab.rds")
dat <- dat[which(dat$chr == bug),]

dat$sample <- paste0("AB study ",dat$mouse.id, " d",dat$day, " ",
                     group = dat$group, sample = dat$mouse.group)

df_ab <- data.frame(sample = dat$sample,
                    group = dat$mouse.group,
                    day = dat$day)

# reseq
#dat <- readRDS("data-raw/reseq.rds")
#dat <- dat[which(dat$chr == bug),]
#dat$snp_id <- paste0(dat$alteration, " ",dat$POS)
#dat$sample <- paste0("reseq study ",dat$mouse.id, " d",dat$day, "", dat$mouse.group)
#df_reseq <- data.frame(id = dat$snp_id, AF = dat$AF, sample = dat$sample)

omm2_metadata <- rbind(df_claudia, df_ab)


usethis::use_data(omm2_metadata, overwrite = TRUE)