# claudia
dat <- readRDS("data-raw/omm_claudia_new.rds")
dat$sample <-  paste0("claudia ", dat$sample)
dat$group <- dat$mouse.group
df_claudia <- data.frame(sample = dat$sample,
                         group = dat$group, 
                         day = NA,
                         study = "Claudia",
                         mouse = dat$mouse.id)

dat <- readRDS("data-raw/omm_ab.rds")

dat$sample <- paste0("AB study ",dat$mouse.id, " d",dat$day, " ",
                     group = dat$group, sample = dat$mouse.group)

df_ab <- data.frame(sample = dat$sample,
                    group = dat$mouse.group,
                    day = dat$day,
                    study = "AB",
                    mouse = dat$mouse.id)
# reseq
dat <- readRDS("data-raw/reseq.rds")
dat$sample <- paste0("reseq study ",dat$mouse.id, " d",dat$day, "", dat$mouse.group)

df_reseq <- data.frame(sample = dat$sample, group = dat$mouse.group,
                       day = NA,
                       study= "Resequencing",
                       mouse = dat$mouse.id)

omm_metadata <- rbind(df_claudia, df_ab, df_reseq)

usethis::use_data(omm_metadata, overwrite = TRUE)