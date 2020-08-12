bug <- "Akkermansia_muciniphila"
# claudia
dat <- readRDS("data-raw/omm_claudia_new.rds")
dat$sample <-  paste0("claudia ", dat$sample)
dat$group <- dat$mouse.group
df_claudia <- data.frame(sample = dat$sample,
                         group = dat$group, 
                         day = NA)

dat <- readRDS("data-raw/omm_ab.rds")

dat$sample <- paste0("AB study ",dat$mouse.id, " d",dat$day, " ",
                     group = dat$group, sample = dat$mouse.group)

df_ab <- data.frame(sample = dat$sample,
                    group = dat$mouse.group,
                    day = dat$day)

# reseq
dat <- readRDS("data-raw/reseq.rds")
dat$sample <- paste0("reseq study ",dat$mouse.id, " d",dat$day, "", dat$mouse.group)

df_reseq <- data.frame(sample = dat$sample, group = dat$mouse.group, day = NA)

omm2_metadata <- rbind(df_claudia, df_ab, df_reseq)
usethis::use_data(omm2_metadata, overwrite = TRUE)