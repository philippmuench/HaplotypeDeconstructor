syn <- as.data.frame(matrix(0,30,10))
syn[1:10, 1:10] <- t(matrix(rep(1:10)/10, 10, 10)) # H1 in all 10 samples for first 10 SNPs with higher AF over time
syn[11:20, 6:8] <- matrix(0.9, 10, 3) # H2
syn[21:30, 1:5] <- matrix(0.1, 10, 6) # H3

#Heatmap(syn)
#dim(syn)
usethis::use_data(syn, overwrite = TRUE)