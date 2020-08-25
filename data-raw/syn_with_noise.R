
syn_with_noise <- as.data.frame(matrix(0,30,10))

syn_with_noise[1:10, 1:10] <- syn_with_noise[1:10, 1:10] +
  t(matrix(rep(1:10)/10, 10, 10)) # H1 in all 10 samples for first 10 SNPs with higher AF over time
syn_with_noise[11:20, 6:8] <- syn_with_noise[11:20, 6:8] + 
  matrix(0.9, 10, 3) # H2
syn_with_noise[21:30, 1:5] <- syn_with_noise[21:30, 1:5] + 
  matrix(0.5, 10, 6) # H3

usethis::use_data(syn_with_noise, overwrite = TRUE)