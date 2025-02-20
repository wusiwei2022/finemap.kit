#' @export
finemap.wakefield = function(data, chrom = NULL, pos.lower = NULL, pos.upper = NULL, pthresh = 0.0001, wake = 0.04, cs = 0.95){
  if(is.null(chrom) & is.null(pos.lower) & is.null(pos.upper)){
    message("Run Fine Mapping on the complete data")}else 
    if(is.numeric(chrom) & is.numeric(pos.lower) & is.numeric(pos.upper)){
      message("Run Fine Mapping on chromosome ", chrom, "starting from ", pos.lower, " to ", pos.upper)}
  else{stop("Please provide all of the chrom, pos.lower, and pos.upper. Otherwise please provide set all of them as NULL to run on the complete data")
       data = data %>% filter(chr = chrom & pos >= pos.lower & pos <= pos.upper)}
  ## Check if the minimum p value < pthresh
  if(min(data$p) > pthresh){
    message(paste0("No variants in this genomic region reach p threshold of ", pthresh))
    return(NULL)}
  ## Calculate Bayes Factors
  data = data %>% mutate(z = beta/se)
  ## Calculate variances
  data = data %>% mutate(var = se^2)
  ## Calculate approximate Bayes Factors using Wakerfield formula
  data = data %>% mutate(abf = 1/(  (1 / (sqrt(1 - (wake/(var + wake)))) )*exp( -((z * z)/2) * (wake/(var + wake)))  ))
  ## Change Bayes factors to zero for those not meeting P threshold
  data = data %>% mutate(abf2 = ifelse(p < pthresh, abf, 0))
  ## Calculate the posterior probability
  abf2sum = sum(data$abf2, na.rm = TRUE)
  data = data %>% mutate(postprob = abf2/abf2sum)
  ## Calculate the 95% Credible Set
  data = data[order(data$postprob, decreasing = TRUE), ]
  data = data %>% mutate(accum.postprob = cumsum(postprob))
  ## Output the result
  out = list()
  out[["Post.Prob"]] = data
  out[["Cred.Set"]] = data[c(1, which(data$accum.postprob < cs) +1), ] %>% select("rsid")
  return(out)
}
