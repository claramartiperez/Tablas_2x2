############################################################
# 1. VALIDACIÓN Y UTILIDADES
############################################################

# Data de prueba
# tab_test <- matrix(
#   c(9, 50,
#     1, 940),
#   nrow = 2,
#   byrow = TRUE,
#   dimnames = list(
#     Test = c("Test +", "Test -"),
#     Real = c("R+", "R-")
#   )
# )
# 
# tab_test
# 
# sens <- 9 / (9 + 1)       # 0.9000
# spec <- 940 / (50 + 940)  # 0.9495
# ppv  <- 9 / (9 + 50)      # 0.1525
# npv  <- 940 / (1 + 940)   # 0.9989
# acc  <- (9 + 940) / 1000  # 0.9490

check_2x2 <- function(tab){
  if(!is.matrix(tab) || any(dim(tab) != c(2,2)))
    stop("Input must be a 2x2 contingency table.")
  if(any(tab < 0))
    stop("Frequencies must be non-negative.")
  invisible(TRUE)
}

apply_correction <- function(tab, correction = c("none","haldane")){
  correction <- match.arg(correction)
  if(correction == "haldane" && any(tab == 0)){
    tab <- tab + 0.5
    attr(tab,"correction") <- "Haldane–Anscombe (0.5)"
  }
  tab
}

get_cells <- function(tab){
  check_2x2(tab)
  list(
    a = tab[1,1],  # TP
    b = tab[1,2],  # FP
    c = tab[2,1],  # FN
    d = tab[2,2],  # TN
    n = sum(tab)
  )
}

############################################################
# 2. BOOTSTRAP GENÉRICO
############################################################

bootstrap_se <- function(tab, stat_fun, B = 1000){
  cc <- get_cells(tab)
  df <- data.frame(
    obs  = c(rep(1,cc$a),rep(1,cc$c),rep(0,cc$b),rep(0,cc$d)),
    pred = c(rep(1,cc$a),rep(0,cc$c),rep(1,cc$b),rep(0,cc$d))
  )
  vals <- replicate(B,{
    dfi <- df[sample.int(nrow(df),replace=TRUE),]
    tab_i <- matrix(
      c(sum(dfi$pred==1 & dfi$obs==1),  # TP
        sum(dfi$pred==1 & dfi$obs==0),  # FP
        sum(dfi$pred==0 & dfi$obs==1),  # FN
        sum(dfi$pred==0 & dfi$obs==0)), # TN
      nrow = 2,
      byrow = TRUE)
    stat_fun(tab_i,SE=FALSE)$statistic
  })
  sd(vals,na.rm=TRUE)
}
############################################################
# 3. MEDIDAS BÁSICAS
############################################################

ind_sens <- function(tab,SE=TRUE,method="analytic",correction="none",B=1000){
  tab <- apply_correction(tab,correction); cc <- get_cells(tab)
  est <- cc$a/(cc$a+cc$c)
  if(!SE) return(list(statistic=est))
  se <- if(method=="analytic") sqrt(est*(1-est)/(cc$a+cc$c)) else bootstrap_se(tab,ind_sens,B)
  list(statistic=est,se=se)
}

ind_spec <- function(tab,SE=TRUE,method="analytic",correction="none",B=1000){
  tab <- apply_correction(tab,correction); cc <- get_cells(tab)
  est <- cc$d/(cc$b+cc$d)
  if(!SE) return(list(statistic=est))
  se <- if(method=="analytic") sqrt(est*(1-est)/(cc$b+cc$d)) else bootstrap_se(tab,ind_spec,B)
  list(statistic=est,se=se)
}

ind_ppv <- function(tab,SE=TRUE,method="analytic",correction="none",B=1000){
  tab <- apply_correction(tab,correction); cc <- get_cells(tab)
  est <- cc$a/(cc$a+cc$b)
  if(!SE) return(list(statistic=est))
  se <- if(method=="analytic") sqrt(est*(1-est)/(cc$a+cc$b)) else bootstrap_se(tab,ind_ppv,B)
  list(statistic=est,se=se)
}

ind_npv <- function(tab,SE=TRUE,method="analytic",correction="none",B=1000){
  tab <- apply_correction(tab,correction); cc <- get_cells(tab)
  est <- cc$d/(cc$c+cc$d)
  if(!SE) return(list(statistic=est))
  se <- if(method=="analytic") sqrt(est*(1-est)/(cc$c+cc$d)) else bootstrap_se(tab,ind_npv,B)
  list(statistic=est,se=se)
}

ind_fpr <- function(tab,SE=TRUE,method="analytic",correction="none",B=1000){
  tab <- apply_correction(tab,correction); cc <- get_cells(tab)
  est <- cc$b/(cc$b+cc$d)
  if(!SE) return(list(statistic=est))
  se <- if(method=="analytic") sqrt(est*(1-est)/(cc$b+cc$d)) else bootstrap_se(tab,ind_fpr,B)
  list(statistic=est,se=se)
}

ind_fnr <- function(tab,SE=TRUE,method="analytic",correction="none",B=1000){
  tab <- apply_correction(tab,correction); cc <- get_cells(tab)
  est <- cc$c/(cc$a+cc$c)
  if(!SE) return(list(statistic=est))
  se <- if(method=="analytic") sqrt(est*(1-est)/(cc$a+cc$c)) else bootstrap_se(tab,ind_fnr,B)
  list(statistic=est,se=se)
}

ind_lr_pos <- function(tab,SE=TRUE,method="bootstrap",correction="none",B=1000){
  est <- ind_sens(tab,FALSE,correction=correction)$statistic /
    (1-ind_spec(tab,FALSE,correction=correction)$statistic)
  if(!SE) return(list(statistic=est))
  list(statistic=est,se=bootstrap_se(tab,ind_lr_pos,B))
}

ind_lr_neg <- function(tab,SE=TRUE,method="bootstrap",correction="none",B=1000){
  est <- (1-ind_sens(tab,FALSE,correction=correction)$statistic) /
    ind_spec(tab,FALSE,correction=correction)$statistic
  if(!SE) return(list(statistic=est))
  list(statistic=est,se=bootstrap_se(tab,ind_lr_neg,B))
}

ind_acc <- function(tab,SE=TRUE,method="analytic",correction="none",B=1000){
  tab <- apply_correction(tab,correction); cc <- get_cells(tab)
  est <- (cc$a+cc$d)/cc$n
  if(!SE) return(list(statistic=est))
  se <- if(method=="analytic") sqrt(est*(1-est)/cc$n) else bootstrap_se(tab,ind_acc,B)
  list(statistic=est,se=se)
}
############################################################
# 4. MÉTRICAS COMPUESTAS
############################################################

ind_f1 <- function(tab,SE=TRUE,method="bootstrap",correction="none",B=1000){
  tab <- apply_correction(tab,correction); cc <- get_cells(tab)
  est <- (2*cc$a)/(2*cc$a+cc$b+cc$c)
  if(!SE) return(list(statistic=est))
  list(statistic=est,se=bootstrap_se(tab,ind_f1,B))
}

ind_bacc <- function(tab,SE=TRUE,method="bootstrap",correction="none",B=1000){
  est <- (ind_sens(tab,FALSE,correction=correction)$statistic +
            ind_spec(tab,FALSE,correction=correction)$statistic)/2
  if(!SE) return(list(statistic=est))
  list(statistic=est,se=bootstrap_se(tab,ind_bacc,B))
}

ind_youden <- function(tab,SE=TRUE,method="bootstrap",correction="none",B=1000){
  est <- ind_sens(tab,FALSE,correction=correction)$statistic +
    ind_spec(tab,FALSE,correction=correction)$statistic - 1
  if(!SE) return(list(statistic=est))
  list(statistic=est,se=bootstrap_se(tab,ind_youden,B))
}

ind_markedness <- function(tab,SE=TRUE,method="bootstrap",correction="none",B=1000){
  est <- ind_ppv(tab,FALSE,correction=correction)$statistic +
    ind_npv(tab,FALSE,correction=correction)$statistic - 1
  if(!SE) return(list(statistic=est))
  list(statistic=est,se=bootstrap_se(tab,ind_markedness,B))
}
############################################################
# 5. MEDIDAS DE ASOCIACIÓN
############################################################

ind_phi <- function(tab, SE = TRUE, method = c("analytic", "bootstrap"),
                    correction = "none", B = 1000){
  method <- match.arg(method)
  tab <- apply_correction(tab, correction); cc <- get_cells(tab)
  
  est <- (cc$a*cc$d - cc$b*cc$c) /
    sqrt((cc$a+cc$b)*(cc$a+cc$c)*(cc$d+cc$b)*(cc$d+cc$c))
  
  if(!SE) return(list(statistic = est))
  
  se <- switch(
    method,
    analytic  = (1 - est^2) / sqrt(cc$n - 3),
    bootstrap = bootstrap_se(tab, ind_phi, B)
  )
  
  list(statistic = est, se = se)
}

ind_mcc <- ind_phi

# Evaluando
# ind_phi(tab_test, method = "analytic")
# ind_phi(tab_test, method = "bootstrap", B = 5000)
# ind_mcc(tab_test, method = "analytic")
# ind_mcc(tab_test, method = "bootstrap", B = 5000)

ind_yule_q <- function(tab,SE=TRUE,method="bootstrap",correction="none",B=1000){
  tab <- apply_correction(tab,correction); cc <- get_cells(tab)
  est <- (cc$a*cc$d-cc$b*cc$c)/(cc$a*cc$d+cc$b*cc$c)
  if(!SE) return(list(statistic=est))
  list(statistic=est,se=bootstrap_se(tab,ind_yule_q,B))
}

ind_yule_y <- function(tab,SE=TRUE,method="bootstrap",correction="none",B=1000){
  tab <- apply_correction(tab,correction); cc <- get_cells(tab)
  est <- (sqrt(cc$a*cc$d)-sqrt(cc$b*cc$c)) /
    (sqrt(cc$a*cc$d)+sqrt(cc$b*cc$c))
  if(!SE) return(list(statistic=est))
  list(statistic=est,se=bootstrap_se(tab,ind_yule_y,B))
}

ind_dor <- function(tab,SE=TRUE,method="bootstrap",correction="haldane",B=1000){
  tab <- apply_correction(tab,correction); cc <- get_cells(tab)
  est <- (cc$a*cc$d)/(cc$b*cc$c)
  if(!SE) return(list(statistic=est))
  list(statistic=est,se=bootstrap_se(tab,ind_dor,B))
}
############################################################
# 6. MEDIDAS DE ACUERDO CORREGIDAS POR AZAR
############################################################

ind_kappa <- function(tab, SE = TRUE, method = c("bootstrap", "analytic"),
                      correction = "none", B = 1000){
  method <- match.arg(method)
  tab <- apply_correction(tab, correction); cc <- get_cells(tab)
  Po <- (cc$a + cc$d) / cc$n
  Pe <- ((cc$a+cc$b)*(cc$a+cc$c) + (cc$c+cc$d)*(cc$b+cc$d)) / cc$n^2
  est <- (Po - Pe) / (1 - Pe)
  
  if(!SE) return(list(statistic = est))
  
  se <- switch(
    method,
    bootstrap = bootstrap_se(tab, ind_kappa, B),
    analytic  = NA_real_
  )
  
  list(statistic = est, se = se)
}

ind_ac1 <- function(tab,SE=TRUE,method="bootstrap",correction="none",B=1000){
  tab <- apply_correction(tab,correction); cc <- get_cells(tab)
  Po <- (cc$a+cc$d)/cc$n
  p1 <- (cc$a+cc$b)/cc$n
  p2 <- (cc$a+cc$c)/cc$n
  p  <- (p1+p2)/2
  Pe <- 2*p*(1-p)
  est <- (Po-Pe)/(1-Pe)
  if(!SE) return(list(statistic=est))
  list(statistic=est,se=bootstrap_se(tab,ind_ac1,B))
}

ind_delta <- function(tab, SE = TRUE, method = "bootstrap",
                      correction = "none", B = 1000){
  tab <- apply_correction(tab, correction)
  cc <- get_cells(tab)
  est <- (cc$a + cc$d - 2 * sqrt(cc$b * cc$c)) / cc$n
  
  if(!SE) return(list(statistic = est))
  
  list(
    statistic = est,
    se = bootstrap_se(tab, ind_delta, B)
  )
}

#######################
# Funciones adicionales
#######################

ind_informedness <- ind_youden

ind_fdr <- function(tab, SE = TRUE, method = "bootstrap", correction = "none", B = 1000){
  est <- 1 - ind_ppv(tab, FALSE, correction = correction)$statistic
  if(!SE) return(list(statistic = est))
  list(statistic = est, se = bootstrap_se(tab, ind_fdr, B))
}

ind_frr <- function(tab, SE = TRUE, method = "bootstrap", correction = "none", B = 1000){
  est <- 1 - ind_npv(tab, FALSE, correction = correction)$statistic
  if(!SE) return(list(statistic = est))
  list(statistic = est, se = bootstrap_se(tab, ind_frr, B))
}

ind_pcui <- function(tab, SE = TRUE, method = "bootstrap", correction = "none", B = 1000){
  est <- ind_sens(tab, FALSE, correction = correction)$statistic *
    ind_ppv(tab, FALSE, correction = correction)$statistic
  if(!SE) return(list(statistic = est))
  list(statistic = est, se = bootstrap_se(tab, ind_pcui, B))
}

ind_ncui <- function(tab, SE = TRUE, method = "bootstrap", correction = "none", B = 1000){
  est <- ind_spec(tab, FALSE, correction = correction)$statistic *
    ind_npv(tab, FALSE, correction = correction)$statistic
  if(!SE) return(list(statistic = est))
  list(statistic = est, se = bootstrap_se(tab, ind_ncui, B))
}

ind_ii <- function(tab, SE = TRUE, method = "bootstrap", correction = "none", B = 1000){
  est <- 2 * ind_acc(tab, FALSE, correction = correction)$statistic - 1
  if(!SE) return(list(statistic = est))
  list(statistic = est, se = bootstrap_se(tab, ind_ii, B))
}

ind_csi <- function(tab, SE = TRUE, method = "bootstrap", correction = "none", B = 1000){
  tab <- apply_correction(tab, correction); cc <- get_cells(tab)
  est <- cc$a / (cc$a + cc$b + cc$c)
  if(!SE) return(list(statistic = est))
  list(statistic = est, se = bootstrap_se(tab, ind_csi, B))
}

############################################################
# 7. FUNCIÓN RESUMEN FINAL
############################################################

ind_summary <- function(tab,SE=TRUE,method="analytic",correction="none",B=1000){
  rbind(
    sens   = unlist(ind_sens(tab,SE,method,correction,B)),
    spec   = unlist(ind_spec(tab,SE,method,correction,B)),
    ppv    = unlist(ind_ppv(tab,SE,method,correction,B)),
    npv    = unlist(ind_npv(tab,SE,method,correction,B)),
    fpr    = unlist(ind_fpr(tab,SE,method,correction,B)),
    fnr    = unlist(ind_fnr(tab,SE,method,correction,B)),
    lr_pos = unlist(ind_lr_pos(tab,SE,"bootstrap",correction,B)),
    lr_neg = unlist(ind_lr_neg(tab,SE,"bootstrap",correction,B)),
    acc    = unlist(ind_acc(tab,SE,method,correction,B)),
    f1     = unlist(ind_f1(tab,SE,"bootstrap",correction,B)),
    bacc   = unlist(ind_bacc(tab,SE,"bootstrap",correction,B)),
    youden = unlist(ind_youden(tab,SE,"bootstrap",correction,B)),
    mark   = unlist(ind_markedness(tab,SE,"bootstrap",correction,B)),
    phi    = unlist(ind_phi(tab,SE,"analytic",correction,B)),
    mcc    = unlist(ind_mcc(tab,SE,"analytic",correction,B)),
    yule_q = unlist(ind_yule_q(tab,SE,"bootstrap",correction,B)),
    yule_y = unlist(ind_yule_y(tab,SE,"bootstrap",correction,B)),
    dor    = unlist(ind_dor(tab,SE,"bootstrap","haldane",B)),
    kappa  = unlist(ind_kappa(tab,SE,"bootstrap",correction,B)),
    ac1    = unlist(ind_ac1(tab,SE,"bootstrap",correction,B)),
    delta = unlist(ind_delta(tab, SE, "bootstrap", correction, B)),
    informedness = unlist(ind_informedness(tab,SE,"bootstrap",correction,B)),
    fdr    = unlist(ind_fdr(tab,SE,"bootstrap",correction,B)),
    frr    = unlist(ind_frr(tab,SE,"bootstrap",correction,B)),
    pcui   = unlist(ind_pcui(tab,SE,"bootstrap",correction,B)),
    ncui   = unlist(ind_ncui(tab,SE,"bootstrap",correction,B)),
    ii     = unlist(ind_ii(tab,SE,"bootstrap",correction,B)),
    csi    = unlist(ind_csi(tab,SE,"bootstrap",correction,B))
  )
}

# round(ind_summary(tab_test),4)
# 
# #### 
# # Validación por caret
# real <- c(rep("pos", 10), rep("neg", 990))
# pred <- c(rep("pos", 9), "neg", rep("pos", 50), rep("neg", 940))
# 
# caret::confusionMatrix(
#   data = factor(pred, levels = c("pos", "neg")),
#   reference = factor(real, levels = c("pos", "neg")),
#   positive = "pos"
# )