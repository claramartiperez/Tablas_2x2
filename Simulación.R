############################################################
# FIGURA 2: MCC vs Accuracy
############################################################

set.seed(123)

n_sim <- 2000

results_acc <- data.frame(
  mcc = numeric(n_sim),
  acc = numeric(n_sim)
)

for(i in 1:n_sim){
  
  # Tabla aleatoria
  tab <- matrix(sample(1:100, 4, replace = TRUE), nrow = 2)
  
  # Calcular métricas
  mcc_val <- ind_mcc(tab, SE = FALSE)$statistic
  acc_val <- ind_acc(tab, SE = FALSE)$statistic
  
  results_acc$mcc[i] <- mcc_val
  results_acc$acc[i] <- acc_val
}

# Plot
plot(results_acc$mcc, results_acc$acc,
     xlab = "MCC",
     ylab = "Accuracy",
     main = "",
     pch = 16)

# Línea de tendencia
abline(lm(acc ~ mcc, data = results_acc), lwd = 2)

############################################################
# FIGURA 3: MCC vs F1
############################################################

set.seed(123)

n_sim <- 2000

results <- data.frame(
  mcc = numeric(n_sim),
  f1  = numeric(n_sim)
)

for(i in 1:n_sim){
  
  # Tabla aleatoria
  tab <- matrix(sample(1:100, 4, replace = TRUE), nrow = 2)
  
  # Calcular métricas
  mcc_val <- ind_mcc(tab, SE = FALSE)$statistic
  f1_val  <- ind_f1(tab, SE = FALSE)$statistic
  
  results$mcc[i] <- mcc_val
  results$f1[i]  <- f1_val
}

# Plot
plot(results$mcc, results$f1,
     xlab = "MCC",
     ylab = "F1-score",
     main = "",
     pch = 16)

abline(lm(f1 ~ mcc, data = results))

############################################################
# FIGURA 4: PPV vs prevalencia
############################################################

set.seed(123)

n_sim <- 2000

results_ppv <- data.frame(
  prevalence = numeric(n_sim),
  ppv = numeric(n_sim)
)

for(i in 1:n_sim){
  
  # Tabla con distintos tamaños
  TP <- sample(1:100,1)
  FN <- sample(1:100,1)
  FP <- sample(1:100,1)
  TN <- sample(1:100,1)
  
  tab <- matrix(c(TP, FN, FP, TN), nrow = 2, byrow = TRUE)
  
  # Prevalencia
  prev <- (TP + FN) / sum(tab)
  
  # PPV
  ppv_val <- ind_ppv(tab, SE = FALSE)$statistic
  
  results_ppv$prevalence[i] <- prev
  results_ppv$ppv[i] <- ppv_val
}

# Plot
plot(results_ppv$prevalence, results_ppv$ppv,
     xlab = "Prevalencia",
     ylab = "PPV",
     main = "",
     pch = 16)

# Línea suavizada
lines(lowess(results_ppv$prevalence, results_ppv$ppv), lwd = 2)

############################################################
# FIGURA 5: PPV vs Prevalencia
############################################################

set.seed(123)

n_sim <- 2000

prev_seq <- seq(0.01, 0.99, length.out = 50)

results_prev <- data.frame(
  prevalence = numeric(),
  ppv = numeric()
)

for(p in prev_seq){
  
  for(i in 1:50){
    
    N <- 200
    
    TP <- round(p * N * runif(1, 0.6, 0.9))
    FN <- round(p * N - TP)
    FP <- round((1 - p) * N * runif(1, 0.1, 0.4))
    TN <- round((1 - p) * N - FP)
    
    tab <- matrix(c(TP, FN, FP, TN), nrow = 2, byrow = TRUE)
    
    ppv_val <- ind_ppv(tab, SE = FALSE)$statistic
    
    results_prev <- rbind(results_prev,
                          data.frame(prevalence = p, ppv = ppv_val))
  }
}

plot(results_prev$prevalence, results_prev$ppv,
     xlab = "Prevalencia",
     ylab = "PPV",
     main = "",
     pch = 16, col = rgb(0,0,1,0.3))

lines(lowess(results_prev$prevalence, results_prev$ppv), lwd = 2)

############################################################
# FIGURA 6: MCC y F1 vs Prevalencia
############################################################

results_metrics <- data.frame(
  prevalence = numeric(),
  mcc = numeric(),
  f1 = numeric()
)

for(p in prev_seq){
  
  for(i in 1:50){
    
    N <- 200
    
    TP <- round(p * N * runif(1, 0.6, 0.9))
    FN <- round(p * N - TP)
    FP <- round((1 - p) * N * runif(1, 0.1, 0.4))
    TN <- round((1 - p) * N - FP)
    
    tab <- matrix(c(TP, FN, FP, TN), nrow = 2, byrow = TRUE)
    
    results_metrics <- rbind(results_metrics,
                             data.frame(
                               prevalence = p,
                               mcc = ind_mcc(tab, SE = FALSE)$statistic,
                               f1  = ind_f1(tab, SE = FALSE)$statistic
                             )
    )
  }
}

plot(results_metrics$prevalence, results_metrics$mcc,
     xlab = "Prevalencia",
     ylab = "Métrica",
     main = "",
     pch = 16, col = rgb(1,0,0,0.3))

points(results_metrics$prevalence, results_metrics$f1,
       pch = 16, col = rgb(0,0,1,0.3))

legend("topright",
       legend = c("MCC", "F1-score"),
       col = c("red","blue"),
       pch = 16)

lines(lowess(results_metrics$prevalence, results_metrics$mcc), col="red", lwd=2)
lines(lowess(results_metrics$prevalence, results_metrics$f1), col="blue", lwd=2)

############################################################
# TABLA 4: COMPARACIÓN DE ERRORES ESTÁNDAR
############################################################

set.seed(123)

# Tabla base 
tab <- matrix(c(50,10,5,35), nrow=2, byrow=TRUE)

# Métricas con método analítico
sens_a <- ind_sens(tab, SE=TRUE, method="analytic")
spec_a <- ind_spec(tab, SE=TRUE, method="analytic")
ppv_a  <- ind_ppv(tab, SE=TRUE, method="analytic")
npv_a  <- ind_npv(tab, SE=TRUE, method="analytic")
acc_a  <- ind_acc(tab, SE=TRUE, method="analytic")

# Métricas con bootstrap
sens_b <- ind_sens(tab, SE=TRUE, method="bootstrap")
spec_b <- ind_spec(tab, SE=TRUE, method="bootstrap")
ppv_b  <- ind_ppv(tab, SE=TRUE, method="bootstrap")
npv_b  <- ind_npv(tab, SE=TRUE, method="bootstrap")
acc_b  <- ind_acc(tab, SE=TRUE, method="bootstrap")
f1_b   <- ind_f1(tab, SE=TRUE)
mcc_b  <- ind_mcc(tab, SE=TRUE)

# Tabla final
tabla_se <- data.frame(
  Metrica = c("Sensibilidad","Especificidad","PPV","NPV","Accuracy","F1-score","MCC"),
  SE_Analitico = c(sens_a$se, spec_a$se, ppv_a$se, npv_a$se, acc_a$se, NA, NA),
  SE_Bootstrap = c(sens_b$se, spec_b$se, ppv_b$se, npv_b$se, acc_b$se, f1_b$se, mcc_b$se)
)

tabla_se
############################################################
# FIGURA 7: COMPARACIÓN DE ERRORES ESTÁNDAR
############################################################

library(tidyr)
library(ggplot2)
## Warning: package 'ggplot2' was built under R version 4.4.3
tabla_long <- pivot_longer(tabla_se,
                           cols = c(SE_Analitico, SE_Bootstrap),
                           names_to = "Metodo",
                           values_to = "SE")

ggplot(tabla_long, aes(x=Metrica, y=SE, fill=Metodo)) +
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() +
  labs(x="", y="Error estándar")
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_bar()`).

############################################################
# TABLA 6: CELDAS NULAS (SIN vs CON CORRECCIÓN)
############################################################

# Tabla con cero
tab_zero <- matrix(c(10,0,0,40), nrow=2, byrow=TRUE)

# Sin corrección
res_sin <- ind_summary(tab_zero, SE=FALSE)

# Con corrección Haldane
res_con <- ind_summary(tab_zero, SE=FALSE, correction="haldane")

# Convertir a data.frame
tabla_zero <- data.frame(
  Metrica = rownames(res_sin),
  Sin_correccion = res_sin[,1],
  Con_correccion = res_con[,1],
  row.names = NULL
)

tabla_zero$Sin_correccion <- round(tabla_zero$Sin_correccion, 3)
tabla_zero$Con_correccion <- round(tabla_zero$Con_correccion, 3)

tabla_zero





