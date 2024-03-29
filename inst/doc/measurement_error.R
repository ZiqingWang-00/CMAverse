## -----------------------------------------------------------------------------
set.seed(1)
expit <- function(x) exp(x)/(1+exp(x))
n <- 10000
C1 <- rnorm(n, mean = 1, sd = 0.5)
C1_error <- C1 + rnorm(n, 0, 0.05)
C2 <- rbinom(n, 1, 0.6)
A <- rbinom(n, 1, expit(0.2 + 0.5*C1 + 0.1*C2))
mc <- matrix(c(0.9,0.1,0.1,0.9), nrow = 2)
A_error <- A
for (j in 1:2) {
  A_error[which(A_error == c(0,1)[j])] <-
    sample(x = c(0,1), size = length(which(A_error == c(0,1)[j])),
           prob = mc[, j], replace = TRUE)
}
M <- rbinom(n, 1, expit(1 + 2*A + 1.5*C1 + 0.8*C2))
Y <- rbinom(n, 1, expit(-3 - 0.4*A - 1.2*M + 0.5*A*M + 0.3*C1 - 0.6*C2))
data <- data.frame(A, A_error, M, Y, C1, C1_error, C2)

## -----------------------------------------------------------------------------
library(CMAverse)
cmdag(outcome = "Y", exposure = "A", mediator = "M",
      basec = c("C1", "C2"), postc = NULL, node = TRUE, text_col = "white")

## ----message=F,warning=F,results='hide'---------------------------------------
res_naive_cont <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                        mediator = "M", basec = c("C1_error", "C2"), EMint = TRUE,
                        mreg = list("logistic"), yreg = "logistic",
                        astar = 0, a = 1, mval = list(1), 
                        estimation = "paramfunc", inference = "delta")

## ----message=F,warning=F------------------------------------------------------
summary(res_naive_cont)

## ----message=F,warning=F,results='hide'---------------------------------------
res_rc_cont <- cmsens(object = res_naive_cont, sens = "me", MEmethod = "rc", 
                      MEvariable = "C1_error", MEvartype = "con", MEerror = 0.05)


## ----message=F,warning=F------------------------------------------------------
summary(res_rc_cont)

## ----message=F,warning=F,results='hide'---------------------------------------
res_simex_cont <- cmsens(object = res_naive_cont, sens = "me", MEmethod = "simex", 
                         MEvariable = "C1_error", MEvartype = "con", MEerror = 0.05)

## ----message=F,warning=F------------------------------------------------------
summary(res_simex_cont)

## ----message=F,warning=F,results='hide'---------------------------------------
res_naive_cat <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A_error",
                       mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                       mreg = list("logistic"), yreg = "logistic",
                       astar = 0, a = 1, mval = list(1), 
                       estimation = "paramfunc", inference = "delta")

## ----message=F,warning=F------------------------------------------------------
summary(res_naive_cat)

## ----message=F,warning=F,results='hide'---------------------------------------
res_simex_cat <- cmsens(object = res_naive_cat, sens = "me", MEmethod = "simex", 
                         MEvariable = "A_error", MEvartype = "cat", MEerror = list(mc))


## ----message=F,warning=F------------------------------------------------------
summary(res_simex_cat)

## ----message=F,warning=F,results='hide'---------------------------------------
res_true <- cmest(data = data, model = "rb", outcome = "Y", exposure = "A",
                       mediator = "M", basec = c("C1", "C2"), EMint = TRUE,
                       mreg = list("logistic"), yreg = "logistic",
                       astar = 0, a = 1, mval = list(1), 
                       estimation = "paramfunc", inference = "delta")

## ----message=F,warning=F------------------------------------------------------
summary(res_true)

