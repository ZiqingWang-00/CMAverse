# updated semicompete script
library(tidyverse)
library(mstate)

#semicompete2 <- function(method = "breslow"){
#
#}

# import sample data
dtSurv1 = read_csv("/Users/apple/Desktop/CAUSAL/project/CMAverse_validation/simulated_dat/raw_dat1.csv")
mstate_dtSurv1 = read_csv("/Users/apple/Desktop/CAUSAL/project/CMAverse_validation/simulated_dat/mstate_dat1.csv")

# add another binary covariate C for generality
dtSurv1$C1 = rbinom(nrow(dtSurv1), 1, 0.3)

# set up function arguments
data = dtSurv1
exposure = "A"
mediator = "M"
outcome = "S"
event = "ind_S"
mediator_event = "ind_M"
basec = c("C", "C1")
basecval = c("C" = "1", "C1" = "0")
astar = "0" # the control value for the exposure. Default is 0
a = "1" # the active value for the exposure. Default is 1.
nboot = 100
total_duration = max(dtSurv$S)
n_eval = 5 # calculate RD on 5 time points from [0, total_duration]
formula_terms = c("A", "M", "C", "C1", "A*M")
method = "breslow"

# function code
## set up relevant objects for later use
s_grid = seq(0, total_duration, length.out=n_eval)
RD_vec = rep(NA, length(s_grid))

## create the formula for multistate modeling if not specified
if (!is.null(formula_terms)){
  # convert formula_terms to multistate format
  for (i in 1:length(formula_terms)){
    if (formula_terms[i] == exposure){
      terms = paste(exposure, c(".1", ".2", ".3"), sep = "", collapse="+")
    }else if (formula_terms[i] == mediator){
      terms = paste(terms, paste(mediator, ".3", sep = ""), sep="+")
    }else if (formula_terms[i] %in% basec){
      j = which(basec == formula_terms[i])
      terms = paste(terms, paste(basec[j], c(".1", ".2", ".3"), sep="", collapse="+"), sep="+")
    }else if (grepl(exposure, formula_terms[i])==T & grepl(mediator, formula_terms[i])==T){ # A.3*M.3
      terms = paste(terms, 
                    paste(paste(exposure, ".3", sep=""), "*", paste(mediator, ".3", sep=""), sep=""),
                    sep="+")
    }
  }
  # create the formula for multistate modeling
  mstate_formula = as.formula(paste("Surv(Tstart, Tstop, status) ~ ", terms, "+ strata(trans)"))
}

## prepare the input data in multistate format
trans = transMat(x=list(c(2, 3), c(3), c()), names=c(exposure, mediator, outcome)) # set up transition matrix
covs_df = c(exposure, mediator, basec) # transition-dependent covariates
mstate_data = msprep(time = c(NA, mediator, outcome), status = c(NA, mediator_event, event),
                     data = data, trans = trans, keep = covs_df)
mstate_data = expand.covs(mstate_data, covs_df, append = TRUE, longnames = FALSE)

## create the multistate joint model object
joint_mod = coxph(mstate_formula, data = mstate_data, method = method)

## set up new data for msfit() in the order of exposure, mediator, and outcome
### newd000: A=0 for all 3 transitions ()
A000 = matrix(rep(as.numeric(astar),9), nrow=3)
colnames(A000) = grep(paste0("\\b", exposure, "\\."), names(mstate_data), value=TRUE, perl=TRUE)
M3000 = matrix(rep(0,3), nrow=3)
colnames(M3000) = grep(paste0("\\b", mediator, "\\.3"), names(mstate_data), value=TRUE, perl=TRUE)
C000 = matrix(rep(0, 3*3*length(basec)), nrow=3)
for (i in 1:length(basec)){
  start = 1+(i-1)*3
  end = 1+(i-1)*3+2
  C000[,start:end] = diag(x=as.numeric(basecval[i]), 3)
}
C_colnames = rep(NA, 3*length(basec))
for (i in 1:length(basec)){
  start = 1+(i-1)*3
  end = 1+(i-1)*3+2
  C_colnames[start:end] = grep(paste0("\\b", basec[i], "\\."), names(mstate_data), value=TRUE, perl=TRUE)
}
newd000 = data.frame(A000, M3000, C000)
colnames(newd000)[(ncol(A000)+ncol(M3000)+1):ncol(newd000)] = C_colnames
newd000$trans = c(1,2,3)
newd000$strata = c(1,2,3)
attr(newd000,"trans") <- trans
class(newd000) <- c("msdata", "data.frame")

### newd010: A=0 for trans1, A=1 for trans2, A=0 for trans3
newd010 = newd000
newd010$A.2[2] = as.numeric(a)
attr(newd010,"trans") <- trans
class(newd010) <- c("msdata", "data.frame")

## for each time point, compute RD
for (j in 1:length(s_grid)){
  print(j)
  s = s_grid[j] # one particular time point to compute the probability of surviving beyond
  
  # use msfit() to get predicted cumulative hazards data frames
  cumhaz000_msfit = msfit(joint_mod, newd000, trans=trans)
  cumhaz010_msfit = msfit(joint_mod, newd010, trans=trans)
  
  # extract cumulative hazards from the msfit objects
  cumhaz000 = cumhaz000_msfit$Haz
  cumhaz010 = cumhaz010_msfit$Haz
  
  # extract transitions that are needed
  # time var the same for all data frames below
  cumhaz000_trans1 = subset(cumhaz000, trans==1)
  cumhaz000_trans2 = subset(cumhaz000, trans==2)
  cumhaz000_trans3 = subset(cumhaz000, trans==3) 
  cumhaz010_trans2 = subset(cumhaz010, trans==2)
  
  # populate the hazard grid
  n_grid = nrow(cumhaz000_trans1)
  haz000_trans1 = rep(NA, n_grid) # grid for hazard \alpha_{01}
  for (i in 1:(nrow(cumhaz000_trans1)-1)){
    haz000_trans1[i]=(cumhaz000_trans1$Haz[i+1]-cumhaz000_trans1$Haz[i])/(cumhaz000_trans1$time[i+1]-cumhaz000_trans1$time[i])
  }
  
  # estimate cumulative hazard up to the end time s
  cumhaz000_trans1_s = cumhaz000_trans1 %>% arrange(abs(time-s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{01}(s|A=0,C=1)
  cumhaz000_trans2_s = cumhaz000_trans2 %>% arrange(abs(time-s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{02}(s|A=0,C=1)
  cumhaz010_trans2_s = cumhaz010_trans2 %>% arrange(abs(time-s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{02}(s|A=1,C=1)
  cumhaz000_trans3_s = cumhaz000_trans3 %>% arrange(abs(time-s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{12}(s|A=0, C=1)
  
  # numerical integration
  # numerical integration
  i = 1
  P_01 = 0
  P_g_01 = 0
  while(cumhaz000_trans1$time[i] < s & i < nrow(cumhaz000_trans1)){
    # time = df$time[i] for M.3 = time for newd
    time = cumhaz000_trans1$time[i] 
    
    # set up newdata as needed (transition 3 now has M)
    ## nameing order: A.1=0, A.2 = 0, A.3=0
    ## C.1 = c(1,0,0) means C=1 for trans1, 0 for trans2, and 0 for trans3
    ## newd001: A.1=0, A.2=0, A.3=1
    newd001 = newd000
    newd001$A.3[3] = as.numeric(a)
    newd001$M.3[3] = time
    attr(newd001, "trans") <- trans
    class(newd001) <- c("msdata", "data.frame")
    
    # create the msfit object, extract cumulative hazards data frame
    cumhaz001_msfit = msfit(joint_mod, newd001, trans=trans)
    cumhaz001 = cumhaz001_msfit$Haz
    # get transition-specific cumulative hazards
    cumhaz001_trans3 = subset(cumhaz001, trans==3)
    cumhaz001_trans3_s = cumhaz001_trans3 %>% arrange(abs(time-s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{12}(s|A=1, C=1)
    
    P_01 = P_01 + (exp(-cumhaz000_trans1$Haz[i]-cumhaz000_trans2$Haz[i])*haz000_trans1[i]*exp(-cumhaz000_trans3_s+cumhaz000_trans3$Haz[i])) * (min(s,cumhaz000_trans1$time[i+1])-cumhaz000_trans1$time[i])
    P_g_01 = P_g_01  + (exp(-cumhaz000_trans1$Haz[i]-cumhaz010_trans2$Haz[i])*haz000_trans1[i]*exp(-cumhaz001_trans3_s+cumhaz001_trans3$Haz[i])) * (min(s,cumhaz000_trans1$time[i+1])-cumhaz000_trans1$time[i])
    
    i=i+1
  }
  
  # P_00 and P_g_00
  P_00 = exp(-cumhaz000_trans1_s-cumhaz000_trans2_s)
  P_g_00 = exp(-cumhaz000_trans1_s-cumhaz010_trans2_s)
  
  # compute RD
  RD_vec[j] = (P_g_00 + P_g_01) - (P_00 + P_01)
  
}

## start bootstrap for loop to compute RD and TE, and SD = TE - RD

# collect results in an output dataframe



