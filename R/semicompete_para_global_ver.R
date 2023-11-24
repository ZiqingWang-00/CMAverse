# updated semicompete script
library(tidyverse)
library(mstate)
library(foreach)
library(doParallel)

# semicompete function code
mstate_formula = function(exposure, mediator, basec, formula_terms){
  # create the formula for multistate modeling if not specified
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
    ## create the formula for multistate modeling
    mstate_formula = as.formula(paste("Surv(Tstart, Tstop, status) ~ ", terms, "+ strata(trans)"))
    return(mstate_formula)
  }
}

make_boot_ind = function(data, nboot){
  ## create nboot bootstrapping indices
  boot_ind = lapply(1:nboot, function(i){
    ind = sample(1:nrow(data), nrow(data), replace=T)
    data.frame(boot_iter = i, boot_ind = I(list(ind)))
  })
  boot_ind_df = do.call(rbind, boot_ind)
  return(boot_ind_df)
}

## function to convert the original survival data to mstate format
make_mstate_dat = function(dat, mediator, outcome, mediator_event, event, trans, covs_df){
  ## convert the bootstrap data to multistate format
  mstate_dat <- msprep(time = c(NA, mediator, outcome), status = c(NA, mediator_event, event),
                         data = dat, trans = trans, keep = covs_df)
  #print(str(mstate_data))  # Add this line for debugging
  mstate_dat <- expand.covs(mstate_dat, covs_df, append = TRUE, longnames = FALSE)
  return(mstate_dat)
}

## function to make data frames for msfit
### fixed newd: same across all bootstrap samples
fixed_newd = function(mstate_data, trans, a, astar, exposure, mediator, basec, basecval){
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
  fixed_newd = list(newd000 = newd000, newd010 = newd010)
  return(fixed_newd)
}

### dynamic newd (newd001): slightly different for each bootstrap sample
#### in cmest(), run this once for the biggest s. For all smaller s, only need to slice
dynamic_newd = function(newd000, time_vec, max_s, a, trans){
  # returns a list of newd data frames for msfit
  up_to_ind = which.max(abs(time_vec - max_s) == min(abs(time_vec - max_s)))
  ## set up newdata as needed (transition 3 now has M)
  ## C.1 = c(1,0,0) means C=1 for trans1, 0 for trans2, and 0 for trans3
  ## newd001: A.1=0, A.2=0, A.3=1
  newd001_list = lapply(time_vec[1:up_to_ind], function(time){
    newd001 = newd000
    newd001$A.3[3] = as.numeric(a)
    newd001$M.3[3] = time
    attr(newd001, "trans") <- trans
    class(newd001) <- c("msdata", "data.frame")
    return(newd001)
  })
  return(newd001_list)
}

# updated function calculates the RD on all s elements in the s_grid for each bootstrap sample i 
s_point_est <- function(i, s_grid, mstate_df, newd000, newd010, mstate_form,
                        exposure, mediator, outcome, basec, mediator_event, event, 
                        trans, method,
                        a){

  #environment(mstate_data) <- environment()
  #e <- environment()
  #print(ls(e))
  #mstate_data = e$mstate_data
  # create the multistate joint model object
  mstate_df <<- mstate_bootlist[[i]]
  joint_mod <- coxph(mstate_form, data = mstate_df, method = method)
  #joint_mod$data <- mstate_df
  #joint_mod2 <- eval.parent(joint_mod)
  print("fitted model")
  
  #print("check y-vector reconstruction")
  #mf <- model.frame(joint_mod2)
  #y <- model.response(mf)
  #y2 <- joint_mod2[["y"]]
  #print(paste("is.null(y2) ==", is.null(y2)))
  #print(dim(y2))
  #print(dim(y))

  # use msfit() to get predicted cumulative hazards data frames
  cumhaz000_msfit = msfit(joint_mod, newd000, trans=trans)
  #cumhaz000_msfit = predict(joint_mod2, newd000)
  cumhaz010_msfit = msfit(joint_mod, newd010, trans=trans)
  print("computed the first two cumhaz dfs")
  # extract cumulative hazards from the msfit objects
  cumhaz000 = cumhaz000_msfit$Haz
  cumhaz010 = cumhaz010_msfit$Haz
  # extract transitions that are needed
  ## time var the same for all data frames below
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

  # get the index to slice the newd001 list
  time_vec = cumhaz000_trans1$time
  max_s = max(s_grid)
  newd001_list = dynamic_newd(newd000, time_vec, max_s, a, trans)
  # numerical integration
  ## use lapply to compute each integrand, corresponding to max_s
  ## for integral list for lower values of s, slice the list
  print("start creating integrand list")
  integrand_list = lapply(newd001_list, function(newd001){
    # create the msfit object, extract cumulative hazards data frame
    cumhaz001_msfit = msfit(joint_mod, newd001, trans=trans)
    cumhaz001 = cumhaz001_msfit$Haz
    # get transition-specific cumulative hazards
    cumhaz001_trans3 = subset(cumhaz001, trans==3)
    cumhaz001_trans3_s = cumhaz001_trans3 %>% arrange(abs(time-max_s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{12}(s|A=1, C=1)
    cumhaz000_trans3_s = cumhaz000_trans3 %>% arrange(abs(time-max_s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{12}(s|A=0, C=1)
    # get the row indices to evaluate the integrand over
    time_ind = which(time_vec == newd001[3,"M.3"])
    #time_ind = min(time_ind, up_to_ind)
    if (time_ind < length(newd001_list)){
      P_01_integrand = (exp(-cumhaz000_trans1$Haz[time_ind]-cumhaz000_trans2$Haz[time_ind])*haz000_trans1[time_ind]*exp(-cumhaz000_trans3_s+cumhaz000_trans3$Haz[time_ind])) * (min(max_s,cumhaz000_trans1$time[time_ind+1])-cumhaz000_trans1$time[time_ind])
      P_g_01_integrand = (exp(-cumhaz000_trans1$Haz[time_ind]-cumhaz010_trans2$Haz[time_ind])*haz000_trans1[time_ind]*exp(-cumhaz001_trans3_s+cumhaz001_trans3$Haz[time_ind])) * (min(max_s,cumhaz000_trans1$time[time_ind+1])-cumhaz000_trans1$time[time_ind])
    }else{
      P_01_integrand = 0
      P_g_01_integrand = 0
    }
    return(c(P_01 = P_01_integrand, P_g_01 = P_g_01_integrand))
  })
  print("finished creating integrand list")

  ## sum up the individual integrands as the estimate for P_01
  RD_vec = rep(NA, length(s_grid))
  for (i in 1:length(RD_vec)){
    s = s_grid[i]
    up_to_ind = which.max(abs(time_vec - s) == min(abs(time_vec - s)))
    sums = colSums(do.call(rbind, integrand_list[1:up_to_ind]))
    P_01 = sums[1]
    P_g_01 = sums[2]
    # estimate cumulative hazard up to the end time s
    cumhaz000_trans1_s = cumhaz000_trans1 %>% arrange(abs(time-s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{01}(s|A=0,C=1)
    cumhaz000_trans2_s = cumhaz000_trans2 %>% arrange(abs(time-s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{02}(s|A=0,C=1)
    cumhaz010_trans2_s = cumhaz010_trans2 %>% arrange(abs(time-s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{02}(s|A=1,C=1)
    cumhaz000_trans3_s = cumhaz000_trans3 %>% arrange(abs(time-s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{12}(s|A=0, C=1)
    # P_00 and P_g_00
    P_00 = exp(-cumhaz000_trans1_s-cumhaz000_trans2_s)
    P_g_00 = exp(-cumhaz000_trans1_s-cumhaz010_trans2_s)
    # compute RD
    RD_vec[i] = (P_g_00 + P_g_01) - (P_00 + P_01)
  }

  return(RD_vec)
}

## TEST RUN
# import sample data
dtSurv1 = read_csv("/Users/apple/Desktop/CAUSAL/project/CMAverse_validation/simulated_dat/raw_dat1.csv")
#mstate_dtSurv1 = read_csv("/Users/apple/Desktop/CAUSAL/project/CMAverse_validation/simulated_dat/mstate_dat1.csv")
# add another binary covariate C for generality
dtSurv1$C1 = rbinom(nrow(dtSurv1), 1, 0.3)
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
formula_terms = c("A", "M", "C", "C1", "A*M")
method = "breslow"
s_grid = seq(0,1,length.out = 5)

boot_ind_df = make_boot_ind(data=data, nboot=nboot)
mstate_form <- mstate_formula(exposure=exposure, mediator=mediator, basec=basec, formula_terms=formula_terms)
mstate_form <- as.formula(mstate_form)
## set up relevant quantities for multistate data prep
trans = transMat(x=list(c(2, 3), c(3), c()), names=c(exposure, mediator, outcome)) # set up transition matrix
covs_df = c(exposure, mediator, basec) # transition-dependent covariates
## extract the time vector for making newd001 list
mstate_data_orig = make_mstate_dat(dat=data, mediator=mediator, outcome=outcome, mediator_event=mediator_event, event=event, trans=trans, covs_df=covs_df)
fixed_newd = fixed_newd(mstate_dat=mstate_data_orig, trans=trans, a=a, astar=astar, exposure=exposure, mediator=mediator, basec=basec, basecval=basecval)
joint_mod_orig = coxph(formula=mstate_form, data = mstate_data_orig, method = method)
cumhaz000_msfit = msfit(joint_mod_orig, fixed_newd[[1]], trans=trans)
cumhaz000 = cumhaz000_msfit$Haz
cumhaz000_trans1 = subset(cumhaz000, trans==1)
time_vec = cumhaz000_trans1$time
newd001_list = dynamic_newd(fixed_newd[[1]], time_vec, max_s=max(s_grid), a, trans)
## create a list of mstate data that corresponds to the bootstrap samples
boot_ind_list = boot_ind_df$boot_ind
mstate_bootlist <- lapply(boot_ind_list, function(ind){
  boot_dat = data[ind,]
  mstate_boot_dat <- make_mstate_dat(dat=boot_dat, mediator, outcome, mediator_event, event, trans, covs_df)
  return(mstate_boot_dat)
})

# run
s_point_est(s_grid=s_grid, mstate_df=mstate_bootlist[[2]], newd000=fixed_newd[[1]], newd010=fixed_newd[[2]], mstate_form=mstate_form,
         exposure=exposure, mediator=mediator, outcome=outcome, basec=basec, mediator_event=mediator_event, event=event,
         trans=trans, method=method,
         a=a)
