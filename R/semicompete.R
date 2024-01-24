# updated semicompete script
## JAN 2024 update
### ymreg = "cox"
### EMint = T or F
### revert back to using foreach instead of boot(); attempt to add progress bar
### avoids using the super-assignment operator <<-
#library(tidyverse)
#library(mstate)
#library(foreach)
#library(doParallel)
#library(doSNOW)
#library(progressr) ## use progressr for procession updates
#library(doFuture)  ## attaches also foreach and future

# semicompete function code
## function to create the ymreg formula
mstate_formula = function(data, exposure, mediator, basec, EMint){
  # create the formula for multistate modeling if not specified
  formula_terms = c(exposure, mediator, basec)
  for (i in 1:length(formula_terms)){
    if (formula_terms[i] == exposure){
      if (class(data[,exposure]) %in% c("numeric", "integer")){
        terms = paste(exposure, c(".1", ".2", ".3"), sep = "", collapse="+")
        # if EMint == T, put the interaction term in the formula as well
        if (EMint){
          terms = paste(terms, 
                        paste(paste(exposure, ".3", sep=""), "*", paste(mediator, ".3", sep=""), sep=""),
                        sep="+")
        }
      }else{ # if the exposure is a factor variable
        nlevel = length(levels(data[,exposure]))
        if (nlevel > 2){
          terms = paste(exposure, rep(seq(1,(nlevel-1)), each=3), ".", seq(1,3), sep="", collapse="+")
          exposure_terms = paste(exposure, rep(seq(1,(nlevel-1)), each=3), ".", seq(1,3), sep="")
        }else {
          terms = paste(exposure, ".", seq(1,3), sep="", collapse="+")
          exposure_terms = paste(exposure, ".", seq(1,3), sep="")
        }
        
        if (EMint){
          exp_trans3_terms = exposure_terms[startsWith(exposure_terms, exposure) & endsWith(exposure_terms, ".3")]
          int_terms = paste(exp_trans3_terms, "*", paste(mediator, ".3", sep=""), sep="")
          terms = paste(terms, paste(int_terms, collapse="+"), sep="+")
        }
      }
    }else if (formula_terms[i] == mediator){
      terms = paste(terms, paste(mediator, ".3", sep = ""), sep="+")
    }else if (formula_terms[i] %in% basec){
      j = which(basec == formula_terms[i])
      if (class(data[,formula_terms[i]]) %in% c("numeric", "integer")){
        terms = paste(terms, paste(basec[j], c(".1", ".2", ".3"), sep="", collapse="+"), sep="+")
      }else{
        nlevel = length(levels(data[,basec[j]]))
        if (nlevel > 2){
          terms = paste(terms, paste(basec[j], rep(seq(1,(nlevel-1)), each=3), ".", seq(1,3), sep="", collapse="+"), sep="+")
        }else {
          terms = paste(terms, paste(basec[j], ".", seq(1,3), sep="", collapse="+"), sep="+")
        }
      }
    }
  }
  ## create the formula for multistate modeling
  mstate_formula = as.formula(paste("Surv(Tstart, Tstop, status) ~ ", terms, "+ strata(trans)"))
  return(mstate_formula)
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
fixed_newd = function(mstate_data, trans, a, astar, exposure, mediator, basec, basecval){
  # convert to data frame
  mstate_data = data.frame(mstate_data)
  # create individual components
  ## exposure block
  if (class(mstate_data[,exposure]) %in% c("numeric", "integer")){
    A000 = matrix(rep(as.numeric(astar),9), nrow=3)
    colnames(A000) = paste(exposure, ".", seq(1,3), sep="")
  } else { # if the exposure is a factor
    nlevel = length(levels(mstate_data[,exposure])) # extract the number of levels (Including ref) in the factor variable
    level_selected = which(levels(mstate_data[,exposure]) == astar) # a = astar (ie, a inactive
    #print(paste("nlevel is ", nlevel))
    A000 = matrix(rep(0,3*3*(nlevel-1)), nrow=3) # (nlevel-1) dummy variables created in mstate formatted data
    if (level_selected > 1){ # level selected is not reference level
      start = 3 * level_selected - 5
      end = start + 2
      A000[,start:end] = diag(1,3)
    } 
    if (nlevel > 2){
      colnames(A000) = paste(exposure, rep(seq(1,(nlevel-1)), each=3), ".", seq(1,3), sep="")
    }else {
      colnames(A000) = paste(exposure, ".", seq(1,3), sep="")
    }
  }
  ## mediator block
  M3000 = matrix(rep(0,3), nrow=3) # time to mediator must be numeric
  colnames(M3000) = paste(mediator, ".", 3, sep="")
  ## covariate block
  C_blocks = list()
  for (i in 1:length(basec)){
    if (class(mstate_data[,basec[i]]) %in% c("numeric", "integer")){
      C_blocks[[i]] = diag(x=as.numeric(basecval[i]), 3)
      colnames(C_blocks[[i]]) = paste(basec[i], ".", seq(1,3), sep="")
    } else { # if variable type is factor
      nlevel = length(levels(mstate_data[,basec[i]])) # extract the number of levels (Including ref) in the factor variable
      level_selected = which(levels(mstate_data[,basec[i]]) == basecval[i])
      temp_mat = matrix(rep(0, 3*3*(nlevel-1)), nrow=3)
      # DEBUG FROM HERE
      if (nlevel > 2){
        if (level_selected > 1){
          start = 3*level_selected-5
          end = start + 2
          temp_mat[,start:end] = diag(1, 3)
        } 
        C_blocks[[i]] = temp_mat
        colnames(C_blocks[[i]]) = paste(basec[i], rep(seq(1,(nlevel-1)), each=3), ".", seq(1,3), sep="")
      } else {
        C_blocks[[i]] = temp_mat
        colnames(C_blocks[[i]]) = paste(basec[i], ".", seq(1,3), sep="")
      }
    }
  }
  C000 = do.call(cbind, C_blocks)
  
  ## combine components
  newd000 = data.frame(cbind(A000, M3000, C000))
  newd000$trans = c(1,2,3)
  newd000$strata = c(1,2,3)
  attr(newd000,"trans") <- trans
  class(newd000) <- c("msdata", "data.frame")

  # create newd010, newd100
  newd010 = newd000
  newd100 = newd000
  if (class(mstate_data[,exposure]) %in% c("numeric", "integer")){
    newd010[2,paste(exposure,".",2, sep="")] = as.numeric(a)
    newd100[1,paste(exposure,".",1, sep="")] = as.numeric(a)
  } else{
    nlevel = length(levels(mstate_data[,exposure]))
    level_active = which(levels(mstate_data[,exposure]) == a) 
    if (level_active > 1){ # if active value is not the reference level
      if (nlevel > 2){
        newd010[2,paste(exposure, (level_active-1), ".", 2, sep="")] = 1
        newd100[1,paste(exposure, (level_active-1), ".", 1, sep="")] = 1
      }else {
        newd010[2,paste(exposure,".",2, sep="")] = 1
        newd100[1,paste(exposure,".",1, sep="")] = 1
      }
    } 
  }
  return(list(newd000=newd000, newd010=newd010, newd100=newd100))
}

### dynamic newd (newd001): slightly different for each bootstrap sample
#### in cmest(), run this once for the biggest s. For all smaller s, only need to slice
dynamic_newd = function(mstate_data, exposure, mediator, newd000, time_vec, max_s, a, trans){
  # returns a list of newd data frames for msfit
  up_to_ind = which.max(abs(time_vec - max_s) == min(abs(time_vec - max_s)))
  ## set up newdata as needed (transition 3 now has M)
  ## C.1 = c(1,0,0) means C=1 for trans1, 0 for trans2, and 0 for trans3
  ## newd001: A.1=0, A.2=0, A.3=1
  newd001_list = lapply(time_vec[1:up_to_ind], function(time){
    newd001 = newd000
    if (class(mstate_data[,exposure]) %in% c("numeric", "integer")){
      newd001[3,paste(exposure, ".", 3, sep="")] = as.numeric(a)
    }else { # if the exposure is a factor
      nlevel = length(levels(mstate_data[,exposure]))
      level_active = which(levels(mstate_data[,exposure]) == a) 
      if (level_active > 1){
        if (nlevel > 2){
          newd001[3,paste(exposure, (level_active-1), ".", 3, sep="")] = 1
        }else{
          newd001[3,paste(exposure, ".", 3, sep="")] = 1
        }
      }
    }
    newd001[3,paste(mediator, ".", 3, sep="")] = time
    attr(newd001, "trans") <- trans
    class(newd001) <- c("msdata", "data.frame")
    return(newd001)
  })
  return(newd001_list)
}

## function to make bootstrap indices
make_boot_ind = function(data, nboot){
  ## create nboot bootstrapping indices
  boot_ind = lapply(1:nboot, function(i){
    ind = sample(1:nrow(data), nrow(data), replace=T)
    data.frame(boot_iter = i, boot_ind = I(list(ind)))
  })
  boot_ind_df = do.call(rbind, boot_ind)
  return(boot_ind_df)
}

# updated function calculates the RD on all s elements in the s_grid for each bootstrap sample i 
## passed in as an argument for boot()
s_point_est <- function(i, mstate_bootlist,
                        s_grid, newd000, newd010, newd100, mstate_form,
                        a, astar,
                        exposure, mediator,
                        #exposure, mediator, outcome, basec, mediator_event, event, covs_df, a
                        trans, bh_method, ymreg="coxph"){
  
  mstate_df <- mstate_bootlist[[i]]
  if (ymreg=="coxph"){
    joint_mod <- coxph(mstate_form, data = mstate_df, method = bh_method)
    joint_mod_call <- getCall(joint_mod)
    joint_mod_call$data <- mstate_df
    joint_mod_call$formula <- mstate_form
    joint_mod2 <- eval.parent(joint_mod_call)
  }
  print("fitted model")
  
  # use msfit() to get predicted cumulative hazards data frames
  cumhaz000_msfit = msfit(joint_mod2, newd000, trans=trans) # for RD and SD
  cumhaz010_msfit = msfit(joint_mod2, newd010, trans=trans) # for RD and SD
  cumhaz100_msfit = msfit(joint_mod2, newd100, trans=trans) # for SD only
  # extract cumulative hazards from the msfit objects
  cumhaz000 = cumhaz000_msfit$Haz
  cumhaz010 = cumhaz010_msfit$Haz
  cumhaz100 = cumhaz100_msfit$Haz
  print(sum(cumhaz100 - cumhaz000))
  # extract transitions that are needed
  ## time var the same for all data frames below; all below dfs have the same # of rows
  cumhaz000_trans1 = subset(cumhaz000, trans==1)
  cumhaz000_trans2 = subset(cumhaz000, trans==2)
  cumhaz010_trans2 = subset(cumhaz010, trans==2)
  cumhaz100_trans1 = subset(cumhaz100, trans==1)
  # populate the hazard grid
  time_vec = cumhaz000_trans1$time
  n_grid = nrow(cumhaz000_trans1)
  haz000_trans1 = rep(NA, n_grid) # grids for hazard \alpha_{01}
  haz100_trans1 = rep(NA, n_grid)
  for (i in 1:(nrow(cumhaz000_trans1)-1)){
    haz000_trans1[i]=(cumhaz000_trans1$Haz[i+1]-cumhaz000_trans1$Haz[i])/(time_vec[i+1]-time_vec[i])
    haz100_trans1[i]=(cumhaz100_trans1$Haz[i+1]-cumhaz100_trans1$Haz[i])/(time_vec[i+1]-time_vec[i])
  }
  
  # get the index to slice the newd001 list
  max_s = max(s_grid)
  newd001_trans3_list = dynamic_newd(mstate_df, exposure, mediator, newd000, time_vec, max_s, a, trans) # active - a=1
  newd000_trans3_list = dynamic_newd(mstate_df, exposure, mediator, newd000, time_vec, max_s, astar, trans) # reference - a=0=astar
  # numerical integration
  ## use lapply to compute each integrand, corresponding to max_s
  ## for integral list for lower values of s, slice the list
  print("start creating integrand list")
  integrand_list = lapply(1:length(newd001_trans3_list), function(i){
    # create the msfit object, extract cumulative hazards data frame
    newd001_trans3 = newd001_trans3_list[[i]]
    newd000_trans3 = newd000_trans3_list[[i]]
    cumhaz001_trans3_msfit = msfit(joint_mod2, newd001_trans3, trans=trans)
    cumhaz000_trans3_msfit = msfit(joint_mod2, newd000_trans3, trans=trans)
    cumhaz001_trans3 = cumhaz001_trans3_msfit$Haz
    cumhaz000_trans3 = cumhaz000_trans3_msfit$Haz
    # get transition-specific cumulative hazards
    cumhaz001_trans3 = subset(cumhaz001_trans3, trans==3)
    cumhaz000_trans3 = subset(cumhaz000_trans3, trans==3)
    cumhaz001_trans3_s = cumhaz001_trans3 %>% arrange(abs(time-max_s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{12}(s|A=1, C=1)
    cumhaz000_trans3_s = cumhaz000_trans3 %>% arrange(abs(time-max_s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{12}(s|A=0, C=1)
    # get the row indices to evaluate the integrand over
    #time_ind = which(time_vec == newd001_trans3[3,"M.3"])
    time_ind = which(time_vec == newd001_trans3[3, paste(mediator, ".", 3, sep="")])
    #time_ind = min(time_ind, up_to_ind)
    if (time_ind < length(newd001_trans3_list)){
      # For RD
      P_01_integrand = (exp(-cumhaz000_trans1$Haz[time_ind]-cumhaz000_trans2$Haz[time_ind])*haz000_trans1[time_ind]*exp(-cumhaz000_trans3_s+cumhaz000_trans3$Haz[time_ind])) * (min(max_s,time_vec[time_ind+1])-time_vec[time_ind])
      P_g_01_integrand = (exp(-cumhaz000_trans1$Haz[time_ind]-cumhaz010_trans2$Haz[time_ind])*haz000_trans1[time_ind]*exp(-cumhaz001_trans3_s+cumhaz001_trans3$Haz[time_ind])) * (min(max_s,time_vec[time_ind+1])-time_vec[time_ind])
      # for SD
      sd_integrand1 = (exp(-cumhaz100_trans1$Haz[time_ind]-cumhaz010_trans2$Haz[time_ind]) * haz100_trans1[time_ind] * exp(-cumhaz001_trans3_s+cumhaz001_trans3$Haz[time_ind])) * (min(max_s,time_vec[time_ind+1])-time_vec[time_ind])
    }else{
      P_01_integrand = 0
      P_g_01_integrand = 0
      sd_integrand1 = 0
    }
    return(c(P_01 = P_01_integrand, P_g_01 = P_g_01_integrand, sd_integrand1 = sd_integrand1))
  })
  print("finished creating integrand list")
  
  ## sum up the individual integrands as the estimate for P_01
  RD_vec = rep(NA, length(s_grid))
  names(RD_vec) = paste("RD", s_grid, sep="")
  SD_vec = rep(NA, length(s_grid))
  names(SD_vec) = paste("SD", s_grid, sep="")
  for (i in 1:length(RD_vec)){ # loop through each value in s_grid
    s = s_grid[i]
    up_to_ind = which.max(abs(time_vec - s) == min(abs(time_vec - s)))
    sums = colSums(do.call(rbind, integrand_list[1:up_to_ind]))
    print(sums)
    P_01 = sums[1]
    P_g_01 = sums[2]
    sd_integral = sums[3]
    # estimate cumulative hazard up to the end time s
    cumhaz000_trans1_s = cumhaz000_trans1 %>% arrange(abs(time-s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{01}(s|A=0,C=1)
    cumhaz000_trans2_s = cumhaz000_trans2 %>% arrange(abs(time-s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{02}(s|A=0,C=1)
    cumhaz010_trans2_s = cumhaz010_trans2 %>% arrange(abs(time-s)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{02}(s|A=1,C=1)
    cumhaz100_trans1_s = cumhaz100_trans1 %>% arrange(abs(time-s)) %>% filter(row_number()==1) %>% pull(Haz)
    # P_00 and P_g_00
    P_00 = exp(-cumhaz000_trans1_s-cumhaz000_trans2_s)
    P_g_00 = exp(-cumhaz000_trans1_s-cumhaz010_trans2_s)
    P_10 = exp(-cumhaz100_trans1_s-cumhaz010_trans2_s)
    # compute RD
    RD_vec[i] = (P_g_00 + P_g_01) - (P_00 + P_01)
    SD_vec[i] = (P_10 + sd_integral) - (P_g_00 + P_g_01)
  }
  
  out_df = data.frame(RD = RD_vec, SD = SD_vec, TD = RD_vec+SD_vec)
  
  return(out_df)
}

