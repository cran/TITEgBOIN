#' @title get_oc_TITE_QuasiBOIN
#'
#' @description Obtain the operating characteristics of the model-assisted design 
#' for single agent trials by simulating trials using Bayesian optimal interval (BOIN) (Yuan et al. 2016)/
#' Generalized Bayesian optimal interval (gBOIN) (Mu et al. 2019)/Time-to-event Bayesian optimal interval (TITEBOIN) 
#' (Lin et al. 2020)/Time-to-event generalized Bayesian optimal interval (TITEgBOIN) designs(Takeda et al. 2022).
#'
#' @param target The target toxicity probability (example: \code{target <- 0.30}) or the target normalized equivalent 
#' toxicity score (ETS) (example: \code{target <- 0.47 / 1.5}).
#' @param prob A vector (Bayesian optimal interval (BOIN) or Time-to-event Bayesian optimal interval (TITEBOIN) design)
#' /matrix (Generalized Bayesian optimal interval (gBOIN) or Time-to-event generalized Bayesian optimal interval (TITEgBOIN) design)
#' containing the true toxicity probabilities of the investigational dose levels.
#' @param score For Generalized Bayesian optimal interval (gBOIN)/Time-to-event generalized Bayesian optimal interval (TITEgBOIN), 
#' a vector containing the relative severity of different toxicity grades in terms of dose limiting toxicity (DLTs) in the 
#' dose-finding procedure. As default, toxicity grades of 0/1,2,3, and 4 are assigned values of 0,0.5,1,1.5. 
#' For Bayesian optimal interval (BOIN)/Time-to-event Bayesian optimal interval (TITEBOIN), "NA" should be assigned.
#' @param TITE For Time-to-event Bayesian optimal interval (TITEBOIN)/Time-to-event generalized Bayesian optimal interval (TITEgBOIN), 
#' "TRUE" should be assigned. For Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN), 
#' "FALSE" should be assigned.
#' @param ncohort The total number of cohorts.
#' @param cohortsize The cohort size.
#' @param maxt For Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval (TITEgBOIN), 
#' the maximum follow-up time. for Bayesian optimal interval (BOIN)/ Generalized Bayesian optimal interval (gBOIN), 
#' if you don't need to get 1, the average trial duration needed for the trial, 2, the standard deviation of average trial duration 
#' needed for the trial. Then "NA" should be assigned; If you need to get 1,the average trial duration needed for the trial, 2, 
#' the standard deviation of average trial duration needed for the trial. Then please specify the accrual rate and the maximum 
#' follow-up time.
#' @param accrual For Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval 
#' (TITEgBOIN), the accrual rate, i.e., the number of patients accrued in 1 unit of time, for Bayesian optimal interval (BOIN)/ 
#' Generalized Bayesian optimal interval (gBOIN), if you don't need to get 1,the average trial duration needed for the trial,
#' 2, the standard deviation of average trial duration needed for the trial. Then "NA" should be assigned;
#' if you need to get 1,the average trial duration needed for the trial, 2, the standard deviation of
#' average trial duration needed for the trial, Then please specify the accrual rate and the maximum follow-up time.
#' @param maxpen  For Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval 
#' (TITEgBOIN), the upper limit of the ratio of pending patients. For Bayesian optimal interval (BOIN)/Generalized Bayesian optimal 
#' interval (gBOIN), "NA" should be assigned.
#' @param alpha1  For Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval 
#' (TITEgBOIN), a number from (0,1) that assume toxicity outcomes occurred with probability alpha1 in the last fraction of alpha2 
#' of the assessment window. The default is \code{alpha1=0.5}. For Bayesian optimal interval (BOIN)/Generalized Bayesian optimal 
#' interval (gBOIN), "NA" should be assigned.
#' @param alpha2  For Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized bayesian optimal interval 
#' (TITEgBOIN), a number from (0,1) that assume toxicity outcomes occurred with probability alpha1 in the last fraction of 
#' alpha2 of the assessment window. The default is \code{alpha2=0.5}. For Bayesian optimal interval (BOIN)/Generalized Bayesian 
#' optimal interval (gBOIN), "NA" should be assigned.
#' @param n.earlystop The early stopping parameter and the decision is to stay. If the number of patients treated at the current 
#' dose reaches n.earlystop, stop the trial and select the maxinum tolerated dose (MTD) based on the observed data.
#' The default value \code{n.earlystop=100} essentially turns off this type of early stopping.
#' @param Neli The sample size cutoff for elimination. The default is \code{Neli=3}.
#' @param startdose The starting dose level for the trial.
#' @param p.saf The lower bound. The default value is \code{p.saf=0.6*target}.
#' @param p.tox The upper bound. The default value is \code{p.tox=1.4*target}.
#' @param cutoff.eli The cutoff to eliminate an overly toxic dose for safety. We recommend the default value of 
#' \code{cutoff.eli=0.95} for general use.
#' @param extrasafe Set \code{extrasafe=TRUE} to impose a more stringent stopping rule.
#' @param offset A small positive number (between 0 and 0.5) to control how strict the stopping rule is when 
#' \code{extrasafe=TRUE}. A larger value leads to a more strict stopping rule. The default value \code{offset=0.05} generally works well.
#' @param ntrial The total number of trials to be simulated.
#' @param seed The seed, The default value is \code{seed = 100}
#' @param titration set \code{titration=TRUE} to perform dose escalation with cohort size = 1 
#' to accelerate dose escalation at the beginning of the trial. The default value \code{titration=FALSE}.
#' @param cap.titration cap the \code{titration} up to dose level, set \code{cap.titration=3} 
#' to cap the \code{titration} up to dose level 3 with cohort size = 1. The default value \code{cap.titration=0}.
#' @details This function generates he operating characteristics of the Bayesian optimal interval (BOIN)/ 
#' Generalized Bayesian optimal interval (gBOIN)/Time-to-event Bayesian optimal interval (TITEBOIN)/ Time-to-event 
#' generalized Bayesian optimal interval (TITEgBOIN) designs for trials by simulating trials under the prespecified 
#' true toxicity probabilities of the investigational doses.
#' @return \code{get_oc_TITE_QuasiBOIN()} returns the operating characteristics of the Bayesian optimal interval (BOIN)/
#' Generalized Bayesian optimal interval (gBOIN)/Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized 
#' Bayesian optimal interval (TITEgBOIN) designs as a data frame,
#' including: (1) the percentage of trials that the maximum tolerated dose (MTD) is correctly selected,
#' (2) the percentage of patients that are correctly allocated to the maximum tolerated dose (MTD),
#' (3) the percentage of overdosing selection,
#' (4) the percentage of overdosing allocation,
#' (5) selection percentage at each dose level,
#' (6) the number of patients treated at each dose level,
#' (7) the percentage of patients treated at each dose level,
#' (8) the number of toxicities observed at each dose level,
#' (9) the average number of toxicities,
#' (10) the average number of patients,
#' (11) the percentage of early stopping without selecting the maxinum tolerated dose (MTD),
#' (12) the average trial duration needed for the trial,
#' (13) the standard deviation of average trial duration needed for the trial,
#' (14) simulation set up data frame, include the target toxicity probability/the normalized target equivalent toxicity score (ETS); the true target toxicity probability/
#' the true normalized equivalent toxicity score (ETS) at each dose level based on prob and score, and lambda_e denotes the lower Bayesian optimal boundary
#' and lambda_d denotes the upper Bayesian optimal boundary.
#' @note We should avoid setting the values of p.saf and p.tox very close to the target.
#' This is because the small sample sizes of typical phase I trials prevent us from
#' differentiating the target toxicity rate from the rates close to it. In addition,
#' in most clinical applications, the target toxicity rate is often a rough guess,
#' and finding a dose level with a toxicity rate reasonably close to the target rate
#' will still be of interest to the investigator. In addition, we recommend setting
#' the value of priortox relatively small, for example, priortox=target/2 to accelerate
#' the escalation procedure.
#' @references
#'1. Liu S. and Yuan, Y. (2015). Bayesian optimal interval designs for phase I clinical trials, Journal of the Royal Statistical Society: Series C , 64, 507-523.
#'2. Yuan, Y., Hess, K. R., Hilsenbeck, S. G., & Gilbert, M. R. (2016). Bayesian optimal interval design: a simple and well-performing design for phase I oncology trials. Clinical Cancer Research, 22(17), 4291-4301.
#'3. Zhou, H., Yuan, Y., & Nie, L. (2018). Accuracy, safety, and reliability of novel phase I trial designs. Clinical Cancer Research, 24(18), 4357-4364.
#'4. Zhou, Y., Lin, R., Kuo, Y. W., Lee, J. J., & Yuan, Y. (2021). BOIN Suite: A Software Platform to Design and Implement Novel Early-Phase Clinical Trials. JCO Clinical Cancer Informatics, 5, 91-101.
#'5. Takeda K, Xia Q, Liu S, Rong A. TITE-gBOIN: Time-to-event Bayesian optimal interval design to accelerate dose-finding accounting for toxicity grades. Pharm Stat. 2022 Mar;21(2):496-506. doi: 10.1002/pst.2182. Epub 2021 Dec 3. PMID: 34862715.
#'6. Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K.E. (2018). Time-to-event Bayesian Optimal Interval Design to Accelerate Phase I Trials. Clinical Cancer Research, 24(20): 4921-4930.
#'7. Rongji Mu, Ying Yuan, Jin Xu, Sumithra J. Mandrekar, Jun Yin, gBOIN: A Unified Model-Assisted Phase I Trial Design Accounting for Toxicity Grades, and Binary or Continuous End Points, Journal of the Royal Statistical Society Series C: Applied Statistics, Volume 68, Issue 2, February 2019, Pages 289–308, https://doi.org/10.1111/rssc.12263.
#'8. Lin R, Yuan Y. Time-to-event model-assisted designs for dose-finding trials with delayed toxicity. Biostatistics. 2020 Oct 1;21(4):807-824. doi: 10.1093/biostatistics/kxz007. PMID: 30984972; PMCID: PMC8559898.
#'9. Hsu C, Pan H, Mu R (2022). _UnifiedDoseFinding: Dose-Finding Methods for Non-Binary Outcomes_. R package version 0.1.9, <https://CRAN.R-project.org/package=UnifiedDoseFinding>.
#' @examples
#'
#' #For Bayesian optimal interval (BOIN) design and Output trial duration as an operating
#' #characteristics
#'get_oc_TITE_QuasiBOIN(target=0.3, score=NA,prob=c(0.25,0.30,0.45,0.49,0.53), TITE=FALSE,
#'                      ncohort=10, cohortsize=3,startdose=1,maxt=28,accrual=10,
#'                      maxpen=NA,alpha1=NA,alpha2=NA,cutoff.eli=0.95, ntrial=10,seed=6)
#'
#'
#' #For Bayesian optimal interval (BOIN) design and not Output trial duration as an operating
#' #characteristics
#' get_oc_TITE_QuasiBOIN(target=0.3, score=NA,prob=c(0.25,0.30,0.45,0.49,0.53), TITE=FALSE,
#'                      ncohort=10, cohortsize=3,startdose=1,maxt=NA,accrual=NA,
#'                      maxpen=NA,alpha1=NA,alpha2=NA,cutoff.eli=0.95, ntrial=10,seed=6)
#'
#'
#' #For Generalized Bayesian optimal interval (gBOIN) design and Output trial duration as an
#' #operating characteristics
#' target<-0.47/1.5
#' prob <- matrix(c(0.83,	0.75,	0.62,	0.51,	0.34,	0.19,
#'                 0.12,	0.15,	0.18,	0.19,	0.16,	0.11,
#'                 0.04,	0.07,	0.11,	0.14,	0.15,	0.11,
#'                 0.01,	0.03,	0.09,	0.16,	0.35,	0.59), ncol = 6, byrow = TRUE)
#' get_oc_TITE_QuasiBOIN(target=target, score=c(0,0.5,1,1.5),prob=prob, TITE=FALSE,ncohort=10,
#'                      cohortsize=3,startdose=1,maxt=28,accrual=10, maxpen=NA,alpha1=NA,
#'                      alpha2=NA,cutoff.eli=0.95, ntrial=10,seed=6)
#'
#'
#' #For Generalized Bayesian optimal interval (gBOIN) design and not Output trial duration as
#' #an operating characteristics
#' target<-0.47/1.5
#' prob <- matrix(c(0.83,	0.75,	0.62,	0.51,	0.34,	0.19,
#'                 0.12,	0.15,	0.18,	0.19,	0.16,	0.11,
#'                 0.04,	0.07,	0.11,	0.14,	0.15,	0.11,
#'                 0.01,	0.03,	0.09,	0.16,	0.35,	0.59), ncol = 6, byrow = TRUE)
#' get_oc_TITE_QuasiBOIN(target=target, score=c(0,0.5,1,1.5),prob=prob, TITE=FALSE,ncohort=10,
#'                      cohortsize=3,startdose=1,maxt=NA,accrual=NA, maxpen=NA,alpha1=NA,
#'                      alpha2=NA,cutoff.eli=0.95, ntrial=10,seed=6)
#'
#'
#' #For Time-to-event bayesian optimal interval (TITEBOIN) design
#' get_oc_TITE_QuasiBOIN(target=0.3, score=NA,prob=c(0.25,0.30,0.45,0.49,0.53), TITE=TRUE,
#'                      ncohort=10, cohortsize=3,startdose=1,maxt=28,accrual=10,
#'                      maxpen=0.5,alpha1=0.5,alpha2=0.5,cutoff.eli=0.95,
#'                      ntrial=10,seed=6)
#'
#'
#' #For Time-to-event generalized bayesian optimal interval (TITEgBOIN) design
#' target<-0.47/1.5
#' prob <- matrix(c(0.83,	0.75,	0.62,	0.51,	0.34,	0.19,
#'                 0.12,	0.15,	0.18,	0.19,	0.16,	0.11,
#'                 0.04,	0.07,	0.11,	0.14,	0.15,	0.11,
#'                 0.01,	0.03,	0.09,	0.16,	0.35,	0.59), ncol = 6, byrow = TRUE)
#' get_oc_TITE_QuasiBOIN(target=target, score=c(0,0.5,1,1.5),prob=prob, TITE=TRUE,ncohort=10,
#'                      cohortsize=3,startdose=1,maxt=28,accrual=10, maxpen=0.5,alpha1=0.5,
#'                      alpha2=0.5,cutoff.eli=0.95, ntrial=10,seed=6)
#' @importFrom stats pbeta qbeta rexp rmultinom runif var
#' @export


get_oc_TITE_QuasiBOIN <- function(target, prob, score=c(0,0.5,1.0,1.5), TITE=TRUE, ncohort, cohortsize, maxt=1, accrual=3, maxpen=0.5,
                                  alpha1=0.5,alpha2=0.5,n.earlystop = 100, Neli=3,startdose = 1,
                                  p.saf = 0.6 * target, p.tox =  1.4 * target, cutoff.eli = 0.95,
                                  extrasafe = FALSE, offset = 0.05, ntrial = 1000, seed=100,titration=FALSE,cap.titration=0)
{

  ############################################
  # Function -- get event time and toxicity  #
  ############################################
  gen.tite<-function(n, prob, alpha1,alpha2, Tobs=1)
  {
    ############ subroutines ############
    weib<-function(n, pi, pihalft)
    {
      ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
      alpha = log(log(1-pi)/log(1-pihalft))/log(1/(1-alpha2));
      lambda = -log(1-pi)/(Tobs^alpha);
      t = (-log(runif(n))/lambda)^(1/alpha);
      return(t);
    }

    ############ end of subroutines ############

    tox = rep(0, n);
    grade = rep(0, n);
    t.tox = rep(0, n);
    if (length(prob)==1){
      pi<-prob
    }else {
      pi <- sum(prob[c(2,3,4)])
    }

    pihalft = (1-alpha1)*pi;  # alpha*100% event in (0, 1/2T)
    t.tox = weib(n, pi, pihalft);
    tox[t.tox<=Tobs]=1;
    t.tox[tox==0]=Tobs;
    ##simulation for grade information*****;
    for (i in 1:n){
      if (tox[i]==0){
        grade[i] <- 0
      }else {
        if (length(prob)==1){
          grade[i]<-1
        }else {
          grade[i] <- ((t(rmultinom(1, 1, prob = c(prob[2],prob[3],prob[4]))))%*%c(score[2],score[3],score[4]))/max(score)
        }
      }
    }
    y = sum(grade);

    return(list(tox=tox, t.tox=t.tox, y=y,grade=grade));
  }


  ### simple error checking
  if(target<0.05) {stop("Error: the target is too low! \n"); return();}
  if(target>0.6)  {stop("Error: the target is too high! \n"); return();}
  if((target-p.saf)<(0.1*target)) {stop("Error: the probability deemed safe cannot be higher than or too close to the target! \n"); return();}
  if((p.tox-target)<(0.1*target)) {stop("Error: the probability deemed toxic cannot be lower than or too close to the target! \n"); return();}
  if(offset>=0.5) {stop("Error: the offset is too large! \n"); return();}
  if(n.earlystop<6) {stop("Warning: the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18 \n"); return();}
  if(is.na(maxpen)){maxpen=0.5;}
  if(maxpen<0 || maxpen>0.65) {stop("Error: the value of maxpen should lie within (0,0.65]!  \n"); return();}


  ################################
  # Set parameters -- Calculate  #
  ################################
  set.seed(seed);
  prior.p = rep(1/3,3)
  prior.p = prior.p/sum(prior.p)
  need_duration<-0
  if((!is.na(maxt)) & (!is.na(accrual)) & (TITE==FALSE)){
    need_duration<-1
  }
  if(is.na(score[1])){score = 1}
  if(is.na(maxt)){maxt = 1}
  if(is.na(accrual)){accrual = 1}
  if(is.na(alpha1)){alpha1 = 0.5}
  if(is.na(alpha2)){alpha2 = 0.5}

  #ndose = ncol(prob);
  #p.true <- apply(prob*score_truth,2,sum)/max(score_truth)
  if(is.vector(prob)){
    ndose=length(prob);
    p.true=prob;
    gdesign=0;
  }else{
    ndose = ncol(prob);
    p.true <- apply(prob*score,2,sum)/max(score)
    gdesign=1;
  }

  npts = ncohort*cohortsize
  lambda1 = log((1 - p.saf)/(1 - target))/log(target *
                                                (1 - p.saf)/(p.saf * (1 - target)))
  lambda2 = log((1 - target)/(1 - p.tox))/log(p.tox * (1 -
                                                         target)/(target * (1 - p.tox)))
  # get true MTD
  MTD.true=which(abs(p.true-target)==min(abs(p.true-target))) # when calculating MTD.true, use score and target

  ###############
  # Simulation  #
  ###############
  Y = matrix(rep(0, ndose * ntrial), ncol = ndose);
  N = matrix(rep(0, ndose * ntrial), ncol = ndose);
  dselect = rep(0, ntrial);
  durationV = rep(0, ntrial);
  npendV = rep(0, ntrial);
  a<-b<-1
  elimi.at = rep(0, ntrial)

  for(trial in 1:ntrial)
  {
    ## trial level
    n.d = rep(0, ndose);  # number of toxicity at each dose
    y.d = rep(0, ndose);  # number of patient at each dose
    earlystop = 0; #indicate if trial stops early
    elimi = rep(0, ndose)
    ## patient level
    y=NULL;  #toxicity indicator for each subject
    dv=NULL;  #dose for each subject
    t.enter=NULL; # time enter the study
    t.event=NULL; # time to event
    ## decision point level
    t.decision = 0; # decision making time
    d = startdose;  # current dose level
    npend = 0;

    fillup=0;  ####newly added######
    t.curr=0;  ####newly added######
    
    ####newly added######
    if(titration==TRUE){
      n.titr=min(cap.titration, length(p.true))
      
      for(i in 1:n.titr){
        
        t.enter = c(t.enter, t.curr); 
        
        if(is.vector(prob)){
          obscohort = gen.tite(1, prob[d], alpha1=alpha1,alpha2=alpha2,Tobs=maxt);
        }else {
          obscohort = gen.tite(1, prob[,d], alpha1=alpha1,alpha2=alpha2,Tobs=maxt);
        }
        
        t.event = c(t.event, obscohort$t.tox);
        
        y = c(y, obscohort$tox);
        dv = c(dv, rep(d, 1));
        
        t.decision = t.enter[length(t.enter)]+maxt;
        
        flag=0;
        
        t.bench=min(t.decision,t.event[length(t.enter)])
        
        while(flag==0){
          flag=1
          
          t.curr = t.curr + rexp(1, rate=accrual)
          
          if(t.bench>t.curr){
            
            flag=0
          }
          
        }
        
        cset = (dv==d);
        ntox.curr = sum(y[cset]);
        n.curr = sum(cset);
        
        ## judge end titration or not;
        if(sum(ntox.curr)>0){
          fillup=1
          break
          
        }else if(d==cap.titration & sum(ntox.curr)==0 & cap.titration==ndose){ ## no DLT in titration phase in this way
          fillup=1
          
          # do not change dose.curr
          break
        }else if(d==cap.titration & sum(ntox.curr)==0 & cap.titration<ndose){
          fillup=0
          d=d+1
          break
          
        }
        
        d=d+1
        
      }
    }

    ####end######;
    
    for(i in 1:ncohort)
    {
      ####newly added######
      if(i !=1){
        fillup=0
      }
      ####end#####;
      
      if(fillup==0){
        
        remains=cohortsize
        
        if((length(dv)+cohortsize)>npts){ ## do not have a full cohort remains;
          remains=npts-(length(dv))
        }
        
      # generate data for the new patient
      for(j in 1:remains)
      {
        if((TITE==TRUE) | (TITE==FALSE & need_duration==1 )){
          if(j==1) { 
            if(i==1){t.enter = c(t.enter, t.curr);
            }else{
              t.enter = c(t.enter, t.decision)
            }
          }else {
            t.enter = c(t.enter, t.enter[length(t.enter)] + rexp(1, rate=accrual))
          }
        }else{
          if(j==1) { 
            if(i==1){t.enter = c(t.enter, t.curr);
            }else{
              t.enter = c(t.enter, t.decision)
            }
          }else {
            t.enter = c(t.enter, t.enter[length(t.enter)])
          }
        }
      }
      if(is.vector(prob)){
        obscohort = gen.tite(remains, prob[d], alpha1=alpha1,alpha2=alpha2,Tobs=maxt);
      }else {
        obscohort = gen.tite(remains, prob[,d], alpha1=alpha1,alpha2=alpha2,Tobs=maxt);
      }

      t.event = c(t.event, obscohort$t.tox);
      y = c(y, obscohort$grade);
      dv = c(dv, rep(d, remains));
      t.decision = t.enter[length(t.enter)];
      }else{
        fillup=0
        # generate data for the new patient
        for(j in 2:cohortsize)
        {
          if((TITE==TRUE) | (TITE==FALSE & need_duration==1 )){
            if(j==2) { 
              if(i==1){t.enter = c(t.enter, t.curr);
              }else{
                t.enter = c(t.enter, t.decision)
              }
            }else {
              t.enter = c(t.enter, t.enter[length(t.enter)] + rexp(1, rate=accrual))
            }
          }else{
            if(j==2) { 
              if(i==1){t.enter = c(t.enter, t.curr);
              }else{
                t.enter = c(t.enter, t.decision)
              }
            }else {
              t.enter = c(t.enter, t.enter[length(t.enter)])
            }
          }
        }
        if(is.vector(prob)){
          obscohort = gen.tite(as.numeric(cohortsize)-1, prob[d], alpha1=alpha1,alpha2=alpha2,Tobs=maxt);
        }else {
          obscohort = gen.tite(as.numeric(cohortsize)-1, prob[,d], alpha1=alpha1,alpha2=alpha2,Tobs=maxt);
        }
        
        t.event = c(t.event, obscohort$t.tox);
        y = c(y, obscohort$grade);
        dv = c(dv, rep(d, as.numeric(cohortsize)-1));
        t.decision = t.enter[length(t.enter)];
        
      }
      
      nobs=-1; pending=1;
      d.curr=d;
      npend = npend-1;
      while(pending==1)
      {
        ######################################################
        # pending=1 means suspend decision for a while.      #
        # As long as there is a need to suspend, wait for    #
        # another patient to come and make decision again    #
        ######################################################
        npend = npend+1;
        pending = 0;
        if(i==ncohort | (TITE==FALSE & need_duration==1)) { t.decision = t.decision + min(max(t.event[dv==d]),maxt);
        }else if(i==ncohort | (TITE==FALSE & need_duration!=1)) {
          # If we don't do TITE version, the next decision time and enroll time will be after the last patient (of the current cohort) finish enrollment
          t.decision = t.decision + max(t.event[dv==d]);
        }else {
          t.decision = t.decision + rexp(1, rate=accrual)
        }
        ######################################################
        # t.decision is both the time to decide dose         #
        # for the next cohort, and also the time to enroll   #
        # the next cohort (1st subject).                     #
        ######################################################

        # determine which observation are observed
        delta = ((t.enter+t.event)<=t.decision); # delta=1: non-pending
        t = pmin(t.event, t.decision-t.enter, maxt);   # follow-up time

        # get dose level summary (number of patients, total toxicity, total pending time)
        ntox.curr1<-rep(0,ndose) # total toxicity
        n.curr1<-rep(0,ndose) # number of patients
        n.pend1<-rep(0,ndose) # number of pending patients
        totalt1<-rep(0,ndose) # total pending time

        for(dd in 1:ndose){
          cset1 = dv==dd;
          delta.curr1 = delta[cset1];
          ntox.curr1[dd] = sum((y[cset1])[delta.curr1==1]);
          n.curr1[dd] = sum(cset1);
          n.pend1[dd] = sum(delta[cset1]==0);
          t.curr1 = t[cset1];
          if(n.pend1[dd]>0){
            totalt11 = t.curr1[delta.curr1==0]
            totalt11 = 3*prior.p[1]*totalt11*(totalt11<=maxt/3)+
              ((prior.p[1]-prior.p[2])*maxt+3*prior.p[2]*totalt11)*(maxt/3<totalt11 & totalt11<=2*maxt/3)+
              ((prior.p[1]+prior.p[2]-2*prior.p[3])*maxt+3*prior.p[3]*totalt11)*(2*maxt/3<totalt11 & totalt11<=maxt)
            totalt1[dd] = sum(totalt11)
          }
        }

        ###dose assignment####
        out1=next_TITE_QuasiBOIN(target=target,n=n.curr1,npend=n.pend1, y=ntox.curr1, ft=totalt1, d=d.curr, maxt=maxt, p.saf = p.saf, p.tox = p.tox,elimination=elimi,
                                 cutoff.eli = cutoff.eli, extrasafe = extrasafe, offset=offset,n.earlystop = n.earlystop,maxpen=maxpen,Neli=Neli,print_d = FALSE)

        pending<-out1$pending
        n.pend<-out1$n.pend
        elimi<-out1$elimi
        d<-out1$d
        earlystop<-out1$earlystop
        stop <- out1$stop
        elimineed<-out1$elimineed
        if(earlystop==1|stop==1){ break;}
        if((length(dv))==npts){ break;}
        if(elimineed==1){next;}
        ###end####

      }
      if(earlystop==1|stop==1){t.decision = max(t.enter)+maxt; break;}
      if((length(dv))==npts){t.decision = max(t.enter)+maxt; break;}

    }

    for(k in 1:ndose){
      y.d[k] = sum(y[dv==k]);
      n.d[k] = sum(dv==k);
    }

    npendV[trial]= npend;
    Y[trial, ] = y.d
    N[trial, ] = n.d
    durationV[trial] = t.decision
    if (earlystop == 1) {
      dselect[trial] = 99
    }
    else {dselect[trial] = select_mtd_TITE_QuasiBOIN(target=target,ntox=y.d, npts=n.d,cutoff.eli=cutoff.eli, extrasafe=extrasafe, offset=offset)}


    if (sum(elimi==1)==0){
      elimi.at[trial] = NA     # check if there is any elimination in this trial
    }else {elimi.at[trial] = min(which(elimi==1))}
  }


  # Outputs
  selpercent=rep(0, ndose); # percent of MTD selection
  for(i in 1:ndose) {
    selpercent[i]=sum(dselect==i)/ntrial*100;
  }
  nptsdose = apply(N,2,mean); # number of selecting each dose
  nptsdosepct=nptsdose*100/sum(nptsdose)
  PCS=selpercent[MTD.true]
  PCA=nptsdosepct[MTD.true]
  POS= (sum(dselect>MTD.true)-sum(dselect==99))/ntrial*100 # percent of overdose selection
  if (MTD.true==ndose){POA=0}else{
    POA= sum(N[,(MTD.true+1):ndose])/sum(N)*100 # percent of overdose allocation
  }
  ntoxdose = apply(Y,2,mean);

  simu.setup <- data.frame(dose=as.character(1:ndose),target=as.character(target), p.true=as.character(round(p.true,4)), ncohort=as.character(ncohort), cohortsize = as.character(cohortsize),
                           startdose = as.character(startdose),lambda_e=as.character(round(lambda1,4)), lambda_d=as.character(round(lambda2,4)), cutoff.eli = as.character(cutoff.eli),
                           extrasafe = as.character(extrasafe), offset = as.character(offset),
                           ntrial = as.character(ntrial) )
if(gdesign==0){
  colnames(simu.setup)<-c("dose","the target toxicity probability","     the true target toxicity probability",
                          " ncohort"," cohortsize"," startdose","lambda_e","lambda_d",
                          " cutoff.eli"," extrasafe"," offset"," ntrial")
}else{
  colnames(simu.setup)<-c("dose","the normalized target ETS","     the true normalized ETS",
                          " ncohort"," cohortsize"," startdose","lambda_e","lambda_d",
                          " cutoff.eli"," extrasafe"," offset"," ntrial")
}


  if(TITE==FALSE & need_duration!=1){
    out=list(PCS=round(PCS,1), PCA=round(PCA,1), POS=round(POS,1), POA=round(POA,1),
             selpercent=selpercent, npatients=nptsdose, nptsdosepct=round(nptsdosepct,4),
             ntox=round(ntoxdose,4), totaltox=round(sum(Y)/ntrial,4), totaln=sum(N)/ntrial,
             percentstop=sum(dselect== 99)/ntrial*100, simu.setup=simu.setup);

    return(out);

  }else{
    out=list(PCS=round(PCS,1), PCA=round(PCA,1), POS=round(POS,1), POA=round(POA,1),
             selpercent=selpercent, npatients=nptsdose, nptsdosepct=round(nptsdosepct,4),
             ntox=round(ntoxdose,4), totaltox=round(sum(Y)/ntrial,4), totaln=sum(N)/ntrial,
             percentstop=sum(dselect== 99)/ntrial*100, duration=round(mean(durationV),4),sdduration=sqrt(var(durationV)),simu.setup=simu.setup);

    return(out);

  }

}



