#' @title next_TITE_QuasiBOIN
#'
#' @description Determine the dose for the next cohort of new patients for single-agent trials using Bayesian optimal interval (BOIN)/
#' Generalized Bayesian optimal interval (gBOIN)/Time-to-event bayesian optimal interval (TITEBOIN)/Time-to-event generalized 
#' Bayesian optimal interval (TITEgBOIN) designs.
#'
#' @param target The target toxicity probability (example: \code{target <- 0.30}) or the target normalized equivalent 
#' toxicity score (ETS) (example: \code{target <- 0.47 / 1.5}).
#' @param n Number of patients treated at each dose level.
#' @param npend For Time-to-event Bayesian optimal interval (TITEBOIN)/Time-to-event generalized Bayesian optimal interval (TITEgBOIN), 
#' the number of pending patients at each dose level.For Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN), 
#' "NA" should be assigned.
#' @param y  Number of patients with dose limiting toxicity (DLT) or the sum of Normalized equivalent toxicity score (ETS).
#' @param ft For Time-to-event Bayesian optimal interval (TITEBOIN)/Time-to-event generalized Bayesian optimal interval (TITEgBOIN), 
#' Total follow-up time for pending patients for toxicity at each dose level (days). For Bayesian optimal interval (BOIN)/
#' Generalized Bayesian optimal interval (gBOIN), "NA" should be assigned.
#' @param d Current dose level.
#' @param maxt For Time-to-event Bayesian optimal interval (TITEBOIN)/Time-to-event generalized Bayesian optimal interval (TITEgBOIN),
#' length of assessment window for toxicity (days). For Bayesian optimal interval (BOIN)/Generalized Bayesian optimal interval (gBOIN), 
#' "NA" should be assigned.
#' @param p.saf The lower bound. The default value is \code{p.saf=0.6*target}.
#' @param p.tox The upper bound. The default value is \code{p.tox=1.4*target}.
#' @param elimination Elimination of each dose (0,1 should be assigned, 0 means the dose is not eliminated,
#' 1 means the dose is eliminated due to over toxic(\code{elimination=NA}, 0 is defaulted for each dose level)).
#' @param cutoff.eli The cutoff to eliminate an overly toxic dose for safety.
#' We recommend the default value of (\code{cutoff.eli=0.95}) for general use.
#' @param extrasafe Set \code{extrasafe=TRUE} to impose a more stringent stopping rule
#' @param offset A small positive number (between 0 and 0.5) to control how strict the
#' stopping rule is when \code{extrasafe=TRUE}. A larger value leads to a more
#' strict stopping rule. The default value \code{offset=0.05} generally works well.
#' @param n.earlystop The early stopping parameter. The default value is \code{n.earlystop=100}.
#' @param maxpen  For Time-to-event Bayesian optimal interval (TITEBOIN)/Time-to-event generalized Bayesian optimal interval 
#' (TITEgBOIN), the upper limit of the ratio of pending patients. For Bayesian optimal interval (BOIN)/Generalized Bayesian optimal
#' interval (gBOIN), "NA" should be assigned.
#' @param print_d Print the additional result or not. The default value is \code{print_d=FALSE}.
#' @param Neli The sample size cutoff for elimination. The default is \code{Neli=3}.
#' @param gdesign For Bayesian optimal interval (BOIN) and Time-to-event bayesian optimal interval (TITEBOIN), "FALSE" should be 
#' assigned. For Generalized Bayesian optimal interval (gBOIN) and Time-to-event generalized bayesian optimal interval (TITEgBOIN), 
#' "TRUE" should be assigned . The default is \code{gdesign=FALSE}.
#'
#' @return \code{next_TITE_QuasiBOIN()} returns the toxicity probability and the recommended dose level for the next cohort
#' including: (1) the lower Bayesian optimal boundary (\code{lambda_e})
#' (2) the upper Bayesian optimal boundary (\code{lambda_d})
#' (3) The number of patients or the effective sampe size (ESS) at each dose level (\code{ESS})
#' (4) The dose limiting toxicity (DLT) rate or mu (the estimated quasi-Bernoulli toxicity probability) at each dose level (\code{mu})
#' (5) the recommended dose level for the next cohort as a numeric value under (\code{d})
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
#' #For Bayesian optimal interval (BOIN) design
#' target<-0.3
#' next_TITE_QuasiBOIN(target=target,n=c(3,3,4,4,4,0),npend=NA, y=c(0,0,1,1,1,0), ft=NA,
#'                     d=5, maxt=NA,p.saf= 0.6 * target, p.tox = 1.4 * target,elimination=NA,
#'                     cutoff.eli = 0.95,extrasafe = FALSE, n.earlystop = 10,
#'                     maxpen=NA,print_d = TRUE,gdesign=FALSE)
#' 
#' 
#' #For Generalized Bayesian optimal interval (gBOIN) design
#' target=0.47/1.5
#' next_TITE_QuasiBOIN(target=target,n=c(3,3,4,4,4,0),npend=NA,
#'                     y=c(0, 0, 0.5/1.5, 1.0/1.5, 1.5/1.5, 0),ft=NA, d=5, maxt=NA,
#'                     p.saf= 0.6 * target, p.tox = 1.4 * target,elimination=NA,
#'                     cutoff.eli = 0.95,extrasafe = FALSE, n.earlystop = 10,
#'                     maxpen=NA,print_d = TRUE,gdesign=TRUE)
#' 
#' 
#' #For Time-to-event bayesian optimal interval (TITEBOIN) design
#' target=0.3
#' next_TITE_QuasiBOIN(target=target,n=c(3,3,4,4,4,0),npend=c(0,0,0,1,2,0), y=c(0,0,1,1,1,0),
#'                    ft=c(0, 0, 0, 14, 28, 0),d=5, maxt=28,p.saf= 0.6 * target,
#'                     p.tox = 1.4 * target,elimination=NA,cutoff.eli = 0.95,
#'                     extrasafe = FALSE, n.earlystop = 10,maxpen=0.5,print_d = TRUE,
#'                     gdesign=FALSE)
#' 
#' 
#' #For Time-to-event generalized bayesian optimal interval (TITEgBOIN) design
#' target=0.47/1.5
#' next_TITE_QuasiBOIN(target=target,n=c(3,3,4,4,4,0),npend=c(0,0,0,1,2,0),
#'                     y=c(0, 0, 0.5/1.5, 1.0/1.5, 1.5/1.5, 0),ft=c(0, 0, 0, 14, 28, 0),
#'                     d=5, maxt=28,p.saf= 0.6 * target, p.tox = 1.4 * target,
#'                     elimination=NA,cutoff.eli = 0.95,extrasafe = FALSE,
#'                     n.earlystop = 10,maxpen=0.5,print_d = TRUE,gdesign=TRUE)
#' @importFrom stats pbeta qbeta rexp rmultinom runif var
#' @export

next_TITE_QuasiBOIN <- function(target, n, npend, y, ft, d, maxt=28, p.saf = 0.6 * target, p.tox =  1.4 * target,elimination=NA,
                                cutoff.eli = 0.95, extrasafe = FALSE, offset=0.05,n.earlystop = 100,maxpen=0.50,Neli=3,print_d = FALSE,
                                gdesign=FALSE)
{


  earlystop <- 0
  stop <- 0
  elimineed <-0
  pending <- 0
  suspend<-0

  if(is.na(maxt)){maxt = 1}
  if(is.na(maxpen)){maxpen=0.5;}

  ndose <- length(y)
  if(is.na(npend[1])){npend = rep(0,ndose)}
  if(is.na(ft[1])){ft = rep(0,ndose)}

  if(is.na(elimination[1])){
    elimi = rep(0,ndose)
  }else{
    elimi=elimination
  }

  ntox.curr1<-y
  n.curr1<-n
  n.pend1<-npend
  totalt1<-ft/maxt

  #effective sampe size
  e_n<-(n.curr1-n.pend1+totalt1)
  #mu DLT rate;
  mu<-ntox.curr1/(n.curr1-n.pend1+totalt1)


  lambda1 = log((1 - p.saf)/(1 - target))/log(target *
                                                (1 - p.saf)/(p.saf * (1 - target)))
  lambda2 = log((1 - target)/(1 - p.tox))/log(p.tox * (1 -

                                                         target)/(target * (1 - p.tox)))

  ###dose elimination rule: check all the doses for elimination####
  for(dd in 1:ndose){
    #if (1-pbeta(target, ntox.curr1[dd]+1, (n.curr1[dd]-n.pend1[dd]+totalt1[dd])-ntox.curr1[dd]+1)>cutoff.eli && n.curr1[dd]>=Neli){
    if (1-pbeta(target, ntox.curr1[dd]+1, (n.curr1[dd])-ntox.curr1[dd]+1)>cutoff.eli && n.curr1[dd]>=Neli){
      elimi[dd:ndose]=1;
      break;
    }
  }

  ####current dose#####
  ntox.curr<-ntox.curr1[d]
  n.curr<-n.curr1[d]
  n.pend<-n.pend1[d]
  totalt<-totalt1[d]

  #check whether extra safey rule should be applied
  if(extrasafe)
  {
    if(d==1){
      #if(1-pbeta(target, ntox.curr+1, (n.curr-n.pend+totalt)-ntox.curr+1)>cutoff.eli-offset && n.curr>=Neli) {
      if(1-pbeta(target, ntox.curr+1, n.curr-ntox.curr+1)>cutoff.eli-offset && n.curr>=Neli) {
        earlystop = 1; elimi[1:ndose]=1;  }
    }
  }

  #############################################################################
  # Why use n instead of effective n for the safety/extra safety rule?        #
  # If there are many pending people when we check elimination                #
  # (it is possible since our logic checks elimination first and              #
  # then check whether we have enough non-pending people), using n            #
  # can reduce the chance of eliminating ??false positive?? overly-toxic dose   #
  # comparing to using effective n, which results in higher probability of    #
  # selecting MTD.                                                            #
  #############################################################################

  # check whether the current dose level should be eliminated
  if(elimi[d]==1) {
    d=which(elimi==1)[1]-1
    if(d==0){earlystop = 1;}
    elimineed=1
  }
  # else{#check whether the current dose is toxic based on observed data
  #   if((ntox.curr/n.curr)>=lambda2){
  #     if(d==1){d=d; if(n.pend>0){pending=1};} else{d=d-1};
  #     elimineed=1
  #   }
  # } # This part of code doesn't change results, delete it.

  #check whether current dose reaches the "adequate size" of n.earlystop
  if(n.curr>=n.earlystop){
    #earlystop = 1; #JZ 01/28/2023 edit: Change "earlystop=1" to "stop=1" since we will count earlystop as early termination%
    stop=1;
  }


  #####################################################################
  # Assign Dose:                                                      #
  # 1. First check if should early stop trial                         #
  # 2. Else, check if the current dose should be eliminated           #
  # 3. Else, check if should suspend enrollment                       #
  # 4. Else, make dose assignment by comparing mu with lambda_1,2     #
  #####################################################################
  stay<-0
  if(elimineed==1){

  }else{
    if(n.pend>n.curr*maxpen) {pending=1;suspend<-1}
    else
    {
      #compute the estimate of toxicity rate
      if(n.pend==0){
        phat=ntox.curr/n.curr
      } else {
        #compute the estimated toxicity rate based on the imputed data
        phat = ntox.curr/(n.curr-n.pend+totalt)
      }


      #make dose assignment decisions
      if(phat<=lambda1){
        #check whether the current dose is the highest
        if(d==ndose){d=d;stay<-1}
        else{
          if(elimi[d+1]==1){d=d; stay=1} else{ d=d+1}
        }}
      else if (phat>=lambda2 ){
        if(d==1){d=d;stay=1; if(n.pend>0){pending=1}} else{d=d-1}
      }
      else {d=d; stay=1}
    }
  }

  if(stop==1 & stay==1){
    stop=1
  }else{stop=0}

  if (print_d == FALSE) {
    ddose=list()
    ddose=list(d=d,pending=pending,n.pend=n.pend,earlystop=earlystop,stop=stop,elimineed=elimineed,elimi=elimi)
    return(ddose)
  }else{
    if (earlystop == 1){
      message("\n")
      message("Early stop because the 1st dose is too toxic\n")
    }else if(earlystop == 0){
      if(suspend==1){
        message("\n")
        message("More than 50% of patients have not finished the assessment at the current dose level, suspend the dose allocation\n")
      }else{
        if(sum(elimi)==0){
          ddose=list()
          ddose=list(lambda_e = lambda1,lambda_d = lambda2, e_n=e_n,mu=mu,d=d)
          names(ddose)[1]<-'lambda_e'
          names(ddose)[2]<-'lambda_d'
          if(gdesign=="TRUE"){
          names(ddose)[3]<-'The ESS (the effective sampe size) at each dose level'
          names(ddose)[4]<-'The estimated quasi-Bernoulli toxicity probability at each dose level'
          }else{
            names(ddose)[3]<-'The number of patients at each dose level'
            names(ddose)[4]<-'The DLT rate at each dose level'
          }
          names(ddose)[5]<-'Next dose level'
          print(ddose)
        }else{
          elimi_note<-paste0(paste0("Eliminate dose level [",min(which(elimi==1))),"] and all the doses above this dose level.")
          ddose=list()
          ddose=list(lambda_e = lambda1,lambda_d = lambda2, e_n=e_n,mu=mu,d=d,elimi=elimi_note)
          names(ddose)[1]<-'lambda_e'
          names(ddose)[2]<-'lambda_d'
          if(gdesign=="TRUE"){
          names(ddose)[3]<-'The ESS (the effective sampe size) at each dose level'
          names(ddose)[4]<-'The estimated quasi-Bernoulli toxicity probability at each dose level'
          }else{
            names(ddose)[3]<-'The number of patients at each dose level'
            names(ddose)[4]<-'The DLT rate at each dose level'
          }
          names(ddose)[5]<-'Next dose level'
          names(ddose)[6]<-'Dose elimination due to too toxic'
          print(ddose)
        }
      }
    }

  }
}








