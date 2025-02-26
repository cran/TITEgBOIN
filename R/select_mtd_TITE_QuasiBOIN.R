#' @title select_mtd_TITE_QuasiBOIN
#'
#' @description Obtain the maximum tolerated dose (MTD) of Bayesian optimal interval (BOIN) (Yuan et al. 2016)/
#' Generalized Bayesian optimal interval (gBOIN) (Mu et al. 2019)/Time-to-event Bayesian optimal interval (TITEBOIN) (Lin et al. 2020)/
#' Time-to-event generalized Bayesian optimal interval (TITEgBOIN) (Takeda et al. 2022) designs.
#'
#' @param target The target toxicity probability (example: \code{target <- 0.30}) or the target normalized equivalent 
#' toxicity score (ETS) (example: \code{target <- 0.47 / 1.5}).
#' @param ntox Number of patients with dose limiting toxicity (DLT) or the sum of normalized equivalent toxicity score (ETS).
#' @param npts The number of patients enrolled at each dose level.
#' @param Neli The sample size cutoff for elimination. The default is \code{Neli=3}.
#' @param cutoff.eli The cutoff to eliminate an overly toxic dose for safety.
#' We recommend the default value of (\code{cutoff.eli=0.95}) for general use.
#' @param extrasafe Set \code{extrasafe=TRUE} to impose a more stringent stopping rule
#' @param offset A small positive number (between 0 and 0.5) to control how strict the
#' stopping rule is when \code{extrasafe=TRUE}. A larger value leads to a more
#' strict stopping rule. The default value \code{offset=0.05} generally works well.
#' @param print Print the additional result or not. The default value is \code{print=FALSE}.
#' @param gdesign For Bayesian optimal interval (BOIN) and Time-to-event Bayesian optimal interval (TITEBOIN), "FALSE" should be assigned.
#' For Generalized Bayesian optimal interval (gBOIN) and Time-to-event generalized Bayesian optimal interval (TITEgBOIN), "TRUE" should be assigned . 
#' The default is \code{gdesign=FALSE}.
#' @return \code{select_mtd_TITE_QuasiBOIN()} returns the selected dose.
#' 
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
#' #For Bayesian optimal interval (BOIN) design/Time-to-event bayesian optimal interval (TITEBOIN)
#' #design
#' target<-0.3
#' y<-c(0,0,1,2,3,0)
#' n<-c(3,3,6,9,9,0)
#' select_mtd_TITE_QuasiBOIN(target=target,ntox=y,npts=n,print=TRUE,gdesign=FALSE)
#' 
#' 
#' #For Generalized Bayesian optimal interval (gBOIN) design/Time-to-event generalized bayesian
#' #optimal interval (TITEgBOIN) design
#' target<-0.47/1.5
#' y<-c(0,0,2/1.5,3.5/1.5,5.5/1.5,0)
#' n<-c(3,3,6,9,9,0)
#' select_mtd_TITE_QuasiBOIN(target=target,ntox=y,npts=n,print=TRUE,gdesign=TRUE)
#' @importFrom stats pbeta qbeta rexp rmultinom runif var
#' @export


select_mtd_TITE_QuasiBOIN <- function(target,ntox, npts, Neli=3, cutoff.eli = 0.95, extrasafe = FALSE, offset = 0.05, print = FALSE,gdesign=FALSE) {
  ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1)
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol)))
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }
  ## determine whether the dose has been eliminated during the trial
  y = ntox
  n = npts
  ndose = length(n)
  elimi = rep(0, ndose)
  for (i in 1:ndose) {
    if (n[i] >= Neli) {
      if (1 - pbeta(target, y[i] + 1, n[i] - y[i] + 1) > cutoff.eli) {
        elimi[i:ndose] = 1
        break
      }
    }
  }

  if (extrasafe) {
    if (n[1] >= Neli) {
      if (1 - pbeta(target, y[1] + 1, n[1] - y[1] + 1) > cutoff.eli - offset) {
        elimi[1:ndose] = 1
      }
    }
  }

  ## no dose should be selected (i.e., selectdose=99) if the first dose is already very toxic or all uneliminated doses are never used to treat patients
  if (elimi[1] == 1 || sum(n[elimi == 0]) == 0) {
    selectdose = 99
  }
  else {
    adm.set = (n != 0) & (elimi == 0)
    adm.index = which(adm.set == T)
    y.adm = y[adm.set]
    n.adm = n[adm.set]

    ## poster mean and variance of toxicity probabilities using beta(0.005, 0.005) as the prior
    phat = (y.adm + 0.005)/(n.adm + 0.01)
    phat.var = (y.adm + 0.005) * (n.adm - y.adm + 0.005)/((n.adm + 0.01)^2 * (n.adm + 0.01 + 1))

    ## perform the isotonic transformation using PAVA
    phat = pava(phat, wt = 1/phat.var)
    phat = phat + (1:length(phat)) * 1e-10  ## break ties by adding an increasingly small number
    selectd = sort(abs(phat - target), index.return = T)$ix[1]  ## select dose closest to the target as the MTD
    selectdose = adm.index[selectd]
  }
  if (print == TRUE) {
    if (selectdose == 99) {
      message("All tested doses are overly toxic. No MTD is selected! \n")
    }
    else {
      message("The MTD is dose level ", selectdose, "\n\n")
    }
    trtd = (n != 0)
    poverdose = pava(1 - pbeta(target, y[trtd] + 0.05, n[trtd] -
                                 y[trtd] + 0.05))
    phat.all = pava((y[trtd] + 0.05)/(n[trtd] + 0.1), wt = 1/((y[trtd] +
                                                                 0.05) * (n[trtd] - y[trtd] + 0.05)/((n[trtd] + 0.1)^2 *
                                                                                                       (n[trtd] + 0.1 + 1))))

    if(gdesign==TRUE){

      message("Dose    Posterior normalized ETS estimate        95%                  \n",
          sep = "")
      message("Level     Estimate                       Credible Interval   Pr(toxicity>",
          target, "|data)\n", sep = "")
      for (i in 1:ndose) {
        if (n[i] > 0) {
          message(" ", i, "        ", formatC(phat.all[i],
                                          digits = 2, format = "f"), "                            (", formatC(qbeta(0.025,
                                                                                                               y[i] + 0.05, n[i] - y[i] + 0.05), digits = 2,
                                                                                                         format = "f"), ", ", formatC(qbeta(0.975, y[i] +
                                                                                                                                              0.05, n[i] - y[i] + 0.05), digits = 2, format = "f"),
              ")            ", formatC(poverdose[i], digits = 2,
                                       format = "f"), "\n")
        }

        else {
          message(" ", i, "        ", "----", "                           (",
              "------------", ")            ", "----", "\n")
        }
      }

    }else{
      message("Dose    Posterior DLT rate        95%                  \n",
          sep = "")
      message("Level     Estimate          Credible Interval   Pr(toxicity>",
          target, "|data)\n", sep = "")
      for (i in 1:ndose) {
        if (n[i] > 0) {
          message(" ", i, "        ", formatC(phat.all[i],
                                          digits = 2, format = "f"), "                (", formatC(qbeta(0.025,
                                                                                                  y[i] + 0.05, n[i] - y[i] + 0.05), digits = 2,
                                                                                            format = "f"), ", ", formatC(qbeta(0.975, y[i] +
                                                                                                                                 0.05, n[i] - y[i] + 0.05), digits = 2, format = "f"),
              ")            ", formatC(poverdose[i], digits = 2,
                                       format = "f"), "\n")
        }

        else {
          message(" ", i, "        ", "----", "               (",
              "------------", ")            ", "----", "\n")
        }
      }
    }
    message("NOTE: no estimate is provided for the doses at which no patient was treated.")
  }
  else {
    return(selectdose)
  }
}








