#'Predict cases using generalized SEIR model
#'
#'This function predicts cases of an outbreak using a generalized SEIR model
#'
#'
#' @param country name of the country. It should be a character string.
#' @param start a start date in mm/dd/yy format. Start date can not be earlier than 01/22/20. Start date can not be later than finish date. If start date is \code{NULL} then start date will be 01/22/20.
#' @param finish a finish date in mm/dd/yy format. Finish date can not be earlier than start date. If finish date is \code{NULL} then finish date will be the latest date at John-Hopkins CSSE system.
#' @param Npop total population of the country
#' @param guess initial guess parameters
#' @param dt the time step. This oversamples time to ensure that the algorithm converges
#' @param f number of days for future predictions
#' @param boot if \code{TRUE} bootstrap will be performed to calculate confidence interval
#' @param dt the time step. This oversamples time to ensure that the algorithm converges
#' @param conf confidence level, default is 0.95.
#' @param seed set a seed for reproducible results.
#' @param repeatNumber number of iteration for bootstrap.
#' @param bootSample number of sample for each bootstrap. if \code{NULL} then the number of sample is 80 percent of the original data.
#' @param type a condidence interval type. If \code{"norm"} it calculates based on normal approximation, if \code{"perc"} it calculates based on percentile approximation,
#'
#'@importFrom stats qnorm sd
#'
#'@export predict_SEIQRDP
#'
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @return a list of predicted and actual cases.
#'
#' @examples
#'\donttest{
#'alpha_guess = 0.45
#'beta_guess = 1
#'LT_guess = 2
#'Q_guess = 0.55
#'lambda_guess = c(0.01,0.01,30)
#'kappa_guess = c(0.01,0.001,30)
#'
#'guess = list(alpha_guess,
#'             beta_guess,
#'             1/LT_guess,
#'             Q_guess,
#'             lambda_guess[1],
#'             lambda_guess[2],
#'             lambda_guess[3],
#'             kappa_guess[1],
#'             kappa_guess[2],
#'             kappa_guess[3])
#'
#'
#'pred = predict_SEIQRDP(country = "Germany", start = "10/15/20", finish = "12/15/20",
#'dt = 1, f = 30, conf = 0.95, Npop = 80000000, guess, boot = FALSE,
#'seed = 123, repeatNumber = 100, bootSample = NULL, type = "norm")
#'
#'predict = pred$pred
#'actual = pred$actual
#'}
#' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
#'
#' @references Peng, L., Yang, W., Zhang, D., Zhuge, C., Hong, L. 2020. “Epidemic analysis of COVID-19 in China by dynamical modeling”, arXiv preprint arXiv:2002.06563.
#' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}

predict_SEIQRDP <- function(country, start, finish, Npop = NULL, guess, dt = 1, f = 0, boot = FALSE, conf = 0.95,
                            seed = 123, repeatNumber = 200, bootSample = NULL, type = "norm"){


  dt = dt

  covidData = getDataCOVID(start = start, finish = finish, country = country)
  Recovered = covidData$tableRecovered
  Deaths = covidData$tableDeaths
  Confirmed = covidData$tableConfirmed

  if(nrow(Recovered) == 1){
    name = Recovered$CountryRegion
  }else{
    name = paste0(Recovered$ProvinceState, " (",Recovered$CountryRegion,")")
  }

  recovered = Recovered[ ,5:ncol(covidData$tableRecovered)]
  deaths = Deaths[ ,5:ncol(covidData$tableDeaths)]
  confirmed = Confirmed[ ,5:ncol(covidData$tableConfirmed)]


  Q0 = confirmed[1]-recovered[1]-deaths[1]
  I0 = 0.35*Q0
  E0 = 0.45*Q0
  R0 = recovered[1]
  D0 = deaths[1]

  Active = confirmed-recovered-deaths
  Active[Active<0] <- 0

  Q=Active
  R=recovered
  D = deaths

  time = seq(as.Date(start, format = "%m/%d/%y"), as.Date(finish, format = "%m/%d/%y"), by = "1 day")


  params = fit_SEIQRDP(Q = Active, R = recovered, D = deaths, Npop = Npop, E0 = E0, I0 = I0,
                       time = time, dt = dt, guess = guess, trace = FALSE, ci = FALSE)


  res = SEIQRDP(alpha = params$alpha1, beta = params$beta1,
                gamma = params$gamma1, delta = params$delta1,
                lambda0 = c(params$lambda01, params$lambda02, params$lambda03),
                kappa0 = c(params$kappa01, params$kappa02, params$kappa03),
                Npop, E0, I0, Q0, R0, D0,lambdaFun = params$lambdaFun,
                kappaFun = params$kappaFun, tstart = start, tfinish = finish,
                dt = dt, f =f)

  pred = cbind.data.frame(recovered_pred = res$recovered, quarantined_pred = res$quarantined,
                          dead_pred = res$dead,  susceptible_pred = res$susceptible, exposed_pred = res$exposed,
                          infectious_pred = res$infectious, insusceptible_pred = res$insusceptible, simTime = res$simTime)


  if(boot){

    if(is.null(bootSample)){

      bootSample = ceiling(ncol(recovered)*0.8)

    }
    set.seed(seed)

    bootQuarantined = bootSusceptible = bootExposed =
      bootInfectious = bootRecovered = bootDead = bootInsusceptible = matrix(NA,ncol(confirmed)+f, repeatNumber)
    colnames(bootQuarantined) = colnames(bootSusceptible) = colnames(bootExposed) =
      colnames(bootInfectious) = colnames(bootRecovered) = colnames(bootDead) =
      colnames(bootInsusceptible) = paste0("Boot",1:repeatNumber)

    for(b in 1:repeatNumber){

      s = sample(1:length(recovered),bootSample)
      ss = sort(s)

      recoveredBoot = recovered[ ,ss]
      deathsBoot = deaths[ ,ss]
      confirmedBoot = confirmed[ ,ss]

      Q0 = confirmedBoot[1]-recoveredBoot[1]-deathsBoot[1]
      I0 = 0.35*Q0
      E0 = 0.45*Q0
      R0 = recoveredBoot[1]
      D0 = deathsBoot[1]

      ActiveBoot = confirmedBoot-recoveredBoot-deathsBoot
      ActiveBoot[ActiveBoot<0] <- 0


      time = as.Date(names(confirmedBoot), format = "%m/%d/%y")

      paramsBoot = fit_SEIQRDP(Q = ActiveBoot, R = recoveredBoot, D = deathsBoot, Npop = Npop,
                               E0 = E0, I0 = I0, time = time, dt = dt, guess = guess,
                               trace = FALSE, ci = FALSE)


      resBoot = SEIQRDP(alpha = params$alpha1, beta = params$beta1,
                        gamma = params$gamma1, delta = params$delta1,
                        lambda0 = c(params$lambda01, params$lambda02, params$lambda03),
                        kappa0 = c(params$kappa01, params$kappa02, params$kappa03),
                        Npop, E0, I0, Q0, R0, D0,lambdaFun = params$lambdaFun,
                        kappaFun = params$kappaFun, tstart = start, tfinish = finish,
                        dt = dt, f =f)


      bootQuarantined[,b] = resBoot$quarantined
      bootSusceptible[,b] = resBoot$susceptible
      bootExposed[,b] = resBoot$exposed
      bootInfectious[,b] = resBoot$infectious
      bootRecovered[,b] = resBoot$recovered
      bootDead[,b] = resBoot$dead
      bootInsusceptible[,b] = resBoot$insusceptible

      print(b)

    }

    preds = c("susceptible", "exposed", "infectious", "quarantined", "recovered", "dead", "insusceptible")

    if(type == "norm"){

      m = res[["susceptible"]]
      se = apply(bootSusceptible, 1, sd)
      bias = m-apply(bootSusceptible, 1, mean)
      pred$susceptible_ll = m - bias - qnorm(((1+conf)/2))*se
      pred$susceptible_ul = m  - bias + qnorm(((1+conf)/2))*se

      m = res[["exposed"]]
      se = apply(bootExposed, 1, sd)
      bias = m-apply(bootExposed, 1, mean)
      pred$exposed_ll = m - bias - qnorm(((1+conf)/2))*se
      pred$exposed_ul = m  - bias + qnorm(((1+conf)/2))*se

      m = res[["infectious"]]
      se = apply(bootInfectious, 1, sd)
      bias = m-apply(bootInfectious, 1, mean)
      pred$infectious_ll = m - bias - qnorm(((1+conf)/2))*se
      pred$infectious_ul = m  - bias + qnorm(((1+conf)/2))*se

      m = res[["quarantined"]]
      se = apply(bootQuarantined, 1, sd)
      bias = m-apply(bootQuarantined, 1, mean)
      pred$quarantined_ll = m - bias - qnorm(((1+conf)/2))*se
      pred$quarantined_ul = m  - bias + qnorm(((1+conf)/2))*se

      m = res[["recovered"]]
      se = apply(bootRecovered, 1, sd)
      bias = m-apply(bootRecovered, 1, mean)
      pred$recovered_ll = m - bias - qnorm(((1+conf)/2))*se
      pred$recovered_ul = m  - bias + qnorm(((1+conf)/2))*se

      m = res[["dead"]]
      se = apply(bootDead, 1, sd)
      bias = m-apply(bootDead, 1, mean)
      pred$dead_ll = m - bias - qnorm(((1+conf)/2))*se
      pred$dead_ul = m  - bias + qnorm(((1+conf)/2))*se

      m = res[["insusceptible"]]
      se = apply(bootInsusceptible, 1, sd, na.rm = TRUE)
      bias = m-apply(bootInsusceptible, 1, mean)
      pred$insusceptible_ll = m - bias - qnorm(((1+conf)/2))*se
      pred$insusceptible_ul = m  - bias + qnorm(((1+conf)/2))*se

    }

    if(type == "perc"){

      alpha <- (1 + c(-conf, conf))/2

      norm.inter <- function (t, alpha)
      {
        t <- t[is.finite(t)]
        R <- length(t)
        rk <- (R + 1) * alpha
        if (!all(rk > 1 & rk < R))
          warning("extreme order statistics used as endpoints")
        k <- trunc(rk)
        inds <- seq_along(k)
        out <- inds
        kvs <- k[k > 0 & k < R]
        tstar <- sort(t, partial = sort(union(c(1, R), c(kvs, kvs +
                                                           1))))
        ints <- (k == rk)
        if (any(ints))
          out[inds[ints]] <- tstar[k[inds[ints]]]
        out[k == 0] <- tstar[1L]
        out[k == R] <- tstar[R]
        not <- function(v) xor(rep(TRUE, length(v)), v)
        temp <- inds[not(ints) & k != 0 & k != R]
        temp1 <- qnorm(alpha[temp])
        temp2 <- qnorm(k[temp]/(R + 1))
        temp3 <- qnorm((k[temp] + 1)/(R + 1))
        tk <- tstar[k[temp]]
        tk1 <- tstar[k[temp] + 1L]
        out[temp] <- tk + (temp1 - temp2)/(temp3 - temp2) * (tk1 -
                                                               tk)
        cbind(round(rk, 2), out)
      }


      for(i in 1:nrow(bootSusceptible)){
        qq <- norm.inter(bootSusceptible[i,], alpha)
        pred$susceptible_ll[i] = qq[1,2]
        pred$susceptible_ul[i] = qq[2,2]
      }


      for(i in 1:nrow(bootExposed)){
        qq <- norm.inter(bootExposed[i,], alpha)
        pred$exposed_ll[i] = qq[1,2]
        pred$exposed_ul[i] = qq[2,2]
      }

      for(i in 1:nrow(bootInfectious)){
        qq <- norm.inter(bootInfectious[i,], alpha)
        pred$infectious_ll[i] = qq[1,2]
        pred$infectious_ul[i] = qq[2,2]
      }

      for(i in 1:nrow(bootQuarantined)){
        qq <- norm.inter(bootQuarantined[i,], alpha)
        pred$quarantined_ll[i] = qq[1,2]
        pred$quarantined_ul[i] = qq[2,2]
      }


      for(i in 1:nrow(bootRecovered)){
        qq <- norm.inter(bootRecovered[i,], alpha)
        pred$recovered_ll[i] = qq[1,2]
        pred$recovered_ul[i] = qq[2,2]
      }

      for(i in 1:nrow(bootDead)){
        qq <- norm.inter(bootDead[i,], alpha)
        pred$dead_ll[i] = qq[1,2]
        pred$dead_ul[i] = qq[2,2]
      }

      for(i in 1:nrow(bootInsusceptible)){
        qq <- norm.inter(bootInsusceptible[i,], alpha)
        pred$insusceptible_ll[i] = qq[1,2]
        pred$insusceptible_ul[i] = qq[2,2]
      }

    }

  }
  actual = cbind.data.frame(active = as.numeric(Q), recovered = as.numeric(R), deaths=as.numeric(D), realTime = res$realTime)

  return(list(pred = pred, actual = actual, params = params, dt = dt))
}
