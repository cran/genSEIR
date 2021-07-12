#' Plots for Epidemic Curves
#'
#' This function creates plots for reported and predicted active, recovered and death cases.
#'
#'
#' @param object a predict_SEIQRDP result.
#' @param reported a logical argument. If \code{TRUE} reported official cases will be added to the plot.
#' @param sep a logical argument. If \code{TRUE} seperate plots will be plotted. If \code{FALSE} one plot with all desired states will be plotted.
#' @param show select one or more desired state. \code{S}: Susceptible, \code{E}: Exposed, \code{I}: Infectious, \code{Q}: Quarantined, \code{R}: Recovered, \code{D}: Dead, \code{P}: Insusceptible
#' @param ci a logical argument. If \code{TRUE} a bootstrap confidence intetval will be added to the plot.
#' @param title an optional title for the plot.
#' @param checkRates if \code{TRUE} compares the fitted and calcualted death and recovered ratios through plots
#' @param ... other plot options
#'
#' @importFrom ggplot2 geom_ribbon ggplot geom_line geom_point aes scale_colour_manual scale_fill_manual xlab ylab ggtitle labs theme
#'
#' @export plot_SEIQRDP
#'
#' @return plots for epidemic curves: active cases, recovered and deaths
#'
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @examples
#' \donttest{
#'  alpha_guess = 0.45
#'  beta_guess = 1
#'  LT_guess = 2
#'  Q_guess = 0.55
#'  lambda_guess = c(0.01,0.01,30)
#'  kappa_guess = c(0.01,0.001,30)
#'
#'  guess = list(alpha_guess,
#'               beta_guess,
#'               1/LT_guess,
#'               Q_guess,
#'               lambda_guess[1],
#'               lambda_guess[2],
#'               lambda_guess[3],
#'               kappa_guess[1],
#'               kappa_guess[2],
#'               kappa_guess[3])
#'
#'  pred = predict_SEIQRDP(country = "Germany", start = "10/15/20", finish = "12/15/20",
#'                         dt = 1, f = 30, conf = 0.95, Npop = 80000000, guess, boot = TRUE,
#'                         seed = 123, repeatNumber = 10, bootSample = NULL, type = "norm")
#'
#'
#'  plot_SEIQRDP(object = pred, sep = FALSE, ci = TRUE, show = c("Q", "R", "D"), checkRates = TRUE)
#'
#'
#' }
#'
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}

plot_SEIQRDP <- function(object, reported = TRUE, sep = FALSE, show =  c("S", "E", "I", "Q", "R", "D", "P"),
                         ci = FALSE, title = NULL, checkRates = FALSE, ...){

  realTime = as.Date(object$actual$realTime)
  simTime = as.Date(object$pred$simTime)

  pred = object$pred
  actual = object$actual




  if(sep){

    if("R" %in% show) {pRsim = ggplot(pred) + geom_line(aes(y = recovered_pred, x=simTime, color = "Recovered"))+
      xlab("Time") + ylab("Recovered Cases")+ggtitle(title)+ theme(legend.position="none")}else{pRsim = NULL}
    if("R" %in% show & ci)  pRsim =  pRsim + geom_ribbon(aes(ymin = recovered_ll, ymax = recovered_ul, x = simTime), alpha = 0.3)
    if("Q" %in% show) {pQsim = ggplot(pred) + geom_line(aes(y = quarantined_pred, x=simTime, color="Active"))+
      xlab("Time") + ylab("Active Cases")+ggtitle(title)+ theme(legend.position="none")}
    if("Q" %in% show & ci) {pQsim =  pQsim + geom_ribbon(aes(ymin = quarantined_ll, ymax = quarantined_ul, x = simTime), alpha = 0.3)}
    if("D" %in% show) {pDsim = ggplot(pred) + geom_line(aes(y = dead_pred, x=simTime, color="Deceased"))+
      xlab("Time") + ylab("Death Cases")+ggtitle(title)+ theme(legend.position="none")}
    if("D" %in% show & ci) pDsim =  pDsim + geom_ribbon(aes(ymin = dead_ll, ymax = dead_ul, x = simTime), alpha = 0.3)else{pDsimCi = NULL}
    if("S" %in% show) {pSsim = ggplot(pred) + geom_line(aes(y = susceptible_pred, x=simTime, color="Susceptible"))+
      xlab("Time") + ylab("Susceptible Cases")+ggtitle(title)+ theme(legend.position="none")}else{pSsim = NULL}
    if("S" %in% show & ci) pSsim =  pSsim + geom_ribbon(aes(ymin = susceptible_ll, ymax = susceptible_ul, x = simTime), alpha = 0.3)
    if("E" %in% show) {pEsim = ggplot(pred) + geom_line(aes(y = exposed_pred, x=simTime, color="Exposed"))+
      xlab("Time") + ylab("Exposed Cases")+ggtitle(title)+ theme(legend.position="none")}else{pEsim = NULL}
    if("E" %in% show & ci) pEsim =  pEsim + geom_ribbon(aes(ymin = exposed_ll, ymax = exposed_ul, x = simTime), alpha = 0.3)
    if("I" %in% show) {pIsim = ggplot(pred) + geom_line(aes(y = infectious_pred, x=simTime, color="Infectious"))+
      xlab("Time") + ylab("Infectious Cases")+ggtitle(title)+ theme(legend.position="none")}else{pIsim = NULL}
    if("I" %in% show & ci) pIsim =  pIsim + geom_ribbon(aes(ymin = infectious_ll, ymax = infectious_ul, x = simTime), alpha = 0.3)
    if("P" %in% show) {pPsim = ggplot(pred) + geom_line(aes(y = insusceptible_pred, x=simTime, color="Insusceptible"))+
      xlab("Time") + ylab("Insusceptible Cases")+ggtitle(title)+ theme(legend.position="none")}else{pPsim = NULL}
    if("P" %in% show & ci) pPsim =  pPsim + geom_ribbon(aes(ymin = insusceptible_ll, ymax = insusceptible_ul, x = simTime), alpha = 0.3)



    if("R" %in% show & reported) {pRact = ggplot(data = actual) + geom_point(aes(y = recovered, x = realTime, fill = "Recovered"), size=2, shape=21, stroke=0)+
      xlab("Time") + ylab("Confirmed Recovered Cases")+ggtitle(title) + theme(legend.position="none")}else{pRact = NULL}
    if("Q" %in% show & reported) {pQact = ggplot(data = actual) + geom_point(aes(y = active, x = realTime, fill = "Active"), size=2, shape=21, stroke=0)+
      xlab("Time") + ylab("Confirmed Active Cases")+ggtitle(title) + theme(legend.position="none")}else{pQact = NULL}
    if("D" %in% show & reported) {pDact = ggplot(data = actual) + geom_point(aes(y = deaths, x = realTime, fill = "Deceased"), size=2,shape=21, stroke=0)+
      xlab("Time") + ylab("Confirmed Death Cases")+ggtitle(title) + theme(legend.position="none")}else{pDact = NULL}

    p = list(pRsim, pQsim, pDsim, pSsim, pEsim, pIsim, pPsim, pRact, pQact, pDact)

  }else{


    if("R" %in% show) {pRsim = geom_line(data = pred,  aes(y = recovered_pred, x=simTime, color = "Recovered"))}else{pRsim = NULL}
    if("R" %in% show & ci)  pRsimCi = geom_ribbon(data = pred,  aes(ymin = recovered_ll, ymax = recovered_ul, x = simTime), alpha = 0.3)else{pRsimCi = NULL}
    if("Q" %in% show) {pQsim = geom_line(data = pred,  aes(y = quarantined_pred, x=simTime, color="Active"))}else{pQsim = NULL}
    if("Q" %in% show & ci) {pQsimCi =  geom_ribbon(data = pred,  aes(ymin = quarantined_ll, ymax = quarantined_ul, x = simTime), alpha = 0.3)}else{pQsimCi = NULL}
    if("D" %in% show) {pDsim = geom_line(data = pred,  aes(y = dead_pred, x=simTime, color="Deceased"))}else{pDsim = NULL}
    if("D" %in% show & ci) pDsimCi = geom_ribbon(data = pred,  aes(ymin = dead_ll, ymax = dead_ul, x = simTime), alpha = 0.3)else{pDsimCi = NULL}
    if("S" %in% show) {pSsim = geom_line(data = pred,  aes(y = susceptible_pred, x=simTime, color="Susceptible"))}else{pSsim = NULL}
    if("S" %in% show & ci) pSsimCi = geom_ribbon(data = pred,  aes(ymin = susceptible_ll, ymax = susceptible_ul, x = simTime), alpha = 0.3)else{pSsimCi = NULL}
    if("E" %in% show) {pEsim = geom_line(data = pred,  aes(y = exposed_pred, x=simTime, color="Exposed"))}else{pEsim = NULL}
    if("E" %in% show & ci) pEsimCi = geom_ribbon(data = pred,  aes(ymin = exposed_ll, ymax = exposed_ul, x = simTime), alpha = 0.3)else{pEsimCi = NULL}
    if("I" %in% show) {pIsim = geom_line(data = pred,  aes(y = infectious_pred, x=simTime, color="Infectious"))}else{pIsim = NULL}
    if("I" %in% show & ci) pIsimCi = geom_ribbon(data = pred,  aes(ymin = infectious_ll, ymax = infectious_ul, x = simTime), alpha = 0.3)else{pIsimCi = NULL}
    if("P" %in% show) {pPsim = geom_line(data = pred,  aes(y = insusceptible_pred, x=simTime, color="Insusceptible"))}else{pPsim = NULL}
    if("P" %in% show & ci) pPsimCi = geom_ribbon(data = pred,  aes(ymin = insusceptible_ll, ymax = insusceptible_ul, x = simTime), alpha = 0.3)else{pPsimCi = NULL}

    if("R" %in% show & reported) {pRact = geom_point(data = actual,  aes(y = recovered, x = realTime, fill = "Recovered"), size=2, shape=21, stroke=0)}else{pRact = NULL}
    if("Q" %in% show & reported) {pQact = geom_point(data = actual,  aes(y = active, x = realTime, fill = "Active"), size=2, shape=21, stroke=0)}else{pQact = NULL}
    if("D" %in% show & reported) {pDact =  geom_point(data = actual,  aes(y = deaths, x = realTime, fill = "Deceased"), size=2,shape=21, stroke=0)}else{pDact = NULL}

    p =  ggplot() + pRsim + pRsimCi + pQsim + pQsimCi + pDsim + pDsimCi + pSsim + pSsimCi + pEsim + pEsimCi + pIsim + pIsimCi +
      pPsim + pPsimCi + pRact + pQact + pDact

    p = p + labs(fill="Reported", color="Fitted") +
      xlab("Time") + ylab("Cases")+ggtitle(title) +
      scale_fill_manual(values = c("Recovered" = "blue", "Active" = "red",
                                   "Deceased" = "black")) +
      scale_colour_manual(values = c("Recovered" = "blue", "Active" = "red",
                                     "Deceased" = "black", "Susceptible" = "cyan",
                                     "Exposed" = "pink", "Infectious" = "magenta",
                                     "Insusceptible" = "darkolivegreen"))

  }

  if(checkRates){

    checkRates = checkRates(object$actual$realTime, object$actual$active, object$actual$recovered, object$actual$deaths,
                            object$params$kappaFun, object$params$lambdaFun, c(object$params$kappa01, object$params$kappa02,
                                                                               object$params$kappa03), c(object$params$lambda01, object$params$lambda02, object$params$lambda03),
                            object$dt)



  }

  return(list(p, checkRates))

}

