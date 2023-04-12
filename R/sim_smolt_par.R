##' \code{sim_par} defines the simulation parameters & control scenarios used by \code{sim_smolt}.
##'
##' The movement process used predominantly in the simulation is
##' selected by the \code{move} argument.  Additional
##' parameters include: \code{temp} is movement temperature-dependent; \code{advect} do ocean
##' currents influence movement; \code{growth} do smolts grow in temp-dependent fashion;
##' \code{start.date};
##' \code{start} (location); \code{coa} centre of attraction (can be NULL);
##' \code{mdir} directional bias (a 2-element vector);
##' \code{rho} concentration of step directions (for wrapped-Cauchy, a 2-element vector); ...
##'
##' @title Control Values for \code{simulate}.
##' @param temp logical
##' @param advect logical
##' @param growth logical
##' @param scenario migration scenarios: sobi or mir
##' @param land keep track of sim rep hitting land (TRUE)
##' @param boundary keep track of sim rep hitting sim boundary (TRUE)
##' @param ... additional simulation control parameters
##' @return Returns a list with components
##'   \item{\code{move}}{the main move process}
##'   \item{\code{temp}}{temperature-dependent movements}
##'   \item{\code{advect}}{ocean-current-dependent movements}
##'   \item{\code{noise}}{add random noise to current u,v values}
##'   \item{\code{growth}}{temperature-dependent growth}
##'   \item{\code{shelf}}{movements constrained to stay on shelf}
##'   \item{\code{migs}}{migration strategy}
##'   \item{\code{land}}{record if smolt gets "stuck" on land, thereby ending simulation}
##'   \item{\code{boundary}}{record if smolt hits simulation boundary, thereby ending simulation}
##'   \item{\code{pars}}{list of additional, required control parameters}
##' @export

sim_smolt_par <-
  function(temp = TRUE,
           advect = TRUE,
           growth = TRUE,
           scenario = c("sobi", "mir"),
           land = FALSE,
           boundary = FALSE,
           ...) {

    scenario <- match.arg(scenario)

    dots <- list(...)

    pars <- list(
      N = 1440,
      start.dt = ISOdatetime(2023,05,25,16,00,00, tz = "UTC"),
      start = c(6912, 1465),
      coa = c(7120, 2350),
      tsr = c(4, 10), # C from Daniels et al. ASF sim modelling
      pN = 0.75,
      nu = 1, # strength of bias to CoA
      r = 0.005, # scaling param for magnitude of rho as fn of dist from CoA
      rho = 0.7, # directional persistence for crw
      uvm = 1, # magnitude of current vectors: if uvm < 1 current strength is down-scaled
      buffer = 5,
      bl = 2, # body-lengths / s
      fl0 = 0.146,
      g = 0.006, # growth in forklength as % per day (re-scaled to hourly in simulation)
      surv = 0.9936, ## daily survival rate
      reten = 0.845^(1/60),
      Dreten = 60,
      pdrf = c(5, -0.02), # = p(0.5) @ 250 m  + < 0.01 @ 500 m   [c(4.865, -0.0139)  (~ consistent w HFX line V9 @ high power)]
      beta = c(-7.5,-7.5) # potential fn params to keep kelts off land
    )

    ## overide default control pars
    pars[names(dots)] <- dots

    list(temp = temp,
         advect = advect,
         growth = growth,
         scenario = scenario,
         land = land,
         boundary = boundary,
         pars = pars)
  }

