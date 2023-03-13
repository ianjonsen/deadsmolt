##' \code{sim_kelt_par} defines the simulation parameters & control scenarios used by \code{sim_kelt}.
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
##' @param scenario migration scenarios: rs = repeat spawner
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

sim_kelt_par <-
  function(temp = TRUE,
           advect = TRUE,
           growth = TRUE,
           scenario = "rs",
           land = FALSE,
           boundary = FALSE,
           ...) {

    dots <- list(...)

    pars <- list(
      N = 1440,
      start.dt = ISOdatetime(2023,05,25,16,00,00, tz = "UTC"),
      start = c(6912, 1465),
      coa = list(c(runif(1,6750,7280), runif(1, 1600,1900)), c(6912, 1465)),
      rd = 58, # repeat spawner reconditioning time at sea in days
      tsr = c(2, 17), # C from Daniels et al. ASF sim modelling
      pN = 0.75,
      nu = c(1, 2),
      rho = c(0.5, 0.8), # directional persistence for brw
      ntries = 1,
      psi = 0.9,
      uvm = 1, # magnitude of current vectors: if uvm < 1 current strength is down-scaled
      buffer = 5,
      al = 0, # wiebull scale parameter for move steps; if 0 then move steps are fixed at b
      bl = 2, # weibull scale parameter
      fl0 = 0.75,
      g = 0.001, # growth in forklength as % per day (re-scaled to hourly in simulation)
      surv = 0.998, ## daily survival rate
      pdrf = c(5, -0.02), # = p(0.5) @ 250 m  + < 0.01 @ 500 m   [c(4.865, -0.0139)  (~ consistent w HFX line V9 @ high power)]
      beta = c(-10, -10) # potential fn params to keep kelts off land
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
