#' @title biased, correlated random walk movement kernel for Miramichi repeat-spawner kelts
#'
#' @description utility function not to be called by user
#'
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rweibull
#' @importFrom raster extract xyFromCell
#' @keywords internal
#'
move_kernel_kelt <- function(data, xy = NULL, mpar, s, ts, i) {

  if(mpar$scenario == "rs") {
    ## biased, correlated random walk toward a Center of Attraction
      if(i <= round(mpar$pars$rd * 24 * 0.9)) {
        if(all(!is.na(mpar$pars$coa[[1]]))) {
          delta <- c(mpar$pars$coa[[1]][1] - xy[1], mpar$pars$coa[[1]][2] - xy[2])
          rho <- tanh(sqrt(delta[1]^2+delta[2]^2) * mpar$pars$r[1] + 0.2) - tanh(0.2)
          psi <- atan2(delta[1], delta[2])
          phi <- atan2(sin(xy[3]) + mpar$pars$nu[1] * sin(psi), cos(xy[3]) + mpar$pars$nu[1] * cos(psi))
          mu <- rwrpcauchy(1, phi, rho)
        } else {
          phi <- atan2(sin(xy[3]), cos(xy[3]))
          mu <- rwrpcauchy(1, phi, mpar$pars$rho)
        }


      } else if (i > round(mpar$pars$rd * 24 * 0.9)) {
        delta <- c(mpar$pars$coa[[2]][1] - xy[1], mpar$pars$coa[[2]][2] - xy[2])
        rho <- tanh(sqrt(delta[1]^2+delta[2]^2) * mpar$pars$r[2] + 0.1) - tanh(0.1)
        psi <- atan2(delta[1], delta[2])
        phi <- atan2(sin(xy[3]) + mpar$pars$nu[2] * sin(psi), cos(xy[3]) + mpar$pars$nu[2] * cos(psi))
        mu <- rwrpcauchy(1, phi, rho)
      }

    ## state 2: kelt outside of preferred T range
    ## kelts should reverse when < ts.min,
    ##   speed up when > ts.max to 'catch up' to preferred water mass
    if (ts < mpar$pars$tsr[1]) {
      ## slow down & head southward (presumably to warmer water mass)
      mu <- runif(1, 120, 240) / 180*pi
    }
    ## speed up heading northward (presumably to colder water mass)
    if (ts > mpar$pars$tsr[2]) {
      mu <- (runif(1, -60, 60) / 180*pi) %% (2*pi)
      s <- s * 1.25
    }

    new.xy <- cbind(xy[1] + sin(mu) * s, xy[2] + cos(mu) * s)

    ## calculate potential fn values
    pv <- c(extract(data$grad[[1]], new.xy)[1],
            extract(data$grad[[2]], new.xy)[1])

    new2.xy <- new.xy + pv * mpar$pars$beta

    ## if provisional new2.xy is on land then find a location in water, 1km from land
    if (!is.na(extract(data$land, rbind(new2.xy)))) {
      cat("finding water")
      cells <-
        extract(
          data$grad[[1]],
          rbind(new.xy),
          buffer = mpar$pars$buffer,
          cellnumbers = TRUE,
          df = TRUE
        )
      idx <- which(cells[, 3] == 0)[1]
      cell.water <- cells[idx, 2]
      new.xy <- xyFromCell(data$land, cell.water) %>% rbind()
    } else {
      new.xy <- new2.xy
    }

  cbind(new.xy[1], new.xy[2], mu, s)

  } else {
    stop("can only simulate kelt repeat-spawner movements")
  }

}
