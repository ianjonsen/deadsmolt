#' @title biased random walk movement kernel for Miramichi repeat-spawner kelts
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
    ## biased random walk toward a Center of Attraction
      if(i <= round(mpar$pars$rd * 24 * 0.9)) {
        if(all(!is.na(mpar$pars$coa[1,]))) {
          delta <- c(mpar$pars$coa[1,1] - xy[1], mpar$pars$coa[1,2] - xy[2])
          psi <- atan2(delta[1], delta[2])
          phi <- atan2(sin(xy[3]) + mpar$pars$nu[1] * sin(psi), cos(xy[3]) + mpar$pars$nu[1] * cos(psi))
        } else {
          phi <- atan2(sin(xy[3]), cos(xy[3]))
        }
        mu <- rwrpcauchy(1, phi, mpar$pars$rho[1])

      } else if (i > round(mpar$pars$rd * 24 * 0.9)) {
        delta <- c(mpar$pars$coa[2,1] - xy[1], mpar$pars$coa[2,2] - xy[2])
        psi <- atan2(delta[1], delta[2])
        phi <- atan2(sin(xy[3]) + mpar$pars$nu[2] * sin(psi), cos(xy[3]) + mpar$pars$nu[2] * cos(psi))
        mu <- rwrpcauchy(1, phi, mpar$pars$rho[2])
      }

    ## state 2: smolt outside of preferred T range
    ## smolts should slow down/reverse when < ts.min,
    ##   speed up when > ts.max to 'catch up' to preferred water mass
    if (ts < mpar$pars$tsr[1]) {
      ## slow down & head southward (presumably to warmer water mass)
      mu <- runif(1, 120, 240) / 180*pi
      s <- s * 0.75
    }
    ## speed up heading northward (presumably to colder water mass)
    if (ts > mpar$pars$tsr[2]) {
      mu <- (runif(1, -60, 60) / 180*pi) %% (2*pi)
      s <- s * 1.25
    }

    new.xy <- cbind(xy[1] + sin(mu) * s, xy[2] + cos(mu) * s)

    if (!is.null(data$land)) {
      ## calculate potential fn values
      pv <- c(extract(data$grad[[1]], new.xy)[1],
              extract(data$grad[[2]], new.xy)[1])

      new2.xy <- new.xy + pv * mpar$pars$beta

      ## if provisional new.xy is on land then try again
      if (!is.na(extract(data$land, rbind(new2.xy)))) {
        pv <- c(extract(data$grad[[1]], new2.xy)[1],
                extract(data$grad[[2]], new2.xy)[1])
        new3.xy <- new.xy + pv * (mpar$pars$beta * 3)
        ## if still on land then move back to water
        if (!is.na(extract(data$land, new3.xy))) {
          ## find all nearby cells within mpar$buffer km & select the first one in water
          cells <-
            extract(
              data$land,
              rbind(new.xy),
              buffer = mpar$pars$buffer,
              cellnumbers = TRUE,
              df = TRUE
            )
          idx <- which(is.na(cells[, 3]))[1]
          cell.water <- cells[idx, 2]
          new.xy <- xyFromCell(data$land, cell.water) %>% rbind()
        } else {
          new.xy <- new3.xy
        }
      } else {
        new.xy <- new2.xy
      }
    }

    cbind(new.xy[1], new.xy[2], mu, s)


  } else {
    stop("can only simulate kelt repeat-spawner movements")
  }

}
