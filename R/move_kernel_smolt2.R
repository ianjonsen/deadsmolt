#' @title biased random walk movement kernel for Miramichi smolts
#'
#' @description utility function not to be called by user
#'
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rweibull
#' @importFrom raster extract xyFromCell
#' @keywords internal
#'
move_kernel_smolt2 <- function(data, xy = NULL, mpar, i, s, ts) {
  switch(mpar$scenario,
         sobi = {
           ## state 1: migration biased toward SoBI
           delta <- c(mpar$pars$coa[1] - xy[1],
                      mpar$pars$coa[2] - xy[2])
           # rho <-
           #   tanh(sqrt(delta[1] ^ 2 + delta[2] ^ 2) * mpar$pars$r + 0.2) - tanh(0.2)
           psi <- atan2(delta[1], delta[2])
           phi <-
             atan2(sin(xy[3]) + mpar$pars$nu[1] * sin(psi),
                   cos(xy[3]) + mpar$pars$nu[1] * cos(psi))
           mu <- rwrpcauchy(1, phi, mpar$pars$rho)

           ## state 2: smolt outside of preferred T range
           ## smolts should slow down/reverse when < ts.min,
           ##   speed up when > ts.max to 'catch up' to preferred water mass
           if (ts < mpar$pars$tsr[1]) {
             ## slow down & head southward (presumably to warmer water mass)
             mu <- runif(1, 120, 240) / 180 * pi
             s <- s * 0.5
           }
           ## speed up heading northward (presumably to colder water mass)
           if (ts > mpar$pars$tsr[2])
             s <- s * 1.5

           new.xy <-
             cbind(xy[1] + sin(mu) * s, xy[2] + cos(mu) * s)
         },
         mir = {
           ## state 1: migration biased toward Miramichi Bay receivers
           delta <- c(mpar$pars$coa[1] - xy[1],
                      mpar$pars$coa[2] - xy[2])
           # rho <-
           #   tanh(sqrt(delta[1] ^ 2 + delta[2] ^ 2) * mpar$pars$r + 0.2) - tanh(0.2)
           psi <- atan2(delta[1], delta[2])
           phi <- atan2(sin(xy[3]) + mpar$pars$nu[1] * sin(psi),
                        cos(xy[3]) + mpar$pars$nu[1] * cos(psi))
           mu <- rwrpcauchy(1, phi, mpar$pars$rho)

           new.xy <-
             cbind(xy[1] + sin(mu) * s, xy[2] + cos(mu) * s)

           if (mpar$scenario == "mir") {
             dir2l <- extract(data$dir, new.xy)

             ## movement to avoid eastern island
             if (new.xy[1] > 6918.75 &
                 new.xy[2] > 1468 & new.xy[2] < 1470.25) {
               ## move westward, parallel to land
               mu <- as.numeric(dir2l - 0.5 * pi)
             }

             ## movement to avoid central island
             if (new.xy[1] > 6916.25 &
                 new.xy[2] > 1471 & new.xy[2] < 1473) {
               ## move eastward, parallel to land
               mu <- as.numeric(dir2l + 0.5 * pi)
             }

             if (new.xy[1] <= 6916.25 &
                 new.xy[2] > 1472.25 & new.xy[2] < 1476.4) {
               ## move westward, parallel to island
               mu <- as.numeric(dir2l - 0.5 * pi)
             }

             ## movement to avoid western island
             if (new.xy[1] > 6912.5 &
                 new.xy[1] < 6914.8 & new.xy[2] > 1477.5) {
               ## move southeastward, parallel to island
               mu <- as.numeric(dir2l + 0.5 * pi)
             }

             ## movement to avoid mainland at far northwest of Bay
             if (new.xy[1] <= 6912.5 & new.xy[2] > 1477) {
               ## move south, parallel to land
               mu <- as.numeric(dir2l - 0.5 * pi)
             }
           }

           new.xy <-
             cbind(new.xy[1] + sin(mu) * s, new.xy[2] + cos(mu) * s)

         })

  ## calculate potential fn values
  pv <- c(extract(data$grad[[1]], new.xy)[1],
          extract(data$grad[[2]], new.xy)[1])

  new2.xy <- new.xy + pv * mpar$pars$beta

  ## if provisional new2.xy is on land then find a location in water, 1km from land
  if (!is.na(extract(data$land, rbind(new2.xy)))) {
    message("finding water")
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

}
