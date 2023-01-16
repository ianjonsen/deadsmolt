#' @title random walk movement kernel for Miramichi smolts
#'
#' @description utility function not to be called by user
#'
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rweibull
#' @importFrom raster extract xyFromCell
#' @export
#'
move_kernel <- function(data, xy = NULL, mpar, i, s, ts, d2l, dir2l) {

    if (d2l > mpar$pars$buffer) {
      switch(mpar$scenario,
             sobi = {
                 ## state 1: migration biased toward SoBI
                 delta <- c(mpar$pars$coa[1] - xy[1], mpar$pars$coa[2] - xy[2])

                 mu <- atan2(delta[1], delta[2])

                 ## strength of bias to CoA is a hyberbolic fn of distance to CoA
                 #rho <- 1 + tanh(mpar$pars$r * sqrt(sum(delta^2)))

                 ## fixed rho gives strength of bias to the CoA
                 phi <- rwrpcauchy(1, mu, mpar$pars$rho[1])


                 ## state 2: smolt outside of preferred T range
                 ## smolts should slow down/reverse when < ts.min,
                 ##   speed up when > ts.max to 'catch up' to preferred water mass
                 if (ts < mpar$pars$tsr[1]) {
                   ## slow down & head southward (presumably to warmer water mass)
                   phi <- runif(1, 120, 240) / 180*pi
                   s <- s * 0.5
                 }
                 ## speed up heading northward (presumably to colder water mass)
                 if (ts > mpar$pars$tsr[2]) s <- s * 1.5
             })

    } else {
      ## if close to land (within mpar$pars$buffer km but farther than 2km)
        if(d2l > 2) {
          if(xy[2] < 1475) {
            ## direct kelt to move directly toward CoA (N end of SoBI) if in
            ##  Miramichi Bay or outside of Bay but further South, _OR_ if in SoBI
            phi <- atan2(mpar$pars$coa[1] - xy[1], mpar$pars$coa[2] - xy[2])
          } else {
            ## move parallel to land
            phi <- as.numeric(dir2l + 0.5 * pi)
            if(all(xy[1] > 7097, xy[2] > 1915)) {
              #reverse parallel-to-land direction to travel is to N along NF coast
              phi <- as.numeric(dir2l - 0.5 * pi)
            }
          }
        } else {
        ## move in opposite direction of land if within 2km
        phi <- as.numeric(dir2l + pi)
        }
    }

    new.xy <- c(xy[1] + sin(phi) * s, xy[2] + cos(phi) * s)

    ## check if new xy close/on land
    new.d2l <- extract(data$land, rbind(new.xy))

    ## if new location on land (0) then adjust so it's in water
    if(!is.na(new.d2l) & new.d2l == 1) {
      ## find all nearby cells within 5 km & select the one farthest from land
      cells <- extract(data$d2land, rbind(new.xy), buffer = 5, cellnumbers = TRUE, df = TRUE)
      cell.max <- cells[cells[, 3] == max(cells[, 3], na.rm = TRUE), 2][1]
      new.xy <- xyFromCell(data$d2land, cell.max) %>% as.vector()

    } else if(is.na(new.d2l)) {
      new.xy <- c(NA,NA)
    }

  cbind(new.xy[1], new.xy[2], s)

}
