#' @title random walk movement kernel for Miramichi smolts
#'
#' @description utility function not to be called by user
#'
#' @importFrom CircStats rwrpcauchy
#' @importFrom stats rweibull
#' @importFrom raster extract xyFromCell
#' @export
#'
moveSmir <- function(data, xy = NULL, mpar, i, s, ts, w, d2l, dir2l) {

    phi <- NULL
    ## regular movement

    if (d2l > mpar$pars$buffer) {
      switch(mpar$scenario,
             sobi = {
                 ## state 1: migration toward SoBI
                 delta <- c(mpar$pars$coa[1] - xy[1], mpar$pars$coa[2] - xy[2])

                 mu <- atan2(delta[1], delta[2])
                 #rho <- 1 + tanh(mpar$pars$r * sqrt(sum(delta^2)))
                 phi <- rwrpcauchy(1, mu, mpar$pars$rho[1])
                 ## state 2: stop migration and remain in current water mass
                 ## NEED TO REVISIT THIS - smolts should slow down/reverse when < ts.min,
                 ##   speed up when > ts.max to 'catch up' to preferred water mass
                 if (any(ts < mpar$pars$ts.min, ts > mpar$pars$ts.max)) s <- 0
             })
    } else {
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
    if(!is.na(new.d2l) & new.d2l == 0) {
      ## find all nearby cells within 5 km & select the one farthest from land
      cells <- extract(data$land, rbind(new.xy), buffer = 5, cellnumbers = TRUE, df = TRUE)
      cell.max <- cells[cells[, 3] == max(cells[, 3], na.rm = TRUE), 2][1]
      new.xy <- xyFromCell(data$land, cell.max) %>% as.vector()

    } else if(is.na(new.d2l)) {
      new.xy <- c(NA,NA)
    }

  cbind(new.xy[1], new.xy[2], s)

}
