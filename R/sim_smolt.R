#' @title simulate smolt migration
#'
#' @description simulates smolt migration
#'
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#'
#' @param id - identifier for simulation run (individual animal)
#' @param data - a list of required data from \code{presim}
#' @param mpar - simulation control parameters supplied as a list, see details
#' @param pb - use progress bar (logical)
#' @importFrom raster extract xyFromCell
#' @importFrom CircStats rwrpcauchy
#' @importFrom dplyr %>% mutate lag
#' @importFrom tibble as_tibble
#' @importFrom stats runif rbinom
#' @importFrom lubridate week yday
#' @importFrom stringr str_split
#' @export

sim_smolt <-
  function(id=1,
           data = NULL,
           mpar = sim_par(),
           pb = TRUE
  ) {


    if (is.null(data))
      stop("Can't find output from sim_setup()\n")
    if (class(data$d2land)[1] != "RasterLayer") stop("d2land must be a RasterLayer")

    N <- mpar$pars$N

    ## define location matrix & initialise start position
    ## ds - active swimming displacements
    ## dl - displacements to deflect away from land
    xy <- matrix(NA, N, 2)
    xy[1,] <- cbind(mpar$pars$start)
    ds <- matrix(NA, N, 3) #x, y, s

    ## define other vectors
    ## reten - tag retained = 1, expelled = 0
    ## surv - smolt alive = 1, dead = 0
    ## d2l - distance to land, used in movement kernel
    ## dir2l - direction to land (rad), used in movement kernel
    ## u - advection in E-W (m/s)
    ## v - advection in N-S (m/s)
    ## ts - water temp (C)
    reten <- dir <- surv <- d2l <- dir2l <- u <- v <- ts <- vector("numeric", N)
    reten[1] <- 1
    surv[1] <- 1

    ## initialise growth params or set fixed params if no growth
    if(mpar$growth) {
      ## s - swim velocity based on forklength and specified speed/bodylength (m/s)
      ## fl - forklength in m
      s <- fl <- vector("numeric", N)
      fl[1] <- mpar$pars$fl0
      s[1] <- fl[1] * mpar$pars$b * 3.6 ## initial swim speed (fl * b body-lengths / s) - in km/h
    } else {
      ## No growth - fixed size/speed params
      fl <- rep(mpar$pars$fl0, N)
      s <- rep(fl * mpar$pars$b * 3.6, N)
    }

    ## iterate movement
    for (i in 2:N) {
      if(i==2 && pb)  tpb <- txtProgressBar(min = 2, max = N, style = 3)
      ## extract Temperature
      ts[i - 1] <- extract(data$ts[[yday(mpar$pars$start.dt + i * 3600)]],
                           rbind(xy[i - 1,])) - 273.15
      if (is.na(ts[i - 1])) {
        ## calc mean Temp within 2 km buffer of location @ time i-1
        ts[i - 1] <-
          extract(data$ts[[yday(mpar$pars$start.dt + i * 3600)]],
                  rbind(xy[i - 1,]),
                  buffer = 4,
                  df = TRUE)[, 2] %>%
          mean(., na.rm = TRUE) - 273.15
      }

      ### Apply Energetics
      if (mpar$growth) {
        if (is.na(ts[i - 1])) {
          cat("\n stopping simulation: NA value for temperature")
          break
        }
        ## calculate growth in current time step based on assumed daily % gr rate
        ##  smolts only growth when in optimal temp range defined by tsr
        if (all(ts[i - 1] > mpar$pars$tsr[1], ts[i - 1] <= mpar$pars$tsr[2])) {
          fl[i] <- fl[i-1] * (1 + mpar$pars$g) ^ (1/24) #re-scales daily rate to hourly
        } else {
          fl[i] <- fl[i-1]
        }

        ## determine size-based average step for current time step
        ## assume avg swim speed of b bodylengths/s
        s[i] <- fl[i] * mpar$pars$b * 3.6 ## forklength * b m/s converted to km/h
      }

      ## calculate distance & direction to land
      d2l[i-1] <- extract(data$d2land, rbind(xy[i-1, ]))
      dir2l[i-1] <- extract(data$dir2land, rbind(xy[i-1, ]))

      ## Movement kernel
      ds[i,] <- move_kernel(data,
                         xy = xy[i-1,],
                         mpar = mpar,
                         i,
                         s = s[i],
                         ts = ts[i-1],
                         d2l = d2l[i-1],
                         dir2l = dir2l[i-1])

      ### Current Advection
      if (mpar$advect) {
        ## determine envt'l forcing
        ## determine advection due to current, convert from m/s to km/h
        u[i] <- extract(data$u[[yday(mpar$pars$start.dt + i * 3600)]],
                        rbind(xy[i - 1, ]), method = "simple") * 3.6 * mpar$par$uvm
        v[i] <- extract(data$v[[yday(mpar$pars$start.dt + i * 3600)]],
                        rbind(xy[i - 1, ]), method = "simple") * 3.6 * mpar$par$uvm

        if (any(is.na(u[i]), is.na(v[i]))) {
          ## calc mean Temp within 2 km buffer of location @ time i-1
          u[i] <- v[i] <- 0
        }

        ## turn off advection in sobi.box b/c too challenging to get smolts through w currents...
      } else if(!mpar$advect | all(xy[1] >= data$sobi.box[1],
                                   xy[1] <= data$sobi.box[2],
                                   xy[2] >= data$sobi.box[3],
                                   xy[2] <= data$sobi.box[4])) {
        u[i] <- v[i] <- 0
      }

      xy[i, 1:2] <- cbind(ds[i, 1] + u[i],
                          ds[i, 2] + v[i])
      ## overwrite s[i] if set to 0 (outside of preferred temp range) in movement kernel
      s[i] <- ds[i,3]

      if(!is.na(extract(data$land, rbind(xy[i, ])))  & any(!is.na(xy[i,]))) {
        mpar$land <- TRUE
        cat("\n stopping simulation: stuck on land")
        break
      }

      if(any(is.na(xy[i, ]))) {
        mpar$boundary <- TRUE
        cat("\n stopping simulation: hit a boundary")
        break
      }

      ## determine survival
      if(!is.na(mpar$pars$surv)) {
        surv[i] <-
          rbinom(1, 1, mpar$pars$surv ^ (1 / 24)) # rescales daily survival to hourly
        if (surv[i] == 0) {
          cat("\n smolt is dead")
          break
        }
      }

      ## determine tag retention
      if(!is.na(mpar$pars$reten) & i <= mpar$pars$Dreten*24) {
        reten[i] <- rbinom(1, 1, mpar$pars$reten ^ (1 / 24))
        if(reten[i] == 0) {
          cat("\n tag expulsion")
          break
        }
      } else if(!is.na(mpar$pars$reten) & i > mpar$pars$Dreten*24) {
        reten[i] <- 1
      }

      if(pb){
        setTxtProgressBar(tpb, i)
        if(i==N) close(tpb)
      }

    }

    N <- ifelse(!is.na(which(is.na(xy[,1]))[1] - 1), which(is.na(xy[,1]))[1] - 1, N)
    X <-
      data.frame(
        x = xy[, 1],
        y = xy[, 2],
        dx = ds[, 1] - lag(xy[, 1]),
        dy = ds[, 2] - lag(xy[, 2]),
        u = u,
        v = v,
        ts = ts,
        d2land = d2l,
        dir2land = dir2l,
        fl = fl,
        s = s,
        surv = surv,
        reten = reten
      )[1:N,]

    if(sum(is.na(X$ts)) == nrow(X)) {
      X <- X %>% select(-ts)
    } else {
      X$ts[nrow(X)] <- ifelse(X$ts[nrow(X)] == 0, NA, X$ts[nrow(X)])
    }

    sim <- X %>% as_tibble()

    ## re records after sim is stopped for being stuck on land, etc...
    if(mpar$land | mpar$boundary) {
      sim <- sim %>%
        filter((!is.na(x) & !is.na(y) & fl != 0 & s != 0) | surv != 1 | reten != 1)
    }

    nsim <- nrow(sim)

    sim <- sim %>%
      mutate(id = id) %>%
      mutate(date = seq(mpar$pars$start.dt, by = 3600, length.out = nsim)) %>%
      select(id, date, everything())

    param <- mpar
    out <- list(sim = sim, params = param)
    class(out) <- "deadsmolt"

    return(out)
  }
