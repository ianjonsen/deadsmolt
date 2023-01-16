#' @title Pre-simulation setup
#'
#' @description load required rasters, receiver locations
#'
#'
#' @author Ian Jonsen \email{ian.jonsen@mq.edu.au}
#'
#' @param config - path to config.R script containing file.paths for required & optional environmental layers (see Details)
#' @param deadsmolt - logical; turns on Miramichi & SoBI receiver locations if TRUE
#' @param esrf - logical; turns off receiver & ESRF polygons if FALSE
#' @importFrom raster raster stack brick projectRaster extract
#' @importFrom sp coordinates<- proj4string<- CRS spTransform SpatialPointsDataFrame spsample
#' @importFrom sf st_as_sf st_sample st_coordinates st_distance
#' @importFrom dplyr select filter rename bind_cols %>% tibble distinct
#' @export
#'
sim_setup <-
  function(config = config, deadsmolt = FALSE, esrf = TRUE) {

    suppressWarnings(source(config, local = TRUE, echo=FALSE))
    if(is.null(prj)) prj <- "+proj=stere +lat_0=90 +lon_0=-100 +k=0.933012425899506 +x_0=4245000 +y_0=5295000 +R=6371229 +units=km +no_defs"

    out <- list(
      land = suppressWarnings(raster(land)),
      d2land = suppressWarnings(raster(d2land)),
      dir2land = suppressWarnings(raster(dir2land))
    )

    out[["u"]] <- suppressWarnings(stack(file.path(riops, "riops_doy_u.grd")))
    out[["v"]] <- suppressWarnings(stack(file.path(riops, "riops_doy_v.grd")))
    out[["ts"]] <- suppressWarnings(stack(file.path(riops, "riops_doy_t.grd")))

    if(deadsmolt) {
    out[["recLocs"]] <- ds_rec
    }

    if(esrf) {
    out[["recLocs"]] <- esrf_rec
    out[["recPoly"]] <- recPoly_sf
    out[["sobi.box"]] <- c(980,1030,1230,1275)
    out[["esrfPoly"]] <- readRDS(file.path(polygon.file, "NLpoly.RDS"))
    }

    out[["prj"]] <- prj

    return(out)
  }
