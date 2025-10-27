#' @export
circ.plot <- function(x, r=NULL, weight=TRUE, ci=0.95,cords="geo", sep=0.05, pch=16, mvlength=0.2, mvlwd=3, arrowlength=0.1, arrowlwd=0.1, shrink=1.2, zero=pi/2, stack=T, rotation="clock", col="gray40", mv=TRUE, mvcol="blue", cicol="black", cilwd=2, main=NULL, dist=0.2, citick=0.05) {
  stopifnot(is.numeric(x))
  x <- circular(deg2rad(x))



  if (cords=="geo") {
    cords=c("E", "N", "W", "S")}
  else{
    cords=c(90,0,270,180)

  }


  plot.circular(x, shrink=shrink, sep=sep, stack=stack, pch=pch, col=col, zero=zero,
                main=main, axes=F, rotation=rotation)
  axis.circular(at=c(0, pi/2, pi, 3*pi/2), labels=cords, zero=NULL, rotation=rotation)
  if (!is.null(r)) {
    arrows.circular(x, r, col=col, lwd=arrowlwd, length=arrowlength, zero=zero,  rotation=rotation)
    if (!isFALSE(mv) & isTRUE(weight)) {
      arrows.circular(circular_mean_rad(x, r), circ_rho_zar(x, r), lwd=mvlwd, length=mvlength, col=mvcol, zero=zero, rotation=rotation)
    }

    else
      if (!isFALSE(mv)){
        {  arrows.circular(circular_mean_rad(x), rho.circular(x), lwd=mvlwd, length=mvlength, col=mvcol, zero=zero, rotation=rotation)
        }
      }
  }

  if (isTRUE(weight)) {
    a <- circ_MMRT(rad2deg(x),r)
    p <- as.numeric(c[[1]])
  }
  else {
    a <- rayleigh.test(circular(x))
    p <- as.numeric(a[[2]])
  }


  if (!is.null(ci)) {
    if (isTRUE(weight)) {
      confint <- confidence_interval_weighted(rad2deg(x),w=r, conf.level = ci, axial=F)
    }
    else {
      confint <- confidence_interval_weighted(rad2deg(x),w=NULL, conf.level = ci, axial=F)

    }
    ci_low <- deg2rad(confint$conf.interval[1]%%360)
    ci_high <- deg2rad(confint$conf.interval[2]%%360)

    angles <- if(ci_low < ci_high) {
      seq(ci_low, ci_high, length=200)
    } else {
      c(seq(ci_low, 2*pi, length=100),
        seq(0, ci_high, length=100))
    }

    if (p<0.05){
      lines.circular(circular(angles), y=rep(dist, length(angles)),
                               col=cicol, lwd=cilwd, zero=zero, axes=F, rotation=rotation)
      # endpoint ticks
      for(j in c(ci_low, ci_high)) {
        lines.circular(circular(rep(j,2)),
                       y = c(dist - citick, dist + citick),
                       col=cicol, lwd=cilwd, zero=zero, rotation=rotation)
      }
    }
    else {
      warning("Uniformity test: p>0.05. No CI produced")
    }


  }
}

#' @export
confidence_angle_weighted <- function (x, conf.level = 0.95, w = NULL, axial = FALSE, na.rm = TRUE)
{
  if (axial) {
    f <- 2
  }
  else {
    f <- 1
  }
  Z_alpha <- z_score(conf.level)
  sde <- circular_sd_error_weighted(x, w, axial, na.rm)
  temp <- Z_alpha * sde
  if (temp > 1)
    temp <- 1
  rad2deg(asin(temp) * f)
}

#' @export
circ_rho_zar <- function(a,r) {
  x<-list()
  y<-list()
  for (i in 1:length(a)) {
    x[i] <- r[i] * cos(a[i])
    y[i] <- r[i] * sin(a[i])
  }
  X <- sum(unlist(x))/length(a)
  Y <- sum(unlist(y))/length(a)
  r_mw <- as.numeric(sqrt(X^2+Y^2))
  return(r_mw)
}

# R vector length now calculates the zar weighted r vector length from circ_rho_zar(), similar to Oriana.
#' @export
circular_sd_error_weighted <- function (x, w = NULL, axial = FALSE, na.rm = TRUE)
{
  if (axial) {
    f <- 2
  }
  else {
    f <- 1
  }
  if (is.null(w)) {
    w <- rep(1, times = length(x))
  }
  data <- cbind(x = x, w = w)
  if (na.rm) {
    data <- data[stats::complete.cases(data), ]
  }
  x <- data[, "x"]
  w <- data[, "w"]
  n <- length(x)
  kappa <- est.kappa(x, w = w, axial = axial)
  x <- (x * f)%%360
  R <- circ_rho_zar(deg2rad(x), r = w)
  1/sqrt(n * R * kappa)
}

#' @export
confidence_interval_weighted <- function (x, conf.level = 0.95, w = NULL, axial = FALSE, na.rm = TRUE)
{
  conf.angle <- confidence_angle_weighted(x, conf.level, w, axial,
                                          na.rm)
  mu <- circular_mean_deg(x, w = w)
  list(mu = mu, conf.angle = conf.angle, conf.interval = c(mu -
                                                             conf.angle, mu + conf.angle))
}

#' @export
deg2rad <- function(deg) {
  stopifnot(is.numeric(deg))
  (rad <- (pi/180)*deg )
}

#' @export
rad2deg <- function(rad) {
  stopifnot(is.numeric(rad))
  (deg <- rad/(pi/180) )
}

#' @export
z_score <- function(conf.level) {
  stats::qnorm(1 - (1 - conf.level) / 2)
}

#' @export
est.kappa <- function(x, w = NULL, bias = FALSE, axial = TRUE) {
  # Default weights
  if (is.null(w)) {
    w <- rep(1, length(x))
  } else {
    w <- as.numeric(w)
  }

  f <- if(axial) 2 else 1
  x <- (x*f) %% 360

  # Remove NA pairs
  keep <- !is.na(x) & !is.na(w)
  x <- x[keep]
  w <- w[keep]

  mean.dir <- circular_mean_deg(x, w = w)
  mean_cos <- mean(cos(deg2rad(x - mean.dir)))
  kappa <- abs(A1inv(mean_cos))

  if (bias) {
    n <- sum(w)
    if (kappa < 2) {
      kappa <- max(kappa - 2 / (n * kappa), 0)
    }
    if (kappa >= 2) {
      kappa <- ((n - 1)^3 * kappa) / (n^3 + n)
    }
  }
  kappa
}

#' @export
mean_SC <- function(x, w = NULL, na.rm = TRUE) {
  stopifnot(any(is.numeric(x)), is.logical(na.rm))

  if (is.null(w)) w <- rep(1, times = length(x))

  if (na.rm) {
    keep <- !is.na(x) & !is.na(w)
    x <- x[keep]
    w <- w[keep]
  }

  x <- deg2rad(x)

  Z <- sum(w)

  sinx <- w * sin(x)
  cosx <- w * cos(x)
  # sumsin <- sum(sinx)
  # sumcos <- sum(cosx)
  # meansin <- sumsin / Z
  # meancos<- sumcos / Z
  # cbind(C = meancos, S = meansin)
  #
  # sums <- colSums(cbind(cosx, sinx))
  sums <- c(sum(cosx), sum(sinx))
  setNames(sums / Z, nm = c("C", "S"))
}

#' @export
circular_mean_rad <- function(x, w = NULL, na.rm = TRUE) {
  x <- rad2deg(x)
  m <- mean_SC(x, w, na.rm)
  meanx_rad <- atan2(m["S"], m["C"])
  unname(meanx_rad)
}

#' @export
circular_mean_deg <- function(x, w = NULL, na.rm = TRUE) {

  m <- mean_SC(x, w, na.rm)
  meanx_rad <- atan2(m["S"], m["C"])
  meanx_deg <- rad2deg(meanx_rad) %% (360)
  unname(meanx_deg)
}

# Written by Fredrik Hanslin March 2025
# Performs the Moore's modified Rayleigh test for weighted vectors.
# Moore, B. R. (1980). A modification of the Rayleigh test for vector data. Biometrika, 67(1), 175-180.
# Input: x = angles in degrees, y = resultant vector lengths

#' @export
circ_MMRT <- function(x, y){
  v <- data.frame(as.numeric(x),as.numeric(y))
  v <- v %>% arrange(y)
  X <- list()
  Y <- list()

  for (i in 1:nrow(v)) {
    X[i] = i*cos(deg2rad(v[i,1]))
    Y[i] = i*sin(deg2rad(v[i,1]))
  }

  sumx <- base::sum(unlist(X))
  sumy = base::sum(unlist(Y))
  R = sqrt(sumx^2 + sumy^2)

  R_moore <- R/nrow(v)^(3/2)
  p = exp(-3*R_moore^2);
  t <- list('p'= p, 'R_moore' = R_moore)
  return(t)
}
