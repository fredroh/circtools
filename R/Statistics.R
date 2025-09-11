circ.plot <- function(x, r=NULL, weight=TRUE, ci=0.95,cords="geo", sep=0.05, pch=16, mvlength=0.2, mvlwd=3, arrowlength=0.1, arrowlwd=0.1, shrink=1.2, zero=pi/2, stack=T, rotation="clock", col="gray40", mv=TRUE, mvcol="blue", cicol="black", cilwd=2, main=NULL, dist=0.2, citick=0.05) {
  x <- NISTdegTOradian(x)
  
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
      arrows.circular(mean_circ(x), circ_rho_zar(x, r), lwd=mvlwd, length=mvlength, col=mvcol, zero=zero, rotation=rotation)
    }
    
    else
      if (!isFALSE(mv)){
        {  arrows.circular(mean_circ(x), rho.circular(x), lwd=mvlwd, length=mvlength, col=mvcol, zero=zero, rotation=rotation)
        }
      }
  }
  
  if (!is.null(ci)) {
    if (isTRUE(weight)) {    
      confint <- confidence_interval_weighted(NISTradianTOdeg(x),w=r, conf.level = ci, axial=F)
    }
    else {
      confint <- confidence_interval_weighted(NISTradianTOdeg(x),w=NULL, conf.level = ci, axial=F)
      
    }
    ci_low <- NISTdegTOradian(confint$conf.interval[1]%%360)
    ci_high <- NISTdegTOradian(confint$conf.interval[2]%%360)
    
    angles <- if(ci_low < ci_high) {
      seq(ci_low, ci_high, length=200)
    } else {
      c(seq(ci_low, 2*pi, length=100),
        seq(0, ci_high, length=100))
    }
    
    lines.circular(circular(angles), y=rep(dist, length(angles)),
                   col=cicol, lwd=cilwd, zero=zero, axes=F, rotation=rotation) 
    # endpoint ticks
    for(j in c(ci_low, ci_high)) {
      lines.circular(circular(rep(j,2)), 
                     y = c(dist - citick, dist + citick),
                     col=cicol, lwd=cilwd, zero=zero, rotation=rotation)
    }
    
  } 
}

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
  NISTradianTOdeg(asin(temp) * f)
}

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
  R <- circ_rho_zar(NISTdegTOradian(x), r = w)
  1/sqrt(n * R * kappa)
}

confidence_interval_weighted <- function (x, conf.level = 0.95, w = NULL, axial = FALSE, na.rm = TRUE) 
{
  conf.angle <- confidence_angle_weighted(x, conf.level, w, axial, 
                                          na.rm)
  mu <- circular_mean(x, w = w, axial = axial, na.rm = na.rm)
  list(mu = mu, conf.angle = conf.angle, conf.interval = c(mu - 
                                                             conf.angle, mu + conf.angle))
}
