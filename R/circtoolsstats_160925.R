#' Plot circular data
#'
#' @description
#' Easy plotting of circular data through the 'circular' package.
#'
#' @param x A vector of angles in degrees.
#' @param r A vector of r-values (vector lengths).
#' @param weight Logical. Default: TRUE, to weight your angles by r-values.
#' @param ci The confidence level to calculate. Default: 0.95. ci=NULL to not plot confidence interval.
#' @param labels Plot your data with geographical or numeric labels. c("geo", "deg", "rad").
#' @param sep Distance between dots.
#' @param pch Point charactter to use.
#' @param mvwidth Width of mean vector.
#' @param mvlwd Line width of mean vector.
#' @param arrowwidth Width of individual vectors
#' @param arrowlwd Line width of individual vectors
#' @param shrink Parameter that controls the size of the plotted circle. Default is 1.2.
#' @param zero Location of Zero. Default = pi/2.
#' @param stack Logical. To stack datapoints or not
#' @param rotation Rotation of the circle. c("clock", "counter"). Default: "clock"
#' @param col Color of the points and individual vectors.
#' @param mv Logical. To plot mean vector or not.
#' @param mvcol Color of mean vector
#' @param cicol Color of confidence interval
#' @param bins Number of bins to partition the circle
#' @param cilwd Linewidth of Confidence interval
#' @param main Plot title
#' @param dist Distance from circle to confidence interval. Default: 0.2
#' @param citick Width of confidence interval edges
#' @returns A circular plot
#' @examples
#' x <- rad2deg(rvonmises(50, pi, 1.5))
#' circ.plot(x)


#' @export
circ.plot <- function(x, r=NULL, weight=TRUE, ci=0.95,labels="geo", sep=0.05, pch=16, mvwidth=0.2,
                      mvlwd=3, arrowwidth=0.1, arrowlwd=0.1, shrink=1.2, zero=pi/2, stack=T,
                      rotation="clock", col="gray40", mv=TRUE, mvcol="blue", cicol="black", bins=120,
                      cilwd=2, main=NULL, dist=0.2, citick=0.05) {
  stopifnot(is.numeric(x))

  if (!is.null(r)) {
    stopifnot(is.numeric(r))
  }
  x <- circular(deg2rad(x))

  if (is.null(r)) {
    weight <- FALSE
  }



  if (labels=="geo") {
    labels=c("E", "N", "W", "S")}
  else if(labels=="deg"){
    labels=c(90,0,270,180)
  }
  else if(labels=="rad") {
    labels = c(expression(paste(frac(pi, 2))),
               0,
               expression(paste(frac(3*pi,2))),
               expression(paste(pi))
               )
  }


  plot.circular(x, shrink=shrink, sep=sep, stack=stack, pch=pch, col=col, zero=zero,
                main=main, axes=F, rotation=rotation, bins=bins)
  axis.circular(at=circular(c(0, pi/2, pi, 3*pi/2)), labels=labels, zero=NULL, rotation=rotation)
  if (!is.null(r)) {
    arrows.circular(x, r, col=col, lwd=arrowlwd, length=arrowwidth, zero=zero,  rotation=rotation)
    if (!isFALSE(mv) & isTRUE(weight)) {
      arrows.circular(circular(circular_mean_rad(x, r)), circ_rho_zar(x, r), lwd=mvlwd, length=mvwidth, col=mvcol, zero=zero, rotation=rotation)
    }

    else
      if (!isFALSE(mv)){
        {  arrows.circular(circular(circular_mean_rad(x)), rho.circular(x), lwd=mvlwd, length=mvwidth, col=mvcol, zero=zero, rotation=rotation)
        }
      }
  }

  if (isTRUE(weight)) {
    c <- circ_MMRT(rad2deg(x),r)
    p <- c[[1]]
  }
  else {
    a <- rayleigh.test(x)
    p <- a[[2]]
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
                     col=cicol, lwd=cilwd, zero=zero, rotation=rotation)
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


#' Calculate a weighted mean vector length
#'
#' @param a A vector of angles in radians
#' @param r A vector of r-values (vector lengths)
#'
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

#' Calculate the weighted confidence intervals
#' @description
#' Calculate the confidence intervals based on the weighted mean vector length.
#'
#' @param x A vector of angles in degrees
#' @param w A vector of r-values (vector lengths)
#' @param conf.level The confidence level to calculate. Default: 0.95
#' @export
confidence_interval_weighted <- function (x, conf.level = 0.95, w = NULL, axial = FALSE, na.rm = TRUE)
{
  conf.angle <- confidence_angle_weighted(x, conf.level, w, axial,
                                          na.rm)
  mu <- circular_mean_deg(x, w = w)
  list(mu = mu, conf.angle = conf.angle, conf.interval = c(mu -
                                                             conf.angle, mu + conf.angle))
}

#' Convert degrees to radians
#' @param deg Vector of angles in degrees
#' @export
deg2rad <- function(deg) {
  stopifnot(is.numeric(deg))
  (rad <- (pi/180)*deg )
}

#' Convert radians to degrees
#' @param deg Vector of angles in radians
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
circular_mean <- function(x, w = NULL, na.rm = TRUE, units=NULL) {

  if (missing(units) || !units %in% c("deg", "rad"))
    stop("`units` must be 'deg' or 'rad'", call. = FALSE)

  if (units=="rad") {
    x <- rad2deg(x)
    m <- mean_SC(x, w, na.rm)
    meanx_rad <- atan2(m["S"], m["C"])
    unname(meanx_rad)
  }
  else {
    m <- mean_SC(x, w, na.rm)
    meanx_rad <- atan2(m["S"], m["C"])
    meanx_deg <- rad2deg(meanx_rad) %% (360)
    unname(meanx_deg)
  }
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


#' Moores modified Rayleigh Test
#' @description
#' Conducts the Moores modified Rayleigh test for weighted data.
#' @param x Angles in degrees
#' @param y r-values (vector lengths)

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

#' Circular data structure
#' @description
#' Calculates mean angles, vector lengths and computes the rayleigh test for all .txt files in a folder.
#' @param channel Which encoder channel your data is stored in (1:4)
#' @param phases How many phases in your experiment
#' @param phaselength The length of each phase. In seconds
#' @param samplingrate Sampling rate set on your encoder. In milliseconds. Default: 200
#' @param all Logical. Calculate statistics for all .txt files, or only the ones that finish your set number of phases.
#' @param names Logical. Include specific info from .txt file name. If TRUE, file name should be named with the following info: "EXP_DATE_ID.txt". Lengths of Exp, date and ID can differ, but must be separated with "_".


#' @export
get.mean.data <- function(channel = NULL, phases=NULL, phaselength=NULL, samplingrate=200, all=TRUE, names=F) {
  # Set the path to your folder
  folder_path <- choose.dir("Choose data directory")

  # List all .txt files in the folder
  txt_files <- list.files(path = folder_path,
                          pattern = "\\.txt$",
                          full.names = TRUE)

  # data input
  samplingrate <- samplingrate  # The sampling rate set in the encoder software
  resolution <- 1000 / samplingrate  # in Hz
  phaselength <- phaselength  # Change if necessary
  phases <- phases  # Change if necessary

  # Pre‐allocate one row per file
  N <- length(txt_files)

  changesdata <- as.data.frame(matrix(NA,
                                      nrow = N,
                                      ncol = phases - 1))
  colnames(changesdata) <- paste0("Change", 1:(phases - 1))

  meandata <- as.data.frame(matrix(NA,
                                   nrow = N,
                                   ncol = phases))
  colnames(meandata) <- paste0("Phase", 1:phases, ".A")

  rdata <- as.data.frame(matrix(NA,
                                nrow = N,
                                ncol = phases))
  colnames(rdata) <- paste0("Phase", 1:phases, ".L")

  rayleighdata <- as.data.frame(matrix(NA,
                                       nrow = N,
                                       ncol = phases))
  colnames(rayleighdata) <- paste0("Rayleigh", 1:phases)

  if (isTRUE(names)) {
    infodata <- as.data.frame(matrix(NA,
                                     nrow = N,
                                     ncol = 5),
                              stringsAsFactors = FALSE)
    colnames(infodata) <- c("ID", "Num", "Exp", "Date", "Time")
  }
  else {
    infodata <- as.data.frame(matrix(NA,
                                     nrow = N,
                                     ncol = 3),
                              stringsAsFactors = FALSE)
    colnames(infodata) <- c("ID", "Date", "Time")

  }



  # Iterate over each .txt file by index k
  for (k in seq_along(txt_files)) {
    file   <- txt_files[k]
    length <- length(readLines(file))
    name   <- basename(file)

    max_vals <- phases * phaselength * resolution


    numphase <- floor(length / (phaselength * resolution))
    if (numphase > phases) {
      numphase <- phases
    }

    # Skip the files that are not complete if 'all' = F
    if ((isFALSE(all) & numphase==phases) | isTRUE(all)){
      info <- file.info(file)
      Date <- format(info$mtime, "%d.%m.%Y")
      Time <- format(info$mtime, "%H:%M")

      ID  <- tail(strsplit(name, "\\.")[[1]])[[1]]
      Exp <- tail(strsplit(name, "_|\\.")[[1]], 4)[1]
      Num <- tail(strsplit(name, "_|\\.|-")[[1]])[[3]]

      if (isTRUE(names)) {
        infos <- c(ID, Num, Exp, Date, Time)
      }
      else {
        infos <- c(ID, Date, Time)
      }


      # Read data
      angles <- read.delim(file, header = FALSE, sep = "")[[channel]] * 3

      if (length(angles) > max_vals) {
        angles <- angles[1:max_vals]
        cat("In file", name,
            ". File exceeds size for the set parameters.",
            "Only first", phases, "phases are included\n")
      }

      # Divide data into phases
      sorted_angles <- as.data.frame(
        matrix(angles,
               nrow = phaselength * resolution,
               ncol = numphase)
      )

      circmeans    <- matrix(nrow = 1, ncol = numphase)
      circr        <- matrix(nrow = 1, ncol = numphase)
      circrayleigh <- matrix(nrow = 1, ncol = numphase)

      for (i in seq_len(numphase)) {
        circmeans[1, i]    <- circular_mean(sorted_angles[, i], unit = "deg") %% 360
        circr[1, i]        <- rho.circular(circular(deg2rad(sorted_angles[, i])))
        circrayleigh[1, i] <- rayleigh.test(circular(deg2rad(sorted_angles[, i])))[[2]]
      }

      changes <- matrix(nrow = 1, ncol = numphase - 1)
      for (i in seq_len(numphase - 1)) {
        delta <- (circmeans[, i + 1] - circmeans[, i]) %% 360
        changes[1, i] <- if (delta > 180) abs(delta - 360) else delta
      }


      # Append individual info  to dataframes
      infodata[k,1:5] <- infos
      meandata[k,1:numphase] <- circmeans
      rdata[k,1:numphase] <- circr
      rayleighdata[k,1:numphase] <- circrayleigh
      changesdata[k,1:(numphase - 1)] <- changes
    }

  }



  # Combine into one final data.frame
  data <- cbind(infodata,
                meandata,
                rdata,
                changesdata,
                rayleighdata)

  if (isFALSE(all)){
    data <- na.omit(data)
  }
  data

}








# From Philipp Berens
circ_axial <- function(alphas, m = 1, dim = NULL) {

  # default 'm' if missing or NULL
  if (missing(m) || is.null(m)) {
    m <- 1
  }

  # choose the first non-singleton dimension if 'dim' not provided
  if (is.null(dim)) {
    d <- dim(alphas)
    if (is.null(d)) {
      dim <- 1
    } else {
      dim <- which(d > 1)[1]
      if (is.na(dim)) {
        dim <- 1
      }
    }
  }

  # compute mean of exp(1i * alphas * m) along the specified dimension
  if (is.null(dim(alphas))) {
    z_bar_m <- mean(exp(1i * alphas * m))
  } else {
    margins  <- seq_along(dim(alphas))[-dim]
    z_bar_m  <- apply(exp(1i * alphas * m), margins, mean)
  }

  # extract resultant length and angle
  r  <- Mod(z_bar_m)
  mu <- Arg(z_bar_m) / m

  list(r = r, mu = mu)
}


## From Philipp Berens
circ_raotest <- function(alpha) {
  # Convert to degrees and sort
  alpha_deg <- sort(alpha * 180 / pi)
  n <- length(alpha_deg)

  # Compute Rao's U statistic
  lambda <- 360 / n
  diffs  <- diff(alpha_deg)
  U_raw  <- sum(abs(diffs - lambda))
  # wrap‐around gap
  wrap   <- 360 - alpha_deg[n] + alpha_deg[1]
  U      <- 0.5 * (U_raw + abs(wrap - lambda))

  # Lookup function for p and critical Uc
  getVal <- function(N, U) {
    # significance levels
    alphas <- c(0.001, 0.01, 0.05, 0.10)
    # Table II from Russell & Levitin (1995)
    tbl <- matrix(c(
      4, 247.32, 221.14, 186.45, 168.02,
      5, 245.19, 211.93, 183.44, 168.66,
      6, 236.81, 206.79, 180.65, 166.30,
      7, 229.46, 202.55, 177.83, 165.05,
      8, 224.41, 198.46, 175.68, 163.56,
      9, 219.52, 195.27, 173.68, 162.36,
      10, 215.44, 192.37, 171.98, 161.23,
      11, 211.87, 189.88, 170.45, 160.24,
      12, 208.69, 187.66, 169.09, 159.33,
      13, 205.87, 185.68, 167.87, 158.50,
      14, 203.33, 183.90, 166.76, 157.75,
      15, 201.04, 182.28, 165.75, 157.06,
      16, 198.96, 180.81, 164.83, 156.43,
      17, 197.05, 179.46, 163.98, 155.84,
      18, 195.29, 178.22, 163.20, 155.29,
      19, 193.67, 177.08, 162.47, 154.78,
      20, 192.17, 176.01, 161.79, 154.31,
      21, 190.78, 175.02, 161.16, 153.86,
      22, 189.47, 174.10, 160.56, 153.44,
      23, 188.25, 173.23, 160.01, 153.05,
      24, 187.11, 172.41, 159.48, 152.68,
      25, 186.03, 171.64, 158.99, 152.32,
      26, 185.01, 170.92, 158.52, 151.99,
      27, 184.05, 170.23, 158.07, 151.67,
      28, 183.14, 169.58, 157.65, 151.37,
      29, 182.28, 168.96, 157.25, 151.08,
      30, 181.45, 168.38, 156.87, 150.80,
      35, 177.88, 165.81, 155.19, 149.59,
      40, 174.99, 163.73, 153.82, 148.60,
      45, 172.58, 162.00, 152.68, 147.76,
      50, 170.54, 160.53, 151.70, 147.05,
      75, 163.60, 155.49, 148.34, 144.56,
      100, 159.45, 152.46, 146.29, 143.03,
      150, 154.51, 148.84, 143.83, 141.18,
      200, 151.56, 146.67, 142.35, 140.06,
      300, 148.06, 144.09, 140.57, 138.71,
      400, 145.96, 142.54, 139.50, 137.89,
      500, 144.54, 141.48, 138.77, 137.33,
      600, 143.48, 140.70, 138.23, 136.91,
      700, 142.66, 140.09, 137.80, 136.59,
      800, 142.00, 139.60, 137.46, 136.33,
      900, 141.45, 139.19, 137.18, 136.11,
      1000, 140.99, 138.84, 136.94, 135.92
    ), ncol = 5, byrow = TRUE)
    # pick first row with table N >= sample size
    ridx <- which(tbl[,1] >= N)[1]
    # find first column where critical < observed U
    cidx <- which(tbl[ridx, -1] < U)[1]

    if (!is.na(cidx)) {
      UC <- tbl[ridx, cidx + 1]
      p  <- alphas[cidx]
    } else {
      UC <- tbl[ridx, 4]  # critical at alpha = 0.05
      p  <- 0.5
    }
    list(p = p, UC = UC)
  }

  val <- getVal(n, U)
  list(p = val$p, U = U, UC = val$UC)
}

mean_axial_weighted <- function(theta, r) {
  if (length(theta) != length(r)) {
    stop("Lengths of 'theta' and 'r' must match")
  }

  # normalize weights
  w <- r / sum(r)

  # double the angles
  theta2 <- 2 * theta

  # weighted vector sum
  x2 <- sum(w * cos(theta2))
  y2 <- sum(w * sin(theta2))

  R2 <- sqrt(x2^2 + y2^2)

  # mean direction in doubled space
  mean2 <- atan2(y2, x2)

  # halve and wrap into [0, pi)
  mean_axial <- (mean2 / 2) %% pi

  list(theta=rad2deg(mean_axial), r=R2)
}

