#' Plot circular data
#'
#' @description
#' Easy plotting of circular data through the 'circular' package.
#'
#' @param x A vector of angles in degrees.
#' @param r A vector of r-values (vector lengths).
#' @param weight Logical. Default: TRUE, to weight your angles by r-values.
#' @param axial Logical. Default: FALSE.
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
#' @param cidist Distance from circle to confidence interval. Default: 0.2
#' @param citick Width of confidence interval edges
#' @param norm Method for normalization of vector length. "sum" = sum(weights)^2 / sum(weights^2). "n" = number of observations
#' @returns A circular plot
#' @examples
#' x <- rad2deg(rvonmises(50, pi, 1.5))
#' circ.plot(x)

#' @export
circ.plot <- function(x, r = NULL, weight = TRUE, ci = 0.95,
                      labels = c("geo", "rad", "deg"), cex = 1, sep = 0.05, pch = 16,
                      mvwidth = 0.2, mvlwd = 3, arrowwidth = 0.1, arrowlwd = 0.1,
                      shrink = 1.2, zero = pi/2, stack = TRUE, rotation = "clock",
                      col = "gray40", mv = TRUE, mvcol = "blue", cicol = "black",
                      bins = 120, cilwd = 2, main = NULL, cidist = 0.3, citick = 0.05,
                      axial = FALSE, units = c("deg", "rad"), norm = c("sum", "n")) {
  stopifnot(is.numeric(x))
  if (!is.null(r)) stopifnot(is.numeric(r))

  units <- match.arg(units)
  labels <- match.arg(labels)
  norm <- match.arg(norm)

  # helpers expected in environment: deg2rad, rad2deg, circ.stats, confidence_interval_weighted, circ.rayleigh

  # Keep x in degrees for stats, produce radian circular object for plotting
  if (units == "deg") {
    x_deg <- as.numeric(x)
    x_rad_circ <- circular(deg2rad(x_deg))
  } else {
    x_rad_circ <- circular(as.numeric(x))
    x_deg <- rad2deg(as.numeric(x))
  }

  if (is.null(r)) weight <- FALSE

  type <- if (axial) "axial" else "angles"

  if (labels == "geo") {
    lab <- c("E", "N", "W", "S")
  } else if (labels == "deg") {
    lab <- c(90, 0, 270, 180)
  } else {
    lab <- c(expression(paste(frac(pi, 2))), 0, expression(paste(frac(3*pi, 2))), expression(paste(pi)))
  }


  plot.circular(x_rad_circ, shrink = shrink, cex = cex, sep = sep, stack = stack, pch = pch,
                col = col, zero = zero, main = main, axes = FALSE, rotation = rotation, bins = bins)
  if (axial)  points.circular(x_rad_circ+pi, shrink = shrink, cex = cex, sep = sep, stack = stack, pch = pch,
                              col = col, zero = zero, main = main, rotation = rotation, bins = bins)
  axis.circular(at = circular(c(0, pi/2, pi, 3*pi/2)), labels = lab, zero = NULL, rotation = rotation)

  if (!is.null(r)) {
    arrows.circular(x_rad_circ, r, col = col, lwd = arrowlwd, length = arrowwidth, zero = zero, rotation = rotation)
    if (axial)
      arrows.circular(x_rad_circ+pi, r, col = col, lwd = arrowlwd, length = arrowwidth, zero = zero, rotation = rotation)
  }

  stats_mv <- if (!is.null(r) & isTRUE(weight)) {
    circ.stats(x = x_deg, w = r, units = "deg", type = type, norm = norm)
  } else {
    circ.stats(x = x_deg, w = NULL, units = "deg", type = type, norm = norm)
  }
  mean_deg_axis <- stats_mv$mean_orientation   # degrees; axis-scale when axial=TRUE
  vec_len <- stats_mv$vector_length

  if (!isFALSE(mv)) {
    if (axial) {
      mean_rad_axis <- circular(deg2rad(mean_deg_axis)) %% (2 * pi)
      arrows.circular(mean_rad_axis, vec_len, lwd = mvlwd, length = mvwidth, col = mvcol,
                      zero = zero, rotation = rotation)
      arrows.circular((mean_rad_axis + pi) %% (2 * pi), vec_len, lwd = mvlwd, length = mvwidth, col = mvcol,
                      zero = zero, rotation = rotation)
    } else {
      mean_rad_dir <- circular(deg2rad(mean_deg_axis %% 360))
      arrows.circular(mean_rad_dir, vec_len, lwd = mvlwd, length = mvwidth, col = mvcol,
                      zero = zero, rotation = rotation)
    }
  }

  p <- circ.rayleigh(x, r, axial = axial, nperm = 0, seed = NULL)$p.value

  if (!is.null(ci) && p < 0.05) {

    # pull back to plain degrees for the CI function
    # x_deg <- rad2deg(as.numeric(x))
    w_arg  <- if (weight) r else NULL

    ci_out <- circ.confint(
      x          = x_deg,
      conf.level = ci,
      w          = w_arg,
      axial      = axial,
      method     = "clt"     # or "ori", same as your standalone call
    )
    endpoints <- ci_out$conf.int      # [lower, upper] in degrees

    # debug print: make sure this matches your standalone call
    message("CI endpoints (deg): [",
            round(endpoints[1],1), ", ",
            round(endpoints[2],1), "]")

    # to radians
    lower_rad <- endpoints[1] * pi/180
    upper_rad <- endpoints[2] * pi/180

    # correct wrap-around: π for axial, 2π for circular
    domain <- if (axial) pi else 2*pi
    if (upper_rad < lower_rad) upper_rad <- upper_rad + domain

    # build smooth arc
    theta.seq <- seq(lower_rad, upper_rad, length.out = 200)

    # plot it outside the circle at radius “cidist”
    lines.circular(
      circular(theta.seq, units = "radians"),
      rep(cidist, length(theta.seq)),
      col      = cicol,
      lwd      = cilwd,
      zero     = zero,
      rotation = rotation
    )

    ends <- c(lower_rad, upper_rad)

    for (theta in ends) {
      lines.circular(circular(rep(theta, 2)),
                     y    = c(cidist - citick, cidist + citick),
                     col  = cicol,
                     lwd  = cilwd,
                     zero = zero,
                     rotation = rotation)
    }



    if (axial) {
      lines.circular(
        circular(theta.seq+pi, units = "radians"),
        rep(cidist, length(theta.seq)),
        col      = cicol,
        lwd      = cilwd,
        zero     = zero,
        rotation = rotation
      )
      for (theta in ends) {
        lines.circular(circular(rep(theta+pi, 2)),
                       y    = c(cidist - citick, cidist + citick),
                       col  = cicol,
                       lwd  = cilwd,
                       zero = zero,
                       rotation = rotation)
      }
    }

  } else {
    warning("Uniformity test: p>0.05. No CI produced")
  }


}


# --- helpers -----------------------------------------------------------------

#' Convert degrees to radians
#'
#' @param deg Angles in degrees
#' @export
deg2rad <- function(deg) {
  stopifnot(is.numeric(deg))
  deg * pi / 180
}

#' Convert radians to degrees
#' @param rad Angles in radians
#' @export
rad2deg <- function(rad) {
  stopifnot(is.numeric(rad))
  rad * 180 / pi
}

z_score <- function(conf.level) {
  stats::qnorm(1 - (1 - conf.level) / 2)
}



#' Calculate mean statistics
#'
#' @description
#' Calculate mean direction and vector length. Function supports individual weighting of angles as well as axial data.
#' When supplied with weights, the function calculates the weighted mean direction and vector length based on an effective sample size of n = sum(w)^2 / sum(w^2).
#' Argument norm="n", will use the number of observations as the sample size.
#'
#' @param angles Angles
#' @param weights Weights (r-vector lengths)
#' @param type Type of angles supplied. Default: "angles"
#' @param norm Method for normalization of vector length. "sum" = sum(weights)^2 / sum(weights^2). "n" = number of observations
#' @param units Unit for your angles. c("deg", "rad")

#' @export
circ.stats <- function(x, w = NULL,
                       type = c("angles", "axial"),
                       norm = c("sum", "n"),
                       units = c("deg", "rad")) {
  type <- match.arg(type)
  norm <- match.arg(norm)
  units <- match.arg(units)

  axial <- if (type=="axial") TRUE else FALSE

  if (is.null(w)) w <- rep(1, length(x))

  # Remove NA pairs up front
  keep <- !is.na(x) & !is.na(w)
  x <- x[keep]
  w <- w[keep]

  n <- length(x)
  if (n < 1) stop("At least one angle is required after removing NAs.")
  if (length(w) != n) stop("Length of w must match length of x after NA removal.")
  if (any(w < 0)) stop("All w must be non-negative.")

  f <- if (type == "axial") 2 else 1

  # Convert to radians and double for axial
  if (units == "rad") {
    radians <- x * f
  } else {
    radians <- deg2rad(x) * f
  }

  X <- sum(w * cos(radians))
  Y <- sum(w * sin(radians))

  R <- sqrt(X^2 + Y^2)
  denom <- if (norm == "sum") (sum(w)^2)/(sum(w^2))  else n
  R_bar <- R / denom

  mean_rad <- atan2(Y, X) / f

  if (type == "axial") {
    mean_rad <- mean_rad %% pi
  } else {
    mean_rad <- mean_rad %% (2 * pi)
  }

  mean_out <- if (units == "deg") rad2deg(mean_rad) else mean_rad

  list(
    type = type,
    mean_orientation = mean_out,
    vector_length = R_bar,
    R_unscaled = R,
    normalization = denom)
}



# --- A1 inverse --------------------------------------------------------------
A1inv <- function(R) {
  if (!is.numeric(R) || any(R < 0) || any(R >= 1)) {
    if (any(R == 1)) return(Inf)
    stop("R must be in [0,1)")
  }
  kappa <- numeric(length(R))
  for (i in seq_along(R)) {
    r <- R[i]
    if (r < 0.53) {
      kappa[i] <- 2 * r + r^3 + (5 * r^5) / 6
    } else if (r < 0.85) {
      kappa[i] <- -0.4 + 1.39 * r + 0.43 / (1 - r)
    } else {
      kappa[i] <- 1 / (r^3 - 4 * r^2 + 3 * r)
    }
  }
  kappa
}

# --- circular mean wrapper using circ.stats ---------------------------------
circular_mean <- function(x, w = NULL, axial = FALSE, units = "deg", na.rm = TRUE) {
  type <- if (axial) "axial" else "angles"
  stats <- circ.stats(x = x, w = w, type = type, norm = "sum", units = units)
  stats$mean_orientation
}

# --- estimate kappa ----------------------------------------------
est.kappa <- function(x, w = NULL, bias = TRUE, axial = FALSE) {
  if (is.null(w)) w <- rep(1, length(x))
  keep <- !is.na(x) & !is.na(w)
  x <- x[keep]
  w <- as.numeric(w[keep])
  if (length(x) == 0) stop("No data after removing NAs")
  if (any(w < 0)) stop("w must be non-negative")

  f <- if (axial) 2 else 1
  type <- if (axial) "axial" else "angles"
  R_bar <- circ.stats(x, w=w, type=type)$vector_length
  kappa <- abs(A1inv(R_bar))

  if (bias) {
    n_eff <- (sum(w)^2)/(sum(w^2))
    if (kappa < 2) {
      kappa <- max(kappa - 2 / (n_eff * kappa), 0)
    } else {
      kappa <- ((n_eff - 1)^3 * kappa) / (n_eff^3 + n_eff)
    }
  }
  kappa
}



#' Calculate circular confidence interval
#'
#' @description
#' Calculates the circular confidence intervals. Accepts weighted, and axial data.
#'
#' @param x Numeric vector of angles (in degrees)
#' @param conf.level Confidence level (e.g. 0.95)
#' @param w Optional weights
#' @param axial Logical. if TRUE treat data as axial (0–180°)
#' @param na.rm Logical. if TRUE remove NA pairs
#' @param method “clt” or “ori”. Calculate the confidence interval based on the central limit theorem ("clt", standard) or the small-angle formula ("ori") that Oriana uses. "ori" does not use the concentration parameter k, and instead assumes that the data are tightly clustered. "ori" can work for low n-sizes that are highly clustered.
#' @param bias Bias correction of the concentration parameter for low sample sizes, usually n<30. Only applicable when method = "clt".

#' @export
circ.confint <- function(
    x,
    w          = NULL,
    conf.level = 0.95,
    axial      = FALSE,
    na.rm      = TRUE,
    method     = c("clt", "ori"),
    bias       = TRUE
) {
  method <- match.arg(method)

  # 1. compute half-width (deg)
  half_width <- confidence.angle(
    x         = x,
    conf.level= conf.level,
    w         = w,
    axial     = axial,
    na.rm     = na.rm,
    method    = method,
    bias      = bias
  )

  # 2. compute mean direction (deg)
  mean_deg <- circular_mean(
    x      = x,
    w      = w,
    axial  = axial,
    units  = "deg",
    na.rm  = na.rm
  )

  # 3. compute endpoints and wrap to [0,360) or [0,180)
  mod_base <- if (axial) 180 else 360
  lower <- (mean_deg - half_width) %% mod_base
  upper <- (mean_deg + half_width) %% mod_base

  list(
    mean       = mean_deg,
    half_width = half_width,
    conf.int   = c(lower = lower, upper = upper)
  )
}


confidence.angle <- function(
    x,
    w          = NULL,
    conf.level = 0.95,
    axial      = FALSE,
    na.rm      = TRUE,
    method     = c("clt", "ori"),
    bias       = TRUE
) {
  method <- match.arg(method)

  # 1) handle weights & NAs
  if (is.null(w)) w <- rep(1, length(x))
  df <- data.frame(x = x, w = w)
  if (na.rm) df <- df[stats::complete.cases(df), ]
  x <- df$x;  w <- df$w
  stopifnot(length(x) > 0, all(w >= 0))

  # 2) get R̄ from circ.stats()
  type  <- if (axial) "axial" else "angles"
  stats <- circ.stats(x = x, w = w, type = type,
                      norm   = "sum", units = "deg")
  R_bar <- stats$vector_length

  # 3) true effective N
  sum_w  <- sum(w)
  sum_w2 <- sum(w^2)
  n_eff  <- if (sum_w2 > 0) sum_w^2 / sum_w2 else 0

  # 4) critical z
  Z_alpha <- z_score(conf.level)

  # 5) compute half‐width in radians
  half_r <- switch(
    method,

    # ORIANA small‐angle formula
    ori = {
      if (n_eff <= 0 || R_bar <= 0) {
        Inf
      } else {
        se_ori <- sqrt((1 - R_bar) / (2 * n_eff * R_bar))
        hw     <- Z_alpha * se_ori
        if (axial) hw <- hw / 2
        hw
      }
    },

    # CLT normal‐approximation (unbounded)
    clt = {
      kappa <- est.kappa(x, w = w, axial = axial, bias=bias)
      if (!is.finite(kappa) || kappa <= 0 || R_bar <= 0 || n_eff <= 0) {
        Inf
      } else {
        sde <- 1 / sqrt(n_eff * R_bar * kappa)
        hw  <- Z_alpha * sde
        if (axial) hw <- hw / 2
        hw
      }
    }
  )

  # 6) convert to degrees
  rad2deg(half_r)
}




#' Rayleigh test
#'
#' @description
#' Conducts a test of directedness based on inputs. If no weights and axial=F, then the regular Rayleigh test is conducted.
#' If supplied with weights, the Moore's modified rayleigh test is conducted.
#' If axial=T, angles are doubled, before conducting test based on weights.

#' @param x       Numeric vector of angles (in degrees)
#' @param w       Optional weights
#' @param axial   Logical. if TRUE treat data as axial (0–180°)

#' @export
circ.rayleigh <- function(x, w = NULL, axial = FALSE, nperm = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.numeric(x)) stop("x must be numeric (degrees).")
  n <- length(x)
  if (n < 2) stop("need at least two observations")

  x <- na.omit(x)
  w <- na.omit(w)

  nperm <- if (!is.null(w)) 10000 else 0

  # normalize angles to [0,360)
  x <- (x %% 360)

  # Working scale: double for axial data
  if (!is.logical(axial) || length(axial) != 1) stop("axial must be a single logical value")
  if (axial) {
    work_deg <- (2 * x) %% 360
  } else {
    work_deg <- x
  }
  ang_rad <- work_deg * pi / 180

  out <- list(n = n, axial = axial)

  # Case 1: No w provided => standard Rayleigh on working scale
  if (is.null(w)) {
    out$method <- "Rayleigh"
    S <- sum(exp(1i * ang_rad))
    R_abs <- Mod(S)
    Rbar <- R_abs / n
    Z <- n * Rbar^2

    out$R <- R_abs
    out$Rbar <- Rbar
    out$statistic <- Z

    # analytic p-value (improved large-sample approximation)
    p_analytic <- exp(-Z) * (1 + (2*Z - Z^2) / (4*n) -
                               (24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4) / (288 * n^2))
    p_analytic <- max(min(p_analytic, 1), 0)
    out$p_analytic <- p_analytic

    # permutation p-value if requested
    if (!is.numeric(nperm) || length(nperm) != 1 || nperm < 0) stop("nperm must be a nonnegative integer")
    nperm <- as.integer(nperm)
    if (nperm > 0) {
      sim_stats <- numeric(nperm)
      for (k in seq_len(nperm)) {
        sim_ang <- runif(n, 0, 2*pi)
        sim_S <- sum(exp(1i * sim_ang))
        sim_R <- Mod(sim_S)
        sim_stats[k] <- (sim_R^2) / n
      }
      p_perm <- (sum(sim_stats >= Z) + 1) / (nperm + 1)
      out$p_permutation <- p_perm
      out$nperm <- nperm
      if (nperm <= 2000) out$perm_stats <- sim_stats
      out$p.value <- p_perm
    } else {
      out$p.value <- p_analytic
    }

    return(out)
  }

  # Case 2: w provided => Moore's modified Rayleigh procedure
  if (length(w) != n) stop("x and w must have same length")
  if (!is.numeric(w)) stop("w must be numeric when provided")

  out$method <- "Moore modified Rayleigh"
  df <- data.frame(angle = work_deg, w = as.numeric(w))
  df <- df[order(df$w), , drop = FALSE]

  i_idx <- seq_len(n)
  ang_rad_ordered <- df$angle * pi / 180
  sumx <- sum(i_idx * cos(ang_rad_ordered))
  sumw <- sum(i_idx * sin(ang_rad_ordered))
  R <- sqrt(sumx^2 + sumw^2)
  R_moore <- R / (n^(3/2))

  out$R <- R
  out$R_moore <- R_moore

  # analytic p-value (Moore)
  p_analytic <- exp(-3 * R_moore^2)
  p_analytic <- max(min(p_analytic, 1), 0)
  out$p_analytic <- p_analytic

  # permutation p-value if requested (permute angles relative to w)
  if (!is.numeric(nperm) || length(nperm) != 1 || nperm < 0) stop("nperm must be a nonnegative integer")
  nperm <- as.integer(nperm)
  if (nperm > 0) {
    perm_stats <- numeric(nperm)
    for (k in seq_len(nperm)) {
      perm_angles <- sample(work_deg, size = n, replace = FALSE)
      ang_rad_p <- perm_angles * pi / 180
      sumx_p <- sum(i_idx * cos(ang_rad_p))
      sumw_p <- sum(i_idx * sin(ang_rad_p))
      R_p <- sqrt(sumx_p^2 + sumw_p^2)
      perm_stats[k] <- (R_p / (n^(3/2)))
    }
    p_perm <- (sum(perm_stats >= R_moore) + 1) / (nperm + 1)
    out$p_permutation <- p_perm
    out$nperm <- nperm
    out$p.value <- p_analytic
  }

  return(out)
}



#' Circular data structure
#' @description
#' Uses .txt outputs from a US Digital Encoder to calculate mean statistics for group analysis.
#' Calculates mean angles, vector lengths and computes the rayleigh test for all .txt files in a folder.
#' @param path Directory path for the folder to analyse. If left empty, a pop up appears to browse for the folder
#' @param channel Which encoder channel your data is stored in (1:4)
#' @param phases How many phases in your experiment
#' @param phaselength The length of each phase. In seconds
#' @param axial List of phases that should be analysed axially. E.g. c(1,2) analyses phase 1 and 2 axially
#' @param samplingrate Sampling rate set on your encoder. In milliseconds. Default: 200
#' @param all Logical. Calculate statistics for all .txt files, or only the ones that finish your set number of phases.
#' @param names Logical. Include specific info from .txt file name. If TRUE, file name should be named with the following info: "EXP_DATE_ID.txt". Lengths of Exp, date and ID can differ, but must be separated with "_".

#' @export
get.mean.data <- function(path = NULL, channel = NULL, phases=NULL, phaselength=NULL, axial=NULL, samplingrate=200,  all=TRUE, names=F) {
  # Set the path to your folder
  if (is.null(path)) {
    path <- choose.dir("Choose data directory")
  }



  # List all .txt files in the folder
  txt_files <- list.files(path = path,
                          pattern = "\\.txt$",
                          full.names = TRUE)

  # data input
  samplingrate <- samplingrate  # The sampling rate set in the encoder software
  resolution <- 1000 / samplingrate  # in Hz
  phaselength <- phaselength  # Change if necessary
  phases <- phases  # Change if necessary

  # Pre‐allocate one row per file
  N <- length(txt_files)

  if (phases > 1) {
    changesdata <- as.data.frame(matrix(NA,
                                        nrow = N,
                                        ncol = phases - 1))
    colnames(changesdata) <- paste0("Change", 1:(phases - 1))

  }
  else changesdata <- NULL

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
  } else {
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
      } else {
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

      # Convert axial data
      circmeans    <- matrix(nrow = 1, ncol = numphase)
      circr        <- matrix(nrow = 1, ncol = numphase)
      circrayleigh <- matrix(nrow = 1, ncol = numphase)


      for (i in seq_len(numphase)) {
        type <- if (i %in% axial) "axial" else "angles"
        stats     <- circ.stats(
          sorted_angles[, i],
          type = type,
          unit = "deg"
        )
        circmeans[i]    <- stats$mean_orientation
        circr[i]        <- stats$vector_length
        circrayleigh[i] <- if (type == "axial") {
          circ.rayleigh(sorted_angles[, i], axial = TRUE)$p.value
        } else {
          circ.rayleigh(sorted_angles[, i])$p.value
        }
      }

      # compute changes between phases
      if (phases>1) {
        changes <- numeric(numphase - 1)
        for (i in seq_len(numphase - 1)) {
          a1  <- circmeans[i]
          a2  <- circmeans[i + 1]
          ax1 <- i     %in% axial
          ax2 <- (i+1) %in% axial

          key <- paste(as.integer(ax1), as.integer(ax2), sep = "_")
          changes[i] <- switch( # Switch between the differenct scenarios "key" for each phase
            key,
            # unimodal → unimodal
            "0_0" = (a2-a1) %% 360,
            # axial → axial
            "1_1" = abs(angleDiff(a1*2, a2*2, 360) / 2),
            # unimodal → axial
            "0_1" = {
              d1 <- angleDiff(a1, a2, 360)
              d2 <- angleDiff(a1, a2 + 180, 360)
              abs(if (abs(d1) <= abs(d2)) d1 else d2)
            },
            # axial → unimodal
            "1_0" = {
              d1 <- angleDiff(a1, a2, 360)
              d2 <- angleDiff(a1 + 180, a2, 360)
              abs(if (abs(d1) <= abs(d2)) d1 else d2)
            }
          )
        }
      }
      else changes <- NA



      # Append individual info  to dataframes

      if (isTRUE(names)) {
        infodata[k,1:5] <- infos
      } else{
        infodata[k,1:3] <- infos
      }
      meandata[k,seq_len(numphase)] <- circmeans
      rdata[k,seq_len(numphase)] <- circr
      rayleighdata[k,seq_len(numphase)] <- circrayleigh
      if (phases > 1 && numphase > 1) {
        changesdata[k,1:(numphase - 1)] <- changes
      }
    }

  }



  # Combine into one final data.frame

  pieces <- list(infodata, meandata, rdata, changesdata, rayleighdata)
  pieces <- Filter(Negate(is.null), pieces)
  data <- do.call(cbind, pieces)


  if (isFALSE(all)){
    data <- na.omit(data)
  }
  data

}


# helper: signed shortest‐path on [0,mod)
angleDiff <- function(a1, a2, mod) {
  d <- (a2 - a1) %% mod
  if (d > mod/2) d <- d - mod
  d
}

