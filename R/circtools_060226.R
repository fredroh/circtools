#' Plot circular data
#'
#' @description
#' Easy plotting of circular data through the 'circular' package.
#'
#' @param x A vector of angles in degrees.
#' @param r A vector of r-values (vector lengths).
#' @param weight Logical. Default: TRUE, to weight your angles by r-values for mean vector.
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
#' @param method “clt” or “ori”. Calculate the confidence interval based on the central limit theorem ("clt", standard) or the small-angle formula ("ori") that Oriana uses. Calculations based on circ_confmean in MATLAB, so not exactly as Oriana as these formulas are unknown..
#' @returns A circular plot
#' @examples
#' x <- rad2deg(rvonmises(50, pi, 1.5))
#' circ.plot(x)

#' @export
circ.plot <- function(x, r = NULL, weight = FALSE, ci = 0.95,
                      labels = c("geo", "rad", "deg"), cex = 1, sep = 0.05, pch = 16,
                      mvwidth = 0.2, mvlwd = 3, arrowwidth = 0.1, arrowlwd = 0.1,
                      shrink = 1.2, zero = pi/2, stack = TRUE, rotation = "clock",
                      col = "gray40", mv = TRUE, mvcol = "blue", cicol = "black",
                      bins = 120, cilwd = 2, main = NULL, cidist = 0.3, citick = 0.05,
                      axial = FALSE, units = c("deg", "rad"), norm = c("sum", "n"),
                      method=c("ori", "clt")) {
  stopifnot(is.numeric(x))
  if (!is.null(r)) stopifnot(is.numeric(r))

  units <- match.arg(units)
  labels <- match.arg(labels)
  norm <- match.arg(norm)
  method <- match.arg(method)

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
      method     = method     # or "ori", same as your standalone call
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
    method     = c("ori", "clt"),
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
  # If return(Inf)
  if (!is.finite(half_width)) {
    stop(
      "Confidence interval could not be computed: data are too dispersed for the chosen method.\n",
      "Zar's requirement for defining a confidence interval was not met.\n",
      "Try: method='clt' or lowering confidence level."
    )
  }

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
    method     = c("ori", "clt"),
    bias       = TRUE
) {
  method <- match.arg(method)

  # handle weights & NAs
  if (is.null(w)) w <- rep(1, length(x))
  df <- data.frame(x = x, w = w)
  if (na.rm) df <- df[stats::complete.cases(df), ]
  x <- df$x;  w <- df$w
  stopifnot(length(x) > 0, all(w >= 0))

  # R_bar from circ.stats()
  type  <- if (axial) "axial" else "angles"
  stats <- circ.stats(x = x, w = w, type = type,
                      norm   = "sum", units = "deg")
  R_bar <- stats$vector_length

  # true effective N
  sum_w  <- sum(w)
  sum_w2 <- sum(w^2)
  n_eff  <- if (sum_w2 > 0) sum_w^2 / sum_w2 else 0

  # critical z
  Z_alpha <- z_score(conf.level)

  # half‐width in radians
  half_r <- switch(
    method,

    ori = {
      n_ori <- length(x)
      R_ori <- stats$R_unscaled         
      Rbar  <- if (n_ori > 0) R_ori / n_ori else 0
      
      # 95% small‑angle is only valid for very high concentration
      if (Rbar > 0.9) {
        # Oriana small-angle SE (radians)
        SE_ori <- sqrt((1 - Rbar) / R_ori)
        hw <- Z_alpha * SE_ori
        if (axial) hw <- hw / 2
        return(hw)                      
      }
      
      # (Zar-style) CI used by Oriana/CircStat MATLAB
      # Use chi-square critical value with 1 df (not F)
      # circ_confmean uses c2 = chi2inv(1 - xi, 1); with xi = 1 - conf.level
      c2 <- stats::qchisq(conf.level, df = 1)
      
      # Requirement (Zar): rbar must exceed sqrt(c2 / (2n)), otherwise no finite CI
      # (CircStat returns NaN; we signal Inf so circ.confint() can fallback if desired)
      thresh <- sqrt(c2 / (2 * n_ori))
      if (!(is.finite(Rbar) && Rbar > thresh)) {
        return(Inf)
      }
      
      # Zar / CircStat intermediate 't' (see circ_confmean.m, eq. 26.24 & 26.25 logic)
      # For rbar < 0.9 branch:
      t_sq <- (2 * n_ori * (2 * R_ori^2 - n_ori * c2)) / (4 * n_ori - c2)
      t_sq <- max(t_sq, 0)            # guard small negatives from round-off
      t    <- sqrt(t_sq)
      
      #  transform to half-width (radians): delta = acos( t / R )
      # (CircStat applies this as the “final transform”)
      inside <- t / R_ori
      inside <- min(max(inside, -1), 1)  # numerical guard for acos
      hw <- acos(inside)
      
      if (axial) hw <- hw / 2
      hw
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

  #convert to degrees
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

  if (length(txt_files) == 0) stop("No .txt files found in path")
  
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
    nLines <- length(readLines(file))
    name   <- basename(file)

    max_vals <- phases * phaselength * resolution


    numphase <- floor(nLines / (phaselength * resolution))
    if (numphase == 0) next   # skip very short files
    if (numphase > phases) {
      numphase <- phases
    }

    # Skip the files that are not complete if 'all' = F
    if ((isFALSE(all) & numphase==phases) | isTRUE(all)){
      info <- file.info(file)
      Date <- format(info$mtime, "%d.%m.%Y")
      Time <- format(info$mtime, "%H:%M", tz="US/Central") ## tz="US/Central" ONLY FOR TEXAS DATA

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


# Test function with AI, to include different phase lengths.
get.mean.data2 <- function(path = NULL, channel = NULL, phases = NULL,
                          phaselength = NULL, blocks = NULL,
                          axial = NULL, samplingrate = 200,
                          all = FALSE, names = FALSE, verbose = FALSE) {
  
  if (is.null(path)) path <- choose.dir("Choose data directory")
  txt_files <- list.files(path = path, pattern = "\\.txt$", full.names = TRUE)
  if (length(txt_files) == 0) stop("No .txt files found in path")
  
  # resolution = samples per second (your original logic)
  resolution <- 1000 / samplingrate
  
  # Build blocks vector: if NULL, use phaselength for all phases
  if (is.null(blocks)) {
    if (is.null(phaselength)) stop("Either phaselength or blocks must be provided")
    blocks <- rep(phaselength, phases)
  } else {
    if (length(blocks) != phases) stop("Length of 'blocks' must equal 'phases'")
  }
  
  # samples per phase (absolute value), and sign indicates first/last
  samples_per_phase <- as.integer(round(abs(blocks) * resolution))
  if (any(samples_per_phase <= 0)) stop("blocks and samplingrate produce zero samples for a phase")
  
  cum_samples_needed <- cumsum(samples_per_phase)
  max_vals <- sum(samples_per_phase)
  
  results_list <- vector("list", length = 0)
  
  for (file in txt_files) {
    name <- basename(file)
    if (verbose) cat("Processing:", name, "\n")
    
    df <- tryCatch(read.delim(file, header = FALSE, sep = "", stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (is.null(df)) {
      if (verbose) cat("  read error, skipping\n")
      next
    }
    
    if (channel > ncol(df) || channel < 1) {
      if (verbose) cat("  channel", channel, "not found (ncol =", ncol(df), "), skipping\n")
      next
    }
    
    angles <- df[[channel]] * 3
    L <- length(angles)
    if (verbose) {
      cat("  samples in file:", L, "samples needed for all phases:", max_vals, "\n")
      cat("  samples per phase:", paste(samples_per_phase, collapse = ","), "\n")
      cat("  blocks (sign indicates first/last):", paste(blocks, collapse = ","), "\n")
    }
    
    # Determine how many full phases fit
    num_full_phases <- sum(cum_samples_needed <= L)
    
    # If no full phase fits but there are some samples, allow a partial first phase
    if (num_full_phases == 0 && L >= 1) {
      numphase_available <- 1
    } else {
      numphase_available <- min(phases, num_full_phases)
    }
    
    # If user requires all full phases, skip files that don't contain all full phases
    if (isFALSE(all) && num_full_phases < phases) {
      if (verbose) cat("  skipping (not all full phases present)\n")
      next
    }
    
    # Trim to maximum allowed samples (avoid accidental extra phases)
    if (L > max_vals) {
      angles <- angles[1:max_vals]
      L <- length(angles)
      if (verbose) cat("  trimmed to max_vals:", max_vals, "\n")
    }
    
    # Compute start/end indices for each phase (based on samples_per_phase)
    starts <- c(1, head(cum_samples_needed, -1) + 1)
    ends   <- cum_samples_needed
    
    circmeans    <- rep(NA_real_, phases)
    circr        <- rep(NA_real_, phases)
    circrayleigh <- rep(NA_real_, phases)
    
    for (i in seq_len(numphase_available)) {
      start_full <- starts[i]
      end_full   <- ends[i]
      desired    <- samples_per_phase[i]
      
      # If blocks[i] > 0 take first desired samples inside the phase window
      if (blocks[i] > 0) {
        start_i <- start_full
        end_i   <- min(start_full + desired - 1, end_full, L)
      } else {
        # blocks[i] < 0 take last desired samples inside the phase window
        start_i <- max(start_full, end_full - desired + 1)
        end_i   <- min(end_full, L)
      }
      
      if (start_i > end_i) {
        # no samples available for this phase -> leave NA
        if (verbose) cat("  phase", i, "has no samples (start_i >", "end_i)\n")
        next
      }
      
      phase_vec <- angles[start_i:end_i]
      
      type <- if (!is.null(axial) && i %in% axial) "axial" else "angles"
      
      stats <- circ.stats(phase_vec, type = type, unit = "deg")
      circmeans[i]    <- stats$mean_orientation
      circr[i]        <- stats$vector_length
      circrayleigh[i] <- if (type == "axial") {
        circ.rayleigh(phase_vec, axial = TRUE)$p.value
      } else {
        circ.rayleigh(phase_vec)$p.value
      }
    }
    
    # compute changes between consecutive available phases
    changes <- rep(NA_real_, phases - 1)
    if (phases > 1 && numphase_available > 1) {
      for (i in seq_len(numphase_available - 1)) {
        a1  <- circmeans[i]
        a2  <- circmeans[i + 1]
        ax1 <- !is.null(axial) && i %in% axial
        ax2 <- !is.null(axial) && (i + 1) %in% axial
        
        key <- paste(as.integer(ax1), as.integer(ax2), sep = "_")
        changes[i] <- switch(
          key,
          "0_0" = (a2 - a1) %% 360,
          "1_1" = abs(angleDiff(a1 * 2, a2 * 2, 360) / 2),
          "0_1" = {
            d1 <- angleDiff(a1, a2, 360)
            d2 <- angleDiff(a1, a2 + 180, 360)
            abs(if (abs(d1) <= abs(d2)) d1 else d2)
          },
          "1_0" = {
            d1 <- angleDiff(a1, a2, 360)
            d2 <- angleDiff(a1 + 180, a2, 360)
            abs(if (abs(d1) <= abs(d2)) d1 else d2)
          }
        )
      }
    }
    
    # Build row info
    info <- file.info(file)
    Date <- format(info$mtime, "%d.%m.%Y")
    Time <- format(info$mtime, "%H:%M")
    ID  <- tail(strsplit(name, "\\.")[[1]])[[1]]
    Exp <- tail(strsplit(name, "_|\\.")[[1]], 4)[1]
    Num <- tail(strsplit(name, "_|\\.|-")[[1]])[[3]]
    
    if (isTRUE(names)) {
      infovec <- c(ID, Num, Exp, Date, Time)
    } else {
      infovec <- c(ID, Date, Time)
    }
    
    row <- c(infovec,
             as.list(circmeans),
             as.list(circr),
             as.list(changes),
             as.list(circrayleigh))
    
    results_list[[length(results_list) + 1]] <- row
  }
  
  # If no files processed, return empty data.frame with appropriate column names
  if (length(results_list) == 0) {
    if (isTRUE(names)) {
      infocols <- c("ID", "Num", "Exp", "Date", "Time")
    } else {
      infocols <- c("ID", "Date", "Time")
    }
    meancols <- paste0("Phase", 1:phases, ".A")
    rcols    <- paste0("Phase", 1:phases, ".L")
    changecols <- if (phases > 1) paste0("Change", 1:(phases - 1)) else character(0)
    raycols  <- paste0("Rayleigh", 1:phases)
    allcols <- c(infocols, meancols, rcols, changecols, raycols)
    return(as.data.frame(matrix(ncol = length(allcols), nrow = 0, dimnames = list(NULL, allcols))))
  }
  
  # Convert list of rows to data.frame
  if (isTRUE(names)) {
    infocols <- c("ID", "Num", "Exp", "Date", "Time")
  } else {
    infocols <- c("ID", "Date", "Time")
  }
  meancols <- paste0("Phase", 1:phases, ".A")
  rcols    <- paste0("Phase", 1:phases, ".L")
  changecols <- if (phases > 1) paste0("Change", 1:(phases - 1)) else character(0)
  raycols  <- paste0("Rayleigh", 1:phases)
  colnames_all <- c(infocols, meancols, rcols, changecols, raycols)
  
  df_out <- as.data.frame(do.call(rbind, lapply(results_list, function(x) {
    xvec <- unlist(x, use.names = FALSE)
    len_needed <- length(colnames_all)
    if (length(xvec) < len_needed) xvec <- c(xvec, rep(NA, len_needed - length(xvec)))
    xvec
  })), stringsAsFactors = FALSE)
  
  colnames(df_out) <- colnames_all
  num_cols <- setdiff(colnames(df_out), infocols)
  df_out[num_cols] <- lapply(df_out[num_cols], function(x) as.numeric(as.character(x)))
  
  if (isFALSE(all)) {
    df_out <- df_out[complete.cases(df_out[meancols]), , drop = FALSE]
  }
  
  rownames(df_out) <- NULL
  df_out
}


