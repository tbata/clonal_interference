## =====================================================================
##  Companion code to "Some notes on modelling selective interference"
##  (SelectiveInterferenceModelSketch.md, section 3)
##
##  Illustrates:
##   Fig 1  Closed form for Pi_i (Step 6) vs Monte Carlo of the leader
##          chain (Step 5), compared with SSWM and classical GL weights.
##   Fig 2  log(Pi_i / mu_i) vs log(s_i) across loci for increasing N:
##          the slope (heuristic exponent Gamma) grows from 1 (SSWM).
##   Fig 3  Fitted Gamma vs N, split into displacement cost (-gamma)
##          and displacement gain (log B), and checked against the
##          weighted average of secant exponents (Steps 9 and 11).
##   Fig 4  Continuum limit, exponential spread of effects (Steps 9-10):
##          analytical Gamma(s) vs numerical log-log slope, its
##          components, and the ratio n/gamma.
##
##  Base R only. Output: selective_interference_figures.pdf
## =====================================================================

set.seed(2016)

## ---------------------------------------------------------------------
## 1. Model: closed form of Step 6, computed on the log scale
## ---------------------------------------------------------------------
## mu, s  : per-locus mutation rates and selection coefficients
## N      : population size;  kappa : establishment constant, p(s)=kappa*s
## withB  : TRUE = leader chain (Step 6); FALSE = classical GL (B = 1)
interference_model <- function(mu, s, N, kappa = 2, withB = TRUE) {
  o     <- order(s); mu <- mu[o]; s <- s[o]          # loci ordered by s
  L     <- length(s)
  lam   <- N * mu * kappa * s                         # establishment rates (Step 1)
  Lam   <- sum(lam)
  Fwin  <- log(N * kappa * s) / s                     # sweep window (Step 2)
  Lplus <- rev(cumsum(rev(lam))) - lam                # sum of lam over k > i
  gam   <- Fwin * Lplus                               # superior contenders (Step 3)
  r     <- ifelse(Lplus > 0, lam * (1 - exp(-gam)) / pmax(Lplus, 1e-300), 0)
  logB  <- if (withB) c(0, cumsum(log1p(r))[-L]) else rep(0, L)  # log prod_{j<i}(1+r_j)
  logPi <- log(lam) - log(Lam) - gam + logB           # Step 6, log scale
  logPi <- logPi - log(sum(exp(logPi - max(logPi)))) - max(logPi)  # renormalise
                                                      # (exact for withB=TRUE; needed for GL)
  data.frame(locus = seq_len(L), mu = mu, s = s, lam = lam, F = Fwin,
             gamma = gam, r = r, logB = logB, logPi = logPi, Pi = exp(logPi))
}

## SSWM weights (Section 1.2)
sswm <- function(mu, s) { w <- mu * s; w / sum(w) }

## ---------------------------------------------------------------------
## 2. Monte Carlo of the leader Markov chain (Step 5)
## ---------------------------------------------------------------------
simulate_chain <- function(m, nrep = 1e5) {
  L <- nrow(m); lam <- m$lam; gam <- m$gamma
  out <- integer(nrep)
  for (k in seq_len(nrep)) {
    j <- sample.int(L, 1, prob = lam)                 # first established mutation
    while (j < L && runif(1) > exp(-gam[j])) {        # displaced?
      sup <- (j + 1):L
      j <- if (length(sup) == 1) sup else sample(sup, 1, prob = lam[sup])
    }
    out[k] <- j
  }
  tabulate(out, L) / nrep
}

## ---------------------------------------------------------------------
## 3. Heuristic exponent: OLS slope and its secant decomposition
## ---------------------------------------------------------------------
fit_Gamma <- function(m) {
  x  <- log(m$s)
  y  <- m$logPi - log(m$mu)                           # log(Pi_i / mu_i)
  sl <- function(v) unname(coef(lm(v ~ x))[2])
  ## Linear decomposition of y (Step 9): y = log s - gamma + log B + const
  c(Gamma      = sl(y),
    cost       = sl(-m$gamma),                        # from (gamma_i - gamma_k)
    gain       = sl(m$logB))                          # from sum log(1 + r_j)
}

## Discrete Yitzhaki identity (Step 11): OLS slope = weighted mean of the
## secant slopes between adjacent loci, weights w_k = dx_k * sum_{i>k}(x_i - xbar)
secant_average <- function(m) {
  x  <- log(m$s); y <- m$logPi - log(m$mu); L <- length(x)
  dx <- diff(x); sec <- diff(y) / dx
  w  <- dx * sapply(1:(L - 1), function(k) sum(x[(k + 1):L] - mean(x)))
  sum(w * sec) / sum(w)
}

## ---------------------------------------------------------------------
## 4. Continuum limit, exponential spread of effects (Steps 9-10)
## ---------------------------------------------------------------------
## u(s) = (U/sigma) exp(-s/sigma): mutation-rate density over s
continuum <- function(N, U, sigma, kappa = 2) {
  ell   <- function(s) N * kappa * s * U / sigma * exp(-s / sigma)
  Lplus <- function(s) N * kappa * U * (s + sigma) * exp(-s / sigma)
  Fw    <- function(s) log(N * kappa * s) / s
  gam   <- function(s) Fw(s) * Lplus(s)
  nfun  <- function(s) s * ell(s) * Fw(s)
  lnB   <- function(s) integrate(function(x) ell(x) * (1 - exp(-gam(x))) / Lplus(x),
                                 lower = 1 / (N * kappa) * 3, upper = s,
                                 rel.tol = 1e-10)$value
  lnW   <- function(s) log(s) - gam(s) + lnB(s)      # weight per unit mutation rate
  Gamma_formula <- function(s) {                      # boxed result, Step 9
    g <- gam(s); n <- nfun(s)
    1 + g * (1 - 1 / log(N * kappa * s)) + n * (1 + (1 - exp(-g)) / g)
  }
  Gamma_numeric <- function(s, h = 1e-4)              # d lnW / d ln s by finite difference
    (lnW(s * exp(h)) - lnW(s * exp(-h))) / (2 * h)
  list(gam = gam, nfun = nfun, Gamma_formula = Gamma_formula,
       Gamma_numeric = Gamma_numeric, N = N)
}

## =====================================================================
##  Figures
## =====================================================================
pdf("selective_interference_figures.pdf", width = 8, height = 6)
cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")

## A shared set of loci: effects spread over a decade, heterogeneous mu
L  <- 10
s  <- sort(exp(runif(L, log(0.005), log(0.05))))
mu <- exp(runif(L, log(1e-9), log(1e-8)))

## ---- Fig 1: closed form vs Monte Carlo, SSWM and GL -------------------
N1 <- 1e6
m1 <- interference_model(mu, s, N1)
gl <- interference_model(mu, s, N1, withB = FALSE)
mc <- simulate_chain(m1, nrep = 5e4)
mat <- rbind(SSWM = sswm(m1$mu, m1$s), GL = gl$Pi, Chain = m1$Pi)
bp <- barplot(mat, beside = TRUE, col = cols[1:3], border = NA,
              names.arg = sprintf("%.3f", m1$s), las = 2, ylim = c(0, max(mat, mc) * 1.15),
              xlab = "locus selection coefficient s_i", ylab = "P(locus supplies next fixation)",
              main = sprintf("Fig 1. Pi_i: closed form (Step 6) vs Monte Carlo, N = %.0e", N1))
points(bp[3, ], mc, pch = 19, cex = 0.8)
legend("topleft", bty = "n",
       legend = c("SSWM  mu*s", "classical GL (B = 1)", "leader chain, closed form", "Monte Carlo of chain"),
       fill = c(cols[1:3], NA), border = NA, pch = c(NA, NA, NA, 19))

## ---- Fig 2: log(Pi/mu) vs log(s) for increasing N ------------------------
Ns <- 10^c(4, 5, 6, 6.5, 7)
ms <- lapply(Ns, function(N) interference_model(mu, s, N))
ys <- lapply(ms, function(m) { y <- m$logPi - log(m$mu); y - y[L] })
fits <- list()
plot(NA, xlim = range(log(s)), ylim = range(unlist(ys)),
     xlab = "log s_i", ylab = "log(Pi_i / mu_i)  (centred at the top locus)",
     main = "Fig 2. Cross-locus slope = heuristic exponent Gamma")
for (k in seq_along(Ns)) {
  m <- ms[[k]]; y <- ys[[k]]
  points(log(m$s), y, pch = 19, col = cols[k])
  lines(log(m$s), y, col = cols[k])
  fits[[k]] <- fit_Gamma(m)
}
abline(a = -log(s[L]), b = 1, lty = 3)               # slope 1 = SSWM, through top locus
legend("bottomright", bty = "n", col = c(cols[seq_along(Ns)], "black"),
       pch = c(rep(19, length(Ns)), NA), lty = c(rep(1, length(Ns)), 3),
       legend = c(sprintf("N = 1e%.1f, fitted Gamma = %.2f", log10(Ns),
                          sapply(fits, `[`, "Gamma")), "slope 1 (SSWM)"))

## ---- Fig 3: fitted Gamma vs N, decomposition, Yitzhaki check -------------
Ngrid <- 10^seq(3, 7.3, length.out = 40)
res <- t(sapply(Ngrid, function(N) {
  m <- interference_model(mu, s, N)
  c(fit_Gamma(m), secant = secant_average(m))
}))
plot(Ngrid, res[, "Gamma"], log = "x", type = "l", lwd = 2,
     ylim = c(0, max(res[, "Gamma"]) * 1.05),
     xlab = "population size N", ylab = "exponent",
     main = "Fig 3. Fitted Gamma = 1 + displacement cost + displacement gain")
lines(Ngrid, 1 + res[, "cost"], col = cols[2], lwd = 2, lty = 2)
lines(Ngrid, res[, "gain"], col = cols[3], lwd = 2, lty = 2)
points(Ngrid[seq(1, 40, 3)], res[seq(1, 40, 3), "secant"], pch = 1)
abline(h = 1, lty = 3)
legend("topleft", bty = "n", lwd = c(2, 2, 2, NA, 1), lty = c(1, 2, 2, NA, 3),
       pch = c(NA, NA, NA, 1, NA), col = c("black", cols[2], cols[3], "black", "black"),
       legend = c("fitted Gamma (OLS slope)", "1 + slope of -gamma (cost)",
                  "slope of log B (gain)", "weighted mean of secant slopes (Step 11)",
                  "SSWM"))

## ---- Fig 4: continuum, exponential spread of effects ---------------------
op <- par(mfrow = c(1, 2))
sig <- 0.01; U <- 1e-7
sgrid <- seq(0.005, 0.08, length.out = 60)
for (N in c(1e7, 1e8)) {
  cm <- continuum(N, U, sig)
  Gf <- sapply(sgrid, cm$Gamma_formula)
  Gn <- sapply(sgrid[seq(1, 60, 4)], cm$Gamma_numeric)
  g  <- sapply(sgrid, cm$gam); n <- sapply(sgrid, cm$nfun)
  plot(sgrid / sig, Gf, type = "l", lwd = 2, log = "y",
       ylim = range(c(1, Gf, 1 + g)),
       xlab = "s / sigma", ylab = "local exponent Gamma(s)",
       main = sprintf("Fig 4. Continuum, N = %.0e", N))
  points(sgrid[seq(1, 60, 4)] / sig, Gn, pch = 1)
  lines(sgrid / sig, 1 + g, col = cols[2], lty = 2, lwd = 2)          # superior contenders only
  lines(sgrid / sig, 1 + g + n, col = cols[3], lty = 2, lwd = 2)      # classical GL (B = 1)
  legend("topright", bty = "n", cex = 0.75, lwd = c(2, NA, 2, 2), lty = c(1, NA, 2, 2),
         pch = c(NA, 1, NA, NA), col = c("black", "black", cols[2], cols[3]),
         legend = c("Gamma(s), boxed formula", "numerical d lnW / d ln s",
                    "1 + gamma (earlier approximation)", "1 + gamma + n (classical GL)"))
}
par(op)

## ratio n/gamma = s^2 / (sigma (s + sigma)), independent of N and U
curve(x^2 / (1 + x), from = 0.5, to = 8, lwd = 2,
      xlab = "s / sigma", ylab = "n(s) / gamma(s)",
      main = "Fig 4b. Similar-effect vs superior contenders (exponential spread)")
abline(h = 1, lty = 3)

dev.off()

## ---------------------------------------------------------------------
## Console summary
## ---------------------------------------------------------------------
cat("\nFig 1 check (N = 1e6): max |closed form - Monte Carlo| =",
    signif(max(abs(m1$Pi - mc)), 3), "\n")
cat("Fig 3 check: max |OLS slope - weighted secant mean| =",
    signif(max(abs(res[, "Gamma"] - res[, "secant"])), 3), "\n")
cm <- continuum(1e7, U, sig)
cat("Fig 4 check (s = 3 sigma, N = 1e7): formula =", signif(cm$Gamma_formula(3 * sig), 6),
    " numerical =", signif(cm$Gamma_numeric(3 * sig), 6), "\n")
cat("Figures written to selective_interference_figures.pdf\n")
