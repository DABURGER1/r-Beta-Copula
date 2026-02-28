library(mvtnorm)
library(runjags)
library(coda)

qRectBeta <- function(p, mu, phi, rho, eps = 1e-12) {
  s1 <- 1 - abs(2*mu - 1)
  s2 <- (mu - 0.5*phi*s1)/(1 - phi*s1)
  alpha <- rho*s2
  beta <- rho*(1 - s2)
  w <- phi*s1
  
  invF <- function(u) {
    uniroot(function(y)
      w*y + (1 - w)*pbeta(y, alpha, beta) - u,
      c(eps, 1 - eps))$root
  }
  
  vapply(p, invF, numeric(1))
}

qRectBeta_vec <- Vectorize(qRectBeta, vectorize.args = c("p", "mu"))

rRecBetaGauss <- function(n,
                          X1, beta1,
                          X2, beta2,
                          phi1, rho1, phi2, rho2,
                          tau) {
  mu1 <- plogis(drop(X1%*%beta1))
  mu2 <- plogis(drop(X2%*%beta2))
  
  theta <- sin(0.5*pi*tau)
  R <- matrix(c(1, theta, theta, 1), 2)
  U <- pnorm(mvtnorm::rmvnorm(n, sigma = R))
  
  data.frame(
    Y1 = qRectBeta_vec(U[, 1], mu1, phi1, rho1),
    Y2 = qRectBeta_vec(U[, 2], mu2, phi2, rho2)
  )
}

rectbeta_gauss_model <- "
  model {
    for (k in 1:K1) { Beta1[k] ~ dnorm(0, 0.0001) }
    for (k in 1:K2) { Beta2[k] ~ dnorm(0, 0.0001) }
  
    phi1 ~ dunif(0, 1)
    rho1 ~ dgamma(0.0001, 0.0001)
    phi2 ~ dunif(0, 1)
    rho2 ~ dgamma(0.0001, 0.0001)
  
    tau ~ dunif(-1, 1)
    theta <- sin(0.5*pi*tau)
  
    for (i in 1:N) {
      mu1[i] <- ilogit(inprod(X1[i, ], Beta1[]))
      mu2[i] <- ilogit(inprod(X2[i, ], Beta2[]))
  
      S1[i] <- 1 - abs(2*mu1[i] - 1)
      S2[i] <- (mu1[i] - 0.5*phi1*S1[i])/(1 - phi1*S1[i])
      a1[i] <- rho1*S2[i]
      b1[i] <- rho1*(1 - S2[i])
      w1[i] <- phi1*S1[i]
  
      T1[i] <- 1 - abs(2*mu2[i] - 1)
      T2[i] <- (mu2[i] - 0.5*phi2*T1[i])/(1 - phi2*T1[i])
      a2[i] <- rho2*T2[i]
      b2[i] <- rho2*(1 - T2[i])
      w2[i] <- phi2*T1[i]
  
      pdf1[i] <- w1[i] + (1 - w1[i])*dbeta(Y1[i], a1[i], b1[i])
      pdf2[i] <- w2[i] + (1 - w2[i])*dbeta(Y2[i], a2[i], b2[i])
      cdf1[i] <- w1[i]*Y1[i] + (1 - w1[i])*pbeta(Y1[i], a1[i], b1[i])
      cdf2[i] <- w2[i]*Y2[i] + (1 - w2[i])*pbeta(Y2[i], a2[i], b2[i])
  
      z1[i] <- qnorm(cdf1[i], 0, 1)
      z2[i] <- qnorm(cdf2[i], 0, 1)
  
      logC[i] <- -0.5*log(1 - theta^2) -
                 (z1[i]^2 - 2*theta*z1[i]*z2[i] + z2[i]^2)/
                 (2*(1 - theta^2)) +
                 0.5*(z1[i]^2 + z2[i]^2)
  
      LL[i] <- log(pdf1[i]) + log(pdf2[i]) + logC[i]
      zeros[i] ~ dpois(BigC - LL[i])
    }
  }
"

run_rectbeta_sim <- function(n,
                             beta1,
                             beta2,
                             phi1,
                             rho1,
                             phi2,
                             rho2,
                             tau,
                             quick,
                             plot_data = FALSE,
                             nchains = 3) {
  x <- runif(n, -1, 1)
  X1 <- cbind(1, x)
  X2 <- X1
  dat <- rRecBetaGauss(n, X1, beta1, X2, beta2, phi1, rho1, phi2, rho2, tau)
  
  if (plot_data) {
    library(ggplot2)
    
    df_long <- data.frame(
      x = rep(x, 2),
      y = c(dat$Y1, dat$Y2),
      margin = factor(rep(c("Y1", "Y2"), each = n),
                      levels = c("Y1", "Y2"))
    )
    
    p <- ggplot(df_long, aes(x = x, y = y, colour = margin)) +
      geom_point(alpha = 0.6, size = 1.2) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4,
                 colour = "black", show.legend = FALSE) +
      labs(x = "Covariate  x",
           y = "Simulated response",
           colour = "Margin") +
      theme_classic()
    
    print(p)
  }
  
  jags_data <- list(
    N = n,
    K1 = ncol(X1), K2 = ncol(X2),
    X1 = X1, X2 = X2,
    Y1 = dat$Y1, Y2 = dat$Y2,
    zeros = rep(0, n),
    BigC = 1e5,
    pi = pi
  )
  
  inits <- replicate(
    nchains,
    list(
      Beta1 = rnorm(ncol(X1)),
      Beta2 = rnorm(ncol(X2)),
      phi1 = runif(1, 0.1, 0.5),
      rho1 = rgamma(1, 2, 2),
      phi2 = runif(1, 0.1, 0.5),
      rho2 = rgamma(1, 2, 2),
      tau = runif(1, -0.3, 0.3),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = sample.int(.Machine$integer.max, 1)
    ),
    simplify = FALSE
  )
  
  pars <- c("Beta1", "Beta2", "phi1", "rho1", "phi2", "rho2", "tau")
  
  fit <- if (quick) {
    run.jags(
      rectbeta_gauss_model, jags_data, inits,
      monitor = pars, n.chains = nchains,
      adapt = 500, burnin = 500, sample = 500, thin = 1,
      method = "parallel",
      modules = "glm",
      factories = "bugs::MNormal sampler off",
      silent.jags = FALSE
    )
  } else {
    autorun.jags(
      rectbeta_gauss_model, jags_data, inits,
      monitor = pars, n.chains = nchains,
      method = "parallel",
      adapt = 1000,
      startburnin = 4000,
      startsample = 4000,
      thin = 4,
      psrf.target = 1.1,
      max.time = "60m",
      modules = "glm",
      factories = "bugs::MNormal sampler off",
      silent.jags = FALSE
    )
  }
  
  mcmc_list <- as.mcmc.list(fit)
  all_draws <- as.matrix(mcmc_list)
  pooled <- as.mcmc(all_draws)
  hpd <- HPDinterval(pooled, prob = 0.95)
  
  summ <- data.frame(
    param = rownames(hpd),
    post_mean = colMeans(all_draws),
    post_median = apply(all_draws, 2, median),
    post_sd = apply(all_draws, 2, sd),
    hpd_lower = hpd[, "lower"],
    hpd_upper = hpd[, "upper"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  truth_vec <- c(
    setNames(beta1, c("Beta1[1]", "Beta1[2]")),
    setNames(beta2, c("Beta2[1]", "Beta2[2]")),
    phi1 = phi1, rho1 = rho1, phi2 = phi2, rho2 = rho2, tau = tau
  )
  summ$truth <- truth_vec[match(summ$param, names(truth_vec))]
  
  list(summary = summ, fit = fit)
}

run_many <- function(scenario_id,
                     iterations = 1,
                     out_dir = ".",
                     quick = TRUE,
                     ...) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  dots <- list(...)
  all_rows <- vector("list", iterations)
  fits <- vector("list", iterations)
  
  for (iter in seq_len(iterations)) {
    cat(sprintf("\033[1;95m[Scenario %d] Rep %d/%d\033[0m\n", scenario_id, iter, iterations))
    
    sim <- run_rectbeta_sim(quick = quick, ...)
    
    input_vals <- data.frame(
      scenario = scenario_id,
      replicate = iter,
      n = dots$n,
      phi1 = dots$phi1, phi2 = dots$phi2,
      rho1 = dots$rho1, rho2 = dots$rho2,
      tau = dots$tau,
      beta1_1 = dots$beta1[1], beta1_2 = dots$beta1[2],
      beta2_1 = dots$beta2[1], beta2_2 = dots$beta2[2]
    )
    
    all_rows[[iter]] <- cbind(
      input_vals[rep(1, nrow(sim$summary)), ],
      sim$summary[, c("param", "truth", "post_mean", "post_median", "post_sd", "hpd_lower", "hpd_upper")]
    )
    
    fits[[iter]] <- sim$fit
  }
  
  res_combined <- do.call(rbind, all_rows)
  ts <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  fname <- sprintf("Results_RecBetaGauss_S%02d_%s.csv", scenario_id, ts)
  write.csv(res_combined, file.path(out_dir, fname), row.names = FALSE)
  message("Saved combined results: ", fname)
  
  invisible(fits)
}