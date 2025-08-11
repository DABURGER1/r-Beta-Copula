library(runjags)
library(coda)
library(bridgesampling)
library(dplyr)
library(ggplot2)
library(grid)
library(stringr)
library(patchwork)
library(mvtnorm)
library(DHARMa)
library(copula)
library(ks)
library(reshape2)
library(gridExtra)
library(tibble)
library(tidyr)
library(xtable)

set.seed(2025)

NumChains <- 4

setwd(
  "C:/Users/a239866/OneDrive - Syneos Health/Divan @ Syneos Health - Linked Files/Research/Multivariate Outlier Proportions"
)
dat <-
  read.csv("Application/Datasets/Azorella-Agrostis.csv",
           stringsAsFactors = FALSE
  )

dat$V16[is.na(dat$V16)] <- 0
dat$NV16[is.na(dat$NV16)] <- 0

dat$Y1 <- (100*dat$DS03/dat$AzoArea03 + 0.1)/100
dat$Y2 <- (100*dat$DS16/dat$AzoArea16 + 0.1)/100

dat$BasePerAgr <- 100*dat$AgrArea03/dat$AzoArea03
dat$PostPerAgr <- 100*dat$AgrArea16/dat$AzoArea16
dat$LBaseAzoArea <- log(dat$AzoArea03)
dat$LPostAzoArea <- log(dat$AzoArea16)
dat$BasePerCO <- 100*(dat$V03 + dat$NV03)/dat$AzoArea03
dat$PostPerCO <- 100*(dat$V16 + dat$NV16)/dat$AzoArea16
dat$Altitude <- factor(trimws(dat$Altitude))
dat$Aspect <- factor(trimws(dat$Aspect))

vars_needed <- c(
  "Y1",
  "Y2",
  "Altitude",
  "Aspect",
  "BasePerAgr",
  "LBaseAzoArea",
  "BasePerCO",
  "PostPerAgr",
  "LPostAzoArea",
  "PostPerCO"
)
dat <- dat[complete.cases(dat[, vars_needed]), ]
cat("Rows:", nrow(dat), "\n")

X1 <- model.matrix(
  ~ Altitude + Aspect +
    BasePerAgr + LBaseAzoArea + BasePerCO,
  data = dat
)
X2 <- model.matrix(
  ~ Altitude + Aspect +
    PostPerAgr + LPostAzoArea + PostPerCO,
  data = dat
)
stopifnot(
  nrow(X1) == nrow(X2),
  all(rownames(X1) == rownames(X2))
)

Data <- data.frame(
  ID = as.numeric(as.factor(dat$Plot)),
  Y1 = dat$Y1,
  Y2 = dat$Y2
)

NObs <- nrow(Data)
NRand <- length(unique(Data$ID))
NFix1 <- ncol(X1)
NFix2 <- ncol(X2)

cat(
  sprintf(
    "Ready: %d obs, %d random levels, %d + %d fixed effects\n",
    NObs,
    NRand,
    NFix1,
    NFix2
  )
)

BigC <- 1e5

JAGSData <- list(
  NObs = NObs,
  NRand = NRand,
  NFix1 = NFix1,
  NFix2 = NFix2,
  ID = Data$ID,
  X1 = X1,
  X2 = X2,
  Y1 = Data$Y1,
  Y2 = Data$Y2,
  zeros = rep(0, NObs),
  BigC = BigC,
  pi = pi
)

InitsList <- replicate(
  NumChains,
  list(
    sigma1 = 0.2,
    sigma2 = 0.2,
    U1 = rep(0, NRand),
    U2 = rep(0, NRand),
    Beta1 = rep(0, NFix1),
    Beta2 = rep(0, NFix2),
    phi1 = 0.3,
    rho1 = 5,
    phi2 = 0.3,
    rho2 = 5,
    tau = 0.3,
    .RNG.name = "base::Wichmann-Hill",
    .RNG.seed = sample.int(1e6, 1)
  ),
  simplify = FALSE
)

############################
### Descriptive analyses ###
############################

PlotLookup <- dat %>%
  mutate(AltAsp = paste(Altitude, Aspect, sep = "/")) %>%
  distinct(Plot, AltAsp) %>%
  arrange(factor(Plot, levels = sort(unique(Plot)))) %>%
  mutate(
    PlotNum = row_number(),
    PlotLab = sprintf("Plot %d: %s", PlotNum, AltAsp)
  )

ScatterDF <- dat %>%
  mutate(AltAsp = paste(Altitude, Aspect, sep = "/")) %>%
  left_join(PlotLookup, by = c("Plot", "AltAsp")) %>%
  mutate(PlotLab = factor(PlotLab, levels = PlotLookup$PlotLab))

OverallScatter <- ggplot(ScatterDF, aes(Y1, Y2)) +
  geom_point(
    colour = "#0072B2",
    alpha = 0.75,
    size = 1.6
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  labs(x = "Dead Stem Cover (2003)", y = "Dead Stem Cover (2016)") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), legend.position = "none")

OverallTop <- ggplot(ScatterDF, aes(Y1)) +
  geom_histogram(
    bins = 40,
    fill = "steelblue1",
    colour = "grey25",
    linewidth = 0.3
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_void()

OverallRight <- ggplot(ScatterDF, aes(Y2)) +
  geom_histogram(
    bins = 40,
    fill = "#F0E442",
    colour = "grey25",
    linewidth = 0.3
  ) +
  coord_flip() +
  scale_x_continuous(limits = c(0, 1)) +
  theme_void()

LayoutOverall <- "
 A##
 BC#
 ###
"

OverallPlot <- OverallTop + OverallScatter + OverallRight +
  plot_layout(
    design = LayoutOverall,
    heights = c(1, 4, 0.2),
    widths = c(4, 1, 0.2)
  )

ggsave(
  "Manuscript/Output/DS_Scatter_Overall.pdf",
  OverallPlot,
  width = 5.12,
  height = 4.48,
  device = cairo_pdf
)

FacetLabels <- PlotLookup$PlotLab

MakePanel <- function(df, lab) {
  scatter <- ggplot(df, aes(Y1, Y2)) +
    geom_point(
      colour = "#0072B2",
      alpha = 0.8,
      size = 1.4
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    labs(x = "Dead Stem Cover (2003)", y = "Dead Stem Cover (2016)") +
    theme_bw(base_size = 8) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_text(size = 6),
      axis.title.y = element_text(size = 6)
    )
  
  top <- ggplot(df, aes(Y1)) +
    geom_histogram(
      bins = 25,
      fill = "steelblue1",
      colour = "grey25",
      alpha = 0.8,
      linewidth = 0.25
    ) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(title = lab) +
    theme_void() +
    theme(plot.title = element_text(
      hjust = 0.5,
      size = 8,
      margin = margin(b = 1)
    ))
  
  right <- ggplot(df, aes(Y2)) +
    geom_histogram(
      bins = 25,
      fill = "#F0E442",
      colour = "grey25",
      alpha = 0.8,
      linewidth = 0.25
    ) +
    coord_flip() +
    scale_x_continuous(limits = c(0, 1)) +
    theme_void()
  
  layout_local <- "
 A##
 BC#
 ###
 "
  
  top + scatter + right +
    plot_layout(
      design = layout_local,
      heights = c(1, 4, 0.2),
      widths = c(4, 1, 0.2)
    )
}

PanelList <- lapply(FacetLabels, function(lbl) {
  MakePanel(filter(ScatterDF, PlotLab == lbl), lbl)
})

FullFacets <- wrap_plots(PanelList, ncol = 4, nrow = 2) &
  theme(plot.margin = margin(2, 2, 2, 2))

ggsave(
  "Manuscript/Output/DS_Scatter_by_Cluster.pdf",
  FullFacets,
  width = 8,
  height = 4,
  device = cairo_pdf
)

#############################################################
############# Rectangular-beta + Gaussian copula ############
#############################################################

RecBetaGaussModel <- "
  model {

    for (k in 1:NFix1) {
      Beta1[k] ~ dnorm(0, 0.0001)
    }
    for (k in 1:NFix2) {
      Beta2[k] ~ dnorm(0, 0.0001)
    }

    sigma1 ~ dt(0, 0.25, 2)T(0, )
    sigma2 ~ dt(0, 0.25, 2)T(0, )
    tau1 <- pow(sigma1, -2)
    tau2 <- pow(sigma2, -2)

    for (j in 1:NRand) {
      U1[j] ~ dnorm(0, tau1)
      U2[j] ~ dnorm(0, tau2)
    }

    phi1 ~ dunif(0, 1)
    rho1 ~ dgamma(0.0001, 0.0001)
    phi2 ~ dunif(0, 1)
    rho2 ~ dgamma(0.0001, 0.0001)

    tau ~ dunif(-1, 1)
    rho_cop <- sin(0.5*pi*tau)

    for (i in 1:NObs) {

      LP1[i] <- inprod(X1[i, ], Beta1[]) + U1[ID[i]]
      LP2[i] <- inprod(X2[i, ], Beta2[]) + U2[ID[i]]
      p1[i] <- ilogit(LP1[i])
      p2[i] <- ilogit(LP2[i])

      S1_1[i] <- 1 - abs(2*p1[i] - 1)
      S2_1[i] <- (p1[i] - 0.5*phi1*S1_1[i])/(1 - phi1*S1_1[i])
      a1[i] <- rho1*S2_1[i]
      b1[i] <- rho1*(1 - S2_1[i])

      pdf1[i] <- phi1*S1_1[i] + (1 - phi1*S1_1[i])*dbeta(Y1[i], a1[i], b1[i])
      cdf1[i] <- phi1*S1_1[i]*Y1[i] + (1 - phi1*S1_1[i])*pbeta(Y1[i], a1[i], b1[i])

      S1_2[i] <- 1 - abs(2*p2[i] - 1)
      S2_2[i] <- (p2[i] - 0.5*phi2*S1_2[i])/(1 - phi2*S1_2[i])
      a2[i] <- rho2*S2_2[i]
      b2[i] <- rho2*(1 - S2_2[i])

      pdf2[i] <- phi2*S1_2[i] + (1 - phi2*S1_2[i])*dbeta(Y2[i], a2[i], b2[i])
      cdf2[i] <- phi2*S1_2[i]*Y2[i] + (1 - phi2*S1_2[i])*pbeta(Y2[i], a2[i], b2[i])

      z1[i] <- qnorm(cdf1[i], 0, 1)
      z2[i] <- qnorm(cdf2[i], 0, 1)

      logCop[i] <- -0.5*log(1 - pow(rho_cop, 2)) -
                   (pow(z1[i], 2) - 2*rho_cop*z1[i]*z2[i] + pow(z2[i], 2))/
                   (2*(1 - pow(rho_cop, 2))) +
                   0.5*(pow(z1[i], 2) + pow(z2[i], 2))

      LL[i] <- log(pdf1[i]) + log(pdf2[i]) + logCop[i]
      zeros[i] ~ dpois(BigC - LL[i])
    }
  }
"

RecBetaGaussFit <- function(Data, X1, X2) {
  params <- c(
    "Beta1",
    "Beta2",
    "sigma1",
    "sigma2",
    "U1",
    "U2",
    "phi1",
    "rho1",
    "phi2",
    "rho2",
    "tau"
  )
  
  run.jags(
    model = RecBetaGaussModel,
    data = JAGSData,
    inits = InitsList,
    monitor = params,
    n.chains = NumChains,
    adapt = 1000,
    burnin = 15000,
    sample = 1000,
    thin = 25,
    method = "parallel",
    modules = "glm",
    factories = "bugs::MNormal sampler off",
    silent.jags = FALSE
  )
}

RecBetaGaussRun <- FALSE
RecBetaGaussFitFile <- "Manuscript/Output/Rect_Beta_Gauss_Fit.rds"

if (!file.exists(RecBetaGaussFitFile)) {
  cat("No RDS found. Running Gaussian-copula model...\n")
  RecBetaGaussFitObj <- RecBetaGaussFit(Data, X1, X2)
  saveRDS(RecBetaGaussFitObj, RecBetaGaussFitFile)
} else {
  if (RecBetaGaussRun) {
    cat("RDS exists but RecBetaGaussRun = TRUE. Re-running...\n")
    RecBetaGaussFitObj <- RecBetaGaussFit(Data, X1, X2)
    saveRDS(RecBetaGaussFitObj, RecBetaGaussFitFile)
  } else {
    cat("Loading saved Gaussian-copula fit...\n")
    RecBetaGaussFitObj <- readRDS(RecBetaGaussFitFile)
  }
}

print(summary(RecBetaGaussFitObj))

RecBetaGaussMCMC <- as.mcmc.list(RecBetaGaussFitObj)
RecBetaGaussDraws <- do.call(rbind, RecBetaGaussMCMC)

cat(
  "Samples:",
  nrow(RecBetaGaussDraws),
  " Params:",
  ncol(RecBetaGaussDraws),
  "\n"
)

###########################
### Marginal likelihood ###
###########################

RecBetaGaussLogCop <- function(z1, z2, rho) {
  one_m_r2 <- 1 - rho^2
  -0.5*log(one_m_r2) -
    (z1^2 - 2*rho*z1*z2 + z2^2)/(2*one_m_r2) +
    0.5*(z1^2 + z2^2)
}

RecBetaGaussLogLik <- function(theta, dList) {
  names(theta) <- colnames(RecBetaGaussDraws)
  
  Beta1 <- theta[grep("^Beta1\\[", names(theta))]
  Beta2 <- theta[grep("^Beta2\\[", names(theta))]
  U1 <- theta[grep("^U1\\[", names(theta))]
  U2 <- theta[grep("^U2\\[", names(theta))]
  
  sigma1 <- theta["sigma1"]
  sigma2 <- theta["sigma2"]
  phi1 <- theta["phi1"]
  phi2 <- theta["phi2"]
  rho1 <- theta["rho1"]
  rho2 <- theta["rho2"]
  tau <- theta["tau"]
  rho_cop <- sin(0.5*base::pi*tau)
  
  ll <- 0
  for (i in seq_len(dList$NObs)) {
    lp1 <- sum(dList$X1[i, ]*Beta1) + U1[dList$ID[i]]
    lp2 <- sum(dList$X2[i, ]*Beta2) + U2[dList$ID[i]]
    p1 <- plogis(lp1)
    p2 <- plogis(lp2)
    
    s1 <- 1 - abs(2*p1 - 1)
    s2 <- (p1 - 0.5*phi1*s1)/(1 - phi1*s1)
    a1 <- rho1*s2
    b1 <- rho1*(1 - s2)
    pdf1 <- phi1*s1 + (1 - phi1*s1)*dbeta(dList$Y1[i], a1, b1)
    cdf1 <-
      phi1*s1*dList$Y1[i] + (1 - phi1*s1)*pbeta(dList$Y1[i], a1, b1)
    
    t1 <- 1 - abs(2*p2 - 1)
    t2 <- (p2 - 0.5*phi2*t1)/(1 - phi2*t1)
    a2 <- rho2*t2
    b2 <- rho2*(1 - t2)
    pdf2 <- phi2*t1 + (1 - phi2*t1)*dbeta(dList$Y2[i], a2, b2)
    cdf2 <-
      phi2*t1*dList$Y2[i] + (1 - phi2*t1)*pbeta(dList$Y2[i], a2, b2)
    
    z1 <- qnorm(cdf1)
    z2 <- qnorm(cdf2)
    log_cop <- RecBetaGaussLogCop(z1, z2, rho_cop)
    
    ll <- ll + log(pdf1) + log(pdf2) + log_cop
  }
  ll
}

RecBetaGaussLogPrior <- function(theta, dList) {
  names(theta) <- colnames(RecBetaGaussDraws)
  
  Beta1 <- theta[grep("^Beta1\\[", names(theta))]
  Beta2 <- theta[grep("^Beta2\\[", names(theta))]
  U1 <- theta[grep("^U1\\[", names(theta))]
  U2 <- theta[grep("^U2\\[", names(theta))]
  
  sigma1 <- theta["sigma1"]
  sigma2 <- theta["sigma2"]
  phi1 <- theta["phi1"]
  phi2 <- theta["phi2"]
  rho1 <- theta["rho1"]
  rho2 <- theta["rho2"]
  tau <- theta["tau"]
  rho_cop <- sin(0.5*base::pi*tau)
  
  lp <- sum(dnorm(Beta1, 0, 100, log = TRUE)) +
    sum(dnorm(Beta2, 0, 100, log = TRUE)) +
    sum(dnorm(U1, 0, sigma1, log = TRUE)) +
    sum(dnorm(U2, 0, sigma2, log = TRUE))
  
  half_t_log <- function(x, df = 2, scale = 2) {
    log(2) + dt(x/scale, df = df, log = TRUE) - log(scale)
  }
  
  lp <- lp + half_t_log(sigma1) + half_t_log(sigma2) +
    dunif(phi1, 0, 1, log = TRUE) + dunif(phi2, 0, 1, log = TRUE) +
    dgamma(rho1, 0.0001, 0.0001, log = TRUE) +
    dgamma(rho2, 0.0001, 0.0001, log = TRUE) +
    dunif(tau, -1, 1, log = TRUE)
  
  lp
}

RecBetaGaussLogPost <- function(theta, data) {
  lp <- RecBetaGaussLogPrior(theta, data)
  if (!is.finite(lp)) {
    return(lp)
  }
  ll <- RecBetaGaussLogLik(theta, data)
  if (!is.finite(ll)) {
    return(ll)
  }
  lp + ll
}

RecBetaGaussCN <- colnames(RecBetaGaussDraws)
p <- length(RecBetaGaussCN)
RecBetaGaussLB <- rep(-Inf, p)
RecBetaGaussUB <- rep(Inf, p)

Bound <- function(name, lo, hi) {
  i <- which(RecBetaGaussCN == name)
  if (length(i) == 1) {
    RecBetaGaussLB[i] <<- lo
    RecBetaGaussUB[i] <<- hi
  }
}

Bound("sigma1", 0, Inf)
Bound("sigma2", 0, Inf)
Bound("phi1", 0, 1)
Bound("phi2", 0, 1)
Bound("rho1", 0, Inf)
Bound("rho2", 0, Inf)
Bound("tau", -1, 1)

names(RecBetaGaussLB) <- names(RecBetaGaussUB) <- RecBetaGaussCN

RecBetaGaussList <- c(
  "RecBetaGaussLogCop",
  "RecBetaGaussLogLik",
  "RecBetaGaussLogPrior",
  "RecBetaGaussLogPost",
  "JAGSData",
  "RecBetaGaussCN",
  "RecBetaGaussDraws"
)

RecBetaGaussBridgeFile <-
  "Manuscript/Output/Rect_Beta_Gauss_Bridge.rds"
RecBetaGaussRunBridge <- TRUE

if (RecBetaGaussRunBridge || !file.exists(RecBetaGaussBridgeFile)) {
  RecBetaGaussBridge <- bridge_sampler(
    method = "warp3",
    samples = RecBetaGaussDraws,
    log_posterior = RecBetaGaussLogPost,
    data = JAGSData,
    lb = RecBetaGaussLB,
    ub = RecBetaGaussUB,
    cores = 7,
    packages = "bridgesampling",
    varlist = RecBetaGaussList,
    silent = FALSE,
    envir = .GlobalEnv
  )
  saveRDS(RecBetaGaussBridge, RecBetaGaussBridgeFile)
} else {
  RecBetaGaussBridge <- readRDS(RecBetaGaussBridgeFile)
}

cat(
  "Gaussian-copula log-marginal likelihood:",
  RecBetaGaussBridge$logml,
  "\n"
)
print(summary(RecBetaGaussBridge))

############################
### Residual diagnostics ###
############################

qRectBeta <- function(p, mu, phi, rho, eps = 1e-12) {
  s1 <- 1 - abs(2*mu - 1)
  s2 <- (mu - 0.5*phi*s1)/(1 - phi*s1)
  alpha <- rho*s2
  beta <- rho*(1 - s2)
  w <- phi*s1
  
  invF <- function(u) {
    uniroot(
      function(y) {
        w*y + (1 - w)*pbeta(y, alpha, beta) - u
      },
      c(eps, 1 - eps)
    )$root
  }
  
  vapply(p, invF, numeric(1))
}

qRectBetaVec <- Vectorize(qRectBeta, vectorize.args = c("p", "mu"))

RecBetaGaussPostPred <- function(fit, Data, X1, X2, ID, S = 1000) {
  mcmcMat <- do.call(rbind, as.mcmc.list(fit))
  draws <- mcmcMat[sample(seq_len(nrow(mcmcMat)), S, FALSE), ]
  
  n <- nrow(Data)
  sim1 <- matrix(NA_real_, n, S)
  sim2 <- matrix(NA_real_, n, S)
  fit1 <- numeric(n)
  fit2 <- numeric(n)
  
  for (s in seq_len(S)) {
    d <- draws[s, ]
    B1 <- d[grep("^Beta1\\[", names(d))]
    B2 <- d[grep("^Beta2\\[", names(d))]
    U1 <- d[grep("^U1\\[", names(d))]
    U2 <- d[grep("^U2\\[", names(d))]
    
    phi1 <- d["phi1"]
    rho1 <- d["rho1"]
    phi2 <- d["phi2"]
    rho2 <- d["rho2"]
    rho_cop <- sin(0.5*pi*d["tau"])
    
    p1 <- plogis(X1 %*% B1 + U1[ID])
    p2 <- plogis(X2 %*% B2 + U2[ID])
    
    R <- matrix(c(1, rho_cop, rho_cop, 1), 2)
    Z <- rmvt(n, sigma = R, df = Inf)
    U <- pt(Z, df = Inf)
    
    sim1[, s] <- qRectBetaVec(U[, 1], p1, phi1, rho1)
    sim2[, s] <- qRectBetaVec(U[, 2], p2, phi2, rho2)
    
    if (s == 1) {
      fit1 <- p1
      fit2 <- p2
    }
  }
  list(
    sim1 = sim1,
    sim2 = sim2,
    fit1 = fit1,
    fit2 = fit2
  )
}

RecBetaGaussDiagnostics <- function(predObj,
                                    Data,
                                    outDir = "Manuscript/Output",
                                    prefix = "DeadStemCover_Gauss") {
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE)
  }
  
  keep1 <- !is.na(Data$Y1)
  keep2 <- !is.na(Data$Y2)
  
  res1 <- createDHARMa(predObj$sim1[keep1, , drop = FALSE],
                       Data$Y1[keep1], predObj$fit1[keep1],
                       integerResponse = FALSE
  )
  res2 <- createDHARMa(predObj$sim2[keep2, , drop = FALSE],
                       Data$Y2[keep2], predObj$fit2[keep2],
                       integerResponse = FALSE
  )
  
  pvals <- tibble(
    Response = c("DeadStem2003", "DeadStem2016"),
    Uniformity = c(
      testUniformity(res1)$p.value,
      testUniformity(res2)$p.value
    ),
    Dispersion = c(
      testDispersion(res1)$p.value,
      testDispersion(res2)$p.value
    ),
    Outliers = c(
      testOutliers(res1)$p.value,
      testOutliers(res2)$p.value
    )
  )
  
  write.csv(pvals,
            file.path(outDir, paste0(prefix, "_DHARMa_pvalues.csv")),
            row.names = FALSE
  )
  
  makeQQ <- function(r, ttl) {
    ggplot(
      data.frame(
        Expected = qunif(ppoints(r)),
        Observed = sort(r)
      ),
      aes(Expected, Observed)
    ) +
      geom_point(size = 2, colour = "blue3") +
      geom_abline(
        slope = 1,
        intercept = 0,
        colour = "darkred",
        linewidth = 0.9
      ) +
      labs(title = ttl, x = "Expected Residual", y = "Observed Residual") +
      theme_classic(base_size = 18) +
      theme(
        axis.line = element_line(colour = "black", linewidth = 0.9),
        plot.title = element_text(hjust = 0.5)
      )
  }
  
  qq2003 <- makeQQ(residuals(res1), "Dead Stem Cover (2003)")
  qq2016 <- makeQQ(residuals(res2), "Dead Stem Cover (2016)")
  
  ggsave(
    file.path(outDir, paste0(prefix, "_DHARMa_QQ.pdf")),
    qq2003 | qq2016,
    width = 14,
    height = 7,
    units = "in"
  )
  
  pdf(file.path(outDir, paste0(prefix, "_DHARMa_4panel.pdf")),
      width = 8,
      height = 10
  )
  oldpar <- par(mfrow = c(4, 2))
  plot(res1)
  plot(res2)
  par(oldpar)
  dev.off()
  
  invisible(list(
    res1 = res1,
    res2 = res2,
    pvals = pvals
  ))
}

RecBetaGaussPredObj <- RecBetaGaussPostPred(
  fit = RecBetaGaussFitObj,
  Data = Data,
  X1 = X1,
  X2 = X2,
  ID = Data$ID,
  S = nrow(do.call(rbind, as.mcmc.list(RecBetaGaussFitObj)))
)

RecBetaGaussDHARMa <- RecBetaGaussDiagnostics(RecBetaGaussPredObj, Data)

################################################
### Goodness of fit and dependence structure ###
################################################

pRectBeta <- function(y, p, phi, rho) {
  star1 <- 1 - abs(2*p - 1)
  w <- phi*star1
  star2 <- (p - 0.5*phi*star1)/(1 - phi*star1)
  alpha <- rho*star2
  beta <- rho*(1 - star2)
  w*y + (1 - w)*pbeta(y, alpha, beta)
}

RecBetaGaussGOF <- function(fit, Data, X1, X2, B = 500) {
  mcmcMat <- do.call(rbind, as.mcmc.list(fit))
  n_draws <- nrow(mcmcMat)
  n <- nrow(Data)
  
  U1_draws <- matrix(NA_real_, n, n_draws)
  U2_draws <- matrix(NA_real_, n, n_draws)
  
  for (s in seq_len(n_draws)) {
    d <- mcmcMat[s, ]
    Beta1 <- d[grep("^Beta1\\[", names(d))]
    Beta2 <- d[grep("^Beta2\\[", names(d))]
    U1re <- d[grep("^U1\\[", names(d))]
    U2re <- d[grep("^U2\\[", names(d))]
    
    p1 <- plogis(X1 %*% Beta1 + U1re[Data$ID])
    p2 <- plogis(X2 %*% Beta2 + U2re[Data$ID])
    
    U1_draws[, s] <- pRectBeta(Data$Y1, p1, d["phi1"], d["rho1"])
    U2_draws[, s] <- pRectBeta(Data$Y2, p2, d["phi2"], d["rho2"])
  }
  
  u_seq <- seq(0.05, 0.95, by = 0.05)
  
  chi_fn <- function(U) {
    num <- sapply(u_seq, function(u) mean((U[, 1] > 1 - u) & (U[, 2] > 1 - u)))
    den <- sapply(u_seq, function(u) mean(U[, 2] > 1 - u))
    z <- num/den
    z[!is.finite(z)] <- NA_real_
    z
  }
  
  K_fn <- function(U) {
    sapply(u_seq, function(u) mean((U[, 1] <= u) & (U[, 2] <= u)))
  }
  
  chi_draws <- matrix(NA_real_, length(u_seq), n_draws)
  K_draws <- matrix(NA_real_, length(u_seq), n_draws)
  for (s in seq_len(n_draws)) {
    Umat <- cbind(U1_draws[, s], U2_draws[, s])
    chi_draws[, s] <- chi_fn(Umat)
    K_draws[, s] <- K_fn(Umat)
  }
  chi_emp <- rowMeans(chi_draws, na.rm = TRUE)
  K_emp <- rowMeans(K_draws, na.rm = TRUE)
  
  boot_id <- sample(seq_len(n_draws), B, replace = TRUE)
  chi_sim <- matrix(NA_real_, length(u_seq), B)
  K_sim <- matrix(NA_real_, length(u_seq), B)
  
  for (b in seq_len(B)) {
    tau_b <- mcmcMat[boot_id[b], "tau"]
    rho_b <- sin(0.5*pi*tau_b)
    Usim <- rCopula(n, normalCopula(rho_b, dim = 2))
    chi_sim[, b] <- chi_fn(Usim)
    K_sim[, b] <- K_fn(Usim)
  }
  
  chi_lo <- apply(chi_sim, 1, quantile, 0.025, na.rm = TRUE)
  chi_hi <- apply(chi_sim, 1, quantile, 0.975, na.rm = TRUE)
  K_lo <- apply(K_sim, 1, quantile, 0.025, na.rm = TRUE)
  K_hi <- apply(K_sim, 1, quantile, 0.975, na.rm = TRUE)
  
  plot_env <- function(emp, lo, hi, ylab) {
    data.frame(u = u_seq, emp = emp, lo = lo, hi = hi) |>
      ggplot(aes(u, emp)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#6393ff", alpha = 0.25) +
      geom_line(colour = "red4", linewidth = 1) +
      scale_x_continuous(
        limits = c(0.05, 0.95),
        breaks = seq(0.05, 0.95, by = 0.1),
        labels = function(x) sprintf("%.2f", x),
        expand = c(0, 0.04)
      ) +
      labs(x = "u", y = ylab) +
      theme_classic(base_size = 15)
  }
  
  
  p_chi <- plot_env(chi_emp, chi_lo, chi_hi, expression(chi(u)))
  p_K <- plot_env(K_emp, K_lo, K_hi, expression(K(u)))
  
  if (!dir.exists("Manuscript/Output")) dir.create("Manuscript/Output", recursive = TRUE)
  
  ggsave("Manuscript/Output/DeadStemCover_Gauss_ChiPlot.pdf",
         p_chi,
         width = 5.5, height = 4.5, device = cairo_pdf
  )
  ggsave("Manuscript/Output/DeadStemCover_Gauss_KPlot.pdf",
         p_K,
         width = 5.5, height = 4.5, device = cairo_pdf
  )
  
  cat("Plots saved:\n - Gauss_ChiPlot.pdf\n - Gauss_KPlot.pdf\n")
}

RecBetaGaussGOF(
  fit  = RecBetaGaussFitObj,
  Data = Data,
  X1   = X1,
  X2   = X2
)

###########################################################
############# Rectangular-beta + Gumbel copula ############
###########################################################

RecBetaGumbelModel <- "
  model {

    for (k in 1:NFix1) {
      Beta1[k] ~ dnorm(0, 0.0001)
    }
    for (k in 1:NFix2) {
      Beta2[k] ~ dnorm(0, 0.0001)
    }

    sigma1 ~ dt(0, 0.25, 2)T(0, )
    sigma2 ~ dt(0, 0.25, 2)T(0, )
    tau1 <- pow(sigma1, -2)
    tau2 <- pow(sigma2, -2)

    for (j in 1:NRand) {
      U1[j] ~ dnorm(0, tau1)
      U2[j] ~ dnorm(0, tau2)
    }

    phi1 ~ dunif(0, 1)
    rho1 ~ dgamma(0.0001, 0.0001)
    phi2 ~ dunif(0, 1)
    rho2 ~ dgamma(0.0001, 0.0001)

    tau ~ dunif(0, 1)
    theta_cop <- 1/(1 - tau)
    inv_theta <- 1 - tau

    for (i in 1:NObs) {

      LP1[i] <- inprod(X1[i, ], Beta1[]) + U1[ID[i]]
      LP2[i] <- inprod(X2[i, ], Beta2[]) + U2[ID[i]]
      p1[i] <- ilogit(LP1[i])
      p2[i] <- ilogit(LP2[i])

      S1_1[i] <- 1 - abs(2*p1[i] - 1)
      S2_1[i] <- (p1[i] - 0.5*phi1*S1_1[i])/(1 - phi1*S1_1[i])
      a1[i] <- rho1*S2_1[i]
      b1[i] <- rho1*(1 - S2_1[i])

      pdf1[i] <- phi1*S1_1[i] + (1 - phi1*S1_1[i])*dbeta(Y1[i], a1[i], b1[i])
      cdf1[i] <- phi1*S1_1[i]*Y1[i] + (1 - phi1*S1_1[i])*pbeta(Y1[i], a1[i], b1[i])

      S1_2[i] <- 1 - abs(2*p2[i] - 1)
      S2_2[i] <- (p2[i] - 0.5*phi2*S1_2[i])/(1 - phi2*S1_2[i])
      a2[i] <- rho2*S2_2[i]
      b2[i] <- rho2*(1 - S2_2[i])

      pdf2[i] <- phi2*S1_2[i] + (1 - phi2*S1_2[i])*dbeta(Y2[i], a2[i], b2[i])
      cdf2[i] <- phi2*S1_2[i]*Y2[i] + (1 - phi2*S1_2[i])*pbeta(Y2[i], a2[i], b2[i])

      t1[i] <- pow(-log(cdf1[i]), theta_cop)
      t2[i] <- pow(-log(cdf2[i]), theta_cop)
      s[i] <- t1[i] + t2[i]

      log_phi2[i] <- log(inv_theta) -
                     pow(s[i], inv_theta) +
                     (inv_theta - 2)*log(s[i]) +
                     log(inv_theta*pow(s[i], inv_theta) - (inv_theta - 1))

      log_phi_inv_u[i] <- log(theta_cop) +
                          (theta_cop - 1)*log(-log(cdf1[i])) -
                          log(cdf1[i])

      log_phi_inv_v[i] <- log(theta_cop) +
                          (theta_cop - 1)*log(-log(cdf2[i])) -
                          log(cdf2[i])

      logCop[i] <- log_phi2[i] + log_phi_inv_u[i] + log_phi_inv_v[i]

      LL[i] <- log(pdf1[i]) + log(pdf2[i]) + logCop[i]
      zeros[i] ~ dpois(BigC - LL[i])
    }
  }
"

RecBetaGumbelFit <- function(Data, X1, X2) {
  params <- c(
    "Beta1",
    "Beta2",
    "sigma1",
    "sigma2",
    "U1",
    "U2",
    "phi1",
    "rho1",
    "phi2",
    "rho2",
    "tau",
    "theta_cop"
  )
  
  run.jags(
    model = RecBetaGumbelModel,
    data = JAGSData,
    inits = InitsList,
    monitor = params,
    n.chains = NumChains,
    adapt = 1000,
    burnin = 15000,
    sample = 1000,
    thin = 25,
    method = "parallel",
    modules = "glm",
    factories = "bugs::MNormal sampler off",
    silent.jags = FALSE
  )
}

RecBetaGumbelRun <- FALSE
RecBetaGumbelFitFile <- "Manuscript/Output/Rect_Beta_Gumbel_Fit.rds"

if (!file.exists(RecBetaGumbelFitFile)) {
  cat("No RDS found. Running model...\n")
  RecBetaGumbelFitObj <- RecBetaGumbelFit(Data, X1, X2)
  saveRDS(RecBetaGumbelFitObj, RecBetaGumbelFitFile)
} else {
  if (RecBetaGumbelRun) {
    cat("RDS exists but RecBetaGumbelRun = TRUE. Re-running...\n")
    RecBetaGumbelFitObj <- RecBetaGumbelFit(Data, X1, X2)
    saveRDS(RecBetaGumbelFitObj, RecBetaGumbelFitFile)
  } else {
    cat("Loading saved fit...\n")
    RecBetaGumbelFitObj <- readRDS(RecBetaGumbelFitFile)
  }
}

print(summary(RecBetaGumbelFitObj))

RecBetaGumbelMCMC <- as.mcmc.list(RecBetaGumbelFitObj)
RecBetaGumbelDraws <- do.call(rbind, RecBetaGumbelMCMC)

cat(
  "Samples:",
  nrow(RecBetaGumbelDraws),
  " Params:",
  ncol(RecBetaGumbelDraws),
  "\n"
)

#################################
### Posterior tail dependence ###
#################################

ggsave(
  filename = "Manuscript/Output/DeadStemCover_Gumbel_Lambda_U.pdf",
  plot = ggplot(data.frame(lambdaU = 2 - 2^(1/if ("theta_cop" %in% colnames(RecBetaGumbelDraws)) {
    RecBetaGumbelDraws[, "theta_cop"]
  } else {
    1/(1 - RecBetaGumbelDraws[, "tau"])
  })), aes(lambdaU)) +
    geom_histogram(
      bins = 35,
      colour = "grey30",
      fill = "#3182BD",
      linewidth = 0.4
    ) +
    labs(
      x = expression(lambda[U]),
      y = "Posterior Frequency"
    ) +
    theme_classic(base_size = 18) +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.8),
      axis.ticks = element_line(colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.margin = margin(6, 6, 6, 6)
    ),
  width = 5.5,
  height = 4.5,
  device = cairo_pdf
)

###########################
### Marginal likelihood ###
###########################

RecBetaGumLogCop <- function(u, v, theta_cop) {
  inv_theta <- 1/theta_cop
  
  t1 <- (-log(u))^theta_cop
  t2 <- (-log(v))^theta_cop
  s <- t1 + t2
  
  log_phi2 <- log(inv_theta) -
    s^inv_theta +
    (inv_theta - 2)*log(s) +
    log(inv_theta*s^inv_theta - (inv_theta - 1))
  
  log_phi_inv_u <- log(theta_cop) +
    (theta_cop - 1)*log(-log(u)) -
    log(u)
  
  log_phi_inv_v <- log(theta_cop) +
    (theta_cop - 1)*log(-log(v)) -
    log(v)
  
  log_cop <- log_phi2 + log_phi_inv_u + log_phi_inv_v
  log_cop
}

RecBetaGumLogLik <- function(theta, dList) {
  names(theta) <- colnames(RecBetaGumbelDraws)
  
  Beta1 <- theta[grep("^Beta1\\[", names(theta))]
  Beta2 <- theta[grep("^Beta2\\[", names(theta))]
  U1 <- theta[grep("^U1\\[", names(theta))]
  U2 <- theta[grep("^U2\\[", names(theta))]
  
  sigma1 <- theta["sigma1"]
  sigma2 <- theta["sigma2"]
  phi1 <- theta["phi1"]
  phi2 <- theta["phi2"]
  rho1 <- theta["rho1"]
  rho2 <- theta["rho2"]
  tau <- theta["tau"]
  theta_cop <- 1/(1 - tau)
  
  ll <- 0
  for (i in seq_len(dList$NObs)) {
    lp1 <- sum(dList$X1[i, ]*Beta1) + U1[dList$ID[i]]
    lp2 <- sum(dList$X2[i, ]*Beta2) + U2[dList$ID[i]]
    p1 <- plogis(lp1)
    p2 <- plogis(lp2)
    
    s1 <- 1 - abs(2*p1 - 1)
    s2 <- (p1 - 0.5*phi1*s1)/(1 - phi1*s1)
    a1 <- rho1*s2
    b1 <- rho1*(1 - s2)
    pdf1 <- phi1*s1 + (1 - phi1*s1)*dbeta(dList$Y1[i], a1, b1)
    cdf1 <-
      phi1*s1*dList$Y1[i] + (1 - phi1*s1)*pbeta(dList$Y1[i], a1, b1)
    
    t1 <- 1 - abs(2*p2 - 1)
    t2 <- (p2 - 0.5*phi2*t1)/(1 - phi2*t1)
    a2 <- rho2*t2
    b2 <- rho2*(1 - t2)
    pdf2 <- phi2*t1 + (1 - phi2*t1)*dbeta(dList$Y2[i], a2, b2)
    cdf2 <-
      phi2*t1*dList$Y2[i] + (1 - phi2*t1)*pbeta(dList$Y2[i], a2, b2)
    
    log_cop <- RecBetaGumLogCop(cdf1, cdf2, theta_cop)
    
    ll <- ll + log(pdf1) + log(pdf2) + log_cop
  }
  ll
}

RecBetaGumLogPrior <- function(theta, dList) {
  names(theta) <- colnames(RecBetaGumbelDraws)
  
  Beta1 <- theta[grep("^Beta1\\[", names(theta))]
  Beta2 <- theta[grep("^Beta2\\[", names(theta))]
  U1 <- theta[grep("^U1\\[", names(theta))]
  U2 <- theta[grep("^U2\\[", names(theta))]
  
  sigma1 <- theta["sigma1"]
  sigma2 <- theta["sigma2"]
  phi1 <- theta["phi1"]
  phi2 <- theta["phi2"]
  rho1 <- theta["rho1"]
  rho2 <- theta["rho2"]
  tau <- theta["tau"]
  
  lp <- 0
  lp <- lp + sum(dnorm(Beta1, 0, 100, log = TRUE))
  lp <- lp + sum(dnorm(Beta2, 0, 100, log = TRUE))
  lp <- lp + sum(dnorm(U1, 0, sigma1, log = TRUE))
  lp <- lp + sum(dnorm(U2, 0, sigma2, log = TRUE))
  
  half_t_log <- function(x, df = 2, scale = 2) {
    log(2) + dt(x/scale, df = df, log = TRUE) - log(scale)
  }
  
  lp <- lp + half_t_log(sigma1) + half_t_log(sigma2)
  lp <-
    lp + dunif(phi1, 0, 1, log = TRUE) + dunif(phi2, 0, 1, log = TRUE)
  lp <-
    lp + dgamma(rho1, 0.0001, 0.0001, log = TRUE) + dgamma(rho2, 0.0001, 0.0001, log = TRUE)
  lp <- lp + dunif(tau, 0, 1, log = TRUE)
  
  lp
}

RecBetaGumLogPost <- function(theta, data) {
  lp <- RecBetaGumLogPrior(theta, data)
  if (!is.finite(lp)) {
    return(lp)
  }
  ll <- RecBetaGumLogLik(theta, data)
  if (!is.finite(ll)) {
    return(ll)
  }
  lp + ll
}

RecBetaGumCN <- colnames(RecBetaGumbelDraws)
p <- length(RecBetaGumCN)
RecBetaGumLB <- rep(-Inf, p)
RecBetaGumUB <- rep(Inf, p)

Bound <- function(name, lo, hi) {
  i <- which(RecBetaGumCN == name)
  if (length(i) == 1) {
    RecBetaGumLB[i] <<- lo
    RecBetaGumUB[i] <<- hi
  }
}

Bound("sigma1", 0, Inf)
Bound("sigma2", 0, Inf)
Bound("phi1", 0, 1)
Bound("phi2", 0, 1)
Bound("rho1", 0, Inf)
Bound("rho2", 0, Inf)
Bound("tau", 0, 1)

names(RecBetaGumLB) <- names(RecBetaGumUB) <- RecBetaGumCN

RecBetaGumList <- c(
  "RecBetaGumLogCop",
  "RecBetaGumLogLik",
  "RecBetaGumLogPrior",
  "RecBetaGumLogPost",
  "JAGSData",
  "RecBetaGumCN",
  "RecBetaGumbelDraws"
)

RecBetaGumBridgeFile <-
  "Manuscript/Output/Rect_Beta_Gumbel_Bridge.rds"
RecBetaGumRunBridge <- TRUE

if (RecBetaGumRunBridge || !file.exists(RecBetaGumBridgeFile)) {
  RecBetaGumBridge <- bridge_sampler(
    method = "warp3",
    samples = RecBetaGumbelDraws,
    log_posterior = RecBetaGumLogPost,
    data = JAGSData,
    lb = RecBetaGumLB,
    ub = RecBetaGumUB,
    cores = 7,
    packages = "bridgesampling",
    varlist = RecBetaGumList,
    silent = FALSE,
    envir = .GlobalEnv
  )
  saveRDS(RecBetaGumBridge, RecBetaGumBridgeFile)
} else {
  RecBetaGumBridge <- readRDS(RecBetaGumBridgeFile)
}

cat(
  "Gumbel copula log-marginal likelihood:",
  RecBetaGumBridge$logml,
  "\n"
)
print(summary(RecBetaGumBridge))

############################
### Residual diagnostics ###
############################

RecBetaGumbelPostPred <- function(fit, Data, X1, X2, ID, S = 1000) {
  mcmcMat <- do.call(rbind, as.mcmc.list(fit))
  draws <- mcmcMat[sample(seq_len(nrow(mcmcMat)), S, FALSE), ]
  
  n <- nrow(Data)
  sim1 <- matrix(NA_real_, n, S)
  sim2 <- matrix(NA_real_, n, S)
  fit1 <- numeric(n)
  fit2 <- numeric(n)
  
  for (s in seq_len(S)) {
    d <- draws[s, ]
    B1 <- d[grep("^Beta1\\[", names(d))]
    B2 <- d[grep("^Beta2\\[", names(d))]
    U1 <- d[grep("^U1\\[", names(d))]
    U2 <- d[grep("^U2\\[", names(d))]
    
    phi1 <- d["phi1"]
    rho1 <- d["rho1"]
    phi2 <- d["phi2"]
    rho2 <- d["rho2"]
    
    theta <- 1/(1 - d["tau"])
    cop <- gumbelCopula(theta, dim = 2)
    
    p1 <- plogis(X1 %*% B1 + U1[ID])
    p2 <- plogis(X2 %*% B2 + U2[ID])
    
    U <- rCopula(n, cop)
    
    sim1[, s] <- qRectBetaVec(U[, 1], p1, phi1, rho1)
    sim2[, s] <- qRectBetaVec(U[, 2], p2, phi2, rho2)
    
    if (s == 1) {
      fit1 <- p1
      fit2 <- p2
    }
  }
  list(
    sim1 = sim1,
    sim2 = sim2,
    fit1 = fit1,
    fit2 = fit2
  )
}

GumbelDiagnostics <- function(predObj,
                              Data,
                              outDir = "Manuscript/Output",
                              prefix = "DeadStemCover_Gumbel") {
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE)
  }
  
  keep1 <- !is.na(Data$Y1)
  keep2 <- !is.na(Data$Y2)
  
  res1 <- createDHARMa(predObj$sim1[keep1, , drop = FALSE],
                       Data$Y1[keep1], predObj$fit1[keep1],
                       integerResponse = FALSE
  )
  res2 <- createDHARMa(predObj$sim2[keep2, , drop = FALSE],
                       Data$Y2[keep2], predObj$fit2[keep2],
                       integerResponse = FALSE
  )
  
  pvals <- tibble(
    Response = c("DeadStem2003", "DeadStem2016"),
    Uniformity = c(
      testUniformity(res1)$p.value,
      testUniformity(res2)$p.value
    ),
    Dispersion = c(
      testDispersion(res1)$p.value,
      testDispersion(res2)$p.value
    ),
    Outliers = c(
      testOutliers(res1)$p.value,
      testOutliers(res2)$p.value
    )
  )
  
  write.csv(pvals,
            file.path(outDir, paste0(prefix, "_DHARMa_pvalues.csv")),
            row.names = FALSE
  )
  
  makeQQ <- function(r, ttl) {
    ggplot(
      data.frame(Expected = qunif(ppoints(r)), Observed = sort(r)),
      aes(Expected, Observed)
    ) +
      geom_point(size = 2, colour = "blue3") +
      geom_abline(
        slope = 1,
        intercept = 0,
        colour = "darkred",
        linewidth = 0.9
      ) +
      labs(title = ttl, x = "Expected Residual", y = "Observed Residual") +
      theme_classic(base_size = 18) +
      theme(
        axis.line = element_line(colour = "black", linewidth = 0.9),
        plot.title = element_text(hjust = 0.5)
      )
  }
  
  qq2003 <- makeQQ(residuals(res1), "Dead Stem Cover (2003)")
  qq2016 <- makeQQ(residuals(res2), "Dead Stem Cover (2016)")
  
  ggsave(
    file.path(outDir, paste0(prefix, "_DHARMa_QQ.pdf")),
    qq2003 | qq2016,
    width = 14,
    height = 7,
    units = "in"
  )
  
  pdf(file.path(outDir, paste0(prefix, "_DHARMa_4panel.pdf")),
      width = 8,
      height = 10
  )
  op <- par(mfrow = c(4, 2))
  plot(res1)
  plot(res2)
  par(op)
  dev.off()
  
  invisible(list(
    res1 = res1,
    res2 = res2,
    pvals = pvals
  ))
}

RecBetaGumbelPredObj <- RecBetaGumbelPostPred(
  fit = RecBetaGumbelFitObj,
  Data = Data,
  X1 = X1,
  X2 = X2,
  ID = Data$ID,
  S = nrow(do.call(rbind, as.mcmc.list(RecBetaGumbelFitObj)))
)

RecBetaGumbelDHARMa <- GumbelDiagnostics(RecBetaGumbelPredObj, Data)

################################################
### Goodness of fit and dependence structure ###
################################################

RecBetaGumbelGOF <- function(fit, Data, X1, X2, B = 500) {
  mcmcMat <- do.call(rbind, as.mcmc.list(fit))
  n_draws <- nrow(mcmcMat)
  n <- nrow(Data)
  
  U1_draws <- matrix(NA_real_, n, n_draws)
  U2_draws <- matrix(NA_real_, n, n_draws)
  
  for (s in seq_len(n_draws)) {
    d <- mcmcMat[s, ]
    Beta1 <- d[grep("^Beta1\\[", names(d))]
    Beta2 <- d[grep("^Beta2\\[", names(d))]
    U1re <- d[grep("^U1\\[", names(d))]
    U2re <- d[grep("^U2\\[", names(d))]
    
    p1 <- plogis(X1 %*% Beta1 + U1re[Data$ID])
    p2 <- plogis(X2 %*% Beta2 + U2re[Data$ID])
    
    U1_draws[, s] <- pRectBeta(Data$Y1, p1, d["phi1"], d["rho1"])
    U2_draws[, s] <- pRectBeta(Data$Y2, p2, d["phi2"], d["rho2"])
  }
  
  u_seq <- seq(0.05, 0.95, by = 0.05)
  
  chi_fn <- function(U) {
    num <- sapply(u_seq, function(u) mean((U[, 1] > 1 - u) & (U[, 2] > 1 - u)))
    den <- sapply(u_seq, function(u) mean(U[, 2] > 1 - u))
    z <- num/den
    z[!is.finite(z)] <- NA_real_
    z
  }
  
  K_fn <- function(U) {
    sapply(u_seq, function(u) mean((U[, 1] <= u) & (U[, 2] <= u)))
  }
  
  chi_draws <- matrix(NA_real_, length(u_seq), n_draws)
  K_draws <- matrix(NA_real_, length(u_seq), n_draws)
  for (s in seq_len(n_draws)) {
    Umat <- cbind(U1_draws[, s], U2_draws[, s])
    chi_draws[, s] <- chi_fn(Umat)
    K_draws[, s] <- K_fn(Umat)
  }
  chi_emp <- rowMeans(chi_draws, na.rm = TRUE)
  K_emp <- rowMeans(K_draws, na.rm = TRUE)
  
  boot_id <- sample(seq_len(n_draws), B, replace = TRUE)
  chi_sim <- matrix(NA_real_, length(u_seq), B)
  K_sim <- matrix(NA_real_, length(u_seq), B)
  
  for (b in seq_len(B)) {
    tau_b <- mcmcMat[boot_id[b], "tau"]
    theta_b <- 1/(1 - tau_b)
    Usim <- rCopula(n, gumbelCopula(theta_b, dim = 2))
    chi_sim[, b] <- chi_fn(Usim)
    K_sim[, b] <- K_fn(Usim)
  }
  
  chi_lo <- apply(chi_sim, 1, quantile, 0.025, na.rm = TRUE)
  chi_hi <- apply(chi_sim, 1, quantile, 0.975, na.rm = TRUE)
  K_lo <- apply(K_sim, 1, quantile, 0.025, na.rm = TRUE)
  K_hi <- apply(K_sim, 1, quantile, 0.975, na.rm = TRUE)
  
  plot_env <- function(emp, lo, hi, ylab) {
    data.frame(u = u_seq, emp = emp, lo = lo, hi = hi) |>
      ggplot(aes(u, emp)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#6393ff", alpha = 0.25) +
      geom_line(colour = "red4", linewidth = 1) +
      scale_x_continuous(
        limits = c(0.05, 0.95),
        breaks = seq(0.05, 0.95, by = 0.1),
        labels = function(x) sprintf("%.2f", x),
        expand = c(0, 0.04)
      ) +
      labs(x = "u", y = ylab) +
      theme_classic(base_size = 15)
  }
  
  p_chi <- plot_env(chi_emp, chi_lo, chi_hi, expression(chi(u)))
  p_K <- plot_env(K_emp, K_lo, K_hi, expression(K(u)))
  
  if (!dir.exists("Manuscript/Output")) dir.create("Manuscript/Output", recursive = TRUE)
  
  ggsave("Manuscript/Output/DeadStemCover_Gumbel_ChiPlot.pdf",
         p_chi,
         width = 5.5, height = 4.5, device = cairo_pdf
  )
  ggsave("Manuscript/Output/DeadStemCover_Gumbel_KPlot.pdf",
         p_K,
         width = 5.5, height = 4.5, device = cairo_pdf
  )
  
  cat("Plots saved:\n - DeadStemCover_Gumbel_ChiPlot.pdf\n - DeadStemCover_Gumbel_KPlot.pdf\n")
}

RecBetaGumbelGOF(
  fit  = RecBetaGumbelMCMC,
  Data = Data,
  X1   = X1,
  X2   = X2
)

############################################################
############# Rectangular-beta + Clayton copula ############
############################################################

RecBetaClayModel <- "
  model {

    for (k in 1:NFix1) {
      Beta1[k] ~ dnorm(0, 0.0001)
    }
    for (k in 1:NFix2) {
      Beta2[k] ~ dnorm(0, 0.0001)
    }

    sigma1 ~ dt(0, 0.25, 2)T(0, )
    sigma2 ~ dt(0, 0.25, 2)T(0, )
    tau1 <- pow(sigma1, -2)
    tau2 <- pow(sigma2, -2)

    for (j in 1:NRand) {
      U1[j] ~ dnorm(0, tau1)
      U2[j] ~ dnorm(0, tau2)
    }

    phi1 ~ dunif(0, 1)
    rho1 ~ dgamma(0.0001, 0.0001)
    phi2 ~ dunif(0, 1)
    rho2 ~ dgamma(0.0001, 0.0001)

    tau ~ dunif(0, 1)
    theta_cop <- 2*tau/(1 - tau)

    for (i in 1:NObs) {

      LP1[i] <- inprod(X1[i, ], Beta1[]) + U1[ID[i]]
      LP2[i] <- inprod(X2[i, ], Beta2[]) + U2[ID[i]]
      p1[i] <- ilogit(LP1[i])
      p2[i] <- ilogit(LP2[i])

      S1_1[i] <- 1 - abs(2*p1[i] - 1)
      S2_1[i] <- (p1[i] - 0.5*phi1*S1_1[i])/(1 - phi1*S1_1[i])
      a1[i] <- rho1*S2_1[i]
      b1[i] <- rho1*(1 - S2_1[i])

      pdf1[i] <- phi1*S1_1[i] + (1 - phi1*S1_1[i])*dbeta(Y1[i], a1[i], b1[i])
      cdf1[i] <- phi1*S1_1[i]*Y1[i] + (1 - phi1*S1_1[i])*pbeta(Y1[i], a1[i], b1[i])

      S1_2[i] <- 1 - abs(2*p2[i] - 1)
      S2_2[i] <- (p2[i] - 0.5*phi2*S1_2[i])/(1 - phi2*S1_2[i])
      a2[i] <- rho2*S2_2[i]
      b2[i] <- rho2*(1 - S2_2[i])

      pdf2[i] <- phi2*S1_2[i] + (1 - phi2*S1_2[i])*dbeta(Y2[i], a2[i], b2[i])
      cdf2[i] <- phi2*S1_2[i]*Y2[i] + (1 - phi2*S1_2[i])*pbeta(Y2[i], a2[i], b2[i])

      A[i] <- pow(cdf1[i], -theta_cop) + pow(cdf2[i], -theta_cop) - 1
      logCop[i] <- log(theta_cop + 1) -
                   (theta_cop + 1)*(log(cdf1[i]) + log(cdf2[i])) -
                   (2 + 1/theta_cop)*log(A[i])

      LL[i] <- log(pdf1[i]) + log(pdf2[i]) + logCop[i]
      zeros[i] ~ dpois(BigC - LL[i])
    }
  }
"

RecBetaClayFit <- function(Data, X1, X2) {
  params <- c(
    "Beta1",
    "Beta2",
    "sigma1",
    "sigma2",
    "U1",
    "U2",
    "phi1",
    "rho1",
    "phi2",
    "rho2",
    "tau",
    "theta_cop"
  )
  
  run.jags(
    model = RecBetaClayModel,
    data = JAGSData,
    inits = InitsList,
    monitor = params,
    n.chains = NumChains,
    adapt = 1000,
    burnin = 15000,
    sample = 1000,
    thin = 25,
    method = "parallel",
    modules = "glm",
    factories = "bugs::MNormal sampler off",
    silent.jags = FALSE
  )
}

RecBetaClayRun <- FALSE
RecBetaClayFitFile <- "Manuscript/Output/Rect_Beta_Clayton_Fit.rds"

if (!file.exists(RecBetaClayFitFile)) {
  cat("No RDS found. Running model...\n")
  RecBetaClayFitObj <- RecBetaClayFit(Data, X1, X2)
  saveRDS(RecBetaClayFitObj, RecBetaClayFitFile)
} else {
  if (RecBetaClayRun) {
    cat("RDS exists but RecBetaClayRun = TRUE. Re-running...\n")
    RecBetaClayFitObj <- RecBetaClayFit(Data, X1, X2)
    saveRDS(RecBetaClayFitObj, RecBetaClayFitFile)
  } else {
    cat("Loading saved fit...\n")
    RecBetaClayFitObj <- readRDS(RecBetaClayFitFile)
  }
}

print(summary(RecBetaClayFitObj))

RecBetaClayMCMC <- as.mcmc.list(RecBetaClayFitObj)
RecBetaClayDraws <- do.call(rbind, RecBetaClayMCMC)

cat(
  "Samples:",
  nrow(RecBetaClayDraws),
  " Params:",
  ncol(RecBetaClayDraws),
  "\n"
)

#################################
### Posterior tail dependence ###
#################################

ggsave(
  filename = "Manuscript/Output/DeadStemCover_Clayton_Lambda_L.pdf",
  plot = ggplot(data.frame(lambdaU = 2^(-1/if ("theta_cop" %in% colnames(RecBetaClayDraws)) {
    RecBetaClayDraws[, "theta_cop"]
  } else {
    1/(1 - RecBetaClayDraws[, "tau"])
  })), aes(lambdaU)) +
    geom_histogram(
      bins = 35,
      colour = "grey30",
      fill = "#6BAED6",
      linewidth = 0.4
    ) +
    labs(
      x = expression(lambda[L]),
      y = "Posterior Frequency"
    ) +
    theme_classic(base_size = 18) +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.8),
      axis.ticks = element_line(colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.margin = margin(6, 6, 6, 6)
    ),
  width = 5.5,
  height = 4.5,
  device = cairo_pdf
)

###########################
### Marginal likelihood ###
###########################

RecBetaClayLogCop <- function(u, v, theta) {
  log(theta + 1) -
    (theta + 1)*(log(u) + log(v)) -
    (2 + 1/theta)*log(u^(-theta) + v^(-theta) - 1)
}

RecBetaClayLogLik <- function(theta, dList) {
  names(theta) <- colnames(RecBetaClayDraws)
  
  Beta1 <- theta[grep("^Beta1\\[", names(theta))]
  Beta2 <- theta[grep("^Beta2\\[", names(theta))]
  U1 <- theta[grep("^U1\\[", names(theta))]
  U2 <- theta[grep("^U2\\[", names(theta))]
  
  sigma1 <- theta["sigma1"]
  sigma2 <- theta["sigma2"]
  phi1 <- theta["phi1"]
  phi2 <- theta["phi2"]
  rho1 <- theta["rho1"]
  rho2 <- theta["rho2"]
  tau <- theta["tau"]
  theta_cop <- 2*tau/(1 - tau)
  
  ll <- 0
  for (i in seq_len(dList$NObs)) {
    lp1 <- sum(dList$X1[i, ]*Beta1) + U1[dList$ID[i]]
    lp2 <- sum(dList$X2[i, ]*Beta2) + U2[dList$ID[i]]
    p1 <- plogis(lp1)
    p2 <- plogis(lp2)
    
    s1 <- 1 - abs(2*p1 - 1)
    s2 <- (p1 - 0.5*phi1*s1)/(1 - phi1*s1)
    a1 <- rho1*s2
    b1 <- rho1*(1 - s2)
    pdf1 <- phi1*s1 + (1 - phi1*s1)*dbeta(dList$Y1[i], a1, b1)
    cdf1 <-
      phi1*s1*dList$Y1[i] + (1 - phi1*s1)*pbeta(dList$Y1[i], a1, b1)
    
    t1 <- 1 - abs(2*p2 - 1)
    t2 <- (p2 - 0.5*phi2*t1)/(1 - phi2*t1)
    a2 <- rho2*t2
    b2 <- rho2*(1 - t2)
    pdf2 <- phi2*t1 + (1 - phi2*t1)*dbeta(dList$Y2[i], a2, b2)
    cdf2 <-
      phi2*t1*dList$Y2[i] + (1 - phi2*t1)*pbeta(dList$Y2[i], a2, b2)
    
    log_cop <- RecBetaClayLogCop(cdf1, cdf2, theta_cop)
    
    ll <- ll + log(pdf1) + log(pdf2) + log_cop
  }
  ll
}

RecBetaClayLogPrior <- function(theta, dList) {
  names(theta) <- colnames(RecBetaClayDraws)
  
  Beta1 <- theta[grep("^Beta1\\[", names(theta))]
  Beta2 <- theta[grep("^Beta2\\[", names(theta))]
  U1 <- theta[grep("^U1\\[", names(theta))]
  U2 <- theta[grep("^U2\\[", names(theta))]
  
  sigma1 <- theta["sigma1"]
  sigma2 <- theta["sigma2"]
  phi1 <- theta["phi1"]
  phi2 <- theta["phi2"]
  rho1 <- theta["rho1"]
  rho2 <- theta["rho2"]
  tau <- theta["tau"]
  
  lp <- 0
  lp <- lp + sum(dnorm(Beta1, 0, 100, log = TRUE))
  lp <- lp + sum(dnorm(Beta2, 0, 100, log = TRUE))
  lp <- lp + sum(dnorm(U1, 0, sigma1, log = TRUE))
  lp <- lp + sum(dnorm(U2, 0, sigma2, log = TRUE))
  
  half_t <-
    function(x) {
      log(2) + dt(x/2, df = 2, log = TRUE) - log(2)
    }
  lp <- lp + half_t(sigma1) + half_t(sigma2)
  
  lp <-
    lp + dunif(phi1, 0, 1, log = TRUE) + dunif(phi2, 0, 1, log = TRUE)
  lp <-
    lp + dgamma(rho1, 0.0001, 0.0001, log = TRUE) + dgamma(rho2, 0.0001, 0.0001, log = TRUE)
  lp <- lp + dunif(tau, 0, 1, log = TRUE)
  
  lp
}

RecBetaClayLogPost <- function(theta, data) {
  lp <- RecBetaClayLogPrior(theta, data)
  if (!is.finite(lp)) {
    return(lp)
  }
  ll <- RecBetaClayLogLik(theta, data)
  if (!is.finite(ll)) {
    return(ll)
  }
  lp + ll
}

RecBetaClayCN <- colnames(RecBetaClayDraws)
p <- length(RecBetaClayCN)
RecBetaClayLB <- rep(-Inf, p)
RecBetaClayUB <- rep(Inf, p)

setBnd <- function(name, lo, hi) {
  i <- which(RecBetaClayCN == name)
  if (length(i) == 1) {
    RecBetaClayLB[i] <<- lo
    RecBetaClayUB[i] <<- hi
  }
}

setBnd("sigma1", 0, Inf)
setBnd("sigma2", 0, Inf)
setBnd("phi1", 0, 1)
setBnd("phi2", 0, 1)
setBnd("rho1", 0, Inf)
setBnd("rho2", 0, Inf)
setBnd("tau", 0, 1)

names(RecBetaClayLB) <- names(RecBetaClayUB) <- RecBetaClayCN

RecBetaClayList <- c(
  "RecBetaClayLogCop",
  "RecBetaClayLogLik",
  "RecBetaClayLogPrior",
  "RecBetaClayLogPost",
  "JAGSData",
  "RecBetaClayCN",
  "RecBetaClayDraws"
)

RecBetaClayBridgeFile <-
  "Manuscript/Output/Rect_Beta_Clayton_Bridge.rds"
RecBetaClayRunBridge <- TRUE

if (RecBetaClayRunBridge || !file.exists(RecBetaClayBridgeFile)) {
  RecBetaClayBridge <- bridge_sampler(
    method = "warp3",
    samples = RecBetaClayDraws,
    log_posterior = RecBetaClayLogPost,
    data = JAGSData,
    lb = RecBetaClayLB,
    ub = RecBetaClayUB,
    cores = 7,
    packages = "bridgesampling",
    varlist = RecBetaClayList,
    silent = FALSE,
    envir = .GlobalEnv
  )
  saveRDS(RecBetaClayBridge, RecBetaClayBridgeFile)
} else {
  RecBetaClayBridge <- readRDS(RecBetaClayBridgeFile)
}

cat(
  "Clayton copula log-marginal likelihood:",
  RecBetaClayBridge$logml,
  "\n"
)
print(summary(RecBetaClayBridge))

############################
### Residual diagnostics ###
############################

RecBetaClayPostPred <- function(fit, Data, X1, X2, ID, S = 1000) {
  mcmcMat <- do.call(rbind, as.mcmc.list(fit))
  draws <- mcmcMat[sample(seq_len(nrow(mcmcMat)), S, FALSE), ]
  
  n <- nrow(Data)
  sim1 <- matrix(NA_real_, n, S)
  sim2 <- matrix(NA_real_, n, S)
  fit1 <- numeric(n)
  fit2 <- numeric(n)
  
  for (s in seq_len(S)) {
    d <- draws[s, ]
    B1 <- d[grep("^Beta1\\[", names(d))]
    B2 <- d[grep("^Beta2\\[", names(d))]
    U1 <- d[grep("^U1\\[", names(d))]
    U2 <- d[grep("^U2\\[", names(d))]
    
    phi1 <- d["phi1"]
    rho1 <- d["rho1"]
    phi2 <- d["phi2"]
    rho2 <- d["rho2"]
    
    theta <- 2*d["tau"]/(1 - d["tau"])
    cop <- claytonCopula(theta, dim = 2)
    
    p1 <- plogis(X1 %*% B1 + U1[ID])
    p2 <- plogis(X2 %*% B2 + U2[ID])
    
    U <- rCopula(n, cop)
    
    sim1[, s] <- qRectBetaVec(U[, 1], p1, phi1, rho1)
    sim2[, s] <- qRectBetaVec(U[, 2], p2, phi2, rho2)
    
    if (s == 1) {
      fit1 <- p1
      fit2 <- p2
    }
  }
  
  list(
    sim1 = sim1,
    sim2 = sim2,
    fit1 = fit1,
    fit2 = fit2
  )
}

ClaytonDiagnostics <- function(predObj,
                               Data,
                               outDir = "Manuscript/Output",
                               prefix = "DeadStemCover_Clayton") {
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE)
  }
  
  keep1 <- !is.na(Data$Y1)
  keep2 <- !is.na(Data$Y2)
  
  res1 <- createDHARMa(predObj$sim1[keep1, , drop = FALSE],
                       Data$Y1[keep1], predObj$fit1[keep1],
                       integerResponse = FALSE
  )
  res2 <- createDHARMa(predObj$sim2[keep2, , drop = FALSE],
                       Data$Y2[keep2], predObj$fit2[keep2],
                       integerResponse = FALSE
  )
  
  pvals <- tibble(
    Response = c("DeadStem2003", "DeadStem2016"),
    Uniformity = c(
      testUniformity(res1)$p.value,
      testUniformity(res2)$p.value
    ),
    Dispersion = c(
      testDispersion(res1)$p.value,
      testDispersion(res2)$p.value
    ),
    Outliers = c(
      testOutliers(res1)$p.value,
      testOutliers(res2)$p.value
    )
  )
  
  write.csv(pvals,
            file.path(outDir, paste0(prefix, "_DHARMa_pvalues.csv")),
            row.names = FALSE
  )
  
  makeQQ <- function(r, ttl) {
    ggplot(
      data.frame(Expected = qunif(ppoints(r)), Observed = sort(r)),
      aes(Expected, Observed)
    ) +
      geom_point(size = 2, colour = "blue3") +
      geom_abline(
        slope = 1,
        intercept = 0,
        colour = "darkred",
        linewidth = 0.9
      ) +
      labs(title = ttl, x = "Expected Residual", y = "Observed Residual") +
      theme_classic(base_size = 18) +
      theme(
        axis.line = element_line(colour = "black", linewidth = 0.9),
        plot.title = element_text(hjust = 0.5)
      )
  }
  
  qq2003 <- makeQQ(residuals(res1), "Dead Stem Cover (2003)")
  qq2016 <- makeQQ(residuals(res2), "Dead Stem Cover (2016)")
  
  ggsave(
    file.path(outDir, paste0(prefix, "_DHARMa_QQ.pdf")),
    qq2003 | qq2016,
    width = 14,
    height = 7,
    units = "in"
  )
  
  pdf(file.path(outDir, paste0(prefix, "_DHARMa_4panel.pdf")),
      width = 8,
      height = 10
  )
  op <- par(mfrow = c(4, 2))
  plot(res1)
  plot(res2)
  par(op)
  dev.off()
  
  invisible(list(
    res1 = res1,
    res2 = res2,
    pvals = pvals
  ))
}

RecBetaClayPredObj <- RecBetaClayPostPred(
  fit = RecBetaClayFitObj,
  Data = Data,
  X1 = X1,
  X2 = X2,
  ID = Data$ID,
  S = nrow(do.call(rbind, as.mcmc.list(RecBetaClayFitObj)))
)

RecBetaClayDHARMa <- ClaytonDiagnostics(RecBetaClayPredObj, Data)

################################################
### Goodness of fit and dependence structure ###
################################################

RecBetaClaytonGOF <- function(fit, Data, X1, X2, B = 500) {
  mcmcMat <- do.call(rbind, as.mcmc.list(fit))
  n_draws <- nrow(mcmcMat)
  n <- nrow(Data)
  
  U1_draws <- matrix(NA_real_, n, n_draws)
  U2_draws <- matrix(NA_real_, n, n_draws)
  
  for (s in seq_len(n_draws)) {
    d <- mcmcMat[s, ]
    Beta1 <- d[grep("^Beta1\\[", names(d))]
    Beta2 <- d[grep("^Beta2\\[", names(d))]
    U1re <- d[grep("^U1\\[", names(d))]
    U2re <- d[grep("^U2\\[", names(d))]
    
    p1 <- plogis(X1 %*% Beta1 + U1re[Data$ID])
    p2 <- plogis(X2 %*% Beta2 + U2re[Data$ID])
    
    U1_draws[, s] <- pRectBeta(Data$Y1, p1, d["phi1"], d["rho1"])
    U2_draws[, s] <- pRectBeta(Data$Y2, p2, d["phi2"], d["rho2"])
  }
  
  u_seq <- seq(0.05, 0.95, by = 0.05)
  
  chiL_fn <- function(U) {
    num <- sapply(u_seq, function(u) mean((U[, 1] <= u) & (U[, 2] <= u)))
    den <- sapply(u_seq, function(u) mean(U[, 2] <= u))
    z <- num/den
    z[!is.finite(z)] <- NA_real_
    z
  }
  
  K_fn <- function(U) {
    sapply(u_seq, function(u) mean((U[, 1] <= u) & (U[, 2] <= u)))
  }
  
  chi_draws <- matrix(NA_real_, length(u_seq), n_draws)
  K_draws <- matrix(NA_real_, length(u_seq), n_draws)
  for (s in seq_len(n_draws)) {
    Umat <- cbind(U1_draws[, s], U2_draws[, s])
    chi_draws[, s] <- chiL_fn(Umat)
    K_draws[, s] <- K_fn(Umat)
  }
  chi_emp <- rowMeans(chi_draws, na.rm = TRUE)
  K_emp <- rowMeans(K_draws, na.rm = TRUE)
  
  boot_id <- sample(seq_len(n_draws), B, replace = TRUE)
  chi_sim <- matrix(NA_real_, length(u_seq), B)
  K_sim <- matrix(NA_real_, length(u_seq), B)
  
  for (b in seq_len(B)) {
    tau_b <- mcmcMat[boot_id[b], "tau"]
    theta_b <- 2*tau_b/(1 - tau_b)
    Usim <- rCopula(n, claytonCopula(theta_b, dim = 2))
    chi_sim[, b] <- chiL_fn(Usim)
    K_sim[, b] <- K_fn(Usim)
  }
  
  chi_lo <- apply(chi_sim, 1, quantile, 0.025, na.rm = TRUE)
  chi_hi <- apply(chi_sim, 1, quantile, 0.975, na.rm = TRUE)
  K_lo <- apply(K_sim, 1, quantile, 0.025, na.rm = TRUE)
  K_hi <- apply(K_sim, 1, quantile, 0.975, na.rm = TRUE)
  
  plot_env <- function(emp, lo, hi, ylab) {
    data.frame(u = u_seq, emp = emp, lo = lo, hi = hi) |>
      ggplot(aes(u, emp)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#6393ff", alpha = 0.25) +
      geom_line(colour = "red4", linewidth = 1) +
      scale_x_continuous(
        limits = c(0.05, 0.95),
        breaks = seq(0.05, 0.95, by = 0.1),
        labels = function(x) sprintf("%.2f", x),
        expand = c(0, 0.04)
      ) +
      labs(x = "u", y = ylab) +
      theme_classic(base_size = 15)
  }
  
  p_chi <- plot_env(chi_emp, chi_lo, chi_hi, expression(chi[L](u)))
  p_K <- plot_env(K_emp, K_lo, K_hi, expression(K(u)))
  
  if (!dir.exists("Manuscript/Output")) dir.create("Manuscript/Output", recursive = TRUE)
  
  ggsave("Manuscript/Output/DeadStemCover_Clayton_ChiPlot.pdf",
         p_chi,
         width = 5.5, height = 4.5, device = cairo_pdf
  )
  ggsave("Manuscript/Output/DeadStemCover_Clayton_KPlot.pdf",
         p_K,
         width = 5.5, height = 4.5, device = cairo_pdf
  )
  
  cat("Plots saved:\n - DeadStemCover_Clayton_ChiPlot.pdf\n - DeadStemCover_Clayton_KPlot.pdf\n")
}

RecBetaClaytonGOF(
  fit = RecBetaClayFitObj,
  Data = Data,
  X1 = X1,
  X2 = X2
)

########################################################
############# Independent rectangular beta #############
########################################################

RecBetaIndepModel <- "
  model {

    for (k in 1:NFix1) {
      Beta1[k] ~ dnorm(0, 0.0001)
    }
    for (k in 1:NFix2) {
      Beta2[k] ~ dnorm(0, 0.0001)
    }

    sigma1 ~ dt(0, 0.25, 2)T(0, )
    sigma2 ~ dt(0, 0.25, 2)T(0, )
    tau1 <- 1/(sigma1*sigma1)
    tau2 <- 1/(sigma2*sigma2)

    for (j in 1:NRand) {
      U1[j] ~ dnorm(0, tau1)
      U2[j] ~ dnorm(0, tau2)
    }

    phi1 ~ dunif(0, 1)
    rho1 ~ dgamma(0.0001, 0.0001)
    phi2 ~ dunif(0, 1)
    rho2 ~ dgamma(0.0001, 0.0001)

    for (i in 1:NObs) {

      LP1[i] <- inprod(X1[i, ], Beta1[]) + U1[ID[i]]
      LP2[i] <- inprod(X2[i, ], Beta2[]) + U2[ID[i]]
      p1[i] <- ilogit(LP1[i])
      p2[i] <- ilogit(LP2[i])

      Star1_1[i] <- 1 - abs(2*p1[i] - 1)
      Star2_1[i] <- (p1[i] - 0.5*phi1*Star1_1[i])/(1 - phi1*Star1_1[i])
      alpha1[i] <- rho1*Star2_1[i]
      beta1[i] <- rho1*(1 - Star2_1[i])
      pdf1[i] <- phi1*Star1_1[i] +
                 (1 - phi1*Star1_1[i])*dbeta(Y1[i], alpha1[i], beta1[i])

      Star1_2[i] <- 1 - abs(2*p2[i] - 1)
      Star2_2[i] <- (p2[i] - 0.5*phi2*Star1_2[i])/(1 - phi2*Star1_2[i])
      alpha2[i] <- rho2*Star2_2[i]
      beta2[i] <- rho2*(1 - Star2_2[i])
      pdf2[i] <- phi2*Star1_2[i] +
                 (1 - phi2*Star1_2[i])*dbeta(Y2[i], alpha2[i], beta2[i])

      LL[i] <- log(pdf1[i]) + log(pdf2[i])
      zeros[i] ~ dpois(BigC - LL[i])
    }
  }
"

RecBetaIndepFit <- function(Data, X1, X2) {
  params <- c(
    "Beta1",
    "Beta2",
    "sigma1",
    "sigma2",
    "U1",
    "U2",
    "phi1",
    "rho1",
    "phi2",
    "rho2"
  )
  
  run.jags(
    model = RecBetaIndepModel,
    data = JAGSData,
    inits = InitsList,
    monitor = params,
    n.chains = NumChains,
    adapt = 1000,
    burnin = 15000,
    sample = 1000,
    thin = 25,
    method = "parallel",
    modules = "glm",
    factories = "bugs::MNormal sampler off",
    silent.jags = FALSE
  )
}

RecBetaIndepRun <- FALSE
RecBetaIndepFitFile <- "Manuscript/Output/Rect_Beta_Indep_Fit.rds"

if (!file.exists(RecBetaIndepFitFile)) {
  cat("No RDS found. Running model...\n")
  RecBetaIndepFitObj <- RecBetaIndepFit(Data, X1, X2)
  saveRDS(RecBetaIndepFitObj, RecBetaIndepFitFile)
} else {
  if (RecBetaIndepRun) {
    cat("RDS exists but RecBetaIndepRun = TRUE. Re-running...\n")
    RecBetaIndepFitObj <- RecBetaIndepFit(Data, X1, X2)
    saveRDS(RecBetaIndepFitObj, RecBetaIndepFitFile)
  } else {
    cat("Loading saved fit...\n")
    RecBetaIndepFitObj <- readRDS(RecBetaIndepFitFile)
  }
}

print(summary(RecBetaIndepFitObj))

RecBetaIndepMCMCList <- as.mcmc.list(RecBetaIndepFitObj)
RecBetaIndepDraws <- do.call(rbind, RecBetaIndepMCMCList)
cat(
  "Samples:",
  nrow(RecBetaIndepDraws),
  " Params:",
  ncol(RecBetaIndepDraws),
  "\n"
)

###############################
### Log marginal likelihood ###
###############################

RecBetaIndepLogLik <- function(theta, dList) {
  names(theta) <- colnames(RecBetaIndepDraws)
  
  Beta1 <- theta[grep("^Beta1\\[", names(theta))]
  Beta2 <- theta[grep("^Beta2\\[", names(theta))]
  U1 <- theta[grep("^U1\\[", names(theta))]
  U2 <- theta[grep("^U2\\[", names(theta))]
  
  sigma1 <- theta["sigma1"]
  sigma2 <- theta["sigma2"]
  phi1 <- theta["phi1"]
  phi2 <- theta["phi2"]
  rho1 <- theta["rho1"]
  rho2 <- theta["rho2"]
  
  ll <- 0
  for (i in seq_len(dList$NObs)) {
    lp1 <- sum(dList$X1[i, ]*Beta1) + U1[dList$ID[i]]
    lp2 <- sum(dList$X2[i, ]*Beta2) + U2[dList$ID[i]]
    p1 <- 1/(1 + exp(-lp1))
    p2 <- 1/(1 + exp(-lp2))
    
    s1 <- 1 - abs(2*p1 - 1)
    s2 <- (p1 - 0.5*phi1*s1)/(1 - phi1*s1)
    a1 <- rho1*s2
    b1 <- rho1*(1 - s2)
    pdf1 <- phi1*s1 + (1 - phi1*s1)*dbeta(dList$Y1[i], a1, b1)
    
    t1 <- 1 - abs(2*p2 - 1)
    t2 <- (p2 - 0.5*phi2*t1)/(1 - phi2*t1)
    a2 <- rho2*t2
    b2 <- rho2*(1 - t2)
    pdf2 <- phi2*t1 + (1 - phi2*t1)*dbeta(dList$Y2[i], a2, b2)
    
    part <- log(pdf1) + log(pdf2)
    ll <- ll + part
  }
  ll
}

RecBetaIndepLogPrior <- function(theta, dList) {
  names(theta) <- colnames(RecBetaIndepDraws)
  
  Beta1 <- theta[grep("^Beta1\\[", names(theta))]
  Beta2 <- theta[grep("^Beta2\\[", names(theta))]
  U1 <- theta[grep("^U1\\[", names(theta))]
  U2 <- theta[grep("^U2\\[", names(theta))]
  
  sigma1 <- theta["sigma1"]
  sigma2 <- theta["sigma2"]
  phi1 <- theta["phi1"]
  phi2 <- theta["phi2"]
  rho1 <- theta["rho1"]
  rho2 <- theta["rho2"]
  
  lp <- 0
  lp <- lp + sum(dnorm(Beta1, 0, 100, log = TRUE))
  lp <- lp + sum(dnorm(Beta2, 0, 100, log = TRUE))
  lp <- lp + sum(dnorm(U1, 0, sigma1, log = TRUE))
  lp <- lp + sum(dnorm(U2, 0, sigma2, log = TRUE))
  
  half_t_log <- function(x, df = 2, scale = 2) {
    log(2) + dt(x/scale, df = df, log = TRUE) - log(scale)
  }
  lp <- lp + half_t_log(sigma1)
  lp <- lp + half_t_log(sigma2)
  
  lp <- lp + dunif(phi1, 0, 1, log = TRUE)
  lp <- lp + dunif(phi2, 0, 1, log = TRUE)
  lp <- lp + dgamma(rho1, 0.0001, 0.0001, log = TRUE)
  lp <- lp + dgamma(rho2, 0.0001, 0.0001, log = TRUE)
  
  lp
}

RecBetaIndepLogPost <- function(theta, dList) {
  lp <- RecBetaIndepLogPrior(theta, dList)
  if (!is.finite(lp)) {
    return(lp)
  }
  ll <- RecBetaIndepLogLik(theta, dList)
  if (!is.finite(ll)) {
    return(ll)
  }
  lp + ll
}

RecBetaIndepNTheta <- ncol(RecBetaIndepDraws)
RecBetaIndepLB <- rep(-Inf, RecBetaIndepNTheta)
RecBetaIndepUB <- rep(Inf, RecBetaIndepNTheta)
RecBetaIndepCN <- colnames(RecBetaIndepDraws)

RecBetaIndepBound <- function(name, lo, hi) {
  i <- which(RecBetaIndepCN == name)
  if (length(i) == 1) {
    RecBetaIndepLB[i] <<- lo
    RecBetaIndepUB[i] <<- hi
  }
}
RecBetaIndepBound("sigma1", 0, Inf)
RecBetaIndepBound("sigma2", 0, Inf)
RecBetaIndepBound("phi1", 0, 1)
RecBetaIndepBound("phi2", 0, 1)
RecBetaIndepBound("rho1", 0, Inf)
RecBetaIndepBound("rho2", 0, Inf)
names(RecBetaIndepLB) <- RecBetaIndepCN
names(RecBetaIndepUB) <- RecBetaIndepCN

RecBetaIndepPostWrap <- function(par, data) {
  names(par) <- RecBetaIndepCN
  RecBetaIndepLogPost(par, data)
}

RecBetaIndepList <- c(
  "RecBetaIndepLogLik",
  "RecBetaIndepLogPrior",
  "RecBetaIndepLogPost",
  "JAGSData",
  "RecBetaIndepCN",
  "RecBetaIndepDraws"
)

RecBetaIndepRunBridge <- TRUE
RecBetaIndepBridgeFile <-
  "Manuscript/Output/Rect_Beta_Indep_Bridge.rds"

if (RecBetaIndepRunBridge || !file.exists(RecBetaIndepBridgeFile)) {
  RecBetaIndepBridge <- bridge_sampler(
    method = "warp3",
    samples = RecBetaIndepDraws,
    log_posterior = RecBetaIndepPostWrap,
    data = JAGSData,
    lb = RecBetaIndepLB,
    ub = RecBetaIndepUB,
    cores = 7,
    packages = c("bridgesampling"),
    varlist = RecBetaIndepList,
    silent = FALSE,
    envir = .GlobalEnv
  )
  saveRDS(RecBetaIndepBridge, RecBetaIndepBridgeFile)
} else {
  RecBetaIndepBridge <- readRDS(RecBetaIndepBridgeFile)
}

cat("log-marginal likelihood:", RecBetaIndepBridge$logml, "\n")
print(summary(RecBetaIndepBridge))

############################
### Residual diagnostics ###
############################

RecBetaIndepPostPred <- function(fit, Data, X1, X2, ID, S = 1000) {
  draws <- do.call(rbind, as.mcmc.list(fit))
  sel <- sample(seq_len(nrow(draws)), S, FALSE)
  
  n <- nrow(Data)
  sim1 <- matrix(NA_real_, n, S)
  sim2 <- matrix(NA_real_, n, S)
  fit1 <- numeric(n)
  fit2 <- numeric(n)
  
  for (s in seq_len(S)) {
    d <- draws[sel[s], ]
    B1 <- d[grep("^Beta1\\[", names(d))]
    B2 <- d[grep("^Beta2\\[", names(d))]
    U1 <- d[grep("^U1\\[", names(d))]
    U2 <- d[grep("^U2\\[", names(d))]
    
    phi1 <- d["phi1"]
    rho1 <- d["rho1"]
    phi2 <- d["phi2"]
    rho2 <- d["rho2"]
    
    p1 <- plogis(X1 %*% B1 + U1[ID])
    p2 <- plogis(X2 %*% B2 + U2[ID])
    
    Umat <- matrix(runif(2*n), ncol = 2)
    
    sim1[, s] <- qRectBetaVec(Umat[, 1], p1, phi1, rho1)
    sim2[, s] <- qRectBetaVec(Umat[, 2], p2, phi2, rho2)
    
    if (s == 1) {
      fit1 <- p1
      fit2 <- p2
    }
  }
  
  list(
    sim1 = sim1,
    sim2 = sim2,
    fit1 = fit1,
    fit2 = fit2
  )
}

IndepDiagnostics <- function(predObj,
                             Data,
                             outDir = "Manuscript/Output",
                             prefix = "DeadStemCover_Indep") {
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE)
  }
  
  keep1 <- !is.na(Data$Y1)
  keep2 <- !is.na(Data$Y2)
  
  res1 <- createDHARMa(predObj$sim1[keep1, , drop = FALSE],
                       Data$Y1[keep1], predObj$fit1[keep1],
                       integerResponse = FALSE
  )
  res2 <- createDHARMa(predObj$sim2[keep2, , drop = FALSE],
                       Data$Y2[keep2], predObj$fit2[keep2],
                       integerResponse = FALSE
  )
  
  pvals <- tibble(
    Response = c("DeadStem2003", "DeadStem2016"),
    Uniformity = c(
      testUniformity(res1)$p.value,
      testUniformity(res2)$p.value
    ),
    Dispersion = c(
      testDispersion(res1)$p.value,
      testDispersion(res2)$p.value
    ),
    Outliers = c(
      testOutliers(res1)$p.value,
      testOutliers(res2)$p.value
    )
  )
  
  write.csv(pvals,
            file.path(outDir, paste0(prefix, "_DHARMa_pvalues.csv")),
            row.names = FALSE
  )
  
  makeQQ <- function(r, ttl) {
    ggplot(
      data.frame(Expected = qunif(ppoints(r)), Observed = sort(r)),
      aes(Expected, Observed)
    ) +
      geom_point(size = 2, colour = "blue3") +
      geom_abline(
        slope = 1,
        intercept = 0,
        colour = "darkred",
        linewidth = 0.9
      ) +
      labs(title = ttl, x = "Expected Residual", y = "Observed Residual") +
      theme_classic(base_size = 18) +
      theme(
        axis.line = element_line(colour = "black", linewidth = 0.9),
        plot.title = element_text(hjust = 0.5)
      )
  }
  
  qq2003 <- makeQQ(residuals(res1), "Dead Stem Cover (2003)")
  qq2016 <- makeQQ(residuals(res2), "Dead Stem Cover (2016)")
  
  ggsave(
    file.path(outDir, paste0(prefix, "_DHARMa_QQ.pdf")),
    qq2003 | qq2016,
    width = 14,
    height = 7,
    units = "in"
  )
  
  pdf(file.path(outDir, paste0(prefix, "_DHARMa_4panel.pdf")),
      width = 8,
      height = 10
  )
  op <- par(mfrow = c(4, 2))
  plot(res1)
  plot(res2)
  par(op)
  dev.off()
  
  invisible(list(
    res1 = res1,
    res2 = res2,
    pvals = pvals
  ))
}

RecBetaIndepPredObj <- RecBetaIndepPostPred(
  fit = RecBetaIndepFitObj,
  Data = Data,
  X1 = X1,
  X2 = X2,
  ID = Data$ID,
  S = nrow(do.call(rbind, as.mcmc.list(RecBetaIndepFitObj)))
)

RecBetaIndepDHARMa <- IndepDiagnostics(RecBetaIndepPredObj, Data)

#########################################################
############# Independent conventional beta #############
#########################################################

BetaIndepModel <- "
  model {

    for (k in 1:NFix1) {
      Beta1[k] ~ dnorm(0, 0.0001)
    }
    for (k in 1:NFix2) {
      Beta2[k] ~ dnorm(0, 0.0001)
    }

    sigma1 ~ dt(0, 0.25, 2)T(0, )
    sigma2 ~ dt(0, 0.25, 2)T(0, )
    tau1 <- pow(sigma1, -2)
    tau2 <- pow(sigma2, -2)

    for (j in 1:NRand) {
      U1[j] ~ dnorm(0, tau1)
      U2[j] ~ dnorm(0, tau2)
    }

    phi1 ~ dgamma(0.0001, 0.0001)
    phi2 ~ dgamma(0.0001, 0.0001)

    for (i in 1:NObs) {
      LP1[i] <- inprod(X1[i, ], Beta1[]) + U1[ID[i]]
      LP2[i] <- inprod(X2[i, ], Beta2[]) + U2[ID[i]]
      mu1[i] <- ilogit(LP1[i])
      mu2[i] <- ilogit(LP2[i])

      alpha1[i] <- mu1[i]*phi1
      beta1[i] <- (1 - mu1[i])*phi1
      alpha2[i] <- mu2[i]*phi2
      beta2[i] <- (1 - mu2[i])*phi2

      pdf1[i] <- dbeta(Y1[i], alpha1[i], beta1[i])
      pdf2[i] <- dbeta(Y2[i], alpha2[i], beta2[i])

      LL[i] <- log(pdf1[i]) + log(pdf2[i])
      zeros[i] ~ dpois(BigC - LL[i])
    }
  }
"

BetaIndepFit <- function(Data, X1, X2) {
  params <-
    c(
      "Beta1",
      "Beta2",
      "sigma1",
      "sigma2",
      "U1",
      "U2",
      "phi1",
      "phi2"
    )
  
  run.jags(
    model = BetaIndepModel,
    data = JAGSData,
    inits = InitsList,
    monitor = params,
    n.chains = NumChains,
    adapt = 1000,
    burnin = 15000,
    sample = 1000,
    thin = 25,
    method = "parallel",
    modules = "glm",
    factories = "bugs::MNormal sampler off",
    silent.jags = FALSE
  )
}

BetaIndepRun <- FALSE
BetaIndepFitFile <- "Manuscript/Output/Plain_Beta_Indep_Fit.rds"

if (!file.exists(BetaIndepFitFile)) {
  cat("No RDS found. Running model...\n")
  BetaIndepFitObj <- BetaIndepFit(Data, X1, X2)
  saveRDS(BetaIndepFitObj, BetaIndepFitFile)
} else {
  if (BetaIndepRun) {
    cat("RDS exists but BetaIndepRun = TRUE. Re-running...\n")
    BetaIndepFitObj <- BetaIndepFit(Data, X1, X2)
    saveRDS(BetaIndepFitObj, BetaIndepFitFile)
  } else {
    cat("Loading saved fit...\n")
    BetaIndepFitObj <- readRDS(BetaIndepFitFile)
  }
}

print(summary(BetaIndepFitObj))

BetaIndepMCMC <- as.mcmc.list(BetaIndepFitObj)
BetaIndepDraws <- do.call(rbind, BetaIndepMCMC)

cat(
  "Samples:",
  nrow(BetaIndepDraws),
  " Params:",
  ncol(BetaIndepDraws),
  "\n"
)

###############################
### Log marginal likelihood ###
###############################

BetaIndepLogLik <- function(theta, dList) {
  names(theta) <- colnames(BetaIndepDraws)
  
  Beta1 <- theta[grep("^Beta1\\[", names(theta))]
  Beta2 <- theta[grep("^Beta2\\[", names(theta))]
  U1 <- theta[grep("^U1\\[", names(theta))]
  U2 <- theta[grep("^U2\\[", names(theta))]
  
  sigma1 <- theta["sigma1"]
  sigma2 <- theta["sigma2"]
  phi1 <- theta["phi1"]
  phi2 <- theta["phi2"]
  
  ll <- 0
  for (i in seq_len(dList$NObs)) {
    lp1 <- sum(dList$X1[i, ]*Beta1) + U1[dList$ID[i]]
    lp2 <- sum(dList$X2[i, ]*Beta2) + U2[dList$ID[i]]
    mu1 <- plogis(lp1)
    mu2 <- plogis(lp2)
    
    a1 <- mu1*phi1
    b1 <- (1 - mu1)*phi1
    a2 <- mu2*phi2
    b2 <- (1 - mu2)*phi2
    
    ll <- ll +
      dbeta(dList$Y1[i], a1, b1, log = TRUE) +
      dbeta(dList$Y2[i], a2, b2, log = TRUE)
  }
  ll
}

BetaIndepLogPrior <- function(theta, dList) {
  names(theta) <- colnames(BetaIndepDraws)
  
  Beta1 <- theta[grep("^Beta1\\[", names(theta))]
  Beta2 <- theta[grep("^Beta2\\[", names(theta))]
  U1 <- theta[grep("^U1\\[", names(theta))]
  U2 <- theta[grep("^U2\\[", names(theta))]
  
  sigma1 <- theta["sigma1"]
  sigma2 <- theta["sigma2"]
  phi1 <- theta["phi1"]
  phi2 <- theta["phi2"]
  
  lp <- 0
  lp <- lp + sum(dnorm(Beta1, 0, 100, log = TRUE))
  lp <- lp + sum(dnorm(Beta2, 0, 100, log = TRUE))
  lp <- lp + sum(dnorm(U1, 0, sigma1, log = TRUE))
  lp <- lp + sum(dnorm(U2, 0, sigma2, log = TRUE))
  
  half_t <-
    function(x) {
      log(2) + dt(x/2, df = 2, log = TRUE) - log(2)
    }
  lp <- lp + half_t(sigma1) + half_t(sigma2)
  
  lp <- lp + dgamma(phi1, 0.0001, 0.0001, log = TRUE)
  lp <- lp + dgamma(phi2, 0.0001, 0.0001, log = TRUE)
  
  lp
}

BetaIndepLogPost <- function(theta, data) {
  lp <- BetaIndepLogPrior(theta, data)
  if (!is.finite(lp)) {
    return(lp)
  }
  ll <- BetaIndepLogLik(theta, data)
  if (!is.finite(ll)) {
    return(ll)
  }
  lp + ll
}

BetaIndepCN <- colnames(BetaIndepDraws)
p <- length(BetaIndepCN)
BetaIndepLB <- rep(-Inf, p)
BetaIndepUB <- rep(Inf, p)

setBnd <- function(name, lo, hi) {
  i <- which(BetaIndepCN == name)
  if (length(i) == 1) {
    BetaIndepLB[i] <<- lo
    BetaIndepUB[i] <<- hi
  }
}

setBnd("sigma1", 0, Inf)
setBnd("sigma2", 0, Inf)
setBnd("phi1", 0, Inf)
setBnd("phi2", 0, Inf)

names(BetaIndepLB) <- names(BetaIndepUB) <- BetaIndepCN

BetaIndepList <- c(
  "BetaIndepLogLik",
  "BetaIndepLogPrior",
  "BetaIndepLogPost",
  "JAGSData",
  "BetaIndepCN",
  "BetaIndepDraws"
)

BetaIndepBridgeFile <-
  "Manuscript/Output/Plain_Beta_Indep_Bridge.rds"
BetaIndepRunBridge <- TRUE

if (BetaIndepRunBridge || !file.exists(BetaIndepBridgeFile)) {
  BetaIndepBridge <- bridge_sampler(
    method = "warp3",
    samples = BetaIndepDraws,
    log_posterior = function(par, data) {
      names(par) <- BetaIndepCN
      BetaIndepLogPost(par, data)
    },
    data = JAGSData,
    lb = BetaIndepLB,
    ub = BetaIndepUB,
    cores = 7,
    packages = "bridgesampling",
    varlist = BetaIndepList,
    silent = FALSE,
    envir = .GlobalEnv
  )
  saveRDS(BetaIndepBridge, BetaIndepBridgeFile)
} else {
  BetaIndepBridge <- readRDS(BetaIndepBridgeFile)
}

cat(
  "Plain-Beta (indep) log-marginal likelihood:",
  BetaIndepBridge$logml,
  "\n"
)
print(summary(BetaIndepBridge))

############################
### Residual diagnostics ###
############################

qConvBeta <- function(u, mu, phi) {
  qbeta(u, mu*phi, (1 - mu)*phi)
}

BetaIndepPostPred <- function(fit, Data, X1, X2, ID, S = 1000) {
  mcmcMat <- do.call(rbind, as.mcmc.list(fit))
  draws <- mcmcMat[sample(seq_len(nrow(mcmcMat)), S, FALSE), ]
  
  n <- nrow(Data)
  sim1 <- matrix(NA_real_, n, S)
  sim2 <- matrix(NA_real_, n, S)
  fit1 <- numeric(n)
  fit2 <- numeric(n)
  
  for (s in seq_len(S)) {
    d <- draws[s, ]
    B1 <- d[grep("^Beta1\\[", names(d))]
    B2 <- d[grep("^Beta2\\[", names(d))]
    U1 <- d[grep("^U1\\[", names(d))]
    U2 <- d[grep("^U2\\[", names(d))]
    
    phi1 <- d["phi1"]
    phi2 <- d["phi2"]
    
    mu1 <- plogis(X1 %*% B1 + U1[ID])
    mu2 <- plogis(X2 %*% B2 + U2[ID])
    
    if ("tau" %in% names(d)) {
      theta <- sin(0.5*pi*d["tau"])
      R <- matrix(c(1, theta, theta, 1), 2)
      Umat <- pnorm(rmvnorm(n, sigma = R))
    } else {
      Umat <- matrix(runif(2*n), ncol = 2)
    }
    
    sim1[, s] <- qConvBeta(Umat[, 1], mu1, phi1)
    sim2[, s] <- qConvBeta(Umat[, 2], mu2, phi2)
    
    if (s == 1) {
      fit1 <- mu1
      fit2 <- mu2
    }
  }
  
  list(
    sim1 = sim1,
    sim2 = sim2,
    fit1 = fit1,
    fit2 = fit2
  )
}

PlainBetaDiagnostics <- function(predObj,
                                 Data,
                                 outDir = "Manuscript/Output",
                                 prefix = "DeadStemCover_PlainBeta") {
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE)
  }
  
  keep1 <- !is.na(Data$Y1)
  keep2 <- !is.na(Data$Y2)
  
  res1 <- createDHARMa(predObj$sim1[keep1, , drop = FALSE],
                       Data$Y1[keep1], predObj$fit1[keep1],
                       integerResponse = FALSE
  )
  res2 <- createDHARMa(predObj$sim2[keep2, , drop = FALSE],
                       Data$Y2[keep2], predObj$fit2[keep2],
                       integerResponse = FALSE
  )
  
  pvals <- tibble(
    Response = c("DeadStem2003", "DeadStem2016"),
    Uniformity = c(
      testUniformity(res1)$p.value,
      testUniformity(res2)$p.value
    ),
    Dispersion = c(
      testDispersion(res1)$p.value,
      testDispersion(res2)$p.value
    ),
    Outliers = c(
      testOutliers(res1)$p.value,
      testOutliers(res2)$p.value
    )
  )
  
  write.csv(pvals,
            file.path(outDir, paste0(prefix, "_DHARMa_pvalues.csv")),
            row.names = FALSE
  )
  
  makeQQ <- function(r, ttl) {
    ggplot(
      data.frame(Expected = qunif(ppoints(r)), Observed = sort(r)),
      aes(Expected, Observed)
    ) +
      geom_point(size = 2, colour = "blue3") +
      geom_abline(
        slope = 1,
        intercept = 0,
        colour = "darkred",
        linewidth = 0.9
      ) +
      labs(title = ttl, x = "Expected Residual", y = "Observed Residual") +
      theme_classic(base_size = 18) +
      theme(
        axis.line = element_line(colour = "black", linewidth = 0.9),
        plot.title = element_text(hjust = 0.5)
      )
  }
  
  qq2003 <- makeQQ(residuals(res1), "Dead Stem Cover (2003)")
  qq2016 <- makeQQ(residuals(res2), "Dead Stem Cover (2016)")
  
  ggsave(
    file.path(outDir, paste0(prefix, "_DHARMa_QQ.pdf")),
    qq2003 | qq2016,
    width = 14,
    height = 7,
    units = "in"
  )
  
  pdf(file.path(outDir, paste0(prefix, "_DHARMa_4panel.pdf")),
      width = 8,
      height = 10
  )
  op <- par(mfrow = c(4, 2))
  plot(res1)
  plot(res2)
  par(op)
  dev.off()
  
  invisible(list(
    res1 = res1,
    res2 = res2,
    pvals = pvals
  ))
}

RecBetaBetaIndepPredObj <- BetaIndepPostPred(
  fit = BetaIndepFitObj,
  Data = Data,
  X1 = X1,
  X2 = X2,
  ID = Data$ID,
  S = nrow(do.call(rbind, as.mcmc.list(BetaIndepFitObj)))
)

RecBetaBetaIndepDHARMa <- PlainBetaDiagnostics(RecBetaBetaIndepPredObj, Data)

##################################################
### Summarize beta estimates and HPD intervals ###
##################################################

ModelNames <- c(
  Beta = "Beta[Indep]",
  Indep = "RectBeta[Indep]",
  Gauss = "RectBeta[Gauss]",
  Gumbel = "RectBeta[Gumbel]",
  Clay = "RectBeta[Clayton]"
)

SummarizeBetas <- function(draws, tag) {
  keep <- grep("^Beta[12]\\[", colnames(draws), value = TRUE)
  hpd <- HPDinterval(as.mcmc(draws[, keep]), prob = 0.95)
  data.frame(
    Model = ModelNames[tag],
    Beta = keep,
    Med = apply(draws[, keep], 2, median),
    Lo = hpd[, "lower"],
    Hi = hpd[, "upper"],
    Endpoint = ifelse(
      grepl("^Beta1", keep),
      "Dead Stem Cover (2003)",
      "Dead Stem Cover (2016)"
    ),
    check.names = FALSE,
    row.names = NULL
  )
}

DataBetaPlot <- rbind(
  SummarizeBetas(BetaIndepDraws, "Beta"),
  SummarizeBetas(RecBetaIndepDraws, "Indep"),
  SummarizeBetas(RecBetaGaussDraws, "Gauss"),
  SummarizeBetas(RecBetaGumbelDraws, "Gumbel"),
  SummarizeBetas(RecBetaClayDraws, "Clay")
)

DataBetaPlot$Model <- factor(
  DataBetaPlot$Model,
  levels = c(
    "Beta[Indep]",
    "RectBeta[Indep]",
    "RectBeta[Gauss]",
    "RectBeta[Gumbel]",
    "RectBeta[Clayton]"
  )
)

DataBetaPlot$CoefLab <- vapply(
  DataBetaPlot$Beta,
  function(s) {
    m <- regexec("^Beta([12])\\[(\\d+)\\]$", s)
    mt <- regmatches(s, m)[[1]]
    j <- mt[2]
    p <- as.integer(mt[3]) - 1
    sprintf("beta[%s%d]", j, p)
  },
  character(1)
)

ThemeBeta <- theme_bw(base_size = 14) +
  theme(
    strip.background = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.8, "cm"),
    legend.margin = margin(t = -5),
    legend.box.margin = margin(t = -10),
    plot.margin = margin(5, 5, 5, 5)
  )

BetaColor <- c(
  "Beta[Indep]" = "#B2182B",
  "RectBeta[Indep]" = "#EF3B2C",
  "RectBeta[Gauss]" = "#08519C",
  "RectBeta[Gumbel]" = "#3182BD",
  "RectBeta[Clayton]" = "#6BAED6"
)

BetaShape <- c(
  "Beta[Indep]" = 15,
  "RectBeta[Indep]" = 15,
  "RectBeta[Gauss]" = 16,
  "RectBeta[Gumbel]" = 16,
  "RectBeta[Clayton]" = 16
)

BetaLinetype <- c(
  "Beta[Indep]" = "solid",
  "RectBeta[Indep]" = "solid",
  "RectBeta[Gauss]" = "dashed",
  "RectBeta[Gumbel]" = "dashed",
  "RectBeta[Clayton]" = "dashed"
)

BetaLegendLabs <- lapply(
  names(BetaColor),
  function(lbl) {
    fixed <- sub("^([A-Za-z]+)\\[", "plain(\"\\1\")[", lbl)
    parse(text = fixed)
  }
)

MakeBetaPlot <- function(ep) {
  ggplot(
    DataBetaPlot[DataBetaPlot$Endpoint == ep, ],
    aes(
      x = Model,
      y = Med,
      ymin = Lo,
      ymax = Hi,
      colour = Model,
      shape = Model
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dotted",
      colour = "grey50",
      linewidth = 0.8
    ) +
    geom_linerange(aes(linetype = Model),
                   linewidth = 0.8,
                   position = position_dodge(width = 0.6)
    ) +
    geom_point(size = 3, position = position_dodge(width = 0.6)) +
    facet_wrap(~CoefLab,
               nrow = 2,
               scales = "free_y",
               labeller = label_parsed
    ) +
    scale_colour_manual(
      values = BetaColor,
      labels = BetaLegendLabs,
      breaks = names(BetaColor)
    ) +
    scale_shape_manual(values = BetaShape, breaks = names(BetaShape)) +
    scale_linetype_manual(values = BetaLinetype, breaks = names(BetaLinetype)) +
    guides(
      colour = guide_legend(
        override.aes = list(
          shape = unname(BetaShape),
          linetype = unname(BetaLinetype),
          colour = unname(BetaColor),
          size = 3
        )
      ),
      shape = "none",
      linetype = "none"
    ) +
    labs(y = "Posterior Est. (95% HPD)", x = NULL) +
    ThemeBeta +
    coord_cartesian(clip = "off")
}

BetaFig2003 <- MakeBetaPlot("Dead Stem Cover (2003)")
BetaFig2016 <- MakeBetaPlot("Dead Stem Cover (2016)")

print(BetaFig2003)
print(BetaFig2016)

ggsave(
  "Manuscript/Output/DeadStemCover_2003_betas.pdf",
  BetaFig2003,
  width = 8,
  height = 5,
  device = cairo_pdf
)
ggsave(
  "Manuscript/Output/DeadStemCover_2016_betas.pdf",
  BetaFig2016,
  width = 8,
  height = 5,
  device = cairo_pdf
)

##################################################
### Model comparison and residual checks table ###
##################################################

pvals_combined <- bind_rows(
  RecBetaBetaIndepDHARMa$pvals %>% mutate(Model = "Beta[Indep]"),
  RecBetaIndepDHARMa$pvals %>% mutate(Model = "RectBeta[Indep]"),
  RecBetaGaussDHARMa$pvals %>% mutate(Model = "RectBeta[Gauss]"),
  RecBetaGumbelDHARMa$pvals %>% mutate(Model = "RectBeta[Gumbel]"),
  RecBetaClayDHARMa$pvals %>% mutate(Model = "RectBeta[Clayton]")
)

format_pval <- function(x) ifelse(x < 0.001, "$<$0.001", sprintf("%.3f", x))

pvals_formatted <- pvals_combined %>%
  pivot_longer(
    cols = c("Uniformity", "Dispersion", "Outliers"),
    names_to = "Test", values_to = "pval"
  ) %>%
  mutate(pval = format_pval(pval)) %>%
  unite("Response_Test", Response, Test) %>%
  pivot_wider(names_from = Response_Test, values_from = pval)

logml_df <- tribble(
  ~Model, ~logml,
  "Beta[Indep]", BetaIndepBridge$logml,
  "RectBeta[Indep]", RecBetaIndepBridge$logml,
  "RectBeta[Gauss]", RecBetaGaussBridge$logml,
  "RectBeta[Gumbel]", RecBetaGumBridge$logml,
  "RectBeta[Clayton]", RecBetaClayBridge$logml
) %>%
  mutate(logml = sprintf("%.3f", logml))

final_df <- left_join(logml_df, pvals_formatted, by = "Model")

label_map <- c(
  "Beta[Indep]" = "$\\text{Beta}_{\\text{Indep}}$",
  "RectBeta[Indep]" = "$\\text{RectBeta}_{\\text{Indep}}$",
  "RectBeta[Gauss]" = "$\\text{RectBeta}_{\\text{Gauss}}$",
  "RectBeta[Gumbel]" = "$\\text{RectBeta}_{\\text{Gumbel}}$",
  "RectBeta[Clayton]" = "$\\text{RectBeta}_{\\text{Clayton}}$"
)

final_df$Model <- label_map[final_df$Model]

final_tab <- xtable(final_df,
                    caption = "Model comparison: log marginal likelihood and DHARMa p-values"
)

print(final_tab,
      file = "Manuscript/Output/DeadStemCover_model_comparison.tex",
      include.rownames = FALSE,
      sanitize.text.function = identity
)

###############################
### Posterior summary table ###
###############################

sums <- list(
  BetaIndep = summary(BetaIndepFitObj),
  RecBetaIndep = summary(RecBetaIndepFitObj),
  RecBetaGauss = summary(RecBetaGaussFitObj),
  RecBetaGumbel = summary(RecBetaGumbelFitObj),
  RecBetaClay = summary(RecBetaClayFitObj)
)

nice_param <- function(p) {
  sapply(p, function(x) {
    if (grepl("^Beta1\\[[1-6]\\]$", x)) {
      i <- as.integer(sub("^Beta1\\[(\\d)\\]$", "\\1", x)) - 1
      sprintf("$\\beta_{1%d}$", i)
    } else if (grepl("^Beta2\\[[1-6]\\]$", x)) {
      i <- as.integer(sub("^Beta2\\[(\\d)\\]$", "\\1", x)) - 1
      sprintf("$\\beta_{2%d}$", i)
    } else if (grepl("^sigma[12]$", x)) {
      j <- sub("^sigma([12])$", "\\1", x)
      sprintf("$\\sigma_{%s}$", j)
    } else if (grepl("^phi[12]$", x)) {
      j <- sub("^phi([12])$", "\\1", x)
      sprintf("$\\phi_{%s}$", j)
    } else if (grepl("^rho[12]$", x)) {
      j <- sub("^rho([12])$", "\\1", x)
      sprintf("$\\rho_{%s}$", j)
    } else if (x == "tau") {
      "$\\tau$"
    } else {
      x
    }
  }, USE.NAMES = FALSE)
}

keep_pat <- "^(Beta[12]\\[[1-6]\\]|sigma[12]|phi[12]|rho[12]|tau)$"

extract_model <- function(sm, model_name) {
  df <- as.data.frame(sm)
  df <- df[grep(keep_pat, rownames(df)), c("Median", "Lower95", "Upper95")]
  
  tibble(
    param = nice_param(rownames(df)),
    !!paste0(model_name, "_Median") := sprintf("%.3f", df$Median),
    !!paste0(model_name, "_Lower") := sprintf("[%.3f; ", df$Lower95),
    !!paste0(model_name, "_Upper") := sprintf("%.3f]", df$Upper95)
  )
}

combo <- Reduce(
  function(x, y) full_join(x, y, by = "param"),
  Map(extract_model, sums, names(sums))
) %>%
  arrange(param)

col_align <- c("l", rep("c", ncol(combo)))
xt <- xtable(combo, align = col_align)

print(xt,
      include.rownames = FALSE,
      sanitize.text.function = identity,
      file = "Manuscript/Output/DeadStemCover_posterior_summaries.tex"
)

RecBetaGaussFitObj$timetaken
RecBetaGumbelFitObj$timetaken
RecBetaClayFitObj$timetaken
RecBetaIndepFitObj$timetaken
BetaIndepFitObj$timetaken