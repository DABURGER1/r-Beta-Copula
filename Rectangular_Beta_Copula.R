library(copula)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(grid)
library(gtable)

setwd("C:/Users/a239866/OneDrive - Syneos Health/Divan @ Syneos Health - Linked Files/Research/Multivariate Outlier Proportions")

mu1 <- mu2 <- 0.50
rho <- 10
phi1_vec <- c(0, 0, 0.45, 0.45)
phi2_vec <- c(0.15, 0.45, 0.15, 0.40)
stopifnot(length(phi1_vec) == length(phi2_vec))

TauTarget <- 0.40
rho_c <- sin(pi*TauTarget/2)
rho_c

nu <- Inf
ng <- 800
prob_seq <- seq(0.10, 0.90, 0.10)

fills <- colorRampPalette(brewer.pal(9, "Blues"))(length(prob_seq) + 1)
band_lab <- c("< 10%", sprintf("%d–%d%%", seq(10, 80, 10), seq(20, 90, 10)), "≥ 90%")

xi <- function(mu) 1 - abs(2*mu - 1)
kfn <- function(mu, φ) (mu - 0.5*φ*xi(mu))/(1 - φ*xi(mu))
dRB <- function(y, mu, φ, ρ) {
  k <- kfn(mu, φ)
  φ*xi(mu) + (1 - φ*xi(mu))*dbeta(y, ρ*k, ρ*(1 - k))
}
pRB <- function(y, mu, φ, ρ) {
  k <- kfn(mu, φ)
  φ*xi(mu)*y + (1 - φ*xi(mu))*pbeta(y, ρ*k, ρ*(1 - k))
}

g <- seq(0.001, 0.999, length.out = ng)
grd_ref <- expand.grid(y1 = g, y2 = g)
tc <- tCopula(rho_c, dim = 2, df = nu)

u0 <- cbind(
  pRB(grd_ref$y1, 0.50, 0, rho),
  pRB(grd_ref$y2, 0.50, 0, rho)
)
z_ref <- dCopula(u0, tc) *
  dRB(grd_ref$y1, 0.50, 0, rho) *
  dRB(grd_ref$y2, 0.50, 0, rho)
z_ref <- matrix(z_ref, nrow = ng)

z_sorted <- sort(as.vector(z_ref), decreasing = TRUE)
cellA <- diff(g)[1]^2
cumA <- cumsum(z_sorted*cellA)
ref_breaks <- sapply(c(0.50, 0.75, 0.90), function(p) z_sorted[min(which(cumA >= p))])
ref_ltypes <- c("solid", "22", "11")

phi_fmt <- function(x, digits = 2) {
  as.numeric(sub("\\.?0+$", "", formatC(x, format = "f", digits = digits)))
}

for (i in seq_along(phi1_vec)) {
  φ1 <- phi1_vec[i]
  φ2 <- phi2_vec[i]
  
  grd <- grd_ref
  u <- cbind(
    pRB(grd$y1, 0.50, φ1, rho),
    pRB(grd$y2, 0.50, φ2, rho)
  )
  grd$d <- dCopula(u, tc) *
    dRB(grd$y1, 0.50, φ1, rho) *
    dRB(grd$y2, 0.50, φ2, rho)
  
  z <- sort(grd$d, decreasing = TRUE)
  cum <- cumsum(z*cellA)
  lev <- sapply(prob_seq, function(p) z[min(which(cum >= p))])
  grd$band <- cut(grd$d,
                  breaks = c(0, sort(lev), Inf),
                  labels = band_lab, include.lowest = TRUE
  )
  
  title_hdr <- bquote(
    atop(
      "r-Beta HDR",
      "("*phi[1] == .(phi_fmt(φ1))*"," ~ phi[2] == .(phi_fmt(φ2))*")"
    )
  )
  
  p_core <- ggplot(grd, aes(y1, y2, fill = band)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(
      data = data.frame(expand.grid(y1 = g, y2 = g), z = as.vector(z_ref)),
      aes(y1, y2,
          z = z,
          linetype = after_stat(factor(..level..,
                                       levels = ref_breaks,
                                       labels = c("50%", "75%", "90%")
          ))
      ),
      breaks = ref_breaks,
      colour = "red",
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = fills, drop = FALSE,
      name = title_hdr,
      guide = guide_legend(
        title.hjust = 0.5,
        keywidth = unit(8, "pt"),
        keyheight = unit(8, "pt")
      )
    ) +
    scale_linetype_manual(
      values = ref_ltypes,
      name = expression("Beta HDR ("*phi[1] == 0*"," ~ phi[2] == 0*")")
    ) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    labs(x = expression(italic(y)[1]), y = expression(italic(y)[2])) +
    theme_minimal(base_size = 15.5) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0)
    )
  
  p_fill_only <- ggplot(grd, aes(y1, y2, fill = band)) +
    geom_raster() +
    scale_fill_manual(
      values = fills, drop = FALSE,
      name = title_hdr,
      guide = guide_legend(
        title.hjust = 0.5,
        keywidth = unit(12, "pt"),
        keyheight = unit(12, "pt"),
        label.theme = element_text(size = 15.5)
      )
    ) +
    theme_minimal(base_size = 15.5) +
    theme(
      legend.position = "right",
      plot.margin = margin(0, 0, 0, 0)
    )
  
  leg_blue <- cowplot::get_legend(p_fill_only)
  
  df_grey <- data.frame(
    x = rep(c(0, 1), 3),
    y = rep(1:3, each = 2),
    type = factor(rep(c("50%", "75%", "90%"), each = 2),
                  levels = c("50%", "75%", "90%")
    )
  )
  p_lty_only <- ggplot(df_grey, aes(x, y, linetype = type, group = type)) +
    geom_line(colour = "red", linewidth = 0.8) +
    scale_linetype_manual(
      values = setNames(ref_ltypes, c("50%", "75%", "90%")),
      name = expression("Beta HDR ("*phi[1] == 0*"," ~ phi[2] == 0*")"),
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        title.position = "top",
        override.aes = list(colour = "red4", size = 1),
        title.theme = element_text(size = 15.5),
        label.theme = element_text(size = 15.5)
      )
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 15.5),
      legend.title = element_text(size = 15.5),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  guide_grey <- gtable_filter(ggplotGrob(p_lty_only), "guide-box", trim = TRUE)
  leg_grey <- cowplot::ggdraw(guide_grey)
  
  row_top <- cowplot::plot_grid(p_core, leg_blue,
                                nrow = 1, rel_widths = c(1, 0.35), align = "v"
  )
  
  p_joint <- cowplot::plot_grid(row_top, leg_grey,
                                ncol = 1, rel_heights = c(0.5, 0.05)
  )
  
  outname_j <- sprintf(
    "Manuscript/Output/rec_beta_copula_phi1_%s_phi2_%s.pdf",
    phi_fmt(φ1), phi_fmt(φ2)
  )
  ggsave(outname_j, p_joint, width = 7, height = 6.5, device = cairo_pdf)
  
  dens1 <- dRB(g, 0.50, φ1, rho)
  dens2 <- dRB(g, 0.50, φ2, rho)
  
  p_marg <- ggplot() +
    geom_line(aes(g, dens1, colour = "f1"), linewidth = 0.9) +
    geom_line(aes(g, dens2, colour = "f2"), linewidth = 0.9) +
    scale_colour_manual(
      values = c(f1 = "#0072B2", f2 = "#D55E00"),
      labels = c(
        expression(italic(f)(italic(y)[1])),
        expression(italic(f)(italic(y)[2]))
      ),
      name = ""
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
    labs(x = expression(italic(y)), y = expression(italic(f)(italic(y)))) +
    theme_minimal(base_size = 15.5) +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.35),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.margin = margin(b = -6)
    )
  
  outname_m <- sprintf(
    "Manuscript/Output/rec_beta_marginals_phi1_%s_phi2_%s.pdf",
    phi_fmt(φ1), phi_fmt(φ2)
  )
  ggsave(outname_m, p_marg, width = 4, height = 3.2, device = cairo_pdf)
  
  message("Saved pair ", i, "/", length(phi1_vec), ": ", outname_j)
}