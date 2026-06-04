#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL, required = FALSE) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    if (required) stop(sprintf("Missing required argument: %s", flag))
    return(default)
  }
  args[[idx + 1]]
}

continuous_h2_03_tsv <- get_arg("--continuous-h2-03", required = TRUE)
binary_h2_03_tsv <- get_arg("--binary-h2-03", required = TRUE)
continuous_h2_05_tsv <- get_arg("--continuous-h2-05", required = TRUE)
binary_h2_05_tsv <- get_arg("--binary-h2-05", required = TRUE)
output_pdf <- get_arg("--output-pdf", required = TRUE)

read_qq <- function(path) {
  read.table(
    path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

qq_title <- function(letter, trait, h2) {
  sprintf("%s. %s trait, h² = %.1f", letter, trait, h2)
}

plot_qq <- function(qq, main, show_legend = FALSE, lim = NULL) {
  if (is.null(lim)) {
    ymax <- max(qq$gcta_upper, qq$fed_upper, qq$expected, na.rm = TRUE)
    lim <- c(0, ymax * 1.03)
  }

  plot(
    qq$expected,
    qq$gcta_mean,
    type = "n",
    xlim = lim,
    ylim = lim,
    xlab = expression(bold(Expected ~ -log[10](P))),
    ylab = expression(bold(Observed ~ -log[10](P))),
    main = main,
    bty = "l",
    las = 1,
    cex.axis = 0.95,
    cex.lab = 1.15,
    cex.main = 1.18,
    font.main = 2,
    font.lab = 2,
    col.axis = "#283235",
    col.lab = "#172224",
    col.main = "#172224"
  )

  grid(nx = NA, ny = NULL, col = "#e6e8e6", lwd = 0.8)
  abline(0, 1, col = "#b23a3a", lwd = 1.8, lty = 2)

  polygon(
    c(qq$expected, rev(qq$expected)),
    c(qq$gcta_lower, rev(qq$gcta_upper)),
    col = adjustcolor("#2b6ca3", alpha.f = 0.16),
    border = NA
  )
  polygon(
    c(qq$expected, rev(qq$expected)),
    c(qq$fed_lower, rev(qq$fed_upper)),
    col = adjustcolor("#3a9b63", alpha.f = 0.16),
    border = NA
  )

  lines(qq$expected, qq$gcta_mean, col = "#2b6ca3", lwd = 2.7)
  lines(qq$expected, qq$fed_mean, col = "#3a9b63", lwd = 2.7)
  box(bty = "l", lwd = 1.2, col = "#283235")

  if (show_legend) {
    legend(
      "topleft",
      legend = c("Centralized", "Federated"),
      col = c("#2b6ca3", "#3a9b63"),
      lwd = 2.7,
      bty = "n",
      cex = 0.95,
      text.col = "#172224"
    )
  }
}

cont_h2_03 <- read_qq(continuous_h2_03_tsv)
bin_h2_03 <- read_qq(binary_h2_03_tsv)
cont_h2_05 <- read_qq(continuous_h2_05_tsv)
bin_h2_05 <- read_qq(binary_h2_05_tsv)

all_ymax <- max(
  cont_h2_03$gcta_upper, cont_h2_03$fed_upper, cont_h2_03$expected,
  bin_h2_03$gcta_upper, bin_h2_03$fed_upper, bin_h2_03$expected,
  cont_h2_05$gcta_upper, cont_h2_05$fed_upper, cont_h2_05$expected,
  bin_h2_05$gcta_upper, bin_h2_05$fed_upper, bin_h2_05$expected,
  na.rm = TRUE
)
common_lim <- c(0, all_ymax * 1.03)

dir.create(dirname(output_pdf), showWarnings = FALSE, recursive = TRUE)

pdf(output_pdf, width = 10.5, height = 9.2, useDingbats = FALSE)
op <- par(no.readonly = TRUE)
on.exit({
  par(op)
  dev.off()
}, add = TRUE)

par(
  mfrow = c(2, 2),
  mar = c(4.7, 5.0, 3.1, 1.2),
  oma = c(0.5, 0.5, 2.5, 0.5),
  family = "sans",
  font.main = 2,
  font.lab = 2,
  xaxs = "i",
  yaxs = "i"
)

plot_qq(cont_h2_03, qq_title("A", "Continuous", 0.3), show_legend = TRUE, lim = common_lim)
plot_qq(bin_h2_03, qq_title("B", "Binary", 0.3), show_legend = FALSE, lim = common_lim)
plot_qq(cont_h2_05, qq_title("C", "Continuous", 0.5), show_legend = FALSE, lim = common_lim)
plot_qq(bin_h2_05, qq_title("D", "Binary", 0.5), show_legend = FALSE, lim = common_lim)

mtext(
  "Null-SNP calibration across heritability settings",
  outer = TRUE,
  side = 3,
  line = 0.7,
  font = 2,
  cex = 1.25,
  col = "#172224"
)
