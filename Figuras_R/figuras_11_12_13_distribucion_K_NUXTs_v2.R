library(mclust)

out_dir <- "/home/kubuntu/Andrea/TFG/trans_log/mclust/NUXTs/v4/Definitivo/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

archivos <- list(
  nupt_steno = "/home/kubuntu/Andrea/TFG/de_nuevo/BASEML_results/nupt_steno/k_values.txt",
  numt_steno = "/home/kubuntu/Andrea/TFG/de_nuevo/BASEML_results/numt_steno/k_values.txt",
  nupt_oleif = "/home/kubuntu/Andrea/TFG/de_nuevo/BASEML_results/nupt_oleif/k_values.txt",
  numt_oleif = "/home/kubuntu/Andrea/TFG/de_nuevo/BASEML_results/numt_oleif/k_values.txt"
)

colores <- list(
  nupt_steno = "#7baf48",
  numt_steno = "#df3f2d",
  nupt_oleif = "#4f7a4d",
  numt_oleif = "#a13f2b"
)

etiquetas <- list(
  nupt_steno = expression("NUPTs " * italic("M. stenopetala")),
  numt_steno = expression("NUMTs " * italic("M. stenopetala")),
  nupt_oleif = expression("NUPTs " * italic("M. oleifera")),
  numt_oleif = expression("NUMTs " * italic("M. oleifera"))
)

leer_k <- function(ruta) {
  dat <- read.table(ruta, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  k   <- suppressWarnings(as.numeric(dat[, 2]))
  k[!is.na(k) & k > 0]
}

todos_k     <- lapply(archivos, leer_k)
todos_log_k <- lapply(todos_k, log10)

# Límites globales
xlim_lin    <- c(0, 0.45)
xlim_log    <- range(unlist(todos_log_k)) + c(-0.05, 0.05) * diff(range(unlist(todos_log_k)))
x_ticks_log <- seq(floor(xlim_log[1] * 2) / 2, ceiling(xlim_log[2] * 2) / 2, by = 0.5)
ylim_lin    <- c(0, 14)
ylim_log    <- c(0, 2.5)

# ============================================================
# FIGURAS INDIVIDUALES
# ============================================================
for (muestra in names(archivos)) {
  k_raw <- todos_k[[muestra]]
  log_k <- todos_log_k[[muestra]]
  n     <- length(k_raw)
  col   <- colores[[muestra]]
  etiq  <- etiquetas[[muestra]]

  # Lineal individual
  png(paste0(out_dir, "lineal_", muestra, ".png"),
      width = 1800, height = 1200, res = 200)
  par(mar = c(5, 5, 4, 2), cex.main = 1.3, cex.lab = 1.1, cex.axis = 1.0)
  hist(k_raw, breaks = 80, freq = FALSE,
       col = col, border = "black",
       xlim = xlim_lin, ylim = ylim_lin,
       main = etiq,
       xlab = "K (distancia evolutiva)", ylab = "Frecuencia",
       xaxt = "n")
  axis(1, at = seq(0, 0.45, by = 0.05))
  mtext(paste0("n = ", format(n, big.mark = ".")), side = 3, adj = 1, cex = 0.85)
  dev.off()

  # Log10 individual
  png(paste0(out_dir, "log10_", muestra, ".png"),
      width = 1800, height = 1200, res = 200)
  par(mar = c(5, 2, 4, 6), cex.main = 1.3, cex.lab = 1.1, cex.axis = 1.0)
  hist(log_k, breaks = 60, freq = FALSE,
       col = col, border = "black",
       xlim = xlim_log, ylim = ylim_log,
       main = etiq,
       xlab = expression(log[10](K)), ylab = "",
       xaxt = "n", yaxt = "n")
  axis(1, at = x_ticks_log)
  axis(4)
  mtext("Frecuencia", side = 4, line = 3.5, cex = 1.1)
  mtext(paste0("n = ", format(n, big.mark = ".")), side = 3, adj = 0, cex = 0.85)
  dev.off()
}

# ============================================================
# PANEL LINEAL 2x2
# ============================================================
png(paste0(out_dir, "panel_lineal_K.png"), width = 3600, height = 2400, res = 300)
par(mfrow = c(2, 2), mar = c(5, 5, 5, 2), cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1)
for (i in seq_along(archivos)) {
  muestra <- names(archivos)[i]
  k_raw   <- todos_k[[muestra]]
  n       <- length(k_raw)
  hist(k_raw, breaks = 80, freq = FALSE,
       col = colores[[muestra]], border = "black",
       xlim = xlim_lin, ylim = ylim_lin,
       main = etiquetas[[muestra]],
       xlab = "K (distancia evolutiva)", ylab = "Frecuencia",
       xaxt = "n")
  axis(1, at = seq(0, 0.45, by = 0.05))
  mtext(paste0("n = ", format(n, big.mark = ".")), side = 3, adj = 1, cex = 1.1)
  mtext(LETTERS[i], side = 3, adj = 0, cex = 1.8, font = 2, line = 1.5)
}
dev.off()

# ============================================================
# PANEL LOG10 2x2
# ============================================================
png(paste0(out_dir, "panel_log10_K.png"), width = 3600, height = 2400, res = 300)
par(mfrow = c(2, 2), mar = c(5, 2, 5, 6), cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1)
for (i in seq_along(archivos)) {
  muestra <- names(archivos)[i]
  log_k   <- todos_log_k[[muestra]]
  n       <- length(log_k)
  hist(log_k, breaks = 60, freq = FALSE,
       col = colores[[muestra]], border = "black",
       xlim = xlim_log, ylim = ylim_log,
       main = etiquetas[[muestra]],
       xlab = expression(log[10](K)), ylab = "",
       xaxt = "n", yaxt = "n")
  axis(1, at = x_ticks_log)
  axis(4)
  mtext("Frecuencia", side = 4, line = 3.5, cex = 1.2)
  mtext(paste0("n = ", format(n, big.mark = ".")), side = 3, adj = 0, cex = 1.1)
  mtext(LETTERS[i], side = 3, adj = 0, cex = 1.8, font = 2, line = 1.5)
}
dev.off()

cat("Figuras guardadas en:", out_dir, "\n")
