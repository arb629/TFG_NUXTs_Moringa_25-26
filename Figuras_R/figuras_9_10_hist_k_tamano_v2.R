# ============================================================
# Correlación K vs tamaño + Histogramas de tamaño (NUPTs y NUMTs)
# Salida: PNG alta calidad
# ============================================================

library(ggplot2)
library(gridExtra)
library(ggtext)

rutas <- list(
  nupt_oleif = "/home/kubuntu/Andrea/TFG/de_nuevo/BASEML_results/nupt_oleif/k_values.txt",
  numt_oleif = "/home/kubuntu/Andrea/TFG/de_nuevo/BASEML_results/numt_oleif/k_values.txt",
  nupt_steno = "/home/kubuntu/Andrea/TFG/de_nuevo/BASEML_results/nupt_steno/k_values.txt",
  numt_steno = "/home/kubuntu/Andrea/TFG/de_nuevo/BASEML_results/numt_steno/k_values.txt"
)

etiquetas <- list(
  nupt_oleif = "NUPTs M. oleifera",
  numt_oleif = "NUMTs M. oleifera",
  nupt_steno = "NUPTs M. stenopetala",
  numt_steno = "NUMTs M. stenopetala"
)

outdir <- "/home/kubuntu/Andrea/TFG/tamano_correlacion"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Paleta de colores por dataset
colores <- list(
  nupt_oleif = "#4f7a4d",
  numt_oleif = "#a13f2b",
  nupt_steno = "#7baf48",
  numt_steno = "#df3f2d"
)

# ── Extraer longitud del nombre ───────────────────────────────
extraer_longitud <- function(nombre) {
  m <- regmatches(nombre, regexpr("_(\\d+)-(\\d+)_", nombre))
  if (length(m) == 0) return(NA)
  coords <- as.numeric(strsplit(gsub("_", "", m), "-")[[1]])
  abs(coords[2] - coords[1])
}

# ── Leer datos ────────────────────────────────────────────────
leer_datos <- function(ruta) {
  dat <- read.table(ruta, header = FALSE, sep = "\t",
                    col.names = c("id", "K"))
  dat <- dat[dat$K > 0, ]
  dat$longitud <- sapply(dat$id, extraer_longitud)
  dat[!is.na(dat$longitud) & dat$longitud > 0, ]
}

# ── Título con solo la especie en cursiva ─────────────────────
titulo_mixto <- function(etiqueta) {
  # Separa "NUPTs" o "NUMTs" de la especie
  partes <- strsplit(etiqueta, " ", fixed = TRUE)[[1]]
  tipo   <- partes[1]  # NUPTs o NUMTs
  especie <- paste(partes[-1], collapse = " ")  # M. oleifera / M. stenopetala
  bquote(.(tipo) ~ italic(.(especie)))
}

# ── 2. Histograma tamaño lineal ───────────────────────────────
hacer_hist_lineal <- function(dat, nombre, etiqueta, color) {

  media      <- round(mean(dat$longitud), 0)
  mediana    <- round(median(dat$longitud), 0)
  media_bp   <- media
  mediana_bp <- mediana

  label_stats <- paste0(
    "n = ", nrow(dat), "<br>",
    "<span style='color:blue'>&#9632;</span> Media = ", media_bp, " bp<br>",
    "<span style='color:black'>&#9632;</span> Mediana = ", mediana_bp, " bp"
  )

  p <- ggplot(dat, aes(x = longitud)) +
    geom_histogram(bins = 80, fill = color, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = media,   color = "blue",    linetype = "dashed", linewidth = 0.7) +
    geom_vline(xintercept = mediana, color = "black", linetype = "dashed", linewidth = 0.7) +
    annotate("richtext",
             x = Inf, y = Inf,
             label = label_stats,
             hjust = 1.05, vjust = 1.05,
             size = 2.4, label.size = 0.5,
             label.padding = unit(0.3, "lines"),
             fill = "white", color = "black") +
    scale_x_continuous(limits = c(0, quantile(dat$longitud, 0.99) * 1.05),
                       expand = c(0.02, 0)) +
    labs(title = titulo_mixto(etiqueta),
         x = "Tamaño (bp)",
         y = "Frecuencia") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(size = 11))

  ggsave(file.path(outdir, paste0("hist_tamano_lineal_", nombre, ".png")),
         p, width = 6, height = 4, dpi = 600)
  p
}

# ── 3. Histograma tamaño log10 ────────────────────────────────
hacer_hist_log10 <- function(dat, nombre, etiqueta, color) {

  media_log   <- mean(log10(dat$longitud))
  mediana_log <- median(log10(dat$longitud))
  media_bp    <- round(10^media_log, 0)
  mediana_bp  <- round(10^mediana_log, 0)

  label_stats <- paste0(
    "n = ", nrow(dat), "<br>",
    "<span style='color:blue'>&#9632;</span> Media = ", media_bp, " bp<br>",
    "<span style='color:black'>&#9632;</span> Mediana = ", mediana_bp, " bp"
  )

  # Calcular breaks comunes para el eje x (excepto nupt_oleif que tiene rango distinto)
  x_breaks <- seq(-2.5, 0, by = 0.5)

  p <- ggplot(dat, aes(x = log10(longitud))) +
    geom_histogram(bins = 50, fill = color, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = media_log,   color = "blue",    linetype = "dashed", linewidth = 0.7) +
    geom_vline(xintercept = mediana_log, color = "black", linetype = "dashed", linewidth = 0.7) +
    annotate("richtext",
             x = Inf, y = Inf,
             label = label_stats,
             hjust = 1.05, vjust = 1.05,
             size = 2.4, label.size = 0.5,
             label.padding = unit(0.3, "lines"),
             fill = "white", color = "black") +
    scale_x_continuous(breaks = seq(-4, 5, by = 0.5),
                   expand = c(0.02, 0)) +
    coord_cartesian(xlim = c(min(log10(dat$longitud)) - 0.1,
                          max(log10(dat$longitud)) + 0.1)) +
    labs(title = titulo_mixto(etiqueta),
         x = expression(log[10](Tamaño~(bp))),
         y = "Frecuencia") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(size = 11))

  ggsave(file.path(outdir, paste0("hist_tamano_log10_", nombre, ".png")),
         p, width = 6, height = 4, dpi = 600)
  p
}

# ── Ejecutar ──────────────────────────────────────────────────
datos <- lapply(rutas, leer_datos)

plots_hist_lin <- mapply(hacer_hist_lineal, dat=datos, nombre=names(rutas),
                         etiqueta=etiquetas, color=colores, SIMPLIFY=FALSE)
plots_hist_log <- mapply(hacer_hist_log10,  dat=datos, nombre=names(rutas),
                         etiqueta=etiquetas, color=colores, SIMPLIFY=FALSE)

# ── Paneles 2x2 ───────────────────────────────────────────────
guardar_panel <- function(plots, archivo) {
  panel <- arrangeGrob(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol=2)
  ggsave(file.path(outdir, archivo), panel, width=12, height=8, dpi=600)
}

guardar_panel(plots_hist_lin, "panel_hist_tamano_lineal.png")
guardar_panel(plots_hist_log, "panel_hist_tamano_log10.png")

cat("\nArchivos guardados en:", outdir, "\n")
cat("=== COMPLETADO ===\n")
