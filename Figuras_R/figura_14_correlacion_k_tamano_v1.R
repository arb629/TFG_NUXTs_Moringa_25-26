# ============================================================
# Correlación K vs tamaño (NUPTs y NUMTs)
# Kendall + Spearman, leyenda en recuadro negro a la derecha
# Salida: PNG
# ============================================================

library(ggplot2)
library(gridExtra)
library(ggtext)

# ── Rutas de entrada ─────────────────────────────────────────
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

# ── Título con solo la especie en cursiva ─────────────────────
titulo_mixto <- function(etiqueta) {
  partes  <- strsplit(etiqueta, " ", fixed = TRUE)[[1]]
  tipo    <- partes[1]
  especie <- paste(partes[-1], collapse = " ")
  paste0(tipo, " <i>", especie, "</i>")
}

# ── Extraer longitud del nombre ───────────────────────────────
extraer_longitud <- function(nombre) {
  m <- regmatches(nombre, regexpr("_(\\d+)-(\\d+)_", nombre))
  if (length(m) == 0) return(NA)
  coords <- as.numeric(strsplit(gsub("_", "", m), "-")[[1]])
  abs(coords[2] - coords[1])
}

# ── Función principal ─────────────────────────────────────────
hacer_scatter <- function(ruta, nombre, etiqueta) {

  cat("\nProcesando:", etiqueta, "\n")

  dat <- read.table(ruta, header = FALSE, sep = "\t",
                    col.names = c("id", "K"))
  dat <- dat[dat$K > 0, ]
  dat$longitud <- sapply(dat$id, extraer_longitud)
  dat <- dat[!is.na(dat$longitud) & dat$longitud > 0, ]

  cat("n =", nrow(dat), "\n")

  # Tests de correlación
  test_k <- cor.test(dat$K, dat$longitud, method = "kendall")
  test_s <- cor.test(dat$K, dat$longitud, method = "spearman", exact = FALSE)

  tau <- round(test_k$estimate, 3)
  p_k <- test_k$p.value
  rho <- round(test_s$estimate, 3)
  p_s <- test_s$p.value

  cat("Kendall tau =", tau, "p =", formatC(p_k, format = "e", digits = 2), "\n")
  cat("Spearman rho =", rho, "p =", formatC(p_s, format = "e", digits = 2), "\n")

  col_punto <- colores[[nombre]]

  # Texto leyenda
  sig_k <- ifelse(p_k < 0.001, "***", ifelse(p_k < 0.01, "**", ifelse(p_k < 0.05, "*", "ns")))
  sig_s <- ifelse(p_s < 0.001, "***", ifelse(p_s < 0.01, "**", ifelse(p_s < 0.05, "*", "ns")))

  label_corr <- paste0(
    "Kendall τ = ", tau, " (", sig_k, ")\n",
    "p = ", formatC(p_k, format = "e", digits = 2), "\n",
    "\n",
    "Spearman ρ = ", rho, " (", sig_s, ")\n",
    "p = ", formatC(p_s, format = "e", digits = 2), "\n",
    "\nn = ", nrow(dat)
  )

  p <- ggplot(dat, aes(x = K, y = log10(longitud))) +
    geom_point(alpha = 0.25, size = 0.7, color = col_punto) +
    geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.7) +
    coord_cartesian(ylim = c(1, NA)) +
    annotate("label",
             x = max(dat$K),
             y = max(log10(dat$longitud)),
             label = label_corr,
             hjust = 1, vjust = 1,
             size = 2.5,
             label.size = 0.6,
             label.padding = unit(0.3, "lines"),
             fill = "white",
             color = "black") +
    labs(title = titulo_mixto(etiqueta),
         x = "K (sustituciones/sitio)",
         y = expression(log[10](Tamaño~(bp)))) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_markdown(size = 11))

  # Guardar individual
  ggsave(file.path(outdir, paste0("correlacion_k_tamano_", nombre, ".png")),
         p, width = 6, height = 4, dpi = 300)

  p
}

# ── Ejecutar ──────────────────────────────────────────────────
plots <- mapply(hacer_scatter,
                ruta     = rutas,
                nombre   = names(rutas),
                etiqueta = etiquetas,
                SIMPLIFY = FALSE)

# ── Panel 2x2 ─────────────────────────────────────────────────
panel <- arrangeGrob(
  plots$nupt_oleif,
  plots$numt_oleif,
  plots$nupt_steno,
  plots$numt_steno,
  ncol = 2
)

ggsave(file.path(outdir, "panel_correlacion_k_tamano.png"),
       panel, width = 12, height = 8, dpi = 300)

cat("\nArchivos guardados en:", outdir, "\n")
cat("=== COMPLETADO ===\n")
