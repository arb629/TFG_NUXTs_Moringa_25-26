library(ggplot2)
library(dplyr)
library(readxl)
library(cowplot)

# Cargar datos
s1 <- read_excel("/home/kubuntu/Andrea/TFG/zhang_2020_porcentajes/Table_S1.xlsx", skip = 1)
s2 <- read_excel("/home/kubuntu/Andrea/TFG/zhang_2020_porcentajes/Table_S2.xlsx", skip = 1)

nupt_props <- as.numeric(s1[[11]]) * 100
nupt_props <- nupt_props[!is.na(nupt_props)]

numt_props <- as.numeric(s2[[11]]) * 100
numt_props <- numt_props[!is.na(numt_props)]

df <- data.frame(
  value = c(nupt_props, numt_props),
  type  = factor(c(rep("NUPTs", length(nupt_props)),
                   rep("NUMTs", length(numt_props))),
                 levels = c("NUPTs", "NUMTs"))
)

nupt_mean_val   <- round(mean(nupt_props), 3)
nupt_median_val <- round(median(nupt_props), 3)
numt_mean_val   <- round(mean(numt_props), 3)
numt_median_val <- round(median(numt_props), 3)

col_nupt    <- "#4f7a4d"
col_numt    <- "#a13f2b"
col_nupt_bg <- "#c8dbc7"
col_numt_bg <- "#e8c4bc"

lbl_mean_nupt <- paste0("Media NUPTs: ",   gsub("\\.", ",", nupt_mean_val),   "%")
lbl_med_nupt  <- paste0("Mediana NUPTs: ", gsub("\\.", ",", nupt_median_val), "%")
lbl_mean_numt <- paste0("Media NUMTs: ",   gsub("\\.", ",", numt_mean_val),   "%")
lbl_med_numt  <- paste0("Mediana NUMTs: ", gsub("\\.", ",", numt_median_val), "%")

stats <- data.frame(
  type     = factor(c("NUPTs","NUPTs","NUMTs","NUMTs"),
                    levels = c("NUPTs","NUMTs")),
  stat_val = c(nupt_mean_val, nupt_median_val,
               numt_mean_val, numt_median_val),
  label    = c(lbl_mean_nupt, lbl_med_nupt,
               lbl_mean_numt, lbl_med_numt),
  ltype    = c("dashed","solid","dashed","solid"),
  color    = c(col_nupt, col_nupt, col_numt, col_numt)
)

color_vals <- setNames(stats$color, stats$label)
ltype_vals <- setNames(stats$ltype, stats$label)

mol <- data.frame(
  type  = factor("NUPTs", levels = c("NUPTs","NUMTs")),
  value = 3.29
)

set.seed(42)

# Plot principal sin leyenda
p_main <- ggplot() +
  geom_violin(data = df %>% filter(type == "NUPTs"),
              aes(x = type, y = value),
              fill = col_nupt_bg, color = "gray60",
              alpha = 0.5, width = 0.7, scale = "width", trim = TRUE) +
  geom_violin(data = df %>% filter(type == "NUMTs"),
              aes(x = type, y = value),
              fill = col_numt_bg, color = "gray60",
              alpha = 0.5, width = 0.7, scale = "width", trim = TRUE) +
  geom_jitter(data = df %>% filter(type == "NUPTs"),
              aes(x = type, y = value),
              color = col_nupt, alpha = 0.4, size = 1.5,
              width = 0.12, height = 0) +
  geom_jitter(data = df %>% filter(type == "NUMTs"),
              aes(x = type, y = value),
              color = col_numt, alpha = 0.4, size = 1.5,
              width = 0.12, height = 0) +
  geom_errorbar(data = stats,
                aes(x = type, y = stat_val,
                    ymin = stat_val, ymax = stat_val,
                    color = label, linetype = label),
                width = 0.65, linewidth = 1.3) +
  geom_point(data = mol,
             aes(x = type, y = value),
             color = col_nupt, size = 5, shape = 16) +
  scale_color_manual(name = NULL, values = color_vals) +
  scale_linetype_manual(name = NULL, values = ltype_vals) +
  scale_y_continuous(
    limits = c(0, 3.5),
    breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5),
    labels = c("0", "0,5", "1,0", "1,5", "2,0", "2,5", "3,0", "3,5")
  ) +
  scale_x_discrete(limits = c("NUPTs", "NUMTs")) +
  labs(x = NULL, y = "Fracción del genoma nuclear (%)") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(size = 12, face = "bold"),
    axis.text.y  = element_text(size = 10),
    axis.title.y = element_text(size = 11, margin = margin(r = 10)),
    legend.position = "none"
  )

# Plot fantasma solo para extraer leyenda de lineas
p_legend_lines <- ggplot(stats,
                         aes(x = type, y = stat_val,
                             color = label, linetype = label)) +
  geom_line() +
  scale_color_manual(name = NULL, values = color_vals) +
  scale_linetype_manual(name = NULL, values = ltype_vals) +
  theme_classic() +
  theme(legend.text = element_text(size = 9),
        legend.key.width = unit(1.5, "cm"))

# Plot fantasma para extraer leyenda de M. oleifera
p_legend_mol <- ggplot(mol, aes(x = type, y = value)) +
  geom_point(aes(shape = "italic(M.~oleifera)~'NUPTs: 3,29%'"),
             color = col_nupt, size = 5) +
  scale_shape_manual(
    name = NULL,
    values = c("italic(M.~oleifera)~'NUPTs: 3,29%'" = 16),
    labels = c("italic(M.~oleifera)~'NUPTs: 3,29%'" =
                 expression(italic("M. oleifera")~"NUPTs: 3,29%"))
  ) +
  theme_classic() +
  theme(legend.text = element_text(size = 9),
        legend.key.width = unit(1.5, "cm"))

# Extraer leyendas
leg_lines <- get_legend(p_legend_lines)
leg_mol   <- get_legend(p_legend_mol)

# Combinar leyendas verticalmente
leg_combined <- plot_grid(leg_lines, leg_mol,
                          ncol = 1, align = "v",
                          rel_heights = c(4, 1))

# Combinar plot + leyenda
final <- plot_grid(p_main, leg_combined,
                   ncol = 2, rel_widths = c(3, 1.2))

ggsave("figura_zhang2020.png", plot = final,
       width = 10, height = 7, dpi = 300, bg = "white")

print("Figura guardada correctamente")
