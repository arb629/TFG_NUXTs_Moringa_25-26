# NUMTs oleifera
x <- read.table("/home/kubuntu/Andrea/TFG/M_oleifera/NUXTs/mito_oleif_nuc_fixed.out",
                header=FALSE, sep="\t")
x$start_norm <- pmin(x$V9, x$V10)
x$end_norm   <- pmax(x$V9, x$V10)
x_ord  <- x[order(x$V2, x$start_norm, x$end_norm, -x$V12), ]
x_filt <- x_ord[!duplicated(x_ord[c("V2", "start_norm", "end_norm")]), ]
cat("NUMTs oleifera antes:", nrow(x), "\n")
cat("NUMTs oleifera despues:", nrow(x_filt), "\n")
write.table(x_filt[, 1:12],
  "/home/kubuntu/Andrea/TFG/M_oleifera/NUXTs/mito_oleif_nuc_fixed_nodup.out",
  sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# NUMTs stenopetala
x2 <- read.table("/home/kubuntu/Andrea/TFG/M_stenopetala/NUXTs/mito_steno_nuc_fixed.out",
                 header=FALSE, sep="\t")
x2$start_norm <- pmin(x2$V9, x2$V10)
x2$end_norm   <- pmax(x2$V9, x2$V10)
x2_ord  <- x2[order(x2$V2, x2$start_norm, x2$end_norm, -x2$V12), ]
x2_filt <- x2_ord[!duplicated(x2_ord[c("V2", "start_norm", "end_norm")]), ]
cat("NUMTs stenopetala antes:", nrow(x2), "\n")
cat("NUMTs stenopetala despues:", nrow(x2_filt), "\n")
write.table(x2_filt[, 1:12],
  "/home/kubuntu/Andrea/TFG/M_stenopetala/NUXTs/mito_steno_nuc_fixed_nodup.out",
  sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
