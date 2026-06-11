
#BiocManager::install("DESeq2")
#BiocManager::install("edgeR")
#BiocManager::install("apeglm")

library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(apeglm)
library(ashr)

setwd("/mnt/data/transcriptomica/iperez/proyectofinal/")  
exp_count <- read.table("exp/featureCounts/gene_counts.txt", header = TRUE, row.names = 1, sep = "\t")
rownames(exp_count) <- gsub("^gene-", "", rownames(exp_count))
sample_cols <- grep("SRR.*\\.bam$", colnames(exp_count), value = TRUE)
exp_count <- exp_count[, sample_cols]  # sobrescribir exp_count (¡ahora sí!)
colnames(exp_count) <- gsub("^.*(SRR[0-9]+).*$", "\\1", colnames(exp_count))

feno <- read.csv("src/deSEQ/metadata.csv", stringsAsFactors = FALSE)

all(colnames(exp_count) %in% feno$sample_id) 

feno <- feno[match(colnames(exp_count), feno$sample_id), ]
rownames(feno) <- feno$sample_id
feno$sample_id <- NULL

feno$treatment <- factor(paste(feno$condition, feno$time, sep = "_"))

dds <- DESeqDataSetFromMatrix(countData = exp_count,
                              colData = feno,
                              design = ~ treatment)
levels(dds$treatment)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Identificar genes de fago y de E. coli
phage_genes <- grep("^T4", rownames(dds), value = TRUE)
host_genes <- grep("^b[0-9]{4}", rownames(dds), value = TRUE)  # los E. coli siguen siendo "b0001", etc.
dds <- DESeq(dds)

dds_phage <- dds[phage_genes, ]
dds_host <- dds[host_genes, ]

contrastes <- list(
  c("treatment", "infected_20min", "infected_1min"),
  c("treatment", "infected_20min", "infected_4min"),
  c("treatment", "infected_20min", "infected_7min"),
  c("treatment", "infected_7min", "infected_1min"),
  c("treatment", "infected_7min", "infected_4min"),
  c("treatment", "infected_4min", "infected_1min")
)
nombres <- c("20min_vs_1min", "20min_vs_4min", "20min_vs_7min",
             "7min_vs_1min", "7min_vs_4min", "4min_vs_1min")


#Resultados para todos los genes
for (i in seq_along(contrastes)) {
  res_temp <- results(dds, contrast = contrastes[[i]])
  assign(paste0("res_", nombres[i]), res_temp)
  write.csv(as.data.frame(res_temp), paste0("resultados_todos_", nombres[i], ".csv"), row.names = TRUE)
}

#Filtros
for (i in seq_along(nombres)) {
  res_obj <- get(paste0("res_", nombres[i]))
  
  #Filtrar fago
  res_phage <- res_obj[phage_genes, ]
  write.csv(as.data.frame(res_phage), paste0("resultados_fago_", nombres[i], ".csv"), row.names = TRUE)
  
  #Filtrar E. coli
  res_host <- res_obj[host_genes, ]
  write.csv(as.data.frame(res_host), paste0("resultados_ecoli_", nombres[i], ".csv"), row.names = TRUE)
}

#conteos normalizados
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(normalized_counts), "conteos_normalizados_todos.csv", row.names = TRUE)# normalization and save
write.csv(as.data.frame(normalized_counts[phage_genes, ]), "conteos_normalizados_fago.csv", row.names = TRUE)
write.csv(as.data.frame(normalized_counts[host_genes, ]), "conteos_normalizados_ecoli.csv", row.names = TRUE)

######### PCA ################
vsd <- vst(dds, blind = FALSE)

# PCA usando treatment (combina condition + time)
pcaData <- plotPCA(vsd, intgroup = "treatment", ntop = 1000, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = treatment)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% var")) +
  ylab(paste0("PC2: ", percentVar[2], "% var")) +
  ggtitle("PCA por condición y tiempo") +
  theme_bw(base_size = 14) +
  theme(legend.title = element_blank())

ggsave("PCA.pdf", pca_plot, width = 7, height = 5)

################# Volcano Plot - FAGO ############
library(ggplot2)
library(ggrepel)
library(dplyr)

#Lista de archivos de fago
archivos_fago <- c("resultados_fago_20min_vs_1min.csv",
                   "resultados_fago_20min_vs_4min.csv",
                   "resultados_fago_20min_vs_7min.csv",
                   "resultados_fago_7min_vs_1min.csv",
                   "resultados_fago_7min_vs_4min.csv",
                   "resultados_fago_4min_vs_1min.csv")

#Función volcano_plot
volcano_plot <- function(archivo, titulo) {
  res <- read.csv(archivo, row.names = 1)
  res_df <- data.frame(
    gene = rownames(res),
    log2FoldChange = res$log2FoldChange,
    padj = res$padj
  )
  res_df$sig <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 2,
                       ifelse(res_df$log2FoldChange > 2, "Up", "Down"), "NS")
  
  niveles <- unique(res_df$sig)
  colores <- c()
  if ("Up" %in% niveles) colores <- c(colores, "Up" = "red")
  if ("Down" %in% niveles) colores <- c(colores, "Down" = "blue")
  if ("NS" %in% niveles) colores <- c(colores, "NS" = "grey50")
  
  top_genes <- res_df[res_df$sig != "NS", ]
  if(nrow(top_genes) > 0) {
    top_genes <- top_genes[order(top_genes$padj), ][1:min(10, nrow(top_genes)), ]
  } else {
    top_genes <- data.frame()
  }
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = colores) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(title = titulo, x = "Log2 Fold Change", y = "-Log10 padj") +
    theme_bw() + theme(legend.title = element_blank())
  
  if(nrow(top_genes) > 0) {
    p <- p + geom_text_repel(data = top_genes, aes(label = gene),
                             size = 3, max.overlaps = 15, seed = 123)
  }
  return(p)
}

#Generar y guardar los volcano plots de fago
for (arch in archivos_fago) {
  titulo <- gsub("resultados_fago_|\\.csv", "", arch)
  titulo <- gsub("_vs_", " vs ", titulo)
  titulo <- paste("Fago:", titulo)
  p <- volcano_plot(arch, titulo)
  nombre_pdf <- paste0("volcano_fago_", gsub("resultados_fago_|\\.csv", "", arch), ".pdf")
  ggsave(nombre_pdf, p, width = 8, height = 6)
}

################# Volcano Plot - E. coli ############
#Lista de archivos de E. coli
archivos_ecoli <- c("resultados_ecoli_20min_vs_1min.csv",
                    "resultados_ecoli_20min_vs_4min.csv",
                    "resultados_ecoli_20min_vs_7min.csv",
                    "resultados_ecoli_7min_vs_1min.csv",
                    "resultados_ecoli_7min_vs_4min.csv",
                    "resultados_ecoli_4min_vs_1min.csv")

volcano_plot <- function(archivo, titulo) {
  res <- read.csv(archivo, row.names = 1)
  res_df <- data.frame(
    gene = rownames(res),
    log2FoldChange = res$log2FoldChange,
    padj = res$padj
  )
  res_df$sig <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 2,
                       ifelse(res_df$log2FoldChange > 2, "Up", "Down"), "NS")
  
  niveles <- unique(res_df$sig)
  colores <- c()
  if ("Up" %in% niveles) colores <- c(colores, "Up" = "red")
  if ("Down" %in% niveles) colores <- c(colores, "Down" = "blue")
  if ("NS" %in% niveles) colores <- c(colores, "NS" = "grey50")
  
  top_genes <- res_df[res_df$sig != "NS", ]
  if(nrow(top_genes) > 0) {
    top_genes <- top_genes[order(top_genes$padj), ][1:min(10, nrow(top_genes)), ]
  } else {
    top_genes <- data.frame()
  }
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = colores) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(title = titulo, x = "Log2 Fold Change", y = "-Log10 padj") +
    theme_bw() + theme(legend.title = element_blank())
  
  if(nrow(top_genes) > 0) {
    p <- p + geom_text_repel(data = top_genes, aes(label = gene),
                             size = 3, max.overlaps = 15, seed = 123)
  }
  return(p)
}

#Generar y guardar los volcano plots de E. coli
for (arch in archivos_ecoli) {
  titulo <- gsub("resultados_ecoli_|\\.csv", "", arch)
  titulo <- gsub("_vs_", " vs ", titulo)
  titulo <- paste("E. coli:", titulo)
  p <- volcano_plot(arch, titulo)
  nombre_pdf <- paste0("volcano_ecoli_", gsub("resultados_ecoli_|\\.csv", "", arch), ".pdf")
  ggsave(nombre_pdf, p, width = 8, height = 6)
}

############Heatmap of DEGs ##########
library(pheatmap)
library(dplyr)
library(tidyr)

#Cargar matrices normalizadas
norm_fago <- read.csv("conteos_normalizados_fago.csv", row.names = 1)
norm_ecoli <- read.csv("conteos_normalizados_ecoli.csv", row.names = 1)

#Ordenar el metadata por tiempo
feno <- read.csv("metadata.csv", row.names = 1)
orden_tiempo <- order(match(feno$time, c("0min","1min","4min","7min","20min")))
feno_ordenado <- feno[orden_tiempo, ]

#Anotación de columnas para heatmaps
annotation_col <- data.frame(
  Tiempo = factor(feno_ordenado$time, levels = c("0min","1min","4min","7min","20min")),
  Condicion = feno_ordenado$condition
)
rownames(annotation_col) <- rownames(feno_ordenado)

generar_heatmaps_top50 <- function(organismo, contrastes, norm_matrix, annotation_col) {
  for (cont in contrastes) {
    archivo <- paste0("resultados_", organismo, "_", cont, ".csv")
    if (!file.exists(archivo)) {
      warning("No se encuentra: ", archivo)
      next
    }
    
    res <- read.csv(archivo, row.names = 1)
    sig <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 2, ]
    if (nrow(sig) == 0) {
      cat("No hay DEGs en", organismo, "-", cont, "\n")
      next
    }
    sig <- sig[order(sig$padj), ]
    top_genes <- rownames(sig)[1:min(50, nrow(sig))]
    
    #Subset de la matriz normalizada y ordenar columnas
    mat <- norm_matrix[top_genes, , drop = FALSE]
    mat <- mat[, rownames(annotation_col)]  
    
    #Escalar por filas (z-score)
    mat_scaled <- t(scale(t(mat)))
    
    titulo <- paste(toupper(substr(organismo,1,1)), substr(organismo,2,nchar(organismo)), 
                    " - top 50 DEGs (", gsub("_vs_", " vs ", cont), ")", sep="")
    
    #Guardar PDF
    pdf(paste0("heatmap_", organismo, "_top50_", cont, ".pdf"), width = 10, height = 8)
    pheatmap(mat_scaled,
             scale = "none",
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             show_rownames = TRUE,
             annotation_col = annotation_col,
             main = titulo,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
             border_color = NA)
    dev.off()
    cat("Guardado: heatmap_", organismo, "_top50_", cont, ".pdf\n", sep="")
  }
}

#Lista de contrastes (nombres de archivos)
contrastes <- c("20min_vs_1min", "20min_vs_4min", "20min_vs_7min",
                "7min_vs_1min", "7min_vs_4min", "4min_vs_1min")

#Heatmaps para fago
generar_heatmaps_top50("fago", contrastes, norm_fago, annotation_col)

#Heatmaps para E. coli
generar_heatmaps_top50("ecoli", contrastes, norm_ecoli, annotation_col)

#Extraer genes
organismos <- c("fago", "ecoli")
contrastes <- c("20min_vs_1min", "20min_vs_4min", "20min_vs_7min",
                "7min_vs_1min", "7min_vs_4min", "4min_vs_1min")

for (org in organismos) {
  for (cont in contrastes) {
    # Construir el nombre del archivo de resultados
    archivo <- paste0("resultados_", org, "_", cont, ".csv")
    
    #Verificar que el archivo existe
    if (!file.exists(archivo)) {
      warning(paste("No se encuentra:", archivo))
      next
    }
    
    #Leer resultados
    res <- read.csv(archivo, row.names = 1)
    
    #Extraer genes up (log2FC > 2 y padj < 0.05)
    genes_up <- rownames(res[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 2, ])
    #Extraer genes down (log2FC < -2 y padj < 0.05)
    genes_down <- rownames(res[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -2, ])
    #Genes totales (up + down)
    genes_all <- unique(c(genes_up, genes_down))
    
    #Guardar archivos de texto
    write.table(genes_up,   paste0("genes_up_", org, "_", cont, ".txt"),   row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(genes_down, paste0("genes_down_", org, "_", cont, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(genes_all,  paste0("genes_all_", org, "_", cont, ".txt"),  row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    cat(org, cont, "- Up:", length(genes_up), "Down:", length(genes_down), "All:", length(genes_all), "\n")
  }
}


# Lista de contrastes y organismos (según tus archivos)
organismos <- c("fago", "ecoli")
contrastes <- c("20min_vs_1min", "20min_vs_4min", "20min_vs_7min",
                "7min_vs_1min", "7min_vs_4min", "4min_vs_1min")

for (org in organismos) {
  for (cont in contrastes) {
    # Construir nombre del archivo de resultados
    archivo <- paste0("resultados_", org, "_", cont, ".csv")
    if (!file.exists(archivo)) {
      cat("No existe:", archivo, "\n")
      next
    }
    
    # Leer resultados
    res <- read.csv(archivo, row.names = 1)
    
    # Asegurar que padj y log2FoldChange no sean NA
    res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
    
    # 1. Genes upregulados (log2FC > 2 y padj < 0.05) ordenados por padj ascendente
    up <- res[res$log2FoldChange > 2 & res$padj < 0.05, ]
    up <- up[order(up$padj), ]
    top20_up <- head(rownames(up), 20)
    
    # 2. Genes downregulados (log2FC < -2 y padj < 0.05) ordenados por padj ascendente
    down <- res[res$log2FoldChange < -2 & res$padj < 0.05, ]
    down <- down[order(down$padj), ]
    top20_down <- head(rownames(down), 20)
    
    # Guardar listas (opcional)
    write.table(top20_up,   paste0("top20_up_", org, "_", cont, ".txt"),   row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(top20_down, paste0("top20_down_", org, "_", cont, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # Mostrar resumen
    cat(org, cont, "- Up:", length(top20_up), "genes top20; Down:", length(top20_down), "genes top20\n")
  }
}

# install.packages("uwot")
# install.packages("RColorBrewer")

#UMAP
#Matriz de expresión VST
vsd <- vst(dds, blind = FALSE)
expr_matrix <- t(assay(vsd))   # muestras filas, genes columnas

feno_umap <- feno[rownames(expr_matrix), ]

sample_names <- rownames(expr_matrix)
treatment_vec <- colData(dds)[sample_names, "treatment"]

#UMAP sobre PCA
pca <- prcomp(expr_matrix, scale. = TRUE)
set.seed(123)
umap_pca <- umap(pca$x, n_neighbors = 5, min_dist = 0.3)

df_pca <- data.frame(
  UMAP1 = umap_pca[,1],
  UMAP2 = umap_pca[,2],
  treatment = treatment_vec
)

#UMAP sobre datos crudos
gene_vars <- apply(expr_matrix, 2, var)
top_genes <- order(gene_vars, decreasing = TRUE)[1:500]
expr_top <- expr_matrix[, top_genes]

set.seed(123)
umap_raw <- umap(expr_top, n_neighbors = 5, min_dist = 0.3)

df_raw <- data.frame(
  UMAP1 = umap_raw[,1],
  UMAP2 = umap_raw[,2],
  treatment = treatment_vec
)

#Gráficos
p1 <- ggplot(df_pca, aes(x = UMAP1, y = UMAP2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = "UMAP sobre PCA")

p2 <- ggplot(df_raw, aes(x = UMAP1, y = UMAP2, color = treatment)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = "UMAP sobre datos crudos (top 500 genes)")

# Mostrar y guardar
print(p1)
print(p2)
ggsave("UMAP_PCA.pdf", p1, width = 7, height = 5)
ggsave("UMAP_raw.pdf", p2, width = 7, height = 5)


#Análisis GO enriquecido para E. coli

library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(org.EcK12.eg.db)

#Descargar features
if (!file.exists("ecoli_feature_table.txt.gz")) {
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_feature_table.txt.gz",
    destfile = "ecoli_feature_table.txt.gz"
  )
}

#Leer y procesar la tabla de conversión
ecoli_genes <- read.delim(gzfile("ecoli_feature_table.txt.gz"),
                          header = TRUE, sep = "\t", comment.char = "#")
colnames(ecoli_genes)[15] <- "symbol"
colnames(ecoli_genes)[16] <- "GeneID"
colnames(ecoli_genes)[17] <- "locus_tag"

conv_table <- ecoli_genes %>%
  filter(gene == "gene") %>%
  dplyr::select(locus_tag, symbol, GeneID) %>%
  filter(!is.na(locus_tag) & locus_tag != "") %>%
  distinct()

#Cargar matriz normalizada de E. coli
norm_ecoli <- read.csv("conteos_normalizados_ecoli.csv", row.names = 1)

universo <- data.frame(locus_tag = rownames(norm_ecoli)) %>%
  left_join(conv_table, by = "locus_tag") %>%
  filter(!is.na(GeneID))
cat("Universo E. coli:", nrow(universo), "genes\n")

# GO por contraste para E. coli - todas las categorías visibles
contrastes_ecoli <- c("20min_vs_1min", "20min_vs_4min", "20min_vs_7min",
                      "7min_vs_1min",  "7min_vs_4min",  "4min_vs_1min")

for (cont in contrastes_ecoli) {
  for (dir in c("up", "down")) {
    archivo <- paste0("genes_", dir, "_ecoli_", cont, ".txt")
    if (!file.exists(archivo)) next
    
    genes <- tryCatch(
      read.table(archivo, stringsAsFactors = FALSE)[,1],
      error = function(e) character(0)
    )
    if (length(genes) < 3) next
    
    genes_m <- data.frame(locus_tag = genes) %>%
      left_join(conv_table, by = "locus_tag") %>%
      filter(!is.na(GeneID))
    
    if (nrow(genes_m) < 3) next
    
    for (ont in c("BP", "MF", "CC")) {
      go_res <- enrichGO(gene          = as.character(genes_m$GeneID),
                         universe      = as.character(universo$GeneID),
                         OrgDb         = org.EcK12.eg.db,
                         keyType       = "ENTREZID",
                         ont           = ont,
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 1,   # traer todos para graficar
                         qvalueCutoff  = 1,
                         minGSSize     = 3,
                         readable      = TRUE)
      
      df <- as.data.frame(go_res)
      if (nrow(df) == 0) {
        cat("Sin términos:", cont, dir, ont, "\n")
        next
      }
      
      df$sig      <- ifelse(df$p.adjust < 0.05, "Significativo (padj<0.05)", "No significativo")
      df$log_padj <- -log10(df$p.adjust + 1e-10)
      
      # Top 15 por -log10(padj)
      df_plot <- df %>%
        arrange(desc(log_padj)) %>%
        head(15)
      
      p <- ggplot(df_plot, aes(x = reorder(Description, log_padj),
                               y = log_padj,
                               fill = sig)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = Count), hjust = -0.1, size = 3) +
        coord_flip() +
        scale_fill_manual(values = c("Significativo (padj<0.05)" = "#1D9E75",
                                     "No significativo"          = "#D3D1C7")) +
        geom_hline(yintercept = -log10(0.05),
                   linetype = "dashed", color = "red", linewidth = 0.5) +
        labs(title    = paste("GO", ont, "— E. coli —", dir, "—",
                              gsub("_vs_", " vs ", cont)),
             subtitle = paste("Total genes mapeados:", nrow(genes_m)),
             x = "", y = "-log10(padj)", fill = "") +
        theme_bw(base_size = 11) +
        theme(legend.position = "bottom") +
        expand_limits(y = max(df_plot$log_padj) * 1.2)
      
      print(p)
      ggsave(paste0("GO_ecoli_", ont, "_", dir, "_", cont, ".pdf"),
             p, width = 8, height = 6)
      cat("Guardado: GO_ecoli", ont, dir, cont, "\n")
    }
  }
}

#GOS para fago

# Descargar feature table del fago T4 (NC_000866.4)
download.file(
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/836/945/GCF_000836945.1_ViralProj14044/GCF_000836945.1_ViralProj14044_feature_table.txt.gz",
  destfile = "t4_feature_table.txt.gz",
  mode = "wb"
)

t4_raw <- read.delim(
  gzfile("t4_feature_table.txt.gz"),
  header = TRUE, sep = "\t", comment.char = "#",
  stringsAsFactors = FALSE
)

colnames(t4_raw)
head(t4_raw, 5)

colnames(t4_raw)[14] <- "description"
colnames(t4_raw)[15] <- "product"   # símbolo
colnames(t4_raw)[16] <- "symbol"
colnames(t4_raw)[17] <- "locus_tag"
colnames(t4_raw)[18] <- "GeneID"

conv_t4_full <- t4_raw[t4_raw$gene == "CDS", ] %>%
  dplyr::select(locus_tag, symbol, GeneID, description) %>%
  filter(!is.na(locus_tag) & locus_tag != "") %>%
  distinct()

head(conv_t4_full, 10)

conv_t4_full <- conv_t4_full %>%
  mutate(categoria = case_when(
    grepl("replicat|polymerase|topoisomerase|helicase|primase|ligase|exonuclease|endonuclease|recombin",
          description, ignore.case = TRUE) ~ "DNA_replication_repair",
    grepl("capsid|head|tail|fiber|baseplate|sheath|spike|connector|portal|vertex|whisker",
          description, ignore.case = TRUE) ~ "Virion_structure",
    grepl("lysis|holin|lysin|endolysin|spanin|lysis inhibitor",
          description, ignore.case = TRUE) ~ "Lysis",
    grepl("transcri|sigma|promoter|RNA pol|anti-sigma|MotB|FmdB|AsiA|MotA",
          description, ignore.case = TRUE) ~ "Transcription_regulation",
    grepl("tRNA|ribosom|translat",
          description, ignore.case = TRUE) ~ "Translation",
    grepl("nucleotide|dNTP|thymidyl|dihydrofolate|ribonucleo|dCTPase|dUTPase",
          description, ignore.case = TRUE) ~ "Nucleotide_metabolism",
    grepl("inject|membrane|host|anti-restrict|modif|glucosyl|hydroxymethyl",
          description, ignore.case = TRUE) ~ "Host_takeover",
    grepl("hypothetical|unknown",
          description, ignore.case = TRUE) ~ "Hypothetical",
    TRUE ~ "Other"
  ))

table(conv_t4_full$categoria)

# Universo = todos los genes del fago en tu experimento
universo_fago <- data.frame(locus_tag = rownames(norm_fago)) %>%
  left_join(conv_t4_full, by = "locus_tag") %>%
  filter(!is.na(categoria))

# Función de enriquecimiento por categoría (Fisher's exact test)
enriquecer_fago <- function(genes_query, universo_df, titulo) {
  
  query_df <- data.frame(locus_tag = genes_query) %>%
    left_join(conv_t4_full, by = "locus_tag") %>%
    filter(!is.na(categoria))
  
  if (nrow(query_df) < 3) {
    cat("Muy pocos genes para:", titulo, "\n")
    return(NULL)
  }
  
  categorias <- unique(universo_df$categoria)
  
  resultados <- lapply(categorias, function(cat) {
    # Tabla de contingencia
    a <- sum(query_df$categoria == cat)           # query en categoría
    b <- nrow(query_df) - a                        # query fuera
    c <- sum(universo_df$categoria == cat) - a     # universo en categoría (sin query)
    d <- nrow(universo_df) - nrow(query_df) - c    # universo fuera
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    ft <- fisher.test(mat, alternative = "greater")
    
    data.frame(
      categoria   = cat,
      genes_query = a,
      genes_total = sum(universo_df$categoria == cat),
      GeneRatio   = paste0(a, "/", nrow(query_df)),
      BgRatio     = paste0(sum(universo_df$categoria == cat), "/", nrow(universo_df)),
      pvalue      = ft$p.value,
      OR          = ft$estimate
    )
  })
  
  res_df <- do.call(rbind, resultados)
  res_df$padj <- p.adjust(res_df$pvalue, method = "BH")
  res_df <- res_df[order(res_df$pvalue), ]
  res_df$sig <- res_df$padj < 0.05
  return(res_df)
}

# Correr para todos los contrastes
contrastes_fago <- c("20min_vs_1min", "20min_vs_4min", "20min_vs_7min",
                     "7min_vs_1min",  "7min_vs_4min",  "4min_vs_1min")

resultados_todos <- list()

for (cont in contrastes_fago) {
  for (dir in c("up", "down")) {
    archivo <- paste0("genes_", dir, "_fago_", cont, ".txt")
    if (!file.exists(archivo)) next
    
    genes <- tryCatch(
      read.table(archivo, stringsAsFactors = FALSE)[,1],
      error = function(e) character(0)
    )
    if (length(genes) < 3) next
    
    res <- enriquecer_fago(genes, universo_fago, paste(cont, dir))
    if (is.null(res)) next
    
    res$contraste <- gsub("_vs_", " vs ", cont)
    res$direccion <- dir
    resultados_todos[[paste(cont, dir)]] <- res
  }
}

# Unir todo
df_todos <- do.call(rbind, resultados_todos)

# Graficar cada contraste por separado con todas las categorías
for (cont in contrastes_fago) {
  for (dir in c("up", "down")) {
    nombre <- paste(cont, dir)
    if (!nombre %in% names(resultados_todos)) next
    
    df_plot <- resultados_todos[[nombre]]
    df_plot$sig <- ifelse(df_plot$padj < 0.05, "Significativo (padj<0.05)", "No significativo")
    df_plot$log_padj <- -log10(df_plot$padj + 1e-10)  # evitar log(0)
    
    p <- ggplot(df_plot, aes(x = reorder(categoria, log_padj), 
                             y = log_padj, 
                             fill = sig)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = paste0(genes_query, "/", genes_total)),
                hjust = -0.1, size = 3) +
      coord_flip() +
      scale_fill_manual(values = c("Significativo (padj<0.05)" = "#1D9E75",
                                   "No significativo" = "#D3D1C7")) +
      geom_hline(yintercept = -log10(0.05), 
                 linetype = "dashed", color = "red", linewidth = 0.5) +
      labs(title = paste("Fago T4 —", dir, "—", gsub("_vs_", " vs ", cont)),
           subtitle = paste("Total genes:", 
                            tryCatch(nrow(read.table(paste0("genes_", dir, "_fago_", cont, ".txt"))),
                                     error = function(e) 0)),
           x = "",
           y = "-log10(padj)",
           fill = "") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom") +
      expand_limits(y = max(df_plot$log_padj) * 1.2)
    
    print(p)
    ggsave(paste0("GO_categorias_fago_", dir, "_", cont, ".pdf"), 
           p, width = 7, height = 5)
    cat("Guardado:", paste0("GO_categorias_fago_", dir, "_", cont, ".pdf"), "\n")
  }
}
