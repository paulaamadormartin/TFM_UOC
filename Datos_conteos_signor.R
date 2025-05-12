library(hipathia)
setwd("/Users/paula/Documents/Bioinformatica/tfm/TFM_UOC/signor2hipathia/results/v1.0_R4.3.1")
pathways_signor=readRDS("2024_25_may_sáb_18_57Vhi_3.0.2_mgi.RDS")
pathways_signor$pathigraphs$hsaSIGXOR0GBM$effector.subgraphs$`P-hsaSIGXOR0GBM-18` %>% plot

#sacar efectores
subG_signor=pathways_signor$pathigraphs$hsaSIGXOR0GBM$effector.subgraphs$`P-hsaSIGXOR0GBM-18`
V(subG_signor)[degree(subG_signor, mode = "out")==0]$genesList
igraph::as_data_frame(subG_signor, what = "vertices") %>% View

#nombre de las rutas
all_path_names=sapply(pathways_signor$pathigraphs, function(x) x$path.name)
all_path_names["hsaSIGXOR0AC"]

#función para sacar el gen efector
get_entrez_effector=function(subG){
  V(subG)[degree(subG, mode = "out")==0]$genesList %>% unlist()  
}

#ejemplo
get_entrez_effector(pathways_signor$pathigraphs$hsaSIGXOR0GBM$effector.subgraphs$`P-hsaSIGXOR0GBM-18`)

#sacamos todos los genes efectores que hay en signor
# get_all_effector_genes= function(pathways_signor) {
#   # Usamos sapply para recorrer las vías y las subvías efectoras
#   all_genes = sapply(pathways_signor$pathigraphs, function(pw) {
#     # Extraemos genes de cada subvía efectora dentro de cada vía
#     sapply(pw$effector.subgraphs, function(subG) {
#       # Llamamos a la función auxiliar para obtener los genes de cada subvía
#       genes= get_entrez_effector(subG)
#       return(genes)
#     })
#   })
#   
#   Convertimos la lista resultante a un vector plano y quitamos duplicados
#   all_genes=unique(unlist(all_genes))
#   
#   return(all_genes)
# }

#all_effector_genes=get_all_effector_genes(pathways_signor)
#length(all_effector_genes)  # Cuántos genes hay
#head(all_effector_genes)    # Primeros genes


# Función para extraer efectores y las rutas
get_effector_gene_pathway_map=function(pathways_signor) {
  result=list()
  
  for (pw_id in names(pathways_signor$pathigraphs)) {
    pw= pathways_signor$pathigraphs[[pw_id]]
    pw_name= pw$path.name
    
    for (subgraph_name in names(pw$effector.subgraphs)) {
      subG = pw$effector.subgraphs[[subgraph_name]]
      genes= get_entrez_effector(subG)
      
      if (length(genes) > 0) {
        df <- data.frame(
          gene = genes,
          pathway_id = pw_id,
          pathway_name = pw_name,
          stringsAsFactors = FALSE
        )
        result[[paste(pw_id, subgraph_name, sep = "_")]] <- df
      }
    }
  }
  
  # Unir todo en un único data frame
  final_df= do.call(rbind, result)
  return(final_df)
}

# Ejecutar
effector_gene_pathway_map= get_effector_gene_pathway_map(pathways_signor)

# Ejemplo de uso
head(effector_gene_pathway_map)


library(org.Hs.eg.db)


# Crear el mapeo Entrez -> Symbol
entrez_ids= unique(effector_gene_pathway_map$gene)
gene_symbols= mapIds(org.Hs.eg.db,
                       keys = entrez_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")
# Añadir columna con símbolo a tu dataframe
effector_gene_pathway_map$symbol= gene_symbols[as.character(effector_gene_pathway_map$gene)]


# Función para buscar información de genes efectores
search_multiple_genes <- function(effector_df, search_terms) {
  
  # Filtrar por símbolo de gen (ignora mayúsculas/minúsculas)
  gene_info <- effector_df[toupper(effector_df$symbol) %in% toupper(search_terms), ]
  
  # Si hay genes encontrados, devolver el dataframe
  if (nrow(gene_info) > 0) {
    return(gene_info)
  } else {
    return("Ninguno de los genes solicitados está en la lista de genes efectores.")
  }
}


# Lista de genes a buscar
genes_to_search= c("STAT3", "SMAD2", "MTOR", "MAPK8", "NFKBIA", "MAPK1", "JUN", "AKT1", "MAP2K7")

# Ejecutar la búsqueda en el data.frame completo
result_multiple= search_multiple_genes(effector_gene_pathway_map, genes_to_search)

# Mostrar los resultados
print(result_multiple)





#importar datos
library(readr)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")


#normalizar TMM
library(edgeR)
#partimos de los datos crudos
raw_counts=read.delim("/Users/paula/Documents/Bioinformatica/tfm/codigo/data/GSE97979_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1, check.names = FALSE)
str(raw_counts)


# Vector con nombres de muestra (en el mismo orden que columnas de tu matriz)
samples=c("GSM2584544", "GSM2584545", "GSM2584546", "GSM2584547",
             "GSM2584548", "GSM2584549", "GSM2584550", "GSM2584551",
             "GSM2584552", "GSM2584553", "GSM2584554", "GSM2584555",
             "GSM2584556", "GSM2584557", "GSM2584558", "GSM2584559",
             "GSM2584560", "GSM2584561", "GSM2584562", "GSM2584563",
             "GSM2584564", "GSM2584565", "GSM2584566", "GSM2584567")

# Tratamiento, tiempo, réplica
treatment=rep(c("4OHT", "BSA", "IFNg", "Ly29002", "TGF1b", "TNFa"), each = 4)
time=rep(c("1h", "1h", "4h", "4h"), 6)
replicate=rep(c(1, 2), times = 12)

# Crear el metadata
metadata=data.frame(sample = samples,
                       treatment = treatment,
                       time = time,
                       replicate = replicate)

# Opcional: combinar columnas para usar como etiqueta en gráficas
metadata$condition=paste(metadata$treatment, metadata$time, sep = "_")


head(raw_counts[, 1:3])
dge=DGEList(counts = raw_counts)
dge=calcNormFactors(dge, method = "TMM")
#se obtienen counts per million (cpm) y se hace el log
datos_TMM=cpm(dge, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE) 

#PCA
pca_model= do_pca(datos_TMM[1:ncol(datos_TMM),])
par(mar = c(5, 4, 4, 8), xpd = TRUE) 
pca_plot(pca_model, metadata$condition, legend = TRUE)

#truncamiento al percentil 99
quantile_99=apply(datos_TMM, 1, function(x) quantile(x, 0.99))
datos_TMM_trunc=pmin(datos_TMM, quantile_99)

#PCA2
pca_model2= do_pca(datos_TMM_trunc[1:ncol(datos_TMM_trunc),])
par(mar = c(5, 4, 4, 8), xpd = TRUE) 
pca_plot(pca_model2, metadata$condition, legend = TRUE)

#normalizacion hipathia
datos_norm=normalize_data(datos_TMM , by_quantiles = TRUE)

pca_model3= do_pca(datos_norm[1:ncol(datos_norm),])
par(mar = c(5, 4, 4, 8), xpd = TRUE) 
pca_plot(pca_model3, metadata$condition, legend = TRUE)

#boxplot muestras 
# Crear 12 colores
colores=rainbow(12)

# Convertir condiciones a factor
condiciones=as.factor(metadata$condition)

# Mapear colores a cada muestra
col_vector=colores[condiciones]

# Graficar boxplot
boxplot(datos_norm,
        las = 2,
        main = "Boxplot por muestra",
        ylab = "Expresión",
        col = col_vector,
        names = metadata$condition,
        cex.axis = 0.7)


results=hipathia(datos_norm, pathways_signor, decompose = FALSE, verbose=FALSE)
results

path_vals=get_paths_data(results, matrix = TRUE)