library("hipathia")
library(tidyr)
library(org.Hs.eg.db)
library(tibble)
library(readr)
library(edgeR)
library(MultiAssayExperiment)
library(SummarizedExperiment)
pathways=hipathia::load_pathways(species = "hsa")


## Nombre de las rutas
#Creamos una función que acceda a cada elemento de la lista pathigraphs (que representa una ruta especifica de señalización por ejemplo,hsa03320) y nos devuelva el nombre de la ruta

all_path_names=sapply(pathways$pathigraphs, function(x) x$path.name)
# sapply aplica una función a cada elemento de una lista
# la función toma un solo argumento x 
# x$path.name accede al nombre de la ruta 


## Explorar una ruta específica: hsa03320

all_path_names["hsa03320"] #para ver el nombre de esta ruta
pathways$pathigraphs$hsa03320$subgraphs %>% length() #Cada subgrafo representa una posible trayectoria completa desde un nodo de entrada hasta uno de salida.
pathways$pathigraphs$hsa03320$effector.subgraphs %>% length() #hay menos nodos efectores


##  Explorar un subgrafo de un efector específico

pathways$pathigraphs$hsa03320$effector.subgraphs$`P-hsa03320-37` %>% class
pathways$pathigraphs$hsa03320$effector.subgraphs$`P-hsa03320-37` %>% plot 
pathways$pathigraphs$hsa03320$effector.subgraphs$`P-hsa03320-37` %>%  V() #devuelve los vértices (nodos) del grafo, que representan genes/proteínas


## Extraer los genes efectores

subG=pathways$pathigraphs$hsa03320$effector.subgraphs$`P-hsa03320-37`#trayectoria desde un nodo de entrada a uno efector
V(subG)[degree(subG, mode = "out")==0]$genesList
#se accede a todos los nodos y se selecciona los nodos que cumple la condición de que sean terminales (no tenga nodos posteriores)
#genesList devuelve los genes asociados al nodo efector de esta subvía


#Se convierte a un dataframe para visualizar todos los nodos (vértices)

igraph::as_data_frame(subG, what = "vertices") %>% View


## Función para extraer los nodos efectores de cada subgrafo

get_entrez_effector=function(subG){
  V(subG)[degree(subG, mode = "out")==0]$genesList %>% unlist()  
}
#ejemplo
get_entrez_effector(pathways$pathigraphs$hsa03460$effector.subgraphs$`P-hsa03460-28`)

#Nos devuelve el ID del gen efector 



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
effector_gene_pathway_map= get_effector_gene_pathway_map(pathways)


# Crear el mapeo Entrez -> Symbol
entrez_ids= unique(effector_gene_pathway_map$gene)
gene_symbols= mapIds(org.Hs.eg.db,
                     keys = entrez_ids,
                     column = "SYMBOL",
                     keytype = "ENTREZID",
                     multiVals = "first")
# Añadir columna con símbolo a tu dataframe
effector_gene_pathway_map$symbol= gene_symbols[as.character(effector_gene_pathway_map$gene)]


# Paso 1: pasar los rownames (que contienen el circuito) a una columna llamada "circuito"
effector_gene_pathway_map = rownames_to_column(effector_gene_pathway_map, var = "circuito")

# Paso 2: convertir columnas de tipo lista a cadenas separadas por ;
effector_gene_pathway_map_flat = data.frame(
  lapply(effector_gene_pathway_map, function(col) {
    if (is.list(col)) {
      sapply(col, function(x) paste(x, collapse = ";"))
    } else {
      col
    }
  }),
  stringsAsFactors = FALSE
)

# Paso 3: guardar en CSV con la columna "circuito" incluida
write.csv(effector_gene_pathway_map_flat, file = "todos_efectores_y_rutas_KEGG.csv", row.names = FALSE)



# Función para buscar información de genes efectores
search_multiple_genes = function(effector_df, search_terms) {
  
  # Filtrar por símbolo de gen (ignora mayúsculas/minúsculas)
  gene_info = effector_df[toupper(effector_df$symbol) %in% toupper(search_terms), ]
  
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

# Convertir rownames a columna explícita
#result_multiple= result_multiple %>% 
  #tibble::rownames_to_column(var = "circuito")  # Asigna un nombre a la primera columna

# Desanidar listas
result_multiple_flat=result_multiple %>% 
  unnest(cols = c(symbol))

# Guardar en CSV, asegurando que se incluya 'circuito'
write.csv(result_multiple_flat, file = "efectores_y_rutas_KEGG.csv", row.names = FALSE)


#importamos los datos de los conteos
#partimos de los datos crudos
raw_counts=read.delim("/Users/paula/Documents/Bioinformatica/tfm/codigo/data/GSE97979_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1, check.names = FALSE)
str(raw_counts)

#normalizamos con TMM
dge=DGEList(counts = raw_counts)
dge=calcNormFactors(dge, method = "TMM")
#se obtienen counts per million (cpm) y se hace el log
datos_TMM=cpm(dge, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE) 

# #PCA
# pca_model= do_pca(datos_TMM[1:ncol(datos_TMM),])
# par(mar = c(5, 4, 4, 8), xpd = TRUE) 
# pca_plot(pca_model, metadata$condition, legend = TRUE)

#truncamiento al percentil 99
quantile_99=apply(datos_TMM, 1, function(x) quantile(x, 0.99))
datos_TMM_trunc=pmin(datos_TMM, quantile_99)

# #PCA2
# pca_model2= do_pca(datos_TMM_trunc[1:ncol(datos_TMM_trunc),])
# par(mar = c(5, 4, 4, 8), xpd = TRUE) 
# pca_plot(pca_model2, metadata$condition, legend = TRUE)

#normalizacion hipathia
datos_norm=normalize_data(datos_TMM_trunc , by_quantiles = TRUE)

# pca_model3= do_pca(datos_norm[1:ncol(datos_norm),])
# par(mar = c(5, 4, 4, 8), xpd = TRUE) 
# pca_plot(pca_model3, metadata$condition, legend = TRUE)

# #boxplot muestras 
# # Crear 12 colores
# colores=rainbow(12)
# 
# # Convertir condiciones a factor
# condiciones=as.factor(metadata$condition)
# 
# # Mapear colores a cada muestra
# col_vector=colores[condiciones]
# 
# # Graficar boxplot
# boxplot(datos_norm,
#         las = 2,
#         main = "Boxplot por muestra",
#         ylab = "Expresión",
#         col = col_vector,
#         names = metadata$condition,
#         cex.axis = 0.7)


results=hipathia(datos_norm, pathways, decompose = FALSE, verbose=FALSE)

# Extraer los datos de expresión (matrices) desde cada experimento
nodes_expr= assay(experiments(results)[["nodes"]])
paths_expr= assay(experiments(results)[["paths"]])

# Convertir a data.frame para exportar (añadir nombres de fila si quieres conservarlos)
nodes_df= as.data.frame(nodes_expr)
nodes_df$Gene= rownames(nodes_expr)

paths_df= as.data.frame(paths_expr)
paths_df$Pathway = rownames(paths_expr)

# Guardar en archivos CSV
write.csv(nodes_df, file = "expresion_nodes_KEGG.csv", row.names = TRUE)
write.csv(paths_df, file = "expresion_paths_KEGG.csv", row.names = TRUE)

