#datos fosfoproteomica
raw_data = read.delim("41467_2017_2391_MOESM5_ESM.txt", sep = "\t", header = FALSE, fill = TRUE)
head(raw_data, 3)

fosfoproteomica = read.delim("41467_2017_2391_MOESM5_ESM.txt", 
                              sep = "\t", 
                              header = TRUE, 
                              skip = 1, 
                              fill = TRUE)

# Guardar como CSV
write.csv(fosfoproteomica, file = "datos_fosfoproteomica.csv", row.names = FALSE)
