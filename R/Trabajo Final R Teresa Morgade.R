##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################


### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)
datos <- read.table("archivo-entrada.txt", header = TRUE, sep = "\t", dec = ".", check.names = FALSE)
head(datos)
tail(datos)
dim(datos)
### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)
datos$Media_Control <- apply(datos[, c("C2-TPM", "C1-TPM", "C3-TPM")], 1, mean)
datos$DE_Control <- apply(datos[, c("C2-TPM", "C1-TPM", "C3-TPM")], 1, sd)
datos$Media_Condicion <- apply(datos[, c("1W5-TPM", "1W6-TPM", "1W5-COV")],1, mean)
datos$DE_Condicion <- apply(datos[, c("1W5-TPM", "1W6-TPM", "1W5-COV")],1, sd)
head(datos)
### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).
par(mfrow = c(1, 2))
hist(datos$Media_Control, col = "skyblue", main = "Media Control", xlab = "TPM")
hist(datos$Media_Condicion, col = "pink", main = "Media Condicion", xlab = "TPM")
par(mfrow = c(1, 1))
### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)
boxplot(list(Control = datos$Media_Control, Condicion = datos$Media_Condicion),
	col = c("skyblue", "pink"),
	main = "Comparación Control vs Condición",
	ylab = "Valores TPM")
legend("topleft", legend = c("Control", "Condición"), fill = c("skyblue", "pink"))
### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.
p_values <- apply(datos, 1, function(fila) {
	control <- as.numeric(fila[c("C2-TPM", "C1-TPM", "C3-TPM")])
	condicion <- as.numeric(fila[c("1W5-TPM", "1W6-TPM", "1W5-COV")])
	res <- t.test(control, condicion)
	return(res$p.value)
})
datos$p_value <- p_values
### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)
genes_sig <- datos[datos$p_value <0.05,]
nrow(genes_sig)
sobreexpresados <- sum(genes_sig$Media_Condicion > genes_sig$Media_Control)
reprimidos <- sum(genes_sig$Media_Condicion < genes_sig$Media_Control)
print(paste("Sobreexpresados:", sobreexpresados))
print(paste("Reprimidos:", reprimidos))
### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)
matriz_heatmap <- as.matrix(genes_sig[, c("C2-TPM", "C1-TPM", "C3-TPM", "1W5-TPM", "1W6-TPM", "1W5-COV")])
rownames(matriz_heatmap) <- genes_sig$Gene.Name 
heatmap(matriz_heatmap, scale = "row", col = heat.colors(256),
	main = "Heatmap de Genes Diferencialmente Expresados")