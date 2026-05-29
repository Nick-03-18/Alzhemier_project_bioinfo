#Enfermedad de Alzheimer de Origen Tardío (LOAD): Análisis preliminar de estructura poblacional y frecuencia genotípica del gen APOE en el mundo
#Autores:J. Nicolas Avila M. & Sofia Lombana F.

#Datos: proyecto 1000 genomes fase 3 cromosoma 19 de 2504 individuos

#cargo las liberias
library(data.table)
library(dplyr)
library(pheatmap)
library(reshape2)
library(ggplot2)
library(viridis)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(scatterpie)
library(ggnewscale)
library(grid)

#leo la tabla consolidada de genotipos y estadisticas
datos <- fread(
  "analysis/APOE_summary.tsv",
  sep = "\t",
  header = TRUE
)

#estructuro la matriz de conteos para poder hacer las chi cuadrados

matriz <- datos %>%
  select(SUPERPOP,E2E2,E2E3,E2E4,E3E3,E3E4,E4E4)
#matriz numerica
matriz_counts <- as.matrix(matriz[, -1])
#nobmres de filas
rownames(matriz_counts) <- matriz$SUPERPOP
#imprimo la matriz para verificar
print(matriz_counts)

#aplico la prueba chi cuadrado global y luego la imprimo
chi <- chisq.test(matriz_counts)
print(chi)

#exportamos las freucencias esperadas
write.table(
  round(chi$expected, 2),
  "analysis/chi_expected.tsv",
  sep = "\t",
  quote = FALSE
)

#guardamos los resudiales estandarizados en stdres

stdres <- chi$stdres

#exportamos a una tabla esos residuales
write.table(
  round(stdres, 2),
  "analysis/chi_stdres.tsv",
  sep = "\t",
  quote = FALSE
)

#heatmap de los residuales estandarizados

#creamos el pdf que quedara con el heatmap de residuales estandarizados
pdf(
  "analysis/APOE_stdres_heatmap.pdf",
  width = 10,
  height = 8
)

#pintamos el heatmap con pheatmap, usamos los residuales estandarizados
pheatmap(
  stdres,
  scale = "none",
  #agrupamos similares
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  #mostrar los valores residuales estandarizados en cada cuadrado
  display_numbers = TRUE,
  #un decimal
  number_format = "%.1f",
  #titulo
  main = "APOE standardized residuals"
)
#cerramos el pdf
dev.off()

###AHORA, hacemos comparaciones pareadas con pruebas chicuadrado pareadas
#nombres de poblaciones (superpoblaciones)
pops <- rownames(matriz_counts)
#creamos un dataframe vacio
resultados <- data.frame()
#recorremos superpoblacion una por una
for(i in 1:(length(pops)-1)){
#comparamos entre poblaciones
  for(j in (i+1):length(pops)){
    #extraemos nombre superpop1
    pop1 <- pops[i]
    #extraemos nombre superpop1
    pop2 <- pops[j]
    #seleccionamos solo los datos de las poblaciones relevantes en la matriz
    submat <- matriz_counts[c(pop1, pop2), ]

    # quitar vacios si hubiera
    submat <- submat[, colSums(submat) > 0]
    #aplicamos el chi cuadrado pareado
    test <- chisq.test(
      submat,
      simulate.p.value = TRUE,
      #simulaciones
      B = 10000
    )

#añadimos resultados al dataframe principal
    resultados <- rbind(
      resultados,
      data.frame(
        Population1 = pop1,
        Population2 = pop2,
        Chi_square = test$statistic,
        df = test$parameter,
        p_value = test$p.value
      )
    )
  }
}

##Luego corregimos por FDR (False Discovery Rate)
resultados$FDR <- p.adjust(
  resultados$p_value,
  method = "fdr"
)
#ordenamos por FDR
resultados <- resultados %>%
  arrange(FDR)
#guardamos en .tsv
write.table(
  resultados,
  "analysis/APOE_pairwise_chi2.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


#Construimos el heatmap pareado
#extraer superpoblaciones
pops <- unique(c(resultados$Population1,resultados$Population2))
#matriz vacia
fullmat <- matrix(1,nrow = length(pops),ncol = length(pops) )
#asignamos nombres a las filas y columnas
rownames(fullmat) <- pops
colnames(fullmat) <- pops
#ahora lo que hacemos es recorrer cada comparación
for(i in 1:nrow(resultados)){
#extraemos las poblaciones y fdr
  p1 <- resultados$Population1[i]
  p2 <- resultados$Population2[i]
  fdr <- resultados$FDR[i]
  #llenamos la matriz vacia que hicimos
  fullmat[p1, p2] <- fdr
  fullmat[p2, p1] <- fdr
}
#corregir que las poblaciones iguales sean 1
diag(fullmat) <- 1
#transformación logaritmica para que los más altos sean más significativos
logmat <- -log10(fullmat)
#quitar problemas de ceros e infinito por la trnasformación, si es infinito deja el valor finito mayor
logmat[is.infinite(logmat)] <- max(logmat[is.finite(logmat)])
#que la matriz se vea mejor
dfplot <- melt(logmat)
#nombres de columnas
colnames(dfplot) <- c("Population1","Population2","logFDR")

#ahora si creamos el grafico heatmap
p <- ggplot(
  dfplot,
  #asignamos las variables
  aes(
    x = Population1,
    y = Population2,
    fill = logFDR
  )
) +
#llenamos las celdas
  geom_tile() +
  #darle colores
  scale_fill_viridis(
    name = "-log10(FDR)"
  ) +
  theme_minimal() +
  theme(
  #rotamos nombres de poblaicones evitar el sobrelape
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  ) +
  #damos titulo
  labs(
    title = "Pairwise APOE differences"
  )
#guardamos la figura
ggsave(
  "analysis/APOE_pairwise_heatmap.pdf",
  p,
  width = 10,
  height = 9
)


#Clustering jerarquico
#pasamos de conteos a freq relativas
freqmat <- sweep(
  matriz_counts,
  1,
  rowSums(matriz_counts),
  FUN = "/"
)
#calculamos las distancias entre las filas osea entre las superpoblaciones con el metodo estandar euclidiano
distmat <- dist(
  freqmat,
  method = "euclidean"
)
#consturimos el dendrograma que agrupa poblaciones segun similaridad
hc <- hclust(
  distmat,
  #indicamos como se agrupo
  #ward=agrupar causando el menor aumento posible en la variablidad interna
  method = "ward.D2"
)
#iniciamos un pdf para guardar la figura
pdf(
  "analysis/APOE_population_clustering.pdf",
  width = 10,
  height = 8
)
#dibujamos sobre el pdf
plot(
  hc,
  main = "APOE population clustering"
)
#cerramos el pdf
dev.off()




#MAPA CON PIE CHARTS DE FRECUENCIAS RELATIVAS

#hacemos un data frame donde a las superpoblaciones se les asigna que paises contienen
superpop_map <- data.frame(
  country = c(
    # AFR
    "Nigeria",
    "Sierra Leone",
    "Gambia",
    "Kenya",
    # EUR
    "United Kingdom",
    "Finland",
    "Spain",
    "Italy",
    # EAS
    "China",
    "Japan",
    "Vietnam",
    # SAS
    "India",
    "Pakistan",
    "Bangladesh",
    "Sri Lanka",
    # AMR
    "Mexico",
    "Peru",
    "Colombia",
    "Puerto Rico",
    "United States"

  ),
  SUPERPOP = c(
#repetimos el nombre de superpoblación segun el numero de paises
    rep("AFR", 4),
    rep("EUR", 4),
    rep("EAS", 3),
    rep("SAS", 4),
    rep("AMR", 5)
  )
)

#cargamos el mapa mundial
world <- ne_countries(scale = "medium",
  #formato espacial sf
  returnclass = "sf")

#unimos el mapa base con el data frame, unimos por el nombre de pais para consolidar superpoblaciones con el formato espacial de cada pais
world2 <- left_join( world, superpop_map, by = c("name" = "country"))


#tabla de posiciones del pie chart por latitud y longitud
pie_positions <- data.frame(
  SUPERPOP = c("AFR","EUR","EAS","SAS","AMR"),lon = c(20,15,110,78,-75), lat = c(5,55,35,22,10))
#unimos geografía con datos genotipicos
pies <- left_join( pie_positions, datos,  by = "SUPERPOP")

#creamos nuevas columnas con las frecuencias relativas
pies <- pies %>%
  mutate(
    E2E2_rel = E2E2 / N * 100,
    E2E3_rel = E2E3 / N * 100,
    E2E4_rel = E2E4 / N * 100,
    E3E3_rel = E3E3 / N * 100,
    E3E4_rel = E3E4 / N * 100,
    E4E4_rel = E4E4 / N * 100
  )

#damos escalas de colores diferentes a superpoblaciones y genotipos
  superpop_colors <- c( AFR = "#D73027", EUR = "#4575B4", EAS = "#1A9850", SAS = "#984EA3", AMR = "#FF7F00")
  geno_colors <- c(E2E2_rel = "#1b9e77", E2E3_rel = "#66c2a5", E2E4_rel = "#fc8d62", E3E3_rel = "#8da0cb",E3E4_rel = "#e78ac3",E4E4_rel = "#d73027")

#creamos el mapa
  pmap <- ggplot() +
    #dibujamos los paises
    geom_sf(
      data = world2,
      aes(fill = SUPERPOP),
      color = "black",
      linewidth = 0.2
    ) +
    #coloreamos segun superpoblacion
    scale_fill_manual(
      values = superpop_colors,
      na.value = "gray92",
      name = "Superpopulation"
    ) +
    #asignamos la escala de color diferente para los genotipos
    new_scale_fill() +
    #dibujamos pie charts
    geom_scatterpie(data = pies,aes( x = lon,y = lat),
    #variables del pie chart
      cols = c(
        "E2E2_rel",
        "E2E3_rel",
        "E2E4_rel",
        "E3E3_rel",
        "E3E4_rel",
        "E4E4_rel"
      ),pie_scale = 2.3,color = "black") +
    #colores de las variables del pie chart
    scale_fill_manual(
      values = geno_colors,
      name = "Genotype"
    ) +
    #etiquetas bonitas con nombre de superpoblación por pie chart
    geom_label(data = pies,aes(x = lon,y = lat + 10,label = SUPERPOP),fontface = "bold",size = 4,fill = "white",color = "black",label.size = 0.4,label.padding = unit(  0.18,"lines"),label.r = unit(0.15,"lines")
    )+
    #coordenadas del mapa, limites
    coord_sf(xlim = c(-180, 180),ylim = c(-60, 85),expand = FALSE) +
    #tema
    theme_minimal() +
    theme(
      #oceano azul
      panel.background = element_rect(
        fill = "aliceblue",
        color = NA
      ),
      #lineas referencia
      panel.grid.major = element_line(color = "white",  linewidth = 0.2),legend.position = "right",axis.text = element_blank(),  axis.title = element_blank(),plot.title = element_text(face = "bold",size = 18),
      plot.subtitle = element_text(size = 12)) +
      #titulo y subtitulo
    labs(title = "Global APOE superpopulation structure",subtitle = "Pie charts show APOE genotype frequencies"
  )
#guardar mapa
ggsave(
  "analysis/APOE_worldmap_pies.pdf",
  pmap,
  width = 16,
  height = 9
)
