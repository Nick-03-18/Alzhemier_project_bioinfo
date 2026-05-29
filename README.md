# Alzhemier_project_bioinfo
**Análisis preliminar de estructura poblacional y frecuencia genotípica del gen APOE en el mundo**
**Autores:J. Nicolas avila M. & Sofia Lombana F.**

# Abstract
Este estudio describe la distribución genotípica de APOE asociada con Alzheimer de inicio tardío en cinco superpoblaciones del Proyecto 1000 Genomes. A partir de 2054 archivos VCF del cromosoma 19, se analizaron los SNPs rs429358 y rs7412 para reconstruir seis genotipos APOE. Las frecuencias fueron comparadas mediante chi-cuadrado, residuos estandarizados, pruebas pareadas y agrupamiento jerárquico. Los resultados evidenciaron una distribución no homogénea entre superpoblaciones, con predominio de ε3/ε3 y mayor diferenciación en África, asociada a genotipos con ε4. Estos patrones resaltan la importancia de considerar ancestralidad, diversidad genética y contexto poblacional en salud pública.

Acceso a los datos:
En este github se presentan los siguientes archivos:
  * Archivos de entrada para todo el código:
      * ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
          * Debido a que el archivo pesaba más de 25 MB no fue posible subirlo por lo que se puede descargar desde la terminal con:                 curl -O            https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
      * ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi (indice)
      * integrated_call_samples_v3.20130502.ALL.panel.txt (panel muestras por población)
  * Archivos de trabajo:
      * E0_iniciando_el_proyecto.sh (Archivo ejecutable con bash que crea el ambiente de trabajo y descarga todos los paquetes y librerias necesaria (R))
      *  E1_manejo_datos_y_estadísticas.sh (Archivo ejecutable desde la terminal con bash que hace el análisis de datos inicial, extrae SNPs, construye los alelos y genotipos, determina las frecuencias genotípicas relativas y absolutas y las estadísticas poblacionales por superpoblación, consolidando todo en una tabla)
      *  E2_chi_cuadrado_global_pareado_clusteringjerarquico_mapa_pie.R (Archivo ejecutable desde la terminal con el comando Rscript donde se realizan las pruebas chi cuadrado global y pareada, se extraen residuales estandarizados, se realizan heatmaps y se realiza el clsutering jerarquizado)
  *Archivos comentados auxilares, particularmente para el archivo manejo_datos_y_estadísticas se comparte una versión complementaria comentada en su totalidad y denominada G1_manejo_datos_y_estadísticas.sh debido a que la versión E1_manejo_datos_y_estadísticas.sh tiene sus comentarios recortados para maximizar el flujo de ejecución y eivtar errores debido a la acción de comentar con "#" entre lineas.
    
