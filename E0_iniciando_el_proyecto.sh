#!/bin/bash
#PREPARANDO TODO
# Creamos un entorno conda para tener un espacio propio del proyecto para evitar conflictos entre versiones ylibrerias.
conda create -n bioinfo python=3.11

# Activamos el entorno
conda activate bioinfo

# Configurar canales de donde se descargan las cosas (repositorios)

conda config --add channels defaults #oficial de conda
conda config --add channels bioconda #de bioinformática
conda config --add channels conda-forge #más paquetes

# Establecer el oreen en que se buscan los paquetes en los repositorios
conda config --set channel_priority strict

# Instalar paquetes y librerias de todo el proyecto
conda install bcftools vcftools htslib  plink wget
conda install -c conda-forge r-base
conda install -c conda-forge r-data.table
conda install -c conda-forge r-ggplot2 r-reshape2 r-pheatmap r-rcolorbrewer r-viridis
conda install -c conda-forge r-sf r-rnaturalearth r-rnaturalearthdata r-scatterpie
conda install -c conda-forge r-ggnewscale
