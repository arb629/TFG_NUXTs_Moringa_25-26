# TFG - ESTUDIO COMPARATIVO DE LA FRACCIÓN DE ADN NUCLEAR DE ORIGEN ORGANELAR ENTRE LAS ESPECIES MORINGA STENOPETALA Y M. OLEIFERA
Este repositorio contiene los scripts y análisis desarrollados para el Trabajo de Fin de Grado (TFG), centrado en la detección, análisis y visualización de secuencias NUMTs y NUPTs en diferentes genomas.

## Descripción
Se incluyen scripts en R y Python para:
- Extracción y procesamiento de datos genómicos
- Eliminación de duplicados
- Análisis estadístico de correlaciones
- Visualización de resultados (figuras y diagramas)

## Objetivo del repositorio
Organizar de forma reproducible todo el pipeline de análisis del TFG, permitiendo su reutilización y validación por terceros.

## Estructura del repositorio
- Extraccion_BASEML/        # Scripts de extracción y procesamiento
- Figuras_R/               # Scripts de generación de figuras en R
- NUXTs_comunes_especificos/ # Análisis comparativo en Python
- filtro_numt_no_duplicaciones.R

## Scripts principales
- correlacion_k_tamano_v1.R → análisis de correlación entre K y tamaño
- hist_k_tamano_v2.R → histogramas de distribución
- diagrama_venn.py → comparación de conjuntos NUMTs/NUPTs
- figura_5_zhang2020.R → replicación de figura de estudio previo
- baseml_*.py → extracción y procesamiento de datos BASEML

## Requisitos
- R (>= 4.0)
- Python (>= 3.8)
- Librerías R: ggplot2, dplyr, etc.
- Librerías Python: pandas, matplotlib, seaborn

## Autora
Andrea Ríos Bravo  
Grado en Biotecnología  
TFG 2025-2026

