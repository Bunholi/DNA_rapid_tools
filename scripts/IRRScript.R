# =======================================================
# Inter-Rater Reliability (IRR) Analysis Pipeline
# =======================================================
# Authors: Marcela Alvarenga
# Institution: BIOPOLIS/CIBIO, UPorto, Portugal
# Date: March 2025
# 
# Description:
# This script calculates Fleiss' Kappa to assess interrater agreement
# for multiple raters evaluating categorical data.
# Three different assessments are analyzed:
#  1. Primers assessment (R = 0.769)
#  2. DNA/eDNA-based Rapid Tools assessment (R = 0.551)
#  3. DNA/eDNA-based Rapid Tools reassessment (R = 0.928)
#
# Methods:
# - Fleiss' Kappa (Fleiss, 1971) for multiple raters
# - Optionally, exact Kappa (Conger, 1980) and category-wise Kappas
#
# Requirements:
# - Install required packages: irr
# - Input data: MatrixIRR.xlsx (contains ratings from multiple raters)
#
# =======================================================

install.packages("irr")
library(irr)

####FLEISS'S KAPPA FOR M RATERS####
#Computes Fleiss Kappa as an index of interrater agreement between m raters (the reviewers) on categorical data 
#Additionally, category-wise Kappas could be computed

##Computation##
#kappam.fleiss(ratings, exact = FALSE, detail = FALSE)
#ratings - n*m matrix or dataframe, n subjects m raters.
#exact - a logical indicating whether the exact Kappa (Conger, 1980) or the Kappa described by Fleiss (1971) should be computed.
#detail - a logical indicating whether category-wise Kappas should be computed

## IRR for primers assessment ## 
matrix_ws<-readxl::read_xlsx("MatrixIRR.xlsx", sheet="primer", col_names = TRUE)
str(matrix_ws)

#The null hypothesis Kappa=0
kappam.fleiss(matrix_ws)
kappam.fleiss(matrix_ws, exact = TRUE)
kappam.fleiss(matrix_ws, detail = TRUE)

## IRR for rapid tools assessment ## 
matrix_s<-readxl::read_xlsx("MatrixIRR.xlsx", sheet="rapid", col_names = TRUE)
str(matrix_s)

#The null hypothesis Kappa=0
kappam.fleiss(matrix_s)
kappam.fleiss(matrix_s, exact = TRUE)
kappam.fleiss(matrix_s, detail = TRUE)

## IRR for reassessment ## R = 0.927
matrix_t<-readxl::read_xlsx("MatrixIRR.xlsx", sheet="reassessment-rapid", col_names = TRUE)
str(matrix_t)

#The null hypothesis Kappa=0
kappam.fleiss(matrix_t)
kappam.fleiss(matrix_t, exact = TRUE)
kappam.fleiss(matrix_t, detail = TRUE)

