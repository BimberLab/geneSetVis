
#packrat::on()
#packrat::snapshot(".", snapshot.sources = FALSE)

options(repos=structure(BiocManager::repositories()))

# library(shiny)
# library(shinydashboard)
# library(shinyWidgets)
# library(flexdashboard)
# library(shinytest)
# library(testthat)
# 
# library(dplyr)
# library(ggplot2)
# library(ggupset)
# library(tidyr)
# library(stringr)
# library(plotly)
# library(DT)
# library(formattable)
# library(stats)
# library(knitr)
# library(rmdformats)
# 
# library(org.Hs.eg.db)
# library(org.Mmu.eg.db)
# library(org.Mm.eg.db)
# library(AnnotationHub)
# library(ReactomePA)
# library(clusterProfiler)
# library(DOSE)
# library(enrichplot)
# library(fgsea)
# library(enrichR)
# 
# library(biomaRt)
# library(STRINGdb)
# library(msigdbr)


wd <- '.'
setwd(wd)


source('../R/fxs.R', local = TRUE)
source('../R/GeneAliasing.R', local = TRUE)
source('modules/stringdb.R', local = TRUE)
source('modules/msigdb.R', local = TRUE)
source('modules/reactome.R', local = TRUE)
source('modules/david.R', local = TRUE)
source('modules/dose.R', local = TRUE)
source('modules/ncg.R', local = TRUE)
source('modules/dgn.R', local = TRUE)
source('modules/enrichr.R', local = TRUE)
source('tabs.R', local = TRUE)
source('ui.R', local = TRUE)
source('server.R', local = TRUE)


demo1 <- "HMOX1	1.0596510 \nRNF167	0.9790608 \nHSPA5	0.7293491 \nCDKN1A	0.7265868 \nFCGR2B	0.6369659 \nPFN1	0.5453499 \nLAPTM5	0.5164539 \nAHNAK	0.5045917 \nFN1	0.4090008 \nS100A10	0.3566574 \nVIM	0.3409602 \nYWHAZ	0.2911121 \nFTH1.1	0.2733286 \nPDIA3	0.2555106 \nATP5MPL	-0.2565952 \nLAMTOR4	-0.2574608 \nSMDT1	-0.2589715 \nCOX5A	-0.2610802 \nMTDH	-0.2619066 \nNDUFA2	-0.2638782 \nCOX6C	-0.2679750 \nCOX8A	-0.2756591 \nNDUFA1	-0.2781574 \nH2AFJ	-0.2827520 \nTOMM7	-0.2955068 \nRPL23	-0.3009606 \nCOX7C	-0.3324625 \nCASP1	-0.3531754 \nRPS21	-0.3921719 \nRPL38	-0.3928734 \nFOS	-0.8496947 \nIGFBP1	-2.2179911 \nPPT1	0.2956121 \nHEXB	0.2665466 \nNINJ1	0.3056079 \nFGL2	0.2589270 \nLDHA	0.2736736 \nCD59	-0.3042252 \nGSN	0.2728750 \nANXA2	0.2990603 \nLGALS3	0.2911058 \nSLC2A3	0.4835044 \nMT-CO2	-0.3797473 \nPLIN2	0.2974303 \nPLAUR	0.2979632 \nPPP1R15A	0.3040476"
demo2 <- "HMOX1, RNF167, HSPA5, CDKN1A, FCGR2B, PFN1, LAPTM5, AHNAK, FN1, S100A10, VIM, YWHAZ, FTH1.1, PDIA3, ATP5MPL, LAMTOR4, SMDT1, COX5A, MTDH, NDUFA2, COX6C, COX8A, NDUFA1, H2AFJ, TOMM7, RPL23, COX7C, CASP1, RPS21, RPL38, FOS, IGFBP1, PPT1, HEXB, NINJ1, FGL2, LDHA, CD59, GSN, ANXA2, LGALS3	, SLC2A3, MT-CO2, PLIN2, PLAUR, PPP1R15A"

