FROM bimberlab/oosap

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    wget

RUN R -e "install.packages(pkgs=c('shiny', 'shinydashboard', 'shinyWidgets', 'shinyjs', 'flexdashboard', 'dplyr', 'ggplot2', 'ggupset', 'tidyr', 'stringr', 'plotly', 'DT', 'formattable', 'stats'), repos='https://cran.rstudio.com/')" \
    && echo -e "local({\noptions(repos = BiocManager::repositories())\n})\n" >> ~/.Rprofile.site \
    && Rscript -e "BiocManager::install(c('org.Hs.eg.db', 'org.Mmu.eg.db', 'org.Mm.eg.db', 'AnnotationHub', 'ReactomePA', 'clusterProfiler', 'DOSE', 'enrichplot', 'fgsea', 'enrichR', 'biomaRt', 'STRINGdb', 'msigdbr'), dependencies=TRUE, ask = FALSE)" 

##avoid caching of below commands
ARG CACHEBUST=3

# copy the app to the image
COPY genesetvis /usr/local/src/genesetvis/

# select port
EXPOSE 3838

# run app
CMD ["R", "-e", "shiny::runApp('/usr/local/src/genesetvis/', host ='0.0.0.0', port = 3838, launch.browser = TRUE)"]
