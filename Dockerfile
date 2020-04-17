FROM bimberlab/oosap

RUN apt-get update && \
apt-get install -y git libxml2-dev libssl-dev ghostscript

# Install dependencies & package
RUN Rscript -e "install.packages(c('devtools', 'BiocManager', 'remotes'), dependencies=TRUE, ask = FALSE)" \
    && echo -e "local({\noptions(repos = BiocManager::repositories())\n})\n" >> ~/.Rprofile.site \
    # && Rscript -e "BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db', ‘org.Mmu.eg.db’), dependencies=TRUE, ask = FALSE)" \
    && Rscript -e "devtools::install_github(repo = 'kolabx/geneSetVis@pkg', ref = 'Dev', dependencies = T, upgrade = 'always')" \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# select port
EXPOSE 3838

# run app
CMD geneSetVis::launchGeneSetVis()
