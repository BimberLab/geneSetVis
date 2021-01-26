FROM ghcr.io/bimberlabinternal/cellmembrane:latest

# install geneSetVis pkg
RUN Rscript -e "devtools::install_github(repo = 'kolabx/geneSetVis@pkg', dependencies = T, upgrade = 'always')" \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# select port
EXPOSE 3838

# run app
CMD ["R", "-e", "geneSetVis::launchGeneSetVis(port=3838, host='0.0.0.0')"]
