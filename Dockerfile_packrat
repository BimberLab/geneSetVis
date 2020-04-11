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

# copy the app to the image
COPY . /usr/local/src/genesetvis/

WORKDIR /usr/local/src/genesetvis/

# Install packrat packages
RUN Rscript -e 'install.packages("packrat"); \
                packrat::restore()'
                
# select port
EXPOSE 3838

# run app
CMD ["R", "-e", "shiny::runApp('/usr/local/src/genesetvis/', host ='0.0.0.0', port = 3838, launch.browser = FALSE)"]
