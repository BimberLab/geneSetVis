FROM rocker/shiny

RUN apt-get update && \
    apt-get install -y git libxml2-dev libssl-dev ghostscript
    
# copy the app to the image
COPY . /srv/shiny-server/genesetvis
 
WORKDIR /srv/shiny-server/genesetvis

# Install packrat packages
#RUN Rscript -e 'install.packages("packrat"); \
#                packrat::restore()'
                
# permissions
RUN sudo chown shiny:shiny -R packrat && \
    sudo chown shiny:shiny .gitignore 

# Serve only the genesetvis app 
RUN sed -i 's/\/srv\/shiny-server/\/srv\/shiny-server\/genesetvis/' /etc/shiny-server/shiny-server.conf

# Increase timeout
RUN sed -i '/location \/ {/a app_init_timeout 180;' /etc/shiny-server/shiny-server.conf
#RUN sed -i '/location \/ {/a http_keepalive_timeout 180;' /etc/shiny-server/shiny-server.conf

# select port
EXPOSE 3838

# run app
CMD /usr/bin/shiny-server.sh
