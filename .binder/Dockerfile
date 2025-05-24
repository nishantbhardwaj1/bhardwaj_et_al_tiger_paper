FROM rocker/binder:4.4.1
LABEL maintainer='Nishant_Bhardwaj'
COPY --chown=${NB_USER} . ${HOME}
USER ${NB_USER}



RUN wget https://github.com/nishantbhardwaj1/bhardwaj_et_al_tiger_paper/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2025-05-24/')); devtools::install_deps()"

RUN rm DESCRIPTION.1; exit 0
