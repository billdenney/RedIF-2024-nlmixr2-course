# Script to support installation of all packages need in the workshop
# based on the instructions provided here:
# 
# https://nlmixr2.org/
#

# Installation of nlmixr2 for R versions incl. dependencies

install.packages("nlmixr2",dependencies = TRUE)


# check if everything is installed
nlmixr2::nlmixr2CheckInstall()


  
# install supportive packages
install.packages(c("xpose.nlmixr2", # Additional goodness of fit plots
                   # baesd on xpose
                   "nlmixr2targets", # Simplify work with the
                   # `targets` package
                   "babelmixr2", # Convert/run from nlmixr2-based
                   # models to NONMEM, Monolix, and
                   # initialize models with PKNCA
                   "nonmem2rx", # Convert from NONMEM to
                   # rxode2/nlmixr2-based models
                   "nlmixr2lib", # a model library and model
                   # modification functions that
                   # complement model piping
                   "nlmixr2rpt" # Automated Microsoft Word and
                   # PowerPoint reporting for nlmixr2
),
repos = c('https://nlmixr2.r-universe.dev',
          'https://cloud.r-project.org'))


remotes::install_github("ggPMXdevelopment/ggPMX") # Goodness of fit plots

# ShinyMixR

remotes::install_github("RichardHooijmaijers/shinyMixR") # Shiny run manager (like Piranha)

# tydiverse

install.packages("tidyverse") # The tidyverse is an opinionated collection of R packages designed for data science. 

# PKNCA

install.packages("PKNCA") # The PKNCA R package is designed to perform all noncompartmental analysis (NCA) calculations for pharmacokinetic (PK) data.