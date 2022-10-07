del ^sra_*.tar.gz
Rscript -e "roxygen2::roxygenize('.')"
Rscript version_update.R
Rcmd build .
Rcmd INSTALL sra_*.tar.gz

