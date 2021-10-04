
# README SETUP ------------------------------------------------------------





# TABLE OF DEPENDENCIES ---------------------------------------------------

# ToDo: get versions and create a pretty table for markdown
# ToDo: make this automated to incorporate future changes


# 1. package 'pacman' can apparently do just that - BUT only after it has been installed
# so, for now stick with the manually curated list of dependencies, aka copy'n paste
install.packages("pacman")
library(pacman)
pacman::p_depends(mbecs, local=TRUE)


pkg_dependencies <- data.frame("Package"=c("bapred","cluster","dplyr","ggplot2","gridExtra","limma",
                                           "lme4","lmerTest","pals","permute","pheatmap","rmarkdown",
                                           "ruv","sva","tibble","tidyr","vegan","methods"),
                               "Version"=NA,
                               "Date"=NA,
                               "Repository"=NA)

for( pkg in 1:length(pkg_dependencies$Package) ) {
  pkg_dependencies$Version[pkg] <- toString(utils::packageVersion(eval(pkg_dependencies$Package[pkg])))
  tmp_description <- utils::packageDescription(eval(pkg_dependencies$Package[pkg]))
  pkg_dependencies$Date[pkg] <- toString(tmp_description["Date"])
  pkg_dependencies$Repository[pkg] <- toString(tmp_description["Repository"])
}


knitr::kable(pkg_dependencies,
             align = 'c',
             caption = "MBECS package dependencies",
             label = 'mbecs_dependencies_table')

#


















