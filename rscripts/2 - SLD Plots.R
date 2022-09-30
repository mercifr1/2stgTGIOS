# ####################################################
#
#   2stgTGIOS
#   Descriptive plots for SLD
#
#   Q1-2022
#   Francois Mercier 
#
# ####################################################


# ===============================================
# REQUIRED SOURCES
# ===============================================
# Package names
pkg2 <- c("ggplot2")

# Install packages not yet installed
installed_packages <- pkg2 %in% rownames(installed.packages())
if(any(installed_packages == FALSE)){ 
  install.packages(pkg2[!installed_packages])
}

# Packages loading
invisible(lapply(pkg2, library, character.only = TRUE))


# ===============================================
# LOADING DATASET
# ===============================================
dta <- readRDS("../datalocal/HorizOSTGI.rds")


# ===============================================
# VISUALISING SUM OF LONGEST DIAMETERS (SLD)
# ===============================================
# Display SLD spaghetti plot
ybreaks <- c(0, 50, 150, 300, 600)
ggplot(dta, aes(x=TKYEAR, y=SLD)) +
  geom_line(aes(group=UID), colour="wheat4", alpha=0.2) +
  geom_point(colour="grey33", alpha=0.3, size=0.9) +
  scale_x_continuous("Year", breaks=0.5*(0:5)) +
  scale_y_continuous("SLD (mm)", breaks=ybreaks) +
  theme_minimal() + theme(panel.grid.minor=element_blank())

# Display SLD spaghetti plot in Male vs Female
ggplot(dta, aes(x=TKYEAR, y=SLD)) +
  geom_line(aes(group=UID), colour="wheat4", alpha=0.2) +
  geom_point(colour="grey33", alpha=0.3, size=0.9) +
  facet_wrap(~SEX, labeller=labeller(SEX = c("1"="Male", "2"="Female"))) +
  scale_x_continuous("Year", breaks=0.5*(0:5)) +
  scale_y_continuous("SLD (mm)", breaks=ybreaks) +
  theme_minimal() + theme(panel.grid.minor=element_blank())
