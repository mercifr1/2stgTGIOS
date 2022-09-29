# ####################################################
#
#   2stgTGIOS
#   Data preparation
#   
#   Q1-2022
#   Francois Mercier 
#
# ####################################################


# ===============================================
# REQUIRED SOURCES
# ===============================================
# Package names
pkg1 <- c("tidyverse", "haven")

# Install packages not yet installed
installed_packages <- pkg1 %in% rownames(installed.packages())
if(any(installed_packages == FALSE)){ 
  install.packages(pkg1[!installed_packages])
}

# Packages loading
invisible(lapply(pkg1, library, character.only = TRUE))


# ===============================================
# LOADING DATASETS
# ===============================================
subj <- haven::read_sas("../datalocal/rdpsubj.sas7bdat")
rcist <- haven::read_sas("../datalocal/rdprcist.sas7bdat") 


# ===============================================
# FILTERING VARIABLES
# ===============================================
keepin0 <- c("UID", "LSVISDY", "SEX", "BEV_SDY", "BEV_EDY", "DEATDI", "LDH1_5", "BL_VEGFN", "VEGF_STR")

# Keeping patients in PerProtocol set (PP_SET==1)
# Keeping rows from patients exposed to BEV for at least 3 weeks (i.e. at least 1st TA visit)
subj0 <- subj %>% 
  filter(PP_SET==1, !is.na(BEV_SDY), !is.na(BEV_EDY), BEV_EDY>3*7) %>%
  mutate(UID=RANDCODE) %>% select(one_of(keepin0))
# dim(subj0)
# [1] 645  9

# Overall survival (OS) values
os0 <- rcist %>% mutate(UID=RANDCODE) %>%
  filter(UID %in% subj0$UID) %>% group_by(UID) %>%
  slice(1) %>% ungroup() %>% mutate(OSYEAR=OSTIM/365.25) %>%
  select(UID, DEATFLAG, OSTIM, OSCEN, OSYEAR)
# length(os0$UID)
# 645

# Sum of longest diameters (SLD) values
# Removing rows where SLD is NA
# Removing rows where patients have PBLCNT=NA, i.e. removing patient with at least one post-baseline TA
tk0 <- rcist %>% mutate(UID=RANDCODE) %>%
  filter(UID %in% subj0$UID) %>% drop_na(STLDI) %>%
  filter(!is.na(PBLCNT)) %>%
  mutate(SLD=ifelse(STLDI==0, 2.5, STLDI*10), 
         BSLD=BL_STLDI*10, TKYEAR=ORDYTRT/365.25) 
# length(unique(tk0$UID))
# 640


# ===============================================
# BUILDING ANALYSIS DATASETS
# ===============================================
os1 <- left_join(subj0, os0, by="UID") %>%
  select(UID, SEX, OSTIM, OSCEN, OSYEAR, LDH1_5, BL_VEGFN, VEGF_STR)

tk1 <- left_join(tk0, subj0, by="UID") %>%
  select(UID, ORDYTRT, TKYEAR, BL_VIS, BSLD, SLD, PBLCNT)

# Joining the two datasets and eliminating patients have LDH!_5=NA
joinostk <- left_join(tk1, os1, by="UID") %>% drop_na(LDH1_5)

# Eliminating individuals with 1 observation
pos_rmv1 <- as.numeric(names(table(joinostk$UID)[which(as.vector(table(joinostk$UID)==1))]))
joinostk1 <- filter(joinostk, !(UID %in% pos_rmv1))

# Eliminating individuals without 'time when the treatment was initiated'
pos_rmv2 <- joinostk1$UID[which(joinostk1$BL_VIS=="1" & joinostk1$TKYEAR>0)]
joinostk2 <- filter(joinostk1, !(UID %in% pos_rmv2))

# Eliminating individuals with last longitudinal measurement greater than survival time
pos_rmv3 <- joinostk2$UID[which(joinostk2$TKYEAR > joinostk2$OSYEAR)]
dta <- filter(joinostk2, !(UID %in% pos_rmv3))
# length(unique(dta$UID))
# 593

# saveRDS(dta, file="../datalocal/HorizOSTGI.rds")
