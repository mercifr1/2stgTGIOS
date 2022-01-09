#' ####################################################
#'   
#'   2stgTGIOS
#'   Data preparation
#'   
#'   Q1-2022
#'   Francois Mercier 
#'   
#' ####################################################

#' sessionInfo()


#' ===============================================
#' Import RDPSUBJ
#' ===============================================================

subj<-haven::read_sas("../data/rdpsubj.sas7bdat")
#' subj %>% count(PART)
#' subj %>% count(FULL_SET)
#' subj %>% count(PP_SET)
#' subj %>% count(SEX)

keepin0<-c("UID", "RND_DY", "LSVISDY", "SEX",
          "BEV_SDY", "BEV_EDY", "WHO_STR",
          "DIED", "DEATHDY", "DEATDI", "LDH1_5", 
          "TERMPREM", "TERMREA", "BL_VEGFN", "VEGF_STR")
subj0<-subj %>% 
  filter(PP_SET==1) %>%
  mutate(UID=RANDCODE) %>%
  select(one_of(keepin0))


#' ===============================================
#' Import RDPRCIST
#' ===============================================================

rcist<-haven::read_sas("../data/rdprcist.sas7bdat") 
rcista<-rcist %>% mutate(UID=RANDCODE)
#' length(unique(rcista$UID))
#' [1] 690


#' Extract time independent variables
#' ===============================================================

#' Identify 206 patients whose reason for progression was non-target lesion
NTLUID<-rcista %>% 
  filter(PRGCAUSE %in% c("2", "3", "7")) %>%
  group_by(UID) %>%
  slice(1)

resp0<-rcista %>%
#  filter(UID %nin% NTLUID$UID) %>% # Remove NTLUID patients
  group_by(UID) %>%
    filter(row_number()==n()) %>%
  ungroup() %>%
  mutate(RESP=ifelse(RESPOND=="1", 1, 0)) %>%
  select(UID, TIMETP, BESTRESP, RESP,
         DEATFLAG, OSTIM, OSCEN, 
         PRGCAUSE, BASG4LIV, BASG5CM,
         PROGDY, LSTASSDY, S_SUB, S_SUBT)

#' length(unique(resp0$UID))
#' 690
#' test1<-resp0 %>% filter(S_SUB==1)
#' test0<-resp0 %>% filter(S_SUB==0 & !is.na(S_SUBT))

#' resp0 %>% count(BESTRESP) # BOR: 0=CR, 1=PR, 3=SD, 4=PD, 9=NE
#' resp0 %>% count(RESP) # RESP: 0=Yes, 1=No
#' resp0 %>% count(BASG4LIV) # >4 liver mets at baseline
#' resp0 %>% count(BASG5CM) # Liver mets>5cm at baseline


#' Extract time dependent variables
#' ===============================================================
tk0<-rcista %>%
  select(UID, ORDYTRT, VISIT, BL_VIS, NECAT,
         BL_STLDI, STLDI, PBLCNT)

#' Build recistb
#' ===============================================================
rcistb<-left_join(resp0, subj0, by="UID") %>%
  left_join(., tk0, by="UID")
#' length(unique(rcistb$UID))

#UIDrcistb<-rcistb %>% group_by(UID) %>% slice(1) %>% ungroup()
#UIDrcistb%>% count(DIED)

#' Flag patients who switch to another treatment before end of Bev trt
switchUID<-rcistb %>%
  filter(S_SUB!=1, S_SUBT<BEV_EDY) %>%
  select(UID) %>% distinct()

keepin1<-c("UID", 
           
          #' Variables from SUBJ0          
          "RND_DY", "LSVISDY",
          "BEV_SDY", "BEV_EDY", "SEX", "WHO_STR",
          "DIED", "DEATHDY", "DEATDI", "LDH1_5", 
          "TERMPREM", "TERMREA", "BL_VEGFN", "VEGF_STR",
          
          #' Variables from RESP0
          "TIMETP", "BESTRESP", "RESP",
          "DEATFLAG", "OSTIM", "OSCEN", 
          "PRGCAUSE", "BASG4LIV", "BASG5CM",
          "PROGDY", "LSTASSDY",
          
          #' Variables from TK0          
          "ORDYTRT", "VISIT", "BL_VIS",
          "BL_STLDI", "STLDI", "PBLCNT")


#' Prepare the analysis dataset:
#' - remove rows where SLD is NA
#' - keep rows where NECAT is NA i.e. TA not impacted by irradiation or resection
#' - keep rows from patients exposed to BEV for at least 3 weeks (i.e. at least 1st TA visit)
#' - remove rows from patients who have switched to alternative therapy before end of Bev trt
#' - remove rows when TA date is after end of Bev trt

HorizTKOS<-rcistb %>%
  drop_na(STLDI) %>%
  filter(is.na(NECAT), 
         BEV_EDY>3*7,
         UID %nin% switchUID$UID) %>%
  select(one_of(keepin1)) %>%
  #' Flag patients with last TA visit 28 dats after end of Bev trt
  mutate(SLD=ifelse(STLDI==0, 2.5, STLDI*10), 
         BSLD=BL_STLDI*10,
         TYEAR=ORDYTRT/365.25) %>%
  filter(TYEAR<(BEV_EDY/365.25))

dim(HorizTKOS)
length(unique(HorizTKOS$UID))

#'summary(HorizTKOS)
#'UID234<-HorizTKOS %>% filter(UID==234)


#' ===============================================
#' Display BSLD istribution
#' ===============================================
xbreaks<-c(1, 3, 10, 30, 100, 300)
bslddf<-HorizTKOS %>%
  filter(BL_VIS=="1") %>%
  mutate(BSLDmed=median(BSLD),
         BSLDq10=quantile(BSLD, probs=.1),
         BSLDq90=quantile(BSLD, probs=.9))
  
ggplot(bslddf, aes(x=BSLD))+
  geom_density(aes(y=..count..))+
  geom_vline(aes(xintercept=BSLDmed), lty=2)+
  geom_vline(aes(xintercept=BSLDq10), lty=2)+
  geom_vline(aes(xintercept=BSLDq90), lty=2)+
  scale_x_continuous("Baseline SLD (mm)", trans="log", breaks=xbreaks)+
  theme_minimal()+
  theme(panel.grid.minor=element_blank())


#' ===============================================
#' Display SLD spaghetti in resp/non-resp patients
#' ===============================================
ybreaks<-c(0, 100, 200, 400)
g0<-ggplot(HorizTKOS, aes(x=TYEAR, y=SLD))+
  geom_line(aes(group=UID), colour="wheat4", alpha=0.2)+
  geom_point(colour="grey33", alpha=0.3, size=0.9)+
  #  facet_wrap(~RESP)+
  scale_x_continuous("Year", breaks=0.5*(0:5))+
  scale_y_continuous("SLD (mm)", breaks=ybreaks)+
  theme_minimal()+
  theme(panel.grid.minor=element_blank())
ggsave(".\\M3-SEGTTP\\Fig001 SLDspagh.png", g0, w=8, h=6)


#' ===============================================
#' Check sequence of key time points
#' ===============================================
set.seed(1708)
idtrain<-sample(HorizTKOS$UID, 70, replace=F)

gdf<-HorizTKOS %>% filter(UID %in% idtrain) %>%
  mutate(DURBEV=(BEV_EDY-BEV_SDY)/365.25, CUID=reorder(UID, DURBEV))

ggplot(gdf, aes(CUID, TYEAR))+
  geom_segment(aes(xend=CUID, yend=BEV_SDY/365.25))+
  geom_segment(aes(xend=CUID, y=BEV_SDY/365.25, yend=BEV_EDY/365.25), colour="orange1")+
  geom_point(colour="steelblue1")+
  geom_point(aes(y=TIMETP/365.25), size=2, colour="orange3")+
  geom_point(aes(y=DEATHDY/365.25), pch=3, colour="steelblue3")+
  scale_y_continuous("Year", breaks=0.5*(0:5))+
  xlab("Patient")+
  coord_flip()+
  theme_minimal()+
  theme(panel.grid.minor=element_blank(),
        axis.text.y=element_blank())


#' ===============================================
#' Individual profiles
#' ===============================================
#' Note: facet_wrap_paginate() is from ggforce
#' -------------------------------------------
#' Retain 30 patients only
pdf(".\\M3-SEGTTP\\Fig000.pdf", w=12, h=8)
for(i in 1:69){
  print(ggplot(HorizTKOS, aes(TYEAR, SLD))+
      geom_vline(aes(xintercept=TIMETP/365.25), alpha=.4, lwd=1, colour="grey70")+
      geom_vline(aes(xintercept=PROGDY/365.25), lty=2, alpha=.4, lwd=1, colour="darkblue")+
      geom_vline(aes(xintercept=DEATHDY/365.25), lwd=2, colour="black")+
      geom_point(alpha=.6, size=3)+
      facet_wrap_paginate(~UID, ncol=3, nrow=3, page=i)+
      scale_color_viridis_d(guide=F)+
      theme_minimal())
}
dev.off()


#' ===============================================
#' Saving
#' ===============================================

#saveRDS(HorizTKOS, file=".\\M3-SEGTTP\\HorizTKOS.rds")



