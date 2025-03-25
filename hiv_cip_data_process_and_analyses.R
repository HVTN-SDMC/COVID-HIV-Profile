library(dplyr)
library(tidyr)
#library(e1071)
library(Hmisc)
library(haven) 
library(broom)
library(mice)  
library(corrplot)
library(cowplot)
library(ggplot2)
library(ggbiplot)
library(factoextra)

###########
## Begin: data process
###########

cov.adj.resp <- function(data, covar, markers)
{
  covar <- paste(covar, collapse = " + ")
  out.dat <- data
  for (y in markers)
  {
    mod.formula <- paste0(y, " ~ ", covar)
    y.obs <- data[[y]]
    fit <- glm(mod.formula, data=data)
    alpha <- summary(fit)$coef[1]
    yh <- y.obs-(predict(fit)-rep(alpha,length(y.obs)))
    out.dat[[y]] <- yh
  }
  return(out.dat)
}

# load the cohort data 
coh.data <- read.csv("../adata/Demographic.csv") %>% dplyr::rename(severity=cov.sevpool,bmi=bmibl)
coh.data$severity <- factor(coh.data$severity,levels=c("Asymptomatic","Symptomatic","Hospitalized"))
coh.data$sex <- factor(coh.data$sex, levels = c("Male","Female"))
coh.data$region <- factor(coh.data$region, levels = c("Peru","USA"))

# load the assays data
immu.data <- read.csv("../adata/COV2_markers_original.csv") %>% dplyr::select(-visitno)
resp.vars <- names(immu.data)[grepl("resp",names(immu.data))]
## de-select response variables
immu.data <- immu.data %>% dplyr::select(-all_of(resp.vars))
## list of all response magnitude variables 
# resp.mag <- names(immu.data)[(grepl("BAMA",names(immu.data)) | grepl("ICS",names(immu.data)) 
#                           | grepl("BCP",names(immu.data)) | grepl("ADCC",names(immu.data)) | grepl("SECABA",names(immu.data))
#                           | grepl("ADCP",names(immu.data)) ) & !grepl("resp",names(immu.data))]
## remove 
# old.names <- c("BAMA.IgG1.RBD","BAMA.IgG1.NP","BAMA.IgG1.NTD","BAMA.IgG1.S6P","BAMA.IgG1.S2P",
#                "BAMA.IgG3.RBD","BAMA.IgG3.NP","BAMA.IgG3.NTD","BAMA.IgG3.S6P","BAMA.IgG3.S2P",
#                "BAMA.IgA.RBD","BAMA.IgA.NP","BAMA.IgA.NTD","BAMA.IgA.S6P","BAMA.IgA.S2P",
#                "BCP.S.IgG.Tot","BCP.S.IgA.Tot","BCP.S.IgM.Tot",
#                "BCP.S.IgG.Mem","BCP.S.IgA.Mem","BCP.S.IgM.Mem",
#                "ICS.S.CD4","ICS.N.CD4","ICS.S.CD8","ICS.N.CD8")
# new.names <- c("IgG1.RBD","IgG1.NP","IgG1.NTD","IgG1.S6P","IgG1.S2P",
#                "IgG3.RBD","IgG3.NP","IgG3.NTD","IgG3.S6P","IgG3.S2P",
#                "IgA.RBD","IgA.NP","IgA.NTD","IgA.S6P","IgA.S2P",
#                "S.IgG.Tot","S.IgA.Tot","S.IgM.Tot",
#                "S.IgG.Mem","S.IgA.Mem","S.IgM.Mem",
#                "S.CD4","N.CD4","S.CD8","N.CD8")
# immu.data <- immu.data %>% rename_with(all_of(old.names), .fn = ~ all_of(new.names))
mag.vars <- names(immu.data)[-1]

## log-transformation for all magnitudes except SECABA, iADCC, tADCC, AEC2
trans.vars <- mag.vars[-which(mag.vars %in% c("SECABA","iADCC","tADCC","ACE2"))]

log_trans <- function(y) y=log(y)
immu.data <- immu.data %>% mutate_at(all_of(trans.vars),log_trans)

## check correlation
cor <- cor(immu.data[,-1],method="spearman",use="pairwise.complete.obs")
## three pairs "MSD.IgG.RBD" and "MSD.IgG.S2P" (r=0.99), "VSV.ID50" and "VSV.ID80" (r=0.99),  "BAMA.IgA.S2P" and "BAMA.IgA.S6P" (r=0.95), 
## had spearman correlation >= 0.95
## so keep "MSD.IgG.RBD", "VSV.ID50", and "BAMA.IgA.S6P", one marker from each pair
immu.data <- immu.data %>% dplyr::select(-MSD.IgG.S2P, -VSV.ID80, -BAMA.IgA.S2P)
mag.vars <- names(immu.data)[-1]

adata <- inner_join(coh.data,immu.data,by="ptid") 
#check missing data skim(adata)
## impute missings:  one missing IgG1.NP, IgG1.S6P, 1 missing S.CD4/S.CD8, 2 N.CD4/N.CD8, and one for all 6 BCP markers, 3 for MSD.RBD
## 3 missing SECABA/ADCC, 20 missing ADCP

ptid <- adata %>% dplyr::select(ptid)

## impute missing data in all markers except NAB data
imp.vars <- names(immu.data)[-1] ## [!grepl("NAB",names(immu.data))][-1]
x.dat <- adata %>% dplyr::select(all_of(imp.vars),shedder,hiv,severity,age,sex,region,covstdy,diab,hyper,cea,bmi,smokever,smokenow)
row.names(x.dat) <- NULL
set.seed(13452)
imp.dat <- mice(x.dat,m=1,mehod='lasso.select.norm')
imp.dat <- complete(imp.dat,action = 1)

## combine 
ds <- bind_cols(ptid,imp.dat)

## output data with missing imputed and scaled
ds.out <- ds %>% mutate_at(mag.vars, scale)
## checking the correlations between markers,  the correlation between MSD.RBD and MSD.S2P is 0.95.
## MSD.S2P is excluded. 

## immune markers covariate-adjusted for days since SARS-CoV-2 diagnosis
covar <- c("covstdy")
out.dat <- cov.adj.resp(ds,covar,mag.vars) 
## normalize all adjusted markers to mean of zero and std of one
out.dat <- out.dat %>% mutate_at(mag.vars,scale)
write.csv(out.dat,"../adata/COV2_markers_normalized.csv",row.names = FALSE)

# perform imputation and covariate-adjustment for the endemic measures
immu.data <- read.csv("../adata/Endemic_markers_original.csv")
## perform imputation
set.seed(13452)
imp.dat <- mice(immu.data[,-1],m=1,mehod='lasso.select.norm')
imp.dat <- complete(imp.dat,action = 1)
## combine with ptid
ds <- bind_cols(immu.data[,1],imp.dat)
names(ds)[1] <- "ptid"
## adjust covariate day since diagnosis
mag.vars <- names(ds)[-1]
dat <- left_join(ds,coh.data,by="ptid")
covar <- c("covstdy")
out.dat <- cov.adj.resp(dat,covar,mag.vars)
## select ptid and endemic measures
out.dat <- out.dat %>% dplyr::select(ptid,all_of(mag.vars))
names(out.dat)[-1] <- paste0("End.",names(out.dat)[-1])
write.csv(out.dat, "../adata/Endemic_markers_normalized.csv",row.names = FALSE)
###########
## End: data process
###########

###############
#### Begin: polar-plots and correlation plots
################
my.polar.plot <- function(data,markers.dat, filename)
{
  markers <- markers.dat$variable_name
  n.markers <- length(markers)
  ## polar-plots  
  ## calculate percentile for each immune markers
  percent <- function(x) rank(x)/length(x)
  ds <- data %>% mutate_at(markers,percent)
  ## calcuate mean percentiles by PVS for each marker
  p.wide <- ds %>% group_by(hiv) %>%
    summarise_at(all_of(markers), mean)
  ## change the mean data from a wide form to a long form 
  p.long <- gather(p.wide,variable_name,Mean,all_of(markers))
  p.long <- left_join(p.long,markers.dat,by="variable_name")
  
  p.long$Feature <- factor(p.long$variable_name,levels=markers,labels=1:n.markers)
  p.long$type <- factor(p.long$type,levels=c("S-Cov2 IgA","S-Cov2 IgG1","S-Cov2 IgG3","Endemic IgA","Endemic IgG1",
                                             "Endemic IgG3","SECABA","MSD IgG","Neutralization","Functional","Total B Cells",
                                             "IgA+ B Cells","IgG+ B Cells","IgM+ B Cells","CD4+ T Cells","CD8+ T Cells"))
  type_color <-levels(p.long$code)
  
  x.label <- c("1. IgA to S-CoV2 RBD\n2. IgA to S-Cov2 NP\n3. IgA to S-CoV2 NTD\n4. IgA to S-CoV2 6P Spike",
               "5. IgG1 to S-CoV2 RBD\n6. IgG1 to S-CoV2 NP\n7. IgG1 to S-CoV2 NTD\n8. IgG1 to S-CoV2 2P Spike\n9. IgG1 to S-CoV2 6P Spike",
               "10. IgG3 to S-CoV2 RBD\n11. IgG3 to S-CoV2 NP\n12. IgG3 to S-CoV2 NTD\n13. IgG3 to S-CoV2 2P Spike\n14. IgG3 to S-CoV2 6P Spike",
               "15. IgA to Endemic 229E-RBD\n16. IgA to Endemic HKU1-RBD\n17. IgA to Endemic NL63-RBD\n18. IgA to Endemic OC43-RBD",
               "19. IgG1 to Endemic 229E-RBD\n20. IgG1 to Endemic HKU1-RBD\n21. IgG1 to Endemic NL63-RBD\n22. IgG1 to Endemic OC43-RBD",
               "23. IgG3 to Endemic 229E-RBD\n24. IgG3 to Endemic HKU1-RBD\n25. IgG3 to Endemic NL63-RBD\n26. IgG3 to Endemic OC43-RBD",
               "27. Spike-expressing cell antibody binding",
               "28. IgG to RBD (MSD)\n29. IgG to NP (MSD)",
               "30. Pseudovirus neutralization\n31. ACE-2 Blocking",
               "32. ADCP\n33. Infected cell ADCC\n34. Transfected cell ADCC",
               "35. %Spike+ of Total B Cells\n36. %RBD+ of Total B Cells",
               "37. %Spike+ IgA+ of Total B Cells\n38. %S+ IgA+ of Memory B Cells\n39. %RBD+ IgA+ of Memory B Cells",
               "40. %Spike+ IgG+ of Total B Cells\n41. %Spike+ IgG+ of Memory B Cells\n42. %RBD+ IgG+ of Memory B Cells",
               "43. %Spike+ IgM+ of Total B Cells\n44. %Spike+ IgM+ of Memory B Cells\n45. %RBD+ IgM+ of Memory B Cells\n46. %Spike+ IgM+ IgD+ of Total B Cells",
               "47. CD4+ T Cells to Spike\n48. CD4+ T Cells to NP\n49. CD4+ T Cells to E&M\n50. CD4+ T Cells to ORF3a6\n51. CD4+ T Cells to ORF7a7b8",
               "52. CD8+ T Cells to Spike\n53. CD8+ T Cells to N\n54. CD8+ T Cells to E&M\n55. CD8+ T Cells to ORF3a6\n56. CD8+ T Cells to ORF7a7b8")
  
  
  angles <- 90-360 * (1:n.markers-0.5)/n.markers
  angles <- ifelse(angles < -90, angles+180, angles)
  hjust <- ifelse( angles < -90, 1, 0)
  
  cxc <- ggplot(p.long,aes(x=Feature,y=Mean,fill=type))+ 
    geom_hline(yintercept = seq(0,0.8,by=0.1), colour = "grey100", size = 0.2) +  
    # geom_vline(xintercept ='Nucleoprotein', data=rr,colour = "grey90")+
    # geom_histogram( position=position_dodge(0.7), stat="identity", color = "white",width = 0.7)+
    geom_bar(stat="identity",position = position_dodge2(preserve = "single"))+
    theme_light()+
    #facet_grid(~shedder, switch = "y")+
    facet_grid(~hiv)+
    theme_minimal() +
    #scale_y_continuous(breaks=seq(0,0.8,by=0.1),labels=seq(0,0.8,by=0.1))+
    #geom_hline(yintercept = 0.8, color = "black")+
    geom_hline(yintercept = 0.5, color = "red")+
    geom_hline(yintercept = 0.1, color = "grey")+
    geom_hline(yintercept = 0.2, color = "grey")+
    geom_hline(yintercept = 0.3, color = "grey")+
    geom_hline(yintercept = 0.4, color = "grey")+
    geom_hline(yintercept = 0.6, color = "grey")+
    geom_hline(yintercept = 0.7, color = "grey")+
    geom_hline(yintercept = 0.8, color = "grey")+
    # Remove legend, axes, text, and tick marks
    # scale_fill_manual(name = 'Markers',label=x.label)+
    # guides(fill=guide_legend(ncol=2))+
    scale_fill_manual(labels=x.label,values=type_color)+
    theme(
      strip.text.x = element_text(size = 12,face="bold"),
      legend.position = "right",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10,face="bold", angle=angles,hjust = hjust),
      #  axis.text.y = element_text(face="bold"),
      # axis.text.x = element_text(color="black", size=6, angle= myAng),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 14, face = "bold"),
      #   panel.border = element_blank(),
      # panel.grid  = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_rect(
        color="white", fill="white", size=1),strip.text = element_text(size=8,face='bold'),
      strip.placement = "outside"
    )+ 
    #  guides(color = guide_legend(ncol = 2))+
    geom_text(x=0.5,y=0.1,label="0.1",size=2,color="grey50")+
    geom_text(x=0.5,y=0.2,label="0.2",size=2,color="grey50")+
    geom_text(x=0.5,y=0.3,label="0.3",size=2,color="grey50")+
    geom_text(x=0.5,y=0.4,label="0.4",size=2,color="grey50")+
    geom_text(x=0.5,y=0.5,label="0.5",size=2,color="grey50")+
    geom_text(x=0.5,y=0.6,label="0.6",size=2,color="grey50")+
    geom_text(x=0.5,y=0.7,label="0.7",size=2,color="grey50")+
    geom_text(x=0.5,y=0.8,label="0.8",size=2,color="grey50")+
  labs(fill="Markers",x = "",
       y = 'Mean Percentile ',
       title = 'Polar-plot of Mean Percentile of Immune Markers by HIV Status')+
    theme(plot.title=element_text(hjust=0.5))
  cxc + coord_polar(start=0) 
  
  out.name <- paste0("../graph/hiv_polar_plot_",filename)
  ggsave(paste0(out.name,".pdf"),height = 8,width=12)
  # ggsave(paste0(out.name,".tiff"),height = 8,width=12)
}

my.polar.by_severity_plot <- function(data,markers.dat,filename)
{
  markers <- markers.dat$variable_name
  n.markers <- length(markers)
  ## polar-plots  
  ## calculate percentile for each immune markers
  percent <- function(x) rank(x)/length(x)
  ds <- data %>% mutate_at(markers,percent)
  ## calcuate mean percentiles by PVS for each marker
  p.wide <- ds %>% group_by(severity,hiv) %>%
    summarise_at(all_of(markers), mean)
  ## change the mean data from a wide form to a long form 
  p.long <- gather(p.wide,variable_name,Mean,all_of(markers))
  p.long <- left_join(p.long,markers.dat,by="variable_name")
  type_color <- levels(p.long$code)
  
  p.long$Feature <- factor(p.long$variable_name,levels=markers,labels=1:n.markers)
  
  x.label <- c("1. IgA to S-CoV2 RBD\n2. IgA to S-Cov2 NP\n3. IgA to S-CoV2 NTD\n4. IgA to S-CoV2 6P Spike",
               "5. IgG1 to S-CoV2 RBD\n6. IgG1 to S-CoV2 NP\n7. IgG1 to S-CoV2 NTD\n8. IgG1 to S-CoV2 2P Spike\n9. IgG1 to S-CoV2 6P Spike",
               "10. IgG3 to S-CoV2 RBD\n11. IgG3 to S-CoV2 NP\n12. IgG3 to S-CoV2 NTD\n13. IgG3 to S-CoV2 2P Spike\n14. IgG3 to S-CoV2 6P Spike",
               "15. IgA to Endemic 229E-RBD\n16. IgA to Endemic HKU1-RBD\n17. IgA to Endemic NL63-RBD\n18. IgA to Endemic OC43-RBD",
               "19. IgG1 to Endemic 229E-RBD\n20. IgG1 to Endemic HKU1-RBD\n21. IgG1 to Endemic NL63-RBD\n22. IgG1 to Endemic OC43-RBD",
               "23. IgG3 to Endemic 229E-RBD\n24. IgG3 to Endemic HKU1-RBD\n25. IgG3 to Endemic NL63-RBD\n26. IgG3 to Endemic OC43-RBD",
               "27. Spike-expressing cell antibody binding",
               "28. IgG to RBD (MSD)\n29. IgG to NP (MSD)",
               "30. Pseudovirus neutralization\n31. ACE-2 Blocking",
               "32. ADCP\n33. Infected cell ADCC\n34. Transfected cell ADCC",
               "35. %Spike+ of Total B Cells\n36. %RBD+ of Total B Cells",
               "37. %Spike+ IgA+ of Total B Cells\n38. %S+ IgA+ of Memory B Cells\n39. %RBD+ IgA+ of Memory B Cells",
               "40. %Spike+ IgG+ of Total B Cells\n41. %Spike+ IgG+ of Memory B Cells\n42. %RBD+ IgG+ of Memory B Cells",
               "43. %Spike+ IgM+ of Total B Cells\n44. %Spike+ IgM+ of Memory B Cells\n45. %RBD+ IgM+ of Memory B Cells\n46. %Spike+ IgM+ IgD+ of Total B Cells",
               "47. CD4+ T Cells to Spike\n48. CD4+ T Cells to NP\n49. CD4+ T Cells to E&M\n50. CD4+ T Cells to ORF3a6\n51. CD4+ T Cells to ORF7a7b8",
               "52. CD8+ T Cells to Spike\n53. CD8+ T Cells to N\n54. CD8+ T Cells to E&M\n55. CD8+ T Cells to ORF3a6\n56. CD8+ T Cells to ORF7a7b8")
  
  angles <- 90-360 * (1:n.markers-0.5)/n.markers
  angles <- ifelse(angles < -90, angles+180, angles)
  hjust <- ifelse( angles < -90, 1, 0)
  
  cxc <- ggplot(p.long,aes(x=Feature,y=Mean,fill=type))+ 
    geom_hline(yintercept = seq(0,0.8,by=0.1), colour = "grey100", size = 0.2) +  
    # geom_vline(xintercept ='Nucleoprotein', data=rr,colour = "grey90")+
    # geom_histogram( position=position_dodge(0.7), stat="identity", color = "white",width = 0.7)+
    geom_bar(stat="identity",position = position_dodge2(preserve = "single"))+
    #theme_light()+
    #facet_grid(~shedder, switch = "y")+
    facet_grid(rows=vars(severity),cols=vars(hiv), 
               switch = "y")+
    theme(strip.text.y.left = element_text(angle = 0))+
    theme_minimal() +
    scale_y_continuous(breaks=seq(0,0.8,by=0.1),
                       labels=seq(0,0.8,by=0.1))+
    #geom_hline(yintercept = 0.8, color = "black")+
    geom_hline(yintercept = 0.5, color = "red")+
    geom_hline(yintercept = 0.1, color = "grey")+
    geom_hline(yintercept = 0.2, color = "grey")+
    geom_hline(yintercept = 0.3, color = "grey")+
    geom_hline(yintercept = 0.4, color = "grey")+
    geom_hline(yintercept = 0.6, color = "grey")+
    geom_hline(yintercept = 0.7, color = "grey")+
    geom_hline(yintercept = 0.8, color = "grey")+
    # Remove legend, axes, text, and tick marks
    # scale_fill_manual(name = 'Markers',label=x.label)+
    # guides(fill=guide_legend(ncol=2))+
    scale_fill_manual(labels=x.label,values=type_color)+
    theme(
      strip.text.x = element_text(size = 12,face="bold"),
      legend.position = "right",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=10,face="bold", angle=angles,hjust = hjust),
      #  axis.text.y = element_text(face="bold"),
      # axis.text.x = element_text(color="black", size=6, angle= myAng),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 14, face = "bold"),
      #   panel.border = element_blank(),
      # panel.grid  = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_rect(
        color="white", fill="white", size=1),strip.text = element_text(size=8,face='bold'),
      strip.placement = "outside"
    )+ 
    #  guides(color = guide_legend(ncol = 2))+
    geom_text(x=0.5,y=0.1,label="0.1",size=2,color="grey50")+
    geom_text(x=0.5,y=0.2,label="0.2",size=2,color="grey50")+
    geom_text(x=0.5,y=0.3,label="0.3",size=2,color="grey50")+
    geom_text(x=0.5,y=0.4,label="0.4",size=2,color="grey50")+
    geom_text(x=0.5,y=0.5,label="0.5",size=2,color="grey50")+
    geom_text(x=0.5,y=0.6,label="0.6",size=2,color="grey50")+
    geom_text(x=0.5,y=0.7,label="0.7",size=2,color="grey50")+
    geom_text(x=0.5,y=0.8,label="0.8",size=2,color="grey50")+
  labs(fill="Markers",
       x = "",
       y = 'Mean Percentile ',
       title = 'Polar-plot of Mean Percentile of Immune Markers by HIV Status stratified by COVID-19 Severity')+
    theme(plot.title=element_text(hjust=0.5))
  cxc + coord_polar(start=0) 
  
  out.name <- paste0("../graph/hiv_polar_plot_",filename)
  ggsave(paste0(out.name,".pdf"),height = 10,width=16)
  # ggsave(paste0(out.name,".tiff"),height = 8,width=12)
}

my.corplot <- function(data,resp.mag,lab)
{
  n.markers <- length(resp.mag)
  ds <- data %>% dplyr::select(all_of(resp.mag)) 
  ## correlation among all ppts
  rr <- cor(ds,method="spearman",use="pairwise.complete.obs")
  test.res <- cor.mtest(ds, method="spearman",exact=TRUE,conf.level = 0.95)
  p <- test.res$p 
  # rr0 <- round(rr0,2)
  # p.char <- ifelse(p<0.001,"<0.001", round(p,3))
  out.dat <- rr
  for (i in 1:n.markers)
  {
    for (j in 1:(i-1))
    {
      out.dat[i,j] <- p[i,j]
    }
    out.dat[i,i] <- ""
  }

  ## correlation among PWOH
  dat.0 <- data %>% filter(hiv=="PWOH") %>% dplyr::select(all_of(resp.mag))
  rr0 <- cor(dat.0,method="spearman",use="pairwise.complete.obs")
  test.res <- cor.mtest(dat.0, method="spearman",exact=TRUE,conf.level = 0.95)
  p <- test.res$p 
  # rr0 <- round(rr0,2)
  # p.char <- ifelse(p<0.001,"<0.001", round(p,3))
  out.dat <- rr0
  for (i in 1:n.markers)
  {
    for (j in 1:(i-1))
    {
      out.dat[i,j] <- p[i,j]
    }
    out.dat[i,i] <- ""
  }
  pdf(paste0("../graph/corrplots.PWOH_All_",lab,".pdf"),width = 12,height=10)
  corrplot.mixed(rr0,p.mat=p,sig.level = c(0.001,0.01,0.05),pch.cex = 0.5,
                 insig = 'label_sig', pch.col = 'white',order='alphabet',
                 number.cex=0.5,tl.cex=0.8,tl.pos="lt",tl.col="black")
  dev.off()
  
  ## correlation among PLWH
  dat.1 <- data %>% filter(hiv=="PLWH") %>% dplyr::select(all_of(resp.mag))
  rr1 <- cor(dat.1,method="spearman",use="pairwise.complete.obs")
  test.res <- cor.mtest(dat.1, method="spearman",exact=TRUE,conf.level = 0.95)
  p <- test.res$p 
  n.markers=dim(p)[1]
  # rr1 <- round(rr1,2)
  # p.char <- ifelse(p<0.001,"<0.001", round(p,3))
  out.dat <- rr1
  for (i in 1:n.markers)
  {
    for (j in 1:(i-1))
    {
      out.dat[i,j] <- p[i,j]
    }
    out.dat[i,i] <- ""
  }
  pdf(paste0("../graph/corrplots.PLWH_All_",lab,".pdf"),width = 12,height=10)
  corrplot.mixed(rr1,p.mat=p,sig.level = c(0.001,0.01,0.05),pch.cex = 0.5,
                 insig = 'label_sig', pch.col = 'white',order='alphabet',
                 number.cex=0.5,tl.cex=0.8,tl.pos="lt",tl.col="black")
  dev.off()
  
  ## correlation by severity
  ## among asymptomatic
  dat.0 <- data %>% filter(severity=="Asymptomatic") %>% dplyr::select(all_of(resp.mag))
  rr0 <- cor(dat.0,method="spearman",use="pairwise.complete.obs")
  test.res <- cor.mtest(dat.0, method="spearman",exact=TRUE,conf.level = 0.95)
  p <- test.res$p 
  n.markers=dim(p)[1]
  # rr0 <- round(rr0,2)
  # p.char <- ifelse(p<0.001,"<0.001", round(p,3))
  out.dat <- rr0
  for (i in 1:n.markers)
  {
    for (j in 1:(i-1))
    {
      out.dat[i,j] <- p[i,j]
    }
    out.dat[i,i] <- ""
  }

  ## correlation by severity
  ## among symptomatic
  dat.0 <- data %>% filter(severity=="Symptomatic") %>% dplyr::select(all_of(resp.mag))
  rr0 <- cor(dat.0,method="spearman",use="pairwise.complete.obs")
  test.res <- cor.mtest(dat.0, method="spearman",exact=TRUE,conf.level = 0.95)
  p <- test.res$p 
  n.markers=dim(p)[1]
  # rr0 <- round(rr0,2)
  # p.char <- ifelse(p<0.001,"<0.001", round(p,3))
  out.dat <- rr0
  for (i in 1:n.markers)
  {
    for (j in 1:(i-1))
    {
      out.dat[i,j] <- p[i,j]
    }
    out.dat[i,i] <- ""
  }

  ## correlation by severity
  ## among hospitalized
  dat.0 <- data %>% filter(severity=="Hospitalized") %>% dplyr::select(all_of(resp.mag))
  rr0 <- cor(dat.0,method="spearman",use="pairwise.complete.obs")
  test.res <- cor.mtest(dat.0, method="spearman",exact=TRUE,conf.level = 0.95)
  p <- test.res$p 
  n.markers=dim(p)[1]
  # rr0 <- round(rr0,2)
  # p.char <- ifelse(p<0.001,"<0.001", round(p,3))
  out.dat <- rr0
  for (i in 1:n.markers)
  {
    for (j in 1:(i-1))
    {
      out.dat[i,j] <- p[i,j]
    }
    out.dat[i,i] <- ""
  }

  ## correlation among PWOH by severity
  dat.0 <- data %>% filter(hiv=="PWOH" & severity=="Asymptomatic") %>% dplyr::select(all_of(resp.mag))
  rr0 <- cor(dat.0,method="spearman",use="pairwise.complete.obs")
  test.res <- cor.mtest(dat.0, method="spearman",exact=TRUE,conf.level = 0.95)
  p <- test.res$p 
  n.markers=dim(p)[1]
  # rr0 <- round(rr0,2)
  # p.char <- ifelse(p<0.001,"<0.001", round(p,3))
  out.dat <- rr0
  for (i in 1:n.markers)
  {
    for (j in 1:(i-1))
    {
      out.dat[i,j] <- p[i,j]
    }
    out.dat[i,i] <- ""
  }
  pdf(paste0("../graph/corrplots.PWOH_Asym_",lab,".pdf"),width = 12,height=10)
  corrplot.mixed(rr0,p.mat=p,sig.level = c(0.001,0.01,0.05),pch.cex = 0.5,
                 insig = 'label_sig', pch.col = 'white',order='alphabet',
                 number.cex=0.5,tl.cex=0.8,tl.pos="lt",tl.col="black")
  dev.off()
  
  
  dat.0 <- data %>% filter(hiv=="PWOH" & severity=="Symptomatic") %>% dplyr::select(all_of(resp.mag))
  rr0 <- cor(dat.0,method="spearman",use="pairwise.complete.obs")
  test.res <- cor.mtest(dat.0, method="spearman",exact=TRUE,conf.level = 0.95)
  p <- test.res$p 
  n.markers=dim(p)[1]
  # rr0 <- round(rr0,2)
  # p.char <- ifelse(p<0.001,"<0.001", round(p,3))
  out.dat <- rr0
  for (i in 1:n.markers)
  {
    for (j in 1:(i-1))
    {
      out.dat[i,j] <- p[i,j]
    }
    out.dat[i,i] <- ""
  }
  pdf(paste0("../graph/corrplots.PWOH_Sym_",lab,".pdf"),width = 12,height=10)
  corrplot.mixed(rr0,p.mat=p,sig.level = c(0.001,0.01,0.05),pch.cex = 0.5,
                 insig = 'label_sig', pch.col = 'white',order='alphabet',
                 number.cex=0.5,tl.cex=0.8,tl.pos="lt",tl.col="black")
  dev.off()
  
  dat.0 <- data %>% filter(hiv=="PWOH" & severity=="Hospitalized") %>% dplyr::select(all_of(resp.mag))
  rr0 <- cor(dat.0,method="spearman",use="pairwise.complete.obs")
  test.res <- cor.mtest(dat.0, method="spearman",exact=TRUE,conf.level = 0.95)
  p <- test.res$p 
  n.markers=dim(p)[1]
  # rr0 <- round(rr0,2)
  # p.char <- ifelse(p<0.001,"<0.001", round(p,3))
  out.dat <- rr0
  for (i in 1:n.markers)
  {
    for (j in 1:(i-1))
    {
      out.dat[i,j] <- p[i,j]
    }
    out.dat[i,i] <- ""
  }
  pdf(paste0("../graph/corrplots.PWOH_Hos_",lab,".pdf"),width = 12,height=10)
  corrplot.mixed(rr0,p.mat=p,sig.level = c(0.001,0.01,0.05),pch.cex = 0.5,
                 insig = 'label_sig', pch.col = 'white',order='alphabet',
                 number.cex=0.5,tl.cex=0.8,tl.pos="lt",tl.col="black")
  dev.off()
  
  ## correlation among PLWH by severity
  dat.1 <- data %>% filter(hiv=="PLWH" & severity=="Asymptomatic") %>% dplyr::select(all_of(resp.mag))
  rr1 <- cor(dat.1,method="spearman",use="pairwise.complete.obs")
  test.res <- cor.mtest(dat.1, method="spearman",exact=TRUE,conf.level = 0.95)
  p <- test.res$p 
  n.markers=dim(p)[1]
  # rr1 <- round(rr1,2)
  # p.char <- ifelse(p<0.001,"<0.001", round(p,3))
  out.dat <- rr1
  for (i in 1:n.markers)
  {
    for (j in 1:(i-1))
    {
      out.dat[i,j] <- p[i,j]
    }
    out.dat[i,i] <- ""
  }
  pdf(paste0("../graph/corrplots.PLWH_Asym_",lab,".pdf"),width = 12,height=10)
  corrplot.mixed(rr1,p.mat=p,sig.level = c(0.001,0.01,0.05),pch.cex = 0.6,
                 insig = 'label_sig', pch.col = 'white',order='alphabet',
                 number.cex=0.5,tl.cex=0.8,tl.pos="lt",tl.col="black")
  dev.off()
  
  dat.1 <- data %>% filter(hiv=="PLWH" & severity=="Symptomatic") %>% dplyr::select(all_of(resp.mag))
  rr1 <- cor(dat.1,method="spearman",use="pairwise.complete.obs")
  test.res <- cor.mtest(dat.1, method="spearman",exact=TRUE,conf.level = 0.95)
  p <- test.res$p 
  n.markers=dim(p)[1]
  # rr1 <- round(rr1,2)
  # p.char <- ifelse(p<0.001,"<0.001", round(p,3))
  out.dat <- rr1
  for (i in 1:n.markers)
  {
    for (j in 1:(i-1))
    {
      out.dat[i,j] <- p[i,j]
    }
    out.dat[i,i] <- ""
  }
  
  pdf(paste0("../graph/corrplots.PLWH_Sym_",lab,".pdf"),width = 12,height=10)
  corrplot.mixed(rr1,p.mat=p,sig.level = c(0.001,0.01,0.05),pch.cex = 0.6,
                 insig = 'label_sig', pch.col = 'white',order='alphabet',
                 number.cex=0.5,tl.cex=0.8,tl.pos="lt",tl.col="black")
  dev.off()
  
  dat.1 <- data %>% filter(hiv=="PLWH" & severity=="Hospitalized") %>% dplyr::select(all_of(resp.mag))
  rr1 <- cor(dat.1,method="spearman",use="pairwise.complete.obs")
  test.res <- cor.mtest(dat.1, method="spearman",exact=TRUE,conf.level = 0.95)
  p <- test.res$p 
  n.markers=dim(p)[1]
  # rr1 <- round(rr1,2)
  # p.char <- ifelse(p<0.001,"<0.001", round(p,3))
  out.dat <- rr1
  for (i in 1:n.markers)
  {
    for (j in 1:(i-1))
    {
      out.dat[i,j] <- p[i,j]
    }
    out.dat[i,i] <- ""
  }
  
  pdf(paste0("../graph/corrplots.PLWH_Hos_",lab,".pdf"),width = 12,height=10)
  corrplot.mixed(rr1,p.mat=p,sig.level = c(0.001,0.01,0.05),pch.cex = 0.6,
                 insig = 'label_sig', pch.col = 'white',order='alphabet',
                 number.cex=0.5,tl.cex=0.8,tl.pos="lt",tl.col="black")
  dev.off()
}

data <- read.csv("../adata/COV2_markers_normalized.csv") # %>% filter(severity != "Asymptomatic")
data$hiv <- factor(data$hiv,levels=c(0,1),labels=c("PWOH","PLWH"))
data$severity <- factor(data$severity,levels = c("Asymptomatic","Symptomatic","Hospitalized"))

end.data <- read.csv("../adata/Endemic_markers_normalized.csv")
data <- left_join(data,end.data, by="ptid")

markers.dat <- read.csv("../adata/immune_markers.csv")
markers.dat$type <- factor(markers.dat$type,levels=c("S-Cov2 IgA","S-Cov2 IgG1","S-Cov2 IgG3","Endemic IgA","Endemic IgG1",
                                                     "Endemic IgG3","SECABA","MSD IgG","Neutralization","Functional","Total B Cells",
                                                     "IgA+ B Cells","IgG+ B Cells","IgM+ B Cells","CD4+ T Cells","CD8+ T Cells"))
markers <- markers.dat$variable_name

markers.dat$code <- factor(markers.dat$code, levels=c("#CC3299","#FF00FF","#9932CD","#000080","#0000FF",
                                                      "#4682B4","#87CEEB","#00FFFF","#7FFFD4","#00FF7F","#ADFF2F",
                                                      "#FFFF00","#FFD700","#FF7F00","#FF2400","#CD0000"))
## polar plots
my.polar.plot(data,markers.dat,"_56markers")
## stratified by severity
my.polar.by_severity_plot(data,markers.dat,"_56markers_severity")


## correlation heatmap
ds <- data %>% select(hiv,severity,all_of(markers))
names(ds)[c(-1,-2)] <- paste0("m",c(paste0("0",1:9),10:length(markers)))

markers <- names(ds)[c(-1,-2)]

my.corplot(ds,markers,"_56markers")
###############
#### End: polar-plots and correlation plots
################

######
## Begin: Resampling analysis to evaluate sample size effect on network correlation analysis
#####

coordinated_measure <- function(ds)
{
  ## correlation coeff 
  rho <- cor(ds,method="spearman",use="pairwise.complete.obs")
  rho <- rho[lower.tri(rho)]
  ## p-values
  test.res <- cor.mtest(ds, method="spearman",exact=TRUE,conf.level = 0.95)
  p.val <- test.res$p 
  p.val <- p.val[lower.tri(p.val)]
  
  ## coordinated measure: # of pairs having a strong correlation such as abs(rho) >= 0.7 and q-value <= 0.1 
  q.val <- p.adjust(p.val,method="BH")
  
  res.pos.coor <- sum(rho>=0.7 & q.val<=0.1 & p.val<=0.05)
  res.neg.coor <- sum(rho<=-0.7 & q.val<=0.1 & p.val<=0.05)
  res.coor <- res.pos.coor + res.neg.coor
  
  pos.coor <- sum(rho>0 & q.val<=0.1 & p.val<=0.05)
  neg.coor <- sum(rho<0 & q.val<=0.1 & p.val<=0.05)
  coor <- pos.coor + neg.coor
  
  res.out <- data.frame(res.coor=res.coor, res.pos.coor=res.pos.coor,res.neg.coor=res.neg.coor,
                        coor=coor, pos.coor=pos.coor,neg.coor=neg.coor)
  
  ## coordinated measure: # of pairs having a non-zero correlation 
  n.markers <- dim(ds)[2]
  Q = qvalue(p=p.val)
  ## estimate proportion of null correlations
  pi0hat_0 = Q$pi0
  pi0hat = Q$pi0.lambda[Q$lambda==0.5]
  Size.1 = (1-pi0hat_0)*n.markers*(n.markers-1)/2
  Size.2 = (1-pi0hat)*n.markers*(n.markers-1)/2
  Ave_degree.1 = (1-pi0hat_0)*(n.markers-1)
  Ave_degree.2 = (1-pi0hat)*(n.markers-1)
  
  out <- data.frame(Size.1=Size.1,Size.2=Size.2,Ave_degree.1=Ave_degree.1,Ave_degree.2=Ave_degree.2)
  
  return(list(res.out,out))
}

subsampling <- function(dat,n.samples)
{
  res.out <- NULL
  out <- NULL
  n.boots <- 10000
  N <- dim(dat)[1]
  for (i in 1:n.boots)
  {
    sub <- sample(1:N)
    tmp.ds <- dat[sub[1:n.samples],]
    one.out <- coordinated_measure(tmp.ds)
    res.out <- bind_rows(res.out,one.out[1])
    out <- bind_rows(out,one.out[2])
    print(i)
  }
  return(list(res.out,out))
}

ds <- data %>% select(hiv,severity,all_of(markers))
names(ds)[c(-1,-2)] <- paste0("m",c(paste0("0",1:9),10:length(markers)))

markers <- names(ds)[c(-1,-2)]

table(ds$severity,ds$hiv)
#               PWOH PLWH
# Asymptomatic   50    9
# Symptomatic    85   16
# Hospitalized   81   18
## Symptomatic 
ds.sym.pwoh <- ds %>% filter(severity=="Symptomatic" & hiv=="PWOH") %>%
  select(all_of(markers))

out.sym <- subsampling(ds.sym.pwoh,16)
save(out.sym,file="../result/subsample_coor_Sym_PWOH.RData")

## Hospitalized 
ds.hos.pwoh <- ds %>% filter(severity=="Hospitalized" & hiv=="PWOH") %>%
  select(all_of(markers))

out.hos <- subsampling(ds.hos.pwoh,18)
save(out.hos,file="../result/subsample_coor_Hos_PWOH.RData")

## Asymptomatic 
ds.asym.pwoh <- ds %>% filter(severity=="Asymptomatic" & hiv=="PWOH") %>%
  select(all_of(markers))

out.asym <- subsampling(ds.asym.pwoh,9)
save(out.asym,file="../result/subsample_coor_Asym_PWOH.RData")

## all ppts
ds.all.pwoh <- ds %>% filter(hiv=="PWOH") %>%
  select(all_of(markers))
out.all <- subsampling(ds.all.pwoh,43)
save(out.all,file="../result/subsample_coor_All_PWOH.RData")

## coordination measures by HIV status and severity
out.1 <- NULL
out.2 <- NULL
## all ppts
one.out <- coordinated_measure(ds.all.pwoh)
one.out[[1]]$Group <- "All, PWOH"
one.out[[2]]$Group <- "All, PWOH"
out.1 <- bind_rows(out.1,one.out[[1]])
out.2 <- bind_rows(out.2,one.out[[2]])
ds.all.plwh <- ds %>% filter(hiv=="PLWH") %>%
  select(all_of(markers))
one.out <- coordinated_measure(ds.all.plwh)
## save the coordinate measures of PLWH for comparison with the coordinate measures of subsampling of PWOH
coor.all.plwh <- one.out[[1]]
size.all.plwh <- one.out[[2]]
one.out[[1]]$Group <- "All, PLWH"
one.out[[2]]$Group <- "All, PLWH"
out.1 <- bind_rows(out.1,one.out[[1]])
out.2 <- bind_rows(out.2,one.out[[2]])

## Asymptomatic
one.out <- coordinated_measure(ds.asym.pwoh)
one.out[[1]]$Group <- "Asymp, PWOH"
one.out[[2]]$Group <- "Asymp, PWOH"
out.1 <- bind_rows(out.1,one.out[[1]])
out.2 <- bind_rows(out.2,one.out[[2]])
ds.asym.plwh <- ds %>% filter(severity=="Asymptomatic" & hiv=="PLWH") %>%
  select(all_of(markers))
one.out <- coordinated_measure(ds.asym.plwh)
## save the coordinate measures of PLWH for comparison with the coordinate measures of subsampling of PWOH
coor.asym.plwh <- one.out[[1]]
size.asym.plwh <- one.out[[2]]
one.out[[1]]$Group <- "Asymp, PLWH"
one.out[[2]]$Group <- "Asymp, PLWH"
out.1 <- bind_rows(out.1,one.out[[1]])
out.2 <- bind_rows(out.2,one.out[[2]])

## Symptomatic
one.out <- coordinated_measure(ds.sym.pwoh)
one.out[[1]]$Group <- "Symp, PWOH"
one.out[[2]]$Group <- "Symp, PWOH"
out.1 <- bind_rows(out.1,one.out[[1]])
out.2 <- bind_rows(out.2,one.out[[2]])
ds.sym.plwh <- ds %>% filter(severity=="Symptomatic" & hiv=="PLWH") %>%
  select(all_of(markers))
one.out <- coordinated_measure(ds.sym.plwh)
## save the coordinate measures of PLWH for comparison with the coordinate measures of subsampling of PWOH
coor.sym.plwh <- one.out[[1]]
size.sym.plwh <- one.out[[2]]
one.out[[1]]$Group <- "Symp, PLWH"
one.out[[2]]$Group <- "Symp, PLWH"
out.1 <- bind_rows(out.1,one.out[[1]])
out.2 <- bind_rows(out.2,one.out[[2]])

## Hospitalized
one.out <- coordinated_measure(ds.hos.pwoh)
one.out[[1]]$Group <- "Hos, PWOH"
one.out[[2]]$Group <- "Hos, PWOH"
out.1 <- bind_rows(out.1,one.out[[1]])
out.2 <- bind_rows(out.2,one.out[[2]])
ds.hos.plwh <- ds %>% filter(severity=="Hospitalized" & hiv=="PLWH") %>%
  select(all_of(markers))
one.out <- coordinated_measure(ds.hos.plwh)
## save the coordinate measures of PLWH for comparison with the coordinate measures of subsampling of PWOH
coor.hos.plwh <- one.out[[1]]
size.hos.plwh <- one.out[[2]]
one.out[[1]]$Group <- "Hos, PLWH"
one.out[[2]]$Group <- "Hos, PLWH"
out.1 <- bind_rows(out.1,one.out[[1]])
out.2 <- bind_rows(out.2,one.out[[2]])

write.csv(out.1,"../result/coordinate_restricted.csv")
write.csv(out.2,"../result/coordinate_new_measure.csv")


## summarize the results
summary_res <- function(dat,lab,coor.plwh)
{
  n.vars <- dim(dat)[2]
  n.boots <- dim(dat)[1]
  tot.pairs <- n.vars*(n.vars-1)/2
  ## summarize the coordination measures
  out.mean <- round(dat %>% summarise_all(mean),0)
  out.std <- dat %>% summarise_all(var)
  out.std <- round(out.std %>% mutate_all(sqrt),0)
  out.mean <- paste0(out.mean," (",out.std,")")
  out.mean <- data.frame(Mean_SD=out.mean)
  
  out.med <- round(dat %>% summarise_all(median),0)
  out.q25 <- round(dat %>% summarize_all(quantile, 0.25),0)
  out.q75 <- round(dat %>% summarize_all(quantile, 0.75),0)
  out.med <- paste0(out.med," (",out.q25,", ",out.q75,")")
  out.med <- data.frame(Med_IQR=out.med)
  
  out.min <- round(dat %>% summarise_all(min),0)
  out.max <- round(dat %>% summarise_all(max),0)
  out.range <- paste0("(",out.min,", ",out.max,")")
  out.range <- data.frame(Range=out.range)
  
  ## proportion of n.boot subsets of PWOH with coordinated measure >= the coordinate measure for PLWH 
  tmp.ds <- data.frame(matrix(0,nrow = n.boots, ncol = n.vars ))
  names(tmp.ds) <- names(dat)
  for (i in 1:n.boots)
  {
    tmp.ds[i,] <- as.numeric(dat[i,]>=coor.plwh)
  }
  
  out.prob <- round(tmp.ds %>% summarise_all(mean),3) %>% pivot_longer(everything(),names_to = "Measure",
                                                                       values_to = "Prob_ge_PLWH")
  
  ## in count
  res.out <- bind_cols(out.mean,out.med) %>% bind_cols(., out.range) %>% bind_cols(.,out.prob) %>%
    select(Measure,Mean_SD,Med_IQR,Range,Prob_ge_PLWH)
  write.csv(res.out,paste0("../result/summary_count",lab,".csv"))
  
  ## histogram of the coordination measures
  
  plt.ds <- dat %>% pivot_longer(everything(), names_to = "Measure",values_to = "num_pairs") %>%
    mutate(prop=num_pairs/tot.pairs)
  
  plwh.ds <- coor.plwh %>% pivot_longer(everything(), names_to = "Measure",values_to = "n_pairs_plwh") %>%
    mutate(prop.plwh=n_pairs_plwh/tot.pairs)
  
  
  plt.out <- ggplot(plt.ds,aes(x=num_pairs)) +
    geom_histogram(fill="blue",color="black") +
    facet_wrap(vars(Measure),nrow=2,scales = "free") +
    geom_vline(aes(xintercept = n_pairs_plwh),plwh.ds,linetype="dotted",color="red",size=1.5)+
    theme(axis.title.x=element_text(),
          axis.title.y=element_text(),
          legend.position="none")+
    theme_bw()
  ggsave(paste0("../result/hist_count_",lab,".pdf"),plt.out,width = 10,height = 8)
  
  
  plt.out <- ggplot(plt.ds,aes(x=prop)) +
    geom_histogram(fill="blue",color="black") +
    facet_wrap(vars(Measure),nrow=2,scales = "free") +
    geom_vline(aes(xintercept = prop.plwh),plwh.ds,linetype="dotted",color="red",size=1.5)+
    theme(axis.title.x=element_text(),
          axis.title.y=element_text(),
          legend.position="none")+
    theme_bw()
  ggsave(paste0("../result/hist_prop_",lab,".pdf"),plt.out,width = 10,height = 8)
}
load("../result/subsample_coor_All_PWOH.RData")
summary_res(dat=out.all[[1]],lab="All_sim_res",coor.plwh=coor.all.plwh)
summary_res(out.all[[2]],lab="All_sim",size.all.plwh)

load("../result/subsample_coor_Asym_PWOH.RData")
summary_res(dat=out.asym[[1]],lab="Asym_sim_res",coor.plwh=coor.asym.plwh)
summary_res(out.asym[[2]],lab="Asym_sim",size.asym.plwh)

load("../result/subsample_coor_Sym_PWOH.RData")
summary_res(dat=out.sym[[1]],lab="Sym_sim_res",coor.plwh=coor.sym.plwh)
summary_res(out.sym[[2]],lab="Sym_sim",size.sym.plwh)

load("../result/subsample_coor_Hos_PWOH.RData")
summary_res(out.hos[[1]],lab="Hos_sim_res",coor.hos.plwh)
summary_res(out.hos[[2]],lab="Hos_sim",size.hos.plwh)

######
## End: Resampling analysis to evaluate sample size effect on network correlation analysis
#####

############
## Begin: network analysis
###########################


############
## End: network analysis
###########################

############
## Begin: classification analysis
###########################
library(randomForest);  library(stringr)
library("randomForest");library(MLeval) ##for AUC
library(caret) ##FOR AUC
library(ROCR)  # function to create performance for ROC
library(cvAUC)  # estimate cvAUC
library(MLeval);library(SuperLearner)
library(e1071)#SVM

tic <- Sys.time(); set.seed(1234);  x_all <- x; y_all <- y
max_sim1 <-200 # 200 #100 #200 # number of repeats
ss_m <- as.list(c("-99999" )); acc_v1 <- matrix(-999, nrow=max_sim1, ncol=4); ss_m_optimal <- as.list(c("-99999" )); ss_nn <-NULL
train_list1 <- as.list(c("-99999" )); test_list1 <- as.list(c("-99999" )); auc_v1 <- NULL; auc_v2 <- NULL; roc_p_data <- NULL; roc_p_data2 <- NULL; sup_auc1 <- NULL; sup_roc1 <- as.list(c(-99999))
for (s in 1:max_sim1)## repeat
{
  train_pt1 <- 0.9 # 0.75 # used 0.9 proportion of trnaing data
  inTrain <- createDataPartition(y, p = train_pt1, list = FALSE)[,1]  #<=======   split data
  x_train <- x[ inTrain, ]; x_all <- x
  x_test  <- x[-inTrain, ]
  
  y_train <- y[ inTrain]; y_all <- y
  y_test  <- y[-inTrain]
  
  # #**********************************************************************************************
  ## recursive selection:
  result_rfe1 <- rfe(x = x_train,
                     y = y_train,
                     sizes = c(1:length(x_train[1,])),
                     rfeControl = control);   ##e.g. sizes = c(1:5, 10, 13) for use 1 to 5 features, 10 features, and 13 features.
  
  # Print the results
  result_rfe1
  train_out_t <- postResample(predict(result_rfe1, x_train), y_train); wt0 <- predict(result_rfe1, x_train); wt0[ order(wt0[,1], wt0[,2]),]
  ss_m_optimal[[s]] <- result_rfe1$optVariables; ss_nn <- rbind(ss_nn, result_rfe1$optsize)
  # Print the selected features //  result_rfe1$optVariables  are the list of omtimal variable names selected// result_rfe1$optsize is optimal number of var selected
  # predictors(result_rfe1)
  ww1 <- predictors(result_rfe1);  ww2 <- postResample(predict(result_rfe1, x_test), y_test);
  leave_out_t <- postResample(predict(result_rfe1, x_test), y_test); wt <- predict(result_rfe1, x_test); wt[ order(wt[,1], wt[,2]),]
  train_list1[[s]] <- wt0
  
  ss_m[[s]] <- ww1;   acc_v1[s,1] <- ww2[1]; acc_v1[s,2] <- ww2[2]; acc_v1[s,3] <- leave_out_t[1]; acc_v1[s,4] <- leave_out_t[2]
  print(c(s, max_sim1))
  test_list1[[s]] <- cbind(y_test,wt)
  
  
  print(c(s, max_sim1)) 
}
toc <- Sys.time(); print(toc-tic) 
ss_m  # slection results

#------------now use top selected markers predict in loop--------------------#

top7m <- c("BCP.IgG.Mem.S", "End.IgG1.OC43.RBD", "BCP.IgG.Tot.S", "MSD.IgG.RBD", "BCP.Tot.S", "ACE2","End.IgG1.NL63.RBD") #new 2024_feb23

set.seed(202403)  
sel_auc_v1 <- NULL; sel_roc<- as.list(NULL); sel_roc_dd<- as.list(NULL);  tic1 <- Sys.time()
for (km in 1:50){
  k <- 7 
  nm_m1 <- names(x_all) %in% top7m; sel_vv1 <- which(nm_m1==TRUE)
  Y <- as.numeric(y_all);  Y <- (Y-1) ## use ALL data
  # X <- x_train  ## use ALL var
  X <- as.data.frame(x_all[,sel_vv1]) ## use ALL data
  
  sl_fit = CV.SuperLearner(Y =Y, X = as.data.frame(X),family = binomial(), cvControl =list(V = 10),  innerCvControl = list(list(V = 10)),
                           parallel = "multicore", SL.library = c("SL.randomForest", "SL.svm","SL.extraTrees"   ), method = method.AUC)
  
  sl.auc <- cv_auc(pred = sl_fit$SL.predict, Y = sl_fit$Y)
  xc <- (str_sub(sl.auc[2], 2, str_length(sl.auc[2])-1)); c_loc1 <- str_locate(xc,","); ci_c1 <- str_sub(xc, 1,c_loc1[1]-1); ci_c2 <- str_sub(xc, c_loc1[2]+1, str_length(xc))
  n_auc1 <- c(as.numeric(sl.auc[1]), as.numeric(ci_c1), as.numeric(ci_c2))
  
  sel_auc_v1 <- rbind(sel_auc_v1, n_auc1)
  
  ## calculate sensitivity and specificity for ROC curve
  pred.y=prediction(sl_fit$SL.predict,Y)    ## use str(pred.y)  to tkae a look
  perf=performance(pred.y,"tpr","fpr") ## use str(perf)  to tkae a look
  sens=unlist(perf@y.values)
  spec_1=unlist(perf@x.values)
  roc.dat=data.frame(sens=sens,spec_1=spec_1); sup_roc1[[s]] <- roc.dat
  
  sel_roc[[km]] <- roc.dat; sel_roc_dd[[km]] <- perf;
  print(km)
} ## km LOOP   2.7hr run
toc1 <- Sys.time(); toc1-tic1    #Aug6:
sel_auc_v1 ## slected UAC for coresponding varib;e used (in k)                  

dev.off() 
tiff(paste("","Final_marker_plot_sm1.tiff", sep=""), units="in", width=5, height=5, res=300)        

m_roc_sen <- NULL;  m_roc_spe <- NULL; mean_roc <- NULL;  mean_roc_y <- NULL  #==========================   ROC  plot 
for (h in 1:length(sel_roc_dd)){
  perf1 <- sel_roc_dd[[h]]; this_roc <- sel_roc[[h]]
  xseq <- seq(0,1, 0.01); approx_d <- approx(this_roc$spec_1, this_roc$sens,  xseq, ties = mean)# spec_1 is X,sens-y; interpolation
  mean_roc_x <- cbind(mean_roc_x, approx_d$x)
  mean_roc_y <- cbind(mean_roc_y, approx_d$y) 
  
  # m_roc_sen <- cbind(m_roc_sen, this_roc$sens );  m_roc_spe <- cbind(m_roc_spe, this_roc$spec_1 ); print(c(length(this_roc$sens), length(this_roc$spec_1)))
  if (h==1){plot(perf1, col="grey78", main= paste( paste(" ", "mean AUC=",  as.character(colMeans(sel_auc_v1)[1]), "(Approx.95% CI: 0.xx-0.xx)\n",boost1)))
    abline(a=0, b= 1) #diagnol line
  } ## <==========   chnage the file name
  if(h>=2){plot(perf1, col="grey78",  add=TRUE)} #perf1: x.name: "False positive rate" x-axis //   y.name: "True positive rate", y-axis
  # 
  #text(0.12, 0.97, paste("AUC=",round(m_auc[1],3), sep="")); text(0.12, 0.9, "95%CI");text(0.12, 0.83, paste( "(",round(m_auc[7],2),"-",round(m_auc[8],2), ")",sep=""))
}
#m_sen <- rowMeans(m_roc_sen); m_spe <- rowMeans(m_roc_spe);  #<===== this may not work since each run has different length
m_sen <- rowMeans(mean_roc_y); m_spe <- rowMeans(mean_roc_x); m_sen <- c(0,m_sen ); m_spe <- c(0,m_spe );  #<===== this may not work since each run has different length
perf1_0mk <- perf1
perf2 <- perf1; perf2@x.values <- as.list(m_sen); perf2@y.values <- as.list(m_spe)
#plot(perf1,col="red",  add=TRUE)
lines(m_spe, m_sen, col="red" ) #need spline in above
#  abline(v=c(0.1,0.05), col=c("blue", "purple"), lty=c(1,2), lwd=c(1, 3))#
#  abline(h=c(0.9,0.95), col=c("blue", "purple"), lty=c(1,2), lwd=c(1, 3))#
dev.off() 

### exhaust search among all combinations ####
top7m <- c("BCP.IgG.Mem.S", "End.IgG1.OC43.RBD", "BCP.IgG.Tot.S", "MSD.IgG.RBD", "BCP.Tot.S", "ACE2","End.IgG1.NL63.RBD") #new 2024_feb23
x_ori <- x_all 

tic <- Sys.time();
comb_list11 <- as.list(NULL); comb_list12 <- as.list(NULL); comb_list13 <- as.list(NULL); comb_list14 <- as.list(NULL); comb_list15 <- as.list(NULL); ccc <- 0
for (mk in 1:length(top7m)){  
  #   for (mk in 1:1){
  all_cb <- combn(length(top7m), mk)
  
  for(all_c in 1:length(all_cb[1,]) ){
    x_all <- x_ori
    this_top <- top7m[all_cb[ ,all_c]]
    
    name_use <- names(x_all) %in% this_top; col_sel <- which(name_use==TRUE); 
    x_all <- x_all[,col_sel]; name_use1 <- names(x_all) #use selected var //slec_nm refer above , 
    if (length(col_sel)==1){ sing_nm <- names(x_ori)[col_sel]}
    temp_auc1 <- NULL; sup_list1 <- as.list(NULL);  auc_II <- NULL; set.seed(56789)
    
    #for(k in 1:50){ 
    for(k in 1:30){  #50, was20
      if (length(col_sel)>=2){
        SP_sl = CV.SuperLearner(Y =(as.numeric(y)-1), X = x_all, family = binomial(), cvControl =list(V = 10), parallel = "multicore", SL.library = c("SL.randomForest", "SL.svm","SL.extraTrees"  ), method = method.AUC)#AUC  not converge/using SL.step.interaction replacing SL.glm.interaction 
      }else{
        SP_sl = CV.SuperLearner(Y =(as.numeric(y)-1), X = as.data.frame(x_all), family = binomial(), cvControl =list(V = 10), parallel = "multicore", SL.library = c("SL.randomForest", "SL.svm","SL.extraTrees"  ), method = method.AUC)#AUC  not converge/using SL.step.interaction replacing SL.glm.interaction 
        
      }
      #"SL.step.interaction" is very slow
      
      
      sl <- SP_sl
      mm2 <- data.frame(cbind(sl$SL.predict, (as.numeric(y )-1))) # each variant you need to re-load "py" and "xx_t"  in above for different variant <==== 2/4
      pred2 <- ROCR::prediction(mm2$X1, mm2$X2)
      perf2 <- performance(pred2, "tpr", "fpr") 
      auc_est1 = ROCR::performance(pred2, measure = "auc", x.measure = "cutoff")@y.values[[1]]; temp_auc1 <- rbind(temp_auc1, auc_est1)
      # plot(perf2, colorize=FALSE, col="red", main=paste("Lambda AUC=",auc_est1)); abline(coef = c(0,1), col="grey88", lwd=3, lty=2 ) 
      sup_list1[[k]] <- SP_sl; print(k)
      sl.auc <- cv_auc(pred = sl$SL.predict, Y = sl$Y); #auc_II <- rbind(auc_II, sl.auc)
      xc <- (str_sub(sl.auc[2], 2, str_length(sl.auc[2])-1)); c_loc1 <- str_locate(xc,","); ci_c1 <- str_sub(xc, 1,c_loc1[1]-1); ci_c2 <- str_sub(xc, c_loc1[2]+1, str_length(xc))
      n_auc1 <- c(as.numeric(sl.auc[1]), as.numeric(ci_c1), as.numeric(ci_c2))
      auc_II <- rbind(auc_II, n_auc1)
    }# k loop
    average1 <- colMeans(auc_II)
    ccc <-  ccc+1; print(c(mk, k,ccc))
    comb_list11[ccc] <- mk; 
    if (length(col_sel)>=2){comb_list12[ccc] <- as.data.frame(name_use1);}else{comb_list12[ccc] <- sing_nm} 
    
    
    comb_list13[ccc] <- average1[1]    #as.numeric(paste(sl.auc[1]))
    comb_list14[ccc] <- average1[2]  
    comb_list15[ccc] <- average1[3]  
    
  }# all_c LOOP
  
} # mk loop
############
## Begin: classification analysis
###########################
