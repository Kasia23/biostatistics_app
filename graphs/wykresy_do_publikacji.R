#³adujemy data1

legend.labs <- list( 
  "PTEN.ABN" = c( 'PTEN.WT', 'PTEN.ABN'),
  "PTEN.MUT" = c('PTEN.WT', 'PTEN.MUT'),
  "PTEN.DEL" = c( 'PTEN.WT', 'PTEN.DEL'),
  "Risk.group..MRD.d15" = c('standard risk group\n (<0,1% blast cells)\n', 'intermediate risk group\n (0.1 to <10%  blast cells)\n', 'high risk group\n (>10% blast cells)'),
  "Risk.group.MRD.d33." = c('standard risk group\n (<0,1% blast cells)\n', 'intermediate risk group\n (0.1 to <10%  blast cells)\n', 'high risk group\n (>10% blast cells)'),
  "Risk.group.MRD.d78" = c('standard risk group\n (<0,1% blast cells)\n', 'intermediate risk group\n (0.1 to <10%  blast cells)'),
  "PTEN.zMUTbezDEL" = c( 'PTEN.WT', 'PTEN.MUT'),
  "PTEN.bezMUTzDEL" = c('PTEN.WT', 'PTEN.DEL'),
  "gLoR_gHoR"= c('gLoR', 'gHoR'),
  "Notch1Fbxw7PTEN_2gr" = c('NOTCH1/FBXW7mut&PTEN.WT', 'NOTCH1/FBXW7/PTEN.WT'),
  "NOTCH1.Sanger" = c('NOTCH1.WT', 'NOTCH1.MUT'),
  "Notch1Fbxw7_2gr" = c('NOTCH1/FBXW7.WT', 'NOTCH1/FBXW7.MUT'),
  "PTEN.mono.vs.biallelic.inactivation" = c("PTEN.WT", "PTENbiallelic", "PTENmonoallelic")
)

######################################## S1 ############################################################
########################################################################################################
#ustawienie parametrów wywo³ania:
czas <- "OS"
a <- -0.025
b <- 0.03
plik <- "S2_OS.png"
plik_tiff <- "S2_OS.tiff"
########################################################################################################
czas <- "EFS"
a <- -0.025
b <- 0.03
plik <- "S2_EFS.png"
plik_tiff <- "S2_EFS.tiff"
########################################################################################################
czas <- "DFS"
a <- -0.03
b <- 0.033
plik <- "S2_DFS.png"
plik_tiff <- "S2_DFS.tiff"
########################################################################################################
czas <- "RFS"
a <- -0.03
b <- 0.033
plik <- "S2_RFS.png"
plik_tiff <- "S2_RFS.tiff"
########################################################################################################


ind <- paste0("ind", czas)
data1$risk.group..SR.IR.HR.
#dodatkowe analizy!
subset= (PTEN.ABN == 0 & risk.group..SR.IR.HR. == 2)

s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~ ALL.IC.BFM.2002.vs..2009")), data = data1)
w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~ ALL.IC.BFM.2002.vs..2009")), data = data1)
pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)

ggsurv <- ggsurvplot(s,
                     # pval=round(pv,2),
                     legend.title = "",
                     #title = "Krzywe przezycia dla protokolow", 
                     xlab = "Years", xlim = c(0,10), break.time.by = 1, ylab = paste0("p", czas),
                     #conf.int = TRUE, 
                     risk.table = TRUE, risk.table.pos="out", 
                     risk.table.height = 0.4, risk.table.y.text = FALSE, 
                     risk.table.col='strata', data=data1, pval = round(pv,3),
                     legend =  c(0.8, 0.2),
                     legend.labs = c( 'ALL IC-BFM 2002', 'ALL IC-BFM 2009'),
                     ggtheme =  theme_survminer())

ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))


#wyci¹gamy dane:

type = 1

n5 = ggsurv$data.survtable[(ggsurv$data.survtable$time == 5) & (ggsurv$data.survtable$ALL.IC.BFM.2002.vs..2009 == type), ]$n.risk
surv_table_5 = ggsurv$data.survplot[(ggsurv$data.survplot$n.risk == n5) & (ggsurv$data.survplot$ALL.IC.BFM.2002.vs..2009 == type), ]
surv1 = round(surv_table_5$surv,2)
st_err1 = round(surv_table_5$std.err,2)

type = 2

n5 = ggsurv$data.survtable[(ggsurv$data.survtable$time == 5) & (ggsurv$data.survtable$ALL.IC.BFM.2002.vs..2009 == type), ]$n.risk
if(n5==0){n5=1}
surv_table_5 = ggsurv$data.survplot[(ggsurv$data.survplot$n.risk == n5) & (ggsurv$data.survplot$ALL.IC.BFM.2002.vs..2009 == type), ]
surv2 = round(surv_table_5$surv,2)
st_err2 = round(surv_table_5$std.err,2)

n1 = summary(s)$table[1,'records']
n2 = summary(s)$table[2,'records']

e1 = summary(s)$table[1,'events']
e2 = summary(s)$table[2,'events']

text1 <- ifelse(is.na(surv1), paste0('(n=',n1, ', ', e1,' events)'), paste0(surv1,', SE=', st_err1,' (n=',n1, ', ', e1,' events)'))
text2 <- ifelse(is.na(surv2), paste0('(n=',n2, ', ', e2,' events)'), paste0(surv2,', SE=', st_err2,' (n=',n2, ', ', e2,' events)'))


ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold")) +
        annotate("text", x=5, y=surv1+a, label=text1, size=5) +
        annotate("text", x=5, y=surv2+b, label=text2, size=5)


ggsurv[1]

png(plik, width = 820, height = 480)
ggsurv[1]
dev.off()
plik_tiff <- "S2_EFS_ABN_gr2.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
ggsurv[1]
dev.off()
######################################### S5, S6 ##################################################
c <- 0
##################################################################################################
czas <- "EFS"
i <- "PTEN.bezMUTzDEL"
subset_protokol = 1
a=0.045
b=0.033
dig <- 2
plik <- "S5_EFS_PTEN.DEL_2002.png"
plik_tiff <- "S5_EFS_PTEN.DEL_2002.tiff"
###################################################################################################
czas <- "EFS"
i <- "PTEN.bezMUTzDEL"
subset_protokol = 2
a=0.045
b=0.033
dig <- 4
y_value <- 0.45
plik <- "S5_EFS_PTEN.DEL_2009.png"
plik_tiff <- "S5_EFS_PTEN.DEL_2009.tiff"
###################################################################################################
czas <- "EFS"
i <- "PTEN.zMUTbezDEL"
subset_protokol = 1
a=0.045
b=-0.033
dig <- 2
plik <- "S5_EFS_PTEN.MUT_2002.png"
plik_tiff <- "S5_EFS_PTEN.MUT_2002.tiff"
###################################################################################################
czas <- "EFS"
i <- "PTEN.zMUTbezDEL"
subset_protokol = 2
a=0.045
b=0.033
dig <- 3
y_value <- 0.37
plik <- "S5_EFS_PTEN.MUT_2009.png"
plik_tiff <- "S5_EFS_PTEN.MUT_2009.tiff"
###################################################################################################
czas <- "EFS"
i <- "PTEN.ABN"
subset_protokol = 1
a=0.045
b=0.033
c <- 2
dig <- 2
plik <- "S6_EFS_PTEN.ABN_2002.png"
plik_tiff <- "S6_EFS_PTEN.ABN_2002.tiff"
###################################################################################################
czas <- "EFS"
i <- "PTEN.ABN"
subset_protokol = 2
a=0.045
b=0.033
c <- 0
dig <- -1
y_value <- 0.455
plik <- "S6_EFS_PTEN.ABN_2009.png"
plik_tiff <- "S6_EFS_PTEN.ABN_2009.tiff"
###################################################################################################
###################################################################################################
ind <- paste0("ind", czas)
s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset = ALL.IC.BFM.2002.vs..2009==subset_protokol, data = data1)

w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset = ALL.IC.BFM.2002.vs..2009==subset_protokol, data = data1)
pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)
p = pv

ggsurv <- ggsurvplot(s,
                     # pval=round(pv,2),
                     legend.title = "",
                     #title = paste0("Protokol 2002: ",i), 
                     xlab = "Years", xlim = c(0,10), break.time.by = 1, ylab = paste0("p", czas),
                     #conf.int = TRUE, 
                     risk.table = TRUE, risk.table.pos="out", 
                     risk.table.height = 0.4, risk.table.y.text = FALSE, 
                     risk.table.col='strata', data=data1,
                     legend =  c(0.8, 0.2), pval = ifelse(dig > 0, paste0('p=',formatC(p, digits=dig, format="f")), 'p<0.0001'),
                     legend.labs = legend.labs[[i]],
                     ggtheme =  theme_survminer())

ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))


#wyci¹gamy dane:

type = 1

n5 = ggsurv$data.survtable[(ggsurv$data.survtable$time == 5), ]$n.risk
surv_table_5 = ggsurv$data.survplot[(ggsurv$data.survplot$n.risk == n5[type]),][type,]
surv1 = round(surv_table_5$surv,2)
st_err1 = round(surv_table_5$std.err,2)

type = 2

n5 = ggsurv$data.survtable[(ggsurv$data.survtable$time == 5), ]$n.risk
surv_table_5 = ggsurv$data.survplot[(ggsurv$data.survplot$n.risk == n5[type]),][type,]
surv2 = round(surv_table_5$surv,2)
st_err2 = round(surv_table_5$std.err,2)

n1 = summary(s)$table[1,'records']
n2 = summary(s)$table[2,'records']

e1 = summary(s)$table[1,'events']
e2 = summary(s)$table[2,'events']

text1 <- ifelse(is.na(surv1), paste0('(n=',n1, ', ', e1,' events)'), paste0(surv1,', SE=', st_err1,' (n=',n1, ', ', e1,' events)'))
text2 <- ifelse(is.na(surv2), paste0('(n=',n2, ', ', e2,' events)'), paste0(surv2,', SE=', st_err2,' (n=',n2, ', ', e2,' events)'))

text1
text2

ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold")) +
  annotate("text", x=5, y=surv1+a, label= text1, size=5) +
  annotate("text", x=3+c, y=ifelse(is.na(surv2), y_value, surv2+b), label= text2, size =5)


ggsurv[1]

png(plik, width = 820, height = 480)
ggsurv[1]
dev.off()


ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
ggsurv[1]
dev.off()
######################################### S3, S4, S7 + dod ############################################
generate_plot <- function(czas, i, a, b, dig=2){
  ind <- paste0("ind", czas)
  
  s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~", i)), data = data1)
  
  w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~", i)), data = data1)
  pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)
  pv
  
  ggsurv <- ggsurvplot(s,
                       # pval=round(pv,2),
                       legend.title = "",
                       #title = paste0("Protokol 2002: ",i), 
                       xlab = "Years", xlim = c(0,10), break.time.by = 1, ylab = paste0("p", czas),
                       #conf.int = TRUE, 
                       risk.table = TRUE, risk.table.pos="out", 
                       risk.table.height = 0.4, risk.table.y.text = FALSE, 
                       risk.table.col='strata', data=data1, pval = ifelse(dig > 0, paste0('p=',formatC(p, digits=dig, format="f")), 'p<0.001'), 
                       #pval = "p < 0.001",
                       legend =  c(0.8, 0.2),
                       legend.labs = legend.labs[[i]],
                       ggtheme =  theme_survminer())
  
  ggsurv$plot <- ggsurv$plot + 
    theme(legend.text = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
          axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
          axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
          axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))
  
  
  #wyci¹gamy dane:
  
  type = 1
  
  n5 = ggsurv$data.survtable[(ggsurv$data.survtable$time == 5), ]$n.risk
  surv_table_5 = ggsurv$data.survplot[(ggsurv$data.survplot$n.risk == n5[type]),][type,]
  surv1 = round(surv_table_5$surv,2)
  st_err1 = round(surv_table_5$std.err,2)
  
  type = 2
  
  n5 = ggsurv$data.survtable[(ggsurv$data.survtable$time == 5), ]$n.risk
  surv_table_5 = ggsurv$data.survplot[(ggsurv$data.survplot$n.risk == n5[type]),][type,]
  surv2 = round(surv_table_5$surv,2)
  st_err2 = round(surv_table_5$std.err,2)
  
  n1 = summary(s)$table[1,'records']
  n2 = summary(s)$table[2,'records']
  
  e1 = summary(s)$table[1,'events']
  e2 = summary(s)$table[2,'events']
  
  text1 = paste0(surv1,', SE=', st_err1,' (n=',n1, ', ', e1,' events)')
  text2 = paste0(surv2,', SE=', st_err2,' (n=',n2, ', ', e2,' events)')
  text1
  text2
  
  ggsurv$plot <- ggsurv$plot + 
    theme(legend.text = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
          axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
          axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
          axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold")) +
    annotate("text", x=5, y=surv1+a, label= text1, size=5) +
    annotate("text", x=5, y=surv2+b, label= text2, size =5)
  
  
  return(ggsurv[1])
}
#######################################################################################################
png(plik, width = 820, height = 480)
ggsurv[1]
dev.off()
#######################################################################################################

czas <- "OS" 
i <- "NOTCH1.Sanger"
a=-0.035
b=0.06
dig <- 2
plik <- "S3_NOTCH1_sanger_OS.png"
plik_tiff <- "S3_NOTCH1_sanger_OS.tiff"
#######################################################################################################
czas <- "EFS" 
i <- "NOTCH1.Sanger"
a=-0.035
b=0.033
plik <- "S3_NOTCH1_sanger_EFS.png"
plik_tiff <- "S3_NOTCH1_sanger_EFS.tiff"
#######################################################################################################
czas <- "DFS" 
i <- "NOTCH1.Sanger"
a=-0.035
b=0.04
plik <- "S3_NOTCH1_sanger_DFS.png"
plik_tiff <- "S3_NOTCH1_sanger_DFS.tiff"
#######################################################################################################
czas <- "RFS" 
i <- "NOTCH1.Sanger"
a=0.053
b=-0.035
plik <- "S3_NOTCH1_sanger_RFS.png"
plik_tiff <- "S3_NOTCH1_sanger_RFS.tiff"
#######################################################################################################
czas <- "OS" 
i <- "Notch1Fbxw7_2gr"
a=-0.035
b=0.06
plik <- "S4_NOTCH1_and_FBXW7_OS.png"
plik_tiff <- "S4_NOTCH1_and_FBXW7_OS.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
generate_plot(czas, i, a, b)
dev.off()
#######################################################################################################
czas <- "EFS" 
i <- "Notch1Fbxw7_2gr" 
a=-0.035
b=0.033
plik <- "S4_NOTCH1_and_FBXW7_EFS.png"
plik_tiff <- "S4_NOTCH1_and_FBXW7_EFS.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
generate_plot(czas, i, a, b)
dev.off()
#######################################################################################################
czas <- "DFS" 
i <- "Notch1Fbxw7_2gr"
a=-0.036
b=0.033
plik <- "S4_NOTCH1_and_FBXW7_DFS.png"
plik_tiff <- "S4_NOTCH1_and_FBXW7_DFS.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
generate_plot(czas, i, a, b)
dev.off()
#######################################################################################################
czas <- "RFS" 
i <- "Notch1Fbxw7_2gr"
a=-0.036
b=0.033
plik <- "S4_NOTCH1_and_FBXW7_RFS.png"
plik_tiff <- "S4_NOTCH1_and_FBXW7_RFS.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
generate_plot(czas, i, a, b)
dev.off()
#######################################################################################################
#######################################################################################################
czas <- "OS" 
i <- "gLoR_gHoR"
a=-0.035
b=0.06
plik <- "S7_gLoR_gHoR_OS.png"
plik_tiff <- "S7_gLoR_gHoR_OS.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
generate_plot(czas, i, a, b)
dev.off()
#######################################################################################################
czas <- "EFS" 
i <- "gLoR_gHoR"
a=-0.035
b=0.033
plik <- "S7_gLoR_gHoR_EFS.png"
plik_tiff <- "S7_gLoR_gHoR_EFS.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
generate_plot(czas, i, a, b)
dev.off()
#######################################################################################################
czas <- "DFS" 
i <- "gLoR_gHoR"
a=-0.036
b=0.033
plik <- "S7_gLoR_gHoR_DFS.png"
plik_tiff <- "S7_gLoR_gHoR_DFS.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
generate_plot(czas, i, a, b)
dev.off()
#######################################################################################################
czas <- "RFS" 
i <- "gLoR_gHoR"
a=-0.036
b=0.033
plik <- "S7_gLoR_gHoR_RFS.png"
plik_tiff <- "S7_gLoR_gHoR_RFS.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
generate_plot(czas, i, a, b)
dev.off()
#######################################################################################################
czas <- "EFS" 
i <- "PTEN.DEL"
a=-0.036
b=0.033
plik <- "PTEN_DEL.png"
plik_tiff <- "PTEN_DEL.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
generate_plot(czas, i, a, b, dig = -1)
dev.off()
#######################################################################################################
czas <- "EFS" 
i <- "PTEN.MUT"
a=-0.036
b=0.033
plik <- "PTEN_MUT.png"
plik_tiff <- "PTEN_MUT.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
generate_plot(czas, i, a, b, dig=-1)
dev.off()
#######################################################################################################

#######################################################################
############################## S8, S9 ####################################

i <- "PTEN.bezMUTzDEL"
czas <- "EFS"
val <- 1
day <- "Risk.group..MRD.d33"
a=-0.036
b=-0.033
plik <- "S8_gr1.tiff"

s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group.MRD.d33.==val, data = data1)
w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group.MRD.d33.==val, data = data1)

###################################################################

i <- "PTEN.ABN"
czas <- "EFS"
val <- 2
day <- "Risk.group..MRD.d15"
a=-0.033
b=0.033
plik <- "S9_d15_group2_2009.tiff"
data1$Risk.group..MRD.d15
s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset= (ALL.IC.BFM.2002.vs..2009 ==val & Risk.group..MRD.d15 == 1), data = data1)
w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset= (ALL.IC.BFM.2002.vs..2009 ==val & Risk.group..MRD.d15 == 1), data = data1)
summary(s)
###################################################################

i <- "PTEN.ABN"
czas <- "EFS"
val <- 2
day <- "Risk.group..MRD.d15"
a=-0.033
b=0.033
plik <- "S9_pppp.tiff"
data1$PTEN.ABN
s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~","ALL.IC.BFM.2002.vs..2009")), subset= (PTEN.ABN == 0), data = data1)
w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~","ALL.IC.BFM.2002.vs..2009")), subset= (PTEN.ABN == 0), data = data1)
summary(s)
###################################################################
pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)

ggsurv <- ggsurvplot(s,
                     # pval=round(pv,2), pval = "p<0.001"
                     legend.title = "",
                     #title = paste0("Protokol 2002: ",i), 
                     xlab = "Years", xlim = c(0,10), break.time.by = 1, ylab = paste0("p", czas),
                     #conf.int = TRUE, 
                     risk.table = TRUE, risk.table.pos="out", 
                     risk.table.height = 0.4, risk.table.y.text = FALSE, 
                     risk.table.col='strata', data=data1, pval=round(pv,2), legend =  c(0.8, 0.2),
                     legend.labs = legend.labs[[i]],
                     ggtheme =  theme_survminer())

ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))


#wyci¹gamy dane:

type = 1

n5 = ggsurv$data.survtable[(ggsurv$data.survtable$time == 5), ]$n.risk
surv_table_5 = ggsurv$data.survplot[(ggsurv$data.survplot$n.risk == n5[type]),][type,]
surv1 = round(surv_table_5$surv,2)
st_err1 = round(surv_table_5$std.err,2)

type = 2

n5 = ggsurv$data.survtable[(ggsurv$data.survtable$time == 5), ]$n.risk
surv_table_5 = ggsurv$data.survplot[(ggsurv$data.survplot$n.risk == n5[type]),][type,]
surv2 = round(surv_table_5$surv,2)
st_err2 = round(surv_table_5$std.err,2)

n1 = summary(s)$table[1,'records']
n2 = summary(s)$table[2,'records']

e1 = summary(s)$table[1,'events']
e2 = summary(s)$table[2,'events']


text1 <- ifelse(is.na(surv1), paste0('(n=',n1, ', ', e1,' events)'), paste0(surv1,', SE=', st_err1,' (n=',n1, ', ', e1,' events)'))
text2 <- ifelse(is.na(surv2), paste0('(n=',n2, ', ', e2,' events)'), paste0(surv2,', SE=', st_err2,' (n=',n2, ', ', e2,' events)'))

text1
text2

y_value <- 0.37
y_value <- 0.53
y_value <- 0.28
y1_value <- 0.93
ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold")) +
  annotate("text", x=ifelse(is.na(surv2),6.3-2,5), y=ifelse(is.na(surv2), y1_value, surv1-a), label= text1, size=5) +
  annotate("text", x=ifelse(is.na(surv2),4-1,5), y=ifelse(is.na(surv2), y_value, surv2+b), label= text2, size =5)

ggsurv[1]

tiff(plik, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
ggsurv[1]
dev.off()


#######################################################################################################
png(plik, width = 820, height = 480)
ggsurv[1]
dev.off()


#######################################################################
############################## S10 ####################################
czas <- "RFS"
i <- "Risk.group.MRD.d78"
a=-0.036
b=0.033
l1 <- 6.5
l2 <- 5
l3 <- 3.7
dig <- 3
plik <- "S10_B_MRD78_RFS.png"
plik_tiff <- "S10_B_MRD78_RFS.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
ggsurv[1]
dev.off()
########################################################
czas <- "RFS"
i <- "Risk.group..MRD.d15"
l1 = l2 = l3 = 5
l4 <- 0.75
b <- 0.033
a <- 0.033
plik <- "S10_A_MRD15_RFS.png"
plik_tiff <- "S10_A_MRD15_RFS.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
ggsurv[1]
dev.off()
########################################################
czas <- "RFS"
i <- "Risk.group.MRD.d33."
a = b = -0.033
l4 <- 0.33
l3 <- 3.6
l1 = l2 = 6
plik <- "S6_MRD33_RFS.png"
plik_tiff <- "S6_MRD33_RFS.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
ggsurv[1]
dev.off()
########################################################
czas <- "EFS"
i <- "PTEN.mono.vs.biallelic.inactivation"
a=-0.033
b=-0.033
l4 <- 0.46 - 0.033
l3 <- 6
dig <- -1
plik <- "4G.png"
plik_tiff <- "4G.tiff"
ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
ggsurv[1]
dev.off()
########################################################
ind <- paste0("ind", czas)

s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~", i)), data = data1)

w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~", i)), data = data1)
pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)


ggsurv <- ggsurvplot(s,
                     # pval=round(pv,2),
                     legend.title = "",
                     #title = paste0("Protokol 2002: ",i), 
                     xlab = "Years", xlim = c(0,10), break.time.by = 1, ylab = paste0("p", czas),
                     #conf.int = TRUE, 
                     risk.table = TRUE, risk.table.pos="out", 
                     risk.table.height = 0.4, risk.table.y.text = FALSE, 
                     risk.table.col='strata', data=data1, pval = ifelse(dig > 0, paste0('p=',formatC(p, digits=dig, format="f")), 'p<0.001'),
                     legend =  c(0.8, 0.2),
                     legend.labs = legend.labs[[i]],
                     ggtheme =  theme_survminer())

ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))


#wyci¹gamy dane:

type = 1

n5 = ggsurv$data.survtable[(ggsurv$data.survtable$time == 5), ]$n.risk
surv_table_5 = ggsurv$data.survplot[(ggsurv$data.survplot$n.risk == n5[type]),][type,]
surv1 = round(surv_table_5$surv,2)
st_err1 = round(surv_table_5$std.err,2)

type = 2

n5 = ggsurv$data.survtable[(ggsurv$data.survtable$time == 5), ]$n.risk
surv_table_5 = ggsurv$data.survplot[(ggsurv$data.survplot$n.risk == n5[type]),][type,]
surv2 = round(surv_table_5$surv,2)
st_err2 = round(surv_table_5$std.err,2)

if (length(n5)==3){
  type = 3
  
  n5 = ggsurv$data.survtable[(ggsurv$data.survtable$time == 5), ]$n.risk
  surv_table_5 = ggsurv$data.survplot[(ggsurv$data.survplot$n.risk == n5[type]),][type,]
  surv3 = round(surv_table_5$surv,2)
  st_err3 = round(surv_table_5$std.err,2)

  n3 = summary(s)$table[type,'records']
  e3 = summary(s)$table[type,'events']
  text3 = paste0(surv3,', SE=', st_err3,' (n=',n3, ', ', e3,' events)')
  #text3 = paste0('(n=',n3, ', ', e3,' events)')
    
}

n1 = summary(s)$table[1,'records']
n2 = summary(s)$table[2,'records']

e1 = summary(s)$table[1,'events']
e2 = summary(s)$table[2,'events']

text1 = paste0(surv1,', SE=', st_err1,' (n=',n1, ', ', e1,' events)')
text2 = paste0(surv2,', SE=', st_err2,' (n=',n2, ', ', e2,' events)')

ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold")) +
  annotate("text", x=l1, y=surv1-a, label= text1, size=5) +
  annotate("text", x=l2, y=surv2+b, label= text2, size =5) +
  annotate("text", x=l3, y=surv3-b, label= text3, size =5)

ggsurv[1]
#######################################################################################################
png(plik, width = 820, height = 480)
ggsurv[1]
dev.off()


#####################################################################################################
# Fig 4 bez strata
czas <- "EFS"
i <- "PTEN.ABN"
a <- -0.033
b <- 0.033
dig <- -1
plik <- "Fig4_KM_ABN.png"
plik_tiff <- "Fig4_KM_ABN.tiff"
########################################################################################################
czas <- "EFS"
i <- "PTEN.bezMUTzDEL"
a <- -0.033
b <- 0.033
plik <- "Fig4_KM_DEL.png"
plik_tiff <- "Fig4_KM_DEL.tiff"
########################################################################################################
czas <- "EFS"
i <- "PTEN.zMUTbezDEL"
a <- -0.033
b <- 0.033
dig <- 1
plik <- "Fig4_KM_MUT.png"
plik_tiff <- "Fig4_KM_MUT.tiff"
########################################################################################################

ind <- paste0("ind", czas)
s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), data = data1)
w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), data = data1)
pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)
p = pv
pv
ggsurv <- ggsurvplot(s,
                     # pval=round(pv,3), pval = "p < 0.001",
                     legend.title = "",
                     #title = "Krzywe przezycia dla protokolow", 
                     xlab = "Years", xlim = c(0,10), break.time.by = 1, ylab = paste0("p", czas),
                     #conf.int = TRUE, 
                     risk.table = TRUE, risk.table.pos="out", 
                     risk.table.height = 0.4, risk.table.y.text = FALSE, 
                     risk.table.col='strata', data=data1,  pval = ifelse(dig > 0, paste0('p=',formatC(p, digits=dig, format="f")), 'p<0.001'),
                     legend =  c(0.8, 0.2),
                     legend.labs = legend.labs[[i]],
                     ggtheme =  theme_survminer())

ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))


type = 1

n5 = ggsurv$data.survtable[(ggsurv$data.survtable$time == 5), ]$n.risk
surv_table_5 = ggsurv$data.survplot[(ggsurv$data.survplot$n.risk == n5[type]),][type,]
surv1 = round(surv_table_5$surv,2)
st_err1 = round(surv_table_5$std.err,2)

type = 2

n5 = ggsurv$data.survtable[(ggsurv$data.survtable$time == 5), ]$n.risk
surv_table_5 = ggsurv$data.survplot[(ggsurv$data.survplot$n.risk == n5[type]),][type,]
surv2 = round(surv_table_5$surv,2)
st_err2 = round(surv_table_5$std.err,2)

n1 = summary(s)$table[1,'records']
n2 = summary(s)$table[2,'records']

e1 = summary(s)$table[1,'events']
e2 = summary(s)$table[2,'events']

text1 = paste0(surv1,', SE=', st_err1,' (n=',n1, ', ', e1,' events)')
text2 = paste0(surv2,', SE=', st_err2,' (n=',n2, ', ', e2,' events)')



ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold")) +
  annotate("text", x=5, y=surv1+a, label=text1, size=5) +
  annotate("text", x=5, y=surv2+b, label=text2, size=5)


ggsurv[1]

ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
ggsurv[1]
dev.off()



png(plik, width = 820, height = 480)
ggsurv[1]
dev.off()
