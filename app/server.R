#data1 <<- readRDS("data1.rds")
#data1 <<- readRDS("C:\\Users\\User\\Desktop\\aplikacja_MRD\\data1.rds")
data1 <<- read.table("data1.txt")

legend.labs <- list( 
                     "PTEN.ABN" = c( 'PTEN.WT', 'PTEN.ABN'),
                     "PTEN.MUT" = c('PTEN.WT', 'PTEN.MUT'),
                     "PTEN.DEL" = c( 'PTEN.WT', 'PTEN.DEL'),
                     "Risk.group..MRD.d15" = c('standard risk group (<0,1% blast cells)', 'intermediate risk group (0.1 to <10%  blast cells)', 'high risk group (>10% blast cells)'),
                     "Risk.group.MRD.d33." = c('standard risk group (<0,1% blast cells)', 'intermediate risk group (0.1 to <10%  blast cells)', 'high risk group (>10% blast cells)'),
                     "Risk.group.MRD.d78" = c('standard risk group (<0,1% blast cells)', 'intermediate risk group (0.1 to <10%  blast cells)'),
                     "PTEN.mutated.only" = c( 'PTEN.WT', 'PTEN.MUT'),
                     "PTEN.deleted.only" = c('PTEN.WT', 'PTEN.DEL'),
                     "gLoR_gHoR"= c('gLoR', 'gHoR'),
                     "Notch1Fbxw7PTEN_2gr" = c('NOTCH1/FBXW7mut&PTEN.WT', 'NOTCH1/FBXW7/PTEN.WT'),
                     "NOTCH1.Sanger" = c('NOTCH1.WT', 'NOTCH1.MUT'),
                     "Notch1Fbxw7_2gr" = c('NOTCH1/FBXW7.WT', 'NOTCH1/FBXW7.MUT'),
                     "PTEN.mono.vs.biallelic.inactivation" = c("PTEN.WT", "PTENbiallelic", "PTENmonoallelic")
)


#Sys.setlocale("LC_ALL","Polish")
library("data.table")
library("xlsx")
library("vcd")
library("gplots")
library("survminer")
library("ggplot2")
library("vcd")
library("cowplot")
library("foreign")
library("rms")
library("mstate")
library("stringi")

colnames(data1)[89] <-  "PTEN.mutated.only"
colnames(data1)[90] <-  "PTEN.deleted.only"

function(input, output, session) {
   
   output$kont <- renderPlot({
      
      i <- input$gen1
      j <- input$cecha1
      
      s = xtabs(~ data1[,i] + data1[,j]) 
      balloonplot(s, 
                  xlab = names(data1[i]), 
                  ylab = names(data1[j]),
                  dotcolor = "gray", main = "Licznosc")
   })
   
   
   output$prop <- renderPrint({ 
      
      i <- input$gen1
      j <- input$cecha1
      
      s = xtabs(~ data1[,i] + data1[,j], drop.unused.levels=TRUE) 
      names(dimnames(s)) <- c(i, j)
      pt <- round(prop.table(s,1)*100,2)
      pt
      
      
   })
   
   
   
   output$szans <- renderPrint({
      
      i <- input$gen1
      j <- input$cecha1
      
      s = xtabs(~ data1[,i] + data1[,j], data=data1) 
      names(dimnames(s)) <- c(i, j)
      pt <- prop.table(s,2)
      st=pt/(1-pt) 
      
      print(st)
      
   })
   
   
   output$tabela <- renderPrint({
  #    i <- 42
  #    j <- 57
      i <- input$gen1
      j <- input$cecha1
      
      t <- table(data1[,i], data1[,j])

      
      #liczba wszystkich pacjentow branych pod uwage
      liczba <- sum(t)
      
      
      #cmh - test warstwowy wzgledem protokolu
      t2 <- table(data1[,i], data1[,j], data1$ALL.IC.BFM.2002.vs..2009)
      mt <- mantelhaen.test(t2)
      pmt <- mt$p.value
      omt <- mt$estimate
      met2 <- mt$method
      
      if( dim(table(data1[[i]])) == 2 & dim(table(data1[[j]])) == 2 ){
         
         #test Fishera:
         ft <- fisher.test(t, alternative = "two.sided", conf.level = 0.95)
         met <- ft$method
         
         #P
         p <- ft$p.value
         
         #OR 
         or <- summary(oddsratio(t, log=FALSE))[1]
      
         #CI
         ci <- paste0(round(confint(oddsratio(t, log=FALSE))[1],2), " - ", round(confint(oddsratio(t, log=FALSE))[2],2))
         
         #CI dla testu C-M-H
         l.CI_cmh <- mt$conf.int[1]
         r.CI_cmh <- mt$conf.int[2]
         ci_cmh <- paste0(round(l.CI_cmh,2), " - ", round(r.CI_cmh,2))
         
         #OR dla protokolu 1:
         table1 <- xtabs(~ data1[,i] + data1[,j], subset = data1$ALL.IC.BFM.2002.vs..2009 == 1)
         or1 <- summary(oddsratio(table1, log=FALSE))[1]
         ci1 <- paste0(round(confint(oddsratio(table1, log=FALSE))[1],2), " - ", round(confint(oddsratio(table1, log=FALSE))[2],2))
         
         #OR dla protokolu 2:
         table2 <- xtabs(~ data1[,i] + data1[,j], subset = data1$ALL.IC.BFM.2002.vs..2009 == 2)
         or2 <- summary(oddsratio(table2, log=FALSE))[1]
         ci2 <- paste0(round(confint(oddsratio(table2, log=FALSE))[1],2), " - ", round(confint(oddsratio(table2, log=FALSE))[2],2))
         
         
         
         wynik <- c(colnames(data1[i]), colnames(data1[j]), liczba, round(or,2), ci, 
                    formatC(p,digits=4, format="f"), met, 
                    formatC(pmt, digits=4,format="f"), round(omt,2), ci_cmh,
                    round(or1,2), ci1, round(or2,2), ci2)
         wynik <- data.frame(wynik)
         row.names(wynik) <- c("gen:", "cecha:", "ilosc:", "OR:", "CI:", 
                               "p.value:", "metoda:", 
                               "p.value (podzial na protokoly):", "OR (podzial na protokoly):", "CI (podzial na protokoly):",
                               "OR dla protokolu 1:", "CI dla protokolu 1:","OR dla protokolu 2:", "CI dla protokolu 2:")
         
      }else{
         
         #test Chi^2
         cht <- chisq.test(t)
         met <- cht$method
         
         
         #P
         p <- cht$p.value
         
         
         wynik <- c(colnames(data1[i]), colnames(data1[j]), liczba, 
                    formatC(p,digits=4, format="f"), met , 
                    formatC(pmt, digits=4, format="f"))
         wynik <- data.frame(wynik)
         
         row.names(wynik) <- c("gen:", "cecha:", "ilosc:", 
                               "p.value:", "metoda:" ,
                               "p.value (podzial na protokoly):")
         
         
      }
      
      print(wynik)
      
   })
   
  
   output$protokoly <- renderPlot({
      i <<-"PTEN.ABN"
      czas <<- "OS"
      val <<- 0
     
      i <- input$gen2
      czas <- input$czas1
      
      #print(czas)
      #print(val)
      #print(data1[which(data1["X.ID"]=='T1633'),]["risk.group..SR.IR.HR."])
      #print(data1[which(data1["X.ID"]=='T1701'),]["risk.group..SR.IR.HR."])
      #data2 <- data1[!(is.na(data1$PTEN.ABN) | is.na(data1$OS)),]
      
      #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~ ALL.IC.BFM.2002.vs..2009"))
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
                 risk.table.col='strata', data=data1, pval = round(pv,2),
                 legend =  c(0.8, 0.2),
                 legend.labs = c( 'ALL IC-BFM 2002', 'ALL IC-BFM 2009'),
                 ggtheme =  theme_survminer())
      
      ggsurv$plot <- ggsurv$plot + 
      theme(legend.text = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
            axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
            axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
            axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))
      
      ggsurv
   
   }, height = 700)
   
   output$pv_protokoly <- renderPrint({
      i <- input$gen2
      czas <- input$czas1
      
      #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~ ALL.IC.BFM.2002.vs..2009"))
      s <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~ ALL.IC.BFM.2002.vs..2009")), data = data1)
      #print(paste0("p.value: ", round(1 - pchisq(s$chisq, length(s$n) - 1),5) ))
      pv = round(1 - pchisq(s$chisq, length(s$n) - 1),5)

      #hazard ratio: https://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
      HR = (s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
      up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      df_hr = data.frame("p.value" = pv, "HR" = round(HR,2), "low95" = round(low95,2), "up95" = round(up95,2))
      print(df_hr)
      
   })
   
   
   output$KM_protokoly <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~ ALL.IC.BFM.2002.vs..2009"))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~ ALL.IC.BFM.2002.vs..2009")), data = data1)
     summary(s)
     
   })
   
   output$bez_strata <- renderPlot({
      
      i <- input$gen2
      czas <- input$czas1
      
      #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
      s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), data = data1)
      
      w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), data = data1)
      pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)
      
      
      ggsurvplot(s, 
                 # pval=round(pv,2),
                 legend.title = "",
                 #title = "Krzywe przezycia bez stratyfikacji", 
                 xlab = "Years", break.time.by = 1, ylab = paste0("p", czas),
                 risk.table = TRUE, risk.table.pos="out", 
                 risk.table.height = 0.4, risk.table.y.text = FALSE, 
                 risk.table.col='strata',
                 data=data1,  legend =  c(0.8, 0.2), legend.labs = legend.labs[[i]], ggtheme =  theme_survminer())
      
   }, height = 700)
   
   
   output$pv_bez_strata <- renderPrint({
      i <- input$gen2
      czas <- input$czas1
      
      #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
      s <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), data = data1)
      #print(paste0("p.value: ", round(1 - pchisq(s$chisq, length(s$n) - 1),5) ))
      
      pv = round(1 - pchisq(s$chisq, length(s$n) - 1),5)
      
      #hazard ratio: https://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
      HR = (s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
      up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      df_hr = data.frame("p.value" = pv, "HR" = round(HR,2), "low95" = round(low95,2), "up95" = round(up95,2))
      print(df_hr)
      
   })
   
   
   output$KM_bez_strata <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), data = data1)
     print(summary(s))
     
   })
   
   
   
   output$protokoly1 <- renderPlot({
      
      i <- input$gen2
      czas <- input$czas1
      
      #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
      s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset = ALL.IC.BFM.2002.vs..2009==1, data = data1)
      
      w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset = ALL.IC.BFM.2002.vs..2009==1, data = data1)
      pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)
      
      ggsurvplot(s, 
                 # pval=round(pv,2),
                 legend.title = "",
                 #title = paste0("Protokol 2002: ",i), 
                 xlab = "Years", xlim = c(0,10), 
                 break.time.by = 1,
                 ylab = paste0("p", czas), 
                 risk.table = TRUE, risk.table.pos="out", 
                 risk.table.height = 0.4, risk.table.y.text = FALSE, 
                 risk.table.col='strata', data=data1,  legend =  c(0.8, 0.2), legend.labs = legend.labs[[i]],  ggtheme =  theme_survminer())
   
   }, height = 700)
   
   
   
   output$pv_protokoly1 <- renderPrint({
      i <- input$gen2
      czas <- input$czas1
      
      #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
      s <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset = ALL.IC.BFM.2002.vs..2009==1, data = data1)
      #print(paste0("p.value: ", round(1 - pchisq(s$chisq, length(s$n) - 1),5)))
      
      pv = round(1 - pchisq(s$chisq, length(s$n) - 1),5)
      
      #hazard ratio: https://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
      HR = (s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
      up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      df_hr = data.frame("p.value" = pv, "HR" = round(HR,2), "low95" = round(low95,2), "up95" = round(up95,2))
      print(df_hr)
      
   })
   
   
   
   output$KM_protokoly1 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset = ALL.IC.BFM.2002.vs..2009==1, data = data1)
     print(summary(s))
     
   })
   
   
   output$protokoly2 <- renderPlot({
      
      i <- input$gen2
      czas <- input$czas1
      
      #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
      s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=ALL.IC.BFM.2002.vs..2009==2, data = data1)
      
      w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=ALL.IC.BFM.2002.vs..2009==2, data = data1)
      pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)

      
      ggsurvplot(s, 
                 # pval=round(pv,2),
                 legend.title = "",
                 #title = paste0("Protokol 2009: ",i), 
                 xlab = "Years", xlim = c(0,10), break.time.by = 1,
                 ylab = paste0("p", czas),
                 risk.table = TRUE, risk.table.pos="out", 
                 risk.table.height = 0.4, risk.table.y.text = FALSE, 
                 risk.table.col='strata', data=data1,  legend =  c(0.8, 0.2), legend.labs = legend.labs[[i]],  ggtheme =  theme_survminer())
      
   }, height = 700)
   
   
   output$pv_protokoly2 <- renderPrint({
      i <- input$gen2
      czas <- input$czas1
      
      #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
      s <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=ALL.IC.BFM.2002.vs..2009==2, data = data1)
      #print(paste0("p.value: ", round(1 - pchisq(s$chisq, length(s$n) - 1),5) ))
      
      pv = round(1 - pchisq(s$chisq, length(s$n) - 1),5)
      
      #hazard ratio: https://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
      HR = (s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
      up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      df_hr = data.frame("p.value" = pv, "HR" = round(HR,2), "low95" = round(low95,2), "up95" = round(up95,2))
      print(df_hr)
      
   })
   
   
   output$KM_protokoly2 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=ALL.IC.BFM.2002.vs..2009==2, data = data1)
     print(summary(s))
     
   })
   
   
   output$grupa_ryzyka1 <- renderPlot({
      
      i <- input$gen2
      czas <- input$czas1
      
      #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)"))
      s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)")), subset=ALL.IC.BFM.2002.vs..2009==1, data = data1)
      
      w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)")), subset=ALL.IC.BFM.2002.vs..2009==1, data = data1)
      pv <- round(1 - pchisq(w$chisq, length(w$n) - 1),5)
  
      
      ggsurvplot(s, 
                 # pval=round(pv,2),
                 legend.title = "",
                 #title = paste0("Protokol 2002: ",i), 
                 xlab = "Years", xlim = c(0,10), break.time.by = 1,
                 ylab = paste0("p", czas),
                 risk.table = TRUE, risk.table.pos="out", 
                 risk.table.height = 0.4, risk.table.y.text = FALSE, 
                 risk.table.col='strata',
                 data=data1,  legend =  c(0.8, 0.2),  ggtheme =  theme_survminer())
      
   }, height = 700)
   
   
   output$pv_grupa_ryzyka1 <- renderPrint({
      i <- input$gen2
      czas <- input$czas1
      
      #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)"))
      s <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)")), subset=ALL.IC.BFM.2002.vs..2009==1, data = data1)
      #print(paste0("p.value: ", round(1 - pchisq(s$chisq, length(s$n) - 1),5) ))

      pv = round(1 - pchisq(s$chisq, length(s$n) - 1),5)
      
      #hazard ratio: https://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
      HR = (s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
      up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      df_hr = data.frame("p.value" = pv, "HR" = round(HR,2), "low95" = round(low95,2), "up95" = round(up95,2))
      print(df_hr)
      
   })
   
   
   output$KM_grupa_ryzyka1 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)"))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)")), subset=ALL.IC.BFM.2002.vs..2009==1, data = data1)
     print(summary(s))
     
   })
   
   output$grupa_ryzyka2 <- renderPlot({
      
      i <- input$gen2
      czas <- input$czas1
      
      #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)"))
      s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)")), subset=ALL.IC.BFM.2002.vs..2009==2, data = data1)
      
      w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)")), subset=ALL.IC.BFM.2002.vs..2009==2, data = data1)
      pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)
      
      
      ggsurvplot(s, 
                 # pval=round(pv,2),
                 legend.title = "",
                 #title = paste0("Protokol 2009: ",i), 
                 xlab = "Years", xlim = c(0,10), break.time.by = 1,
                 ylab = paste0("p", czas), 
                 risk.table = TRUE, risk.table.pos="out", 
                 risk.table.height = 0.4, risk.table.y.text = FALSE, 
                 risk.table.col='strata', data=data1, legend =  c(0.8, 0.2), ggtheme =  theme_survminer())
      
   }, height = 700)
   
   
   output$pv_grupa_ryzyka2 <- renderPrint({
      i <- input$gen2
      czas <- input$czas1
      
      #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)"))
      s <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)")), subset=ALL.IC.BFM.2002.vs..2009==2, data = data1)
      #paste0("p.value: ", round(1 - pchisq(s$chisq, length(s$n) - 1),5) )
      pv = round(1 - pchisq(s$chisq, length(s$n) - 1),5)
      
      #hazard ratio: https://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
      HR = (s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
      up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      df_hr = data.frame("p.value" = pv, "HR" = round(HR,2), "low95" = round(low95,2), "up95" = round(up95,2))
      print(df_hr)
      
   })
   
   output$KM_grupa_ryzyka2 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)"))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(risk.group..SR.IR.HR.)")), subset=ALL.IC.BFM.2002.vs..2009==2, data = data1)
     print(summary(s))
     
   })
   
   
   output$wszystko <- renderPlot({
     i <- input$gen2
     czas <- input$czas1 
      #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(ALL.IC.BFM.2002.vs..2009) + strata(risk.group..SR.IR.HR.)"))
      s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(ALL.IC.BFM.2002.vs..2009) + strata(risk.group..SR.IR.HR.)")), data = data1)
      
      w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(ALL.IC.BFM.2002.vs..2009) + strata(risk.group..SR.IR.HR.)")), data = data1)
      pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)
      
      
      ggsurvplot(s, 
                 # pval=round(pv,2),
                 legend.title = "",
                 #title = i, ylab = paste0("p", czas),
                 xlab = "Years", xlim = c(0,10), break.time.by = 1, 
                 risk.table = TRUE, risk.table.pos="out", 
                 risk.table.height = 0.4, risk.table.y.text = FALSE,
                 data=data1,  legend =  c(0.8, 0.2),  ggtheme =  theme_survminer())
      
   }, height = 700)
   
   
   output$MRD_15 <- renderPlot({
     i <- input$gen2
     czas <- input$czas1
     val <- input$val
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group..MRD.d15==val, data = data1)
     
     w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group..MRD.d15==val, data = data1)
     pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)
     
     
     ggsurvplot(s,
                # pval=round(pv,2),
                legend.title = "",
                #title = paste0( i, ' d15'), ylab = paste0("p", czas),
                risk.table = TRUE, risk.table.pos="out", 
                risk.table.height = 0.4, risk.table.y.text = FALSE,
                xlab = "Years", xlim = c(0,10), break.time.by = 1,
                data=data1,  legend =  c(0.8, 0.2), legend.labs = legend.labs[[i]],  ggtheme =  theme_survminer())
     
   }, height = 700)
   

   output$MRD_33 <- renderPlot({
     i <- input$gen2
     czas <- input$czas1
     val <- input$val
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group.MRD.d33.==val, data = data1)
     
     w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group.MRD.d33.==val, data = data1)
     pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)
     
     
     ggsurvplot(s,
                # pval=round(pv,2),
                legend.title = "", 
                #title = paste0( i, ' d33'), ylab = paste0("p", czas),
                xlab = "Years", xlim = c(0,10), break.time.by = 1,
                risk.table = TRUE, risk.table.pos="out", 
                risk.table.height = 0.4, risk.table.y.text = FALSE,
                data=data1,  legend =  c(0.8, 0.2), legend.labs = legend.labs[[i]],  ggtheme =  theme_survminer())
     
   }, height = 700)
      
   
   output$MRD_78 <- renderPlot({
     i <- input$gen2
     czas <- input$czas1
     val <- input$val
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group.MRD.d78==val, data = data1)
     
     w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group.MRD.d78==val, data = data1)
     pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)
     
     ggsurvplot(s, 
                # pval=round(pv,2),
                legend.title = "",
                #title = paste0( i, ' d78'), 
                xlab = "Years", xlim = c(0,10), break.time.by = 1, 
                risk.table = TRUE, risk.table.pos="out", 
                risk.table.height = 0.4, risk.table.y.text = FALSE,
                ylab = paste0("p", czas),
                data=data1,  legend =  c(0.8, 0.2), legend.labs = legend.labs[[i]],  ggtheme =  theme_survminer())
     
   }, height = 700)
   
   
   
   output$pv_wszystko <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(ALL.IC.BFM.2002.vs..2009) + strata(risk.group..SR.IR.HR.)"))
      s <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(ALL.IC.BFM.2002.vs..2009) + strata(risk.group..SR.IR.HR.)")), data = data1)
      #print(paste0("p.value: ", round(1 - pchisq(s$chisq, length(s$n) - 1),5) ))
      
      pv = round(1 - pchisq(s$chisq, length(s$n) - 1),5)
      
      #hazard ratio: https://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
      HR = (s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
      up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
      df_hr = data.frame("p.value" = pv, "HR" = round(HR,2), "low95" = round(low95,2), "up95" = round(up95,2))
      print(df_hr)
      
   })
   
   
   
   output$risk_group1 <- renderPlot({
     
     i <- input$gen2
     czas <- input$czas1
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset = risk.group..SR.IR.HR.==1, data = data1)
     
     w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset = risk.group..SR.IR.HR.==1, data = data1)
     pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)
     
     ggsurvplot(s, 
                # pval=round(pv,2),
                legend.title = "",
                #title = paste0("Protokol 2002: ",i), 
                xlab = "Years", xlim = c(0,10), 
                break.time.by = 1,
                ylab = paste0("p", czas), 
                risk.table = TRUE, risk.table.pos="out", 
                risk.table.height = 0.4, risk.table.y.text = FALSE, 
                risk.table.col='strata', data=data1,  legend =  c(0.8, 0.2), legend.labs = legend.labs[[i]],  ggtheme =  theme_survminer())
     
   }, height = 700)
   
   
   
   output$pv_risk_group1 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset = risk.group..SR.IR.HR.==1, data = data1)
     #print(paste0("p.value: ", round(1 - pchisq(s$chisq, length(s$n) - 1),5)))
     
     pv = round(1 - pchisq(s$chisq, length(s$n) - 1),5)
     
     HR = (s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
     up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
     low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
     df_hr = data.frame("p.value" = pv, "HR" = round(HR,2), "low95" = round(low95,2), "up95" = round(up95,2))
     print(df_hr)
     
   })
   
   
   
   output$KM_risk_group1 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset = risk.group..SR.IR.HR.==1, data = data1)
     print(summary(s))
     
   })
   
   
   output$risk_group2 <- renderPlot({
     
     i <- input$gen2
     czas <- input$czas1
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=risk.group..SR.IR.HR.==2, data = data1)
     
     w <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=risk.group..SR.IR.HR.==2, data = data1)
     pv = round(1 - pchisq(w$chisq, length(w$n) - 1),5)
     
     
     ggsurvplot(s, 
                # pval=round(pv,2),
                legend.title = "",
                #title = paste0("Protokol 2009: ",i), 
                xlab = "Years", xlim = c(0,10), break.time.by = 1,
                ylab = paste0("p", czas),
                risk.table = TRUE, risk.table.pos="out", 
                risk.table.height = 0.4, risk.table.y.text = FALSE, 
                risk.table.col='strata', data=data1,  legend =  c(0.8, 0.2), legend.labs = legend.labs[[i]],  ggtheme =  theme_survminer())
     
   }, height = 700)
   
   
   output$pv_risk_group2 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=risk.group..SR.IR.HR.==2, data = data1)
     #print(paste0("p.value: ", round(1 - pchisq(s$chisq, length(s$n) - 1),5) ))
     
     pv = round(1 - pchisq(s$chisq, length(s$n) - 1),5)
     
     #hazard ratio: https://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
     HR = (s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
     up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
     low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
     df_hr = data.frame("p.value" = pv, "HR" = round(HR,2), "low95" = round(low95,2), "up95" = round(up95,2))
     print(df_hr)
     
   })
   
   
   output$KM_risk_group2 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=risk.group..SR.IR.HR.==2, data = data1)
     print(summary(s))
     
   })
   
   output$pv_MRD_15 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     val <- input$val
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group..MRD.d15==val, data = data1)
     #print(paste0("p.value: ", round(1 - pchisq(s$chisq, length(s$n) - 1),5) ))
     
     pv = round(1 - pchisq(s$chisq, length(s$n) - 1),5)
     
     #hazard ratio: https://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
     HR = (s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
     up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
     low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
     df_hr = data.frame("p.value" = pv, "HR" = round(HR,2), "low95" = round(low95,2), "up95" = round(up95,2))
     print(df_hr)
     
   })
   
   
   output$KM_MRD_15 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     val <- input$val
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group..MRD.d15==val, data = data1)
     print(summary(s))
     
     
   })
   
   
   output$pv_MRD_33 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     val <- input$val
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group.MRD.d33.==val, data = data1)
     #print(paste0("p.value: ", round(1 - pchisq(s$chisq, length(s$n) - 1),5) ))
     
     pv = round(1 - pchisq(s$chisq, length(s$n) - 1),5)
     
     #hazard ratio: https://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
     HR = (s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
     up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
     low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
     df_hr = data.frame("p.value" = pv, "HR" = round(HR,2), "low95" = round(low95,2), "up95" = round(up95,2))
     print(df_hr)
     
   })
   
   
   output$KM_MRD_33 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     val <- input$val
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group.MRD.d33.==val, data = data1)
     print(summary(s))
     
   })
   
   
   output$pv_MRD_78 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     val <- input$val
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survdiff(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group.MRD.d78==val, data = data1)
     #print(paste0("p.value: ", round(1 - pchisq(s$chisq, length(s$n) - 1),5) ))
     
     pv = round(1 - pchisq(s$chisq, length(s$n) - 1),5)
     
     #hazard ratio: https://stats.stackexchange.com/questions/124489/how-to-calculate-the-hr-and-95ci-using-the-log-rank-test-in-r
     HR = (s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
     up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
     low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
     df_hr = data.frame("p.value" = pv, "HR" = round(HR,2), "low95" = round(low95,2), "up95" = round(up95,2))
     print(df_hr)
     
   })
   
   
   output$KM_MRD_78 <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     val <- input$val
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i)), subset=Risk.group.MRD.d78==val, data = data1)
     print(summary(s))
     
   })
   
   
   output$KM_wszystko <- renderPrint({
     i <- input$gen2
     czas <- input$czas1
     
     #f <- as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(ALL.IC.BFM.2002.vs..2009) + strata(risk.group..SR.IR.HR.)"))
     s <- survfit(as.formula(paste0("Surv(",czas,"/365, ind",czas,") ~",i,"+ strata(ALL.IC.BFM.2002.vs..2009) + strata(risk.group..SR.IR.HR.)")), data = data1)
     print(summary(s))
     
   })
   
   
   output$coxPH_stop <- renderTable({
      i <- input$gen3
      czas <- input$czas2
     
      c <- coxph(as.formula(paste0("Surv(",czas,"/365, ind" , czas , ") ~" , i , "+ ALL.IC.BFM.2002.vs..2009 + risk.group..SR.IR.HR.")), data = data1)
      data.table("zmienna" = rownames(summary(c)$coefficients), round(summary(c)$coefficients,2))[,c(1,3,6)]
   })
   
   
   output$coxPH <- renderTable({
      czas <- input$czas2
      
      #vv <- names(data1[,c(20:23,25:30,32,34:53, 73:74)])
      #vv <- names(data1[,c(20:53, 73:74)])
      #vv <- names(data1[,c(20:30,32,31,33:53, 73:74)])
      vv <- c("Age.BFM","leukocytoza..BFM.","sex","PTEN.ABN", "PRED.resp", "Risk.group..MRD.d15")
      
      
      #f <- as.formula(stri_paste("Surv(OS, indOS)" , stri_paste(vv, collapse = " + "), sep= " ~ "))
      #f <-   as.formula("Surv(OS, indOS) ~   
      #   NOTCH1.Sanger.domeny + NOTCH1.Sanger + FBXW7.Sanger + PTEN.Sanger + 
      #   PTEN.MLPA.del + PTEN.MLPA + PTEN.ABN + Notch1.FBXW7 + Notch1.FBXW7.PTEN... + 
      #   WT1 + STAT5B + SIL.TAL.MLPA + SIL.TAL + LEF1.MLPA + LEF1 + 
      #   CASP8AP2.MLPA + CASP8AP2 + MYB.AMPL + EZH2.MLPA + EZH2 + 
      #   MLLT3.MLPA + MLLT3 + CDKN2A.MLPA + CDKN2A + CDKN2B.MLPA" )
   
      c <- coxph(as.formula(stri_paste(paste0("Surv(",czas,"/365, ind",czas,")"), stri_paste(vv, collapse = " + "), sep= " ~ ")), data = data1)
      data.table("zmienna" = rownames(summary(c)$coefficients), round(summary(c)$coefficients,3))[,c(1,3,4,6)]
   
      })
   
   
   output$coxPH_all <- renderTable({
     czas <- input$czas2

     vv <- c("ALL.IC.BFM.2002.vs..2009",
             "Risk.group..MRD.d15",
             "Risk.group.MRD.d33.",
             "Risk.group.MRD.d78",
             "NOTCH1.Sanger",
             "FBXW7.Sanger",
             "PTEN.ABN",
             "WT1",
             "STAT5B",
             "SIL.TAL",
             "LEF1",
             "CASP8AP2",
             "MYB.AMPL",
             "EZH2",
             "MLLT3",
             "CDKN2A",
             "CDKN2B",
             "NUP214.ABL.ampl",
             "LMO1.ampl",
             "LMO2.RAG2",
             "NF1",
             "SUZ12",
             "PTPN2",
             "PHF6",
             "Age.BFM",
             "sex",
             "leukocytoza..BFM.",
             "CNS.status",
             "risk.group..SR.IR.HR.",
             "PRED.resp",
             "BM.d15",
             "X.BM.d33"
             )
     
     c <- coxph(as.formula(stri_paste(paste0("Surv(",czas,"/365, ind",czas,")"), stri_paste(vv, collapse = " + "), sep= " ~ ")), data = data1)
     data.table("zmienna" = rownames(summary(c)$coefficients), round(summary(c)$coefficients,3))[,c(1,3,4,6)]
     
   })
   
   
   output$pCIR1 <- renderPlot({
      i <- input$gen4
      czas <- input$czas3
      ind <- paste0("ind",czas)
   
      fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action=na.omit, subset = data1$ALL.IC.BFM.2002.vs..2009 == 1)
      ggcompetingrisks(fit,  multiple_panels = FALSE, title = paste0("Cumulative incidence functions - ",i," protokol 2002"))
      }, height = 400)
   
   output$pCIR2 <- renderPlot({
      i <- input$gen4
      czas <- input$czas3
      ind <- paste0("ind",czas)
      
      fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action=na.omit, subset = data1$ALL.IC.BFM.2002.vs..2009 == 2)
      ggcompetingrisks(fit,  multiple_panels = FALSE, title = paste0("Cumulative incidence functions - ",i," protokol 2009"))
   }, height = 400)  
   
   
   output$pCIR <- renderPlot({
     i <- input$gen4
     czas <- input$czas3
     ind <- paste0("ind",czas)
     
     fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action = na.omit)
     ggcompetingrisks(fit,  multiple_panels = FALSE, title = paste0("Cumulative incidence functions - ",i))
   }, height = 400)  
   
   
   output$pv_pCIR <- renderTable({
     i <- input$gen4
     czas <- input$czas3
     ind <- paste0("ind",czas)
     
     fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action = na.omit)
     
     gray_pval <- fit$Tests[2] # p-values for comparing the subdistribution for each cause across groups (if the number of groups is >1).
     data.table("Gray test p-value: " = formatC(gray_pval, digits=4,format="f")) 
     
   })
   
   output$pv_pCIR1 <- renderTable({
     i <- input$gen4
     czas <- input$czas3
     ind <- paste0("ind",czas)
     
     fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action=na.omit, subset = data1$ALL.IC.BFM.2002.vs..2009 == 1)
     
     gray_pval <- fit$Tests[2] # p-values for comparing the subdistribution for each cause across groups (if the number of groups is >1).
     data.table("Gray test p-value: " = formatC(gray_pval, digits=4,format="f") )
     
   })
   
   output$pv_pCIR2 <- renderTable({
     i <- input$gen4
     czas <- input$czas3
     ind <- paste0("ind",czas)
     
     fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action=na.omit, subset = data1$ALL.IC.BFM.2002.vs..2009 == 2)
     gray_pval <- fit$Tests[2] # p-values for comparing the subdistribution for each cause across groups (if the number of groups is >1).
     data.table("Gray test p-value: " = formatC(gray_pval, digits=4,format="f") )
     
   })
   
   
   output$pv_pCIR_table <- renderPrint({
     i <- input$gen4
     czas <- input$czas3
     ind <- paste0("ind",czas)
     
     fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action=na.omit, subset = data1$ALL.IC.BFM.2002.vs..2009 == 1)
     print(fit)
   })
   
   output$pv_pCIR1_table <- renderPrint({
     i <- input$gen4
     czas <- input$czas3
     ind <- paste0("ind",czas)
     
     fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action=na.omit, subset = data1$ALL.IC.BFM.2002.vs..2009 == 1)
     print(fit)
   })
   
   output$pv_pCIR2_table <- renderPrint({
     i <- input$gen4
     czas <- input$czas3
     ind <- paste0("ind",czas)
     
     fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action=na.omit, subset = data1$ALL.IC.BFM.2002.vs..2009 == 2)
     print(fit)
   })
   
   output$zalozenia <- renderTable({
     
     
     data.table("Gray test p-value: " = formatC(gray_pval, digits=4,format="f") )
     
   })
   
   
   output$pCIR1_multi <- renderPlot({
     i <- input$gen4
     czas <- input$czas3
     ind <- paste0("ind",czas)
     
     fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action=na.omit, subset = data1$ALL.IC.BFM.2002.vs..2009 == 1)
     ggcompetingrisks(fit,  multiple_panels = FALSE, title = paste0("Cumulative incidence functions - ",i," protokol 2002"))
   }, height = 400)
   
   output$pCIR2_multi <- renderPlot({
     i <- input$gen4
     czas <- input$czas3
     ind <- paste0("ind",czas)
     
     fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action=na.omit, subset = data1$ALL.IC.BFM.2002.vs..2009 == 2)
     ggcompetingrisks(fit,  multiple_panels = FALSE, title = paste0("Cumulative incidence functions - ",i," protokol 2009"))
   }, height = 400)  
   
   
   output$pCIR_multi <- renderPlot({
     i <- input$gen4
     czas <- input$czas3
     ind <- paste0("ind",czas)
     
     fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action = na.omit)
     ggcompetingrisks(fit,  multiple_panels = FALSE, title = paste0("Cumulative incidence functions - ",i))
   }, height = 400)  
   
   
   output$gray_pCIR <- renderTable({
     i <- input$gen4
     czas <- input$czas3
     ind <- paste0("ind",czas)
     
     fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action = na.omit)
     
     gray_pval <- fit$Tests[2] # p-values for comparing the subdistribution for each cause across groups (if the number of groups is >1).
     data.table("Gray test p-value: " = formatC(gray_pval, digits=4,format="f") )
     
   })
   
   output$gray_pCIR1 <- renderTable({
     i <- input$gen4
     czas <- input$czas3
     ind <- paste0("ind",czas)
     
     fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action=na.omit, subset = data1$ALL.IC.BFM.2002.vs..2009 == 1)
     
     gray_pval <- fit$Tests[2] # p-values for comparing the subdistribution for each cause across groups (if the number of groups is >1).
     data.table("Gray test p-value: " = formatC(gray_pval, digits=4,format="f") )
     
   })
   
   output$gray_pCIR2 <- renderTable({
     i <- input$gen4
     czas <- input$czas3
     ind <- paste0("ind",czas)
     
     fit <- cmprsk::cuminc(ftime = data1[,czas]/365, fstatus = data1[,ind], group = data1[,i], na.action=na.omit, subset = data1$ALL.IC.BFM.2002.vs..2009 == 2)
     gray_pval <- fit$Tests[2] # p-values for comparing the subdistribution for each cause across groups (if the number of groups is >1).
     data.table("Gray test p-value: " = formatC(gray_pval, digits=4,format="f") )
     
   })
   
   
   output$mediana <- renderTable({
     
     data2 <- data1[data1[,"ALL.IC.BFM.2002.vs..2009"]==1,]
     data9 <- data1[data1[,"ALL.IC.BFM.2002.vs..2009"]==2,]
     
     stat2002 <- na.omit((data2[,"Last.follow.up..date."] - data2[,"Date.of.Dg"]))
     stat2009 <- na.omit((data9[,"Last.follow.up..date."] - data9[,"Date.of.Dg"]))
   
     max2002 <- max( stat2002 )     
     min2002 <- min( stat2002 )    
     median2002 <- median( stat2002 )     
     mean2002 <- mean( stat2002 )    
     
     max2009 <- max( stat2009 )    
     min2009 <- min( stat2009 )    
     median2009 <- median( stat2009 )     
     mean2009 <- mean( stat2009 ) 
     
     #hist(na.omit((data1["Last.follow.up..date."] - data1["Date.of.Dg"])[[1]])/365)
     
     data.table(" " = c("min", "max", "median", "mean"), 
                "2002" = c(min2002, max2002, median2002, mean2002) , 
                "2009" = c(min2009, max2009, median2009, mean2009))
   })
   
}
