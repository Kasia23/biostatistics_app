# czasy:
colnames(data1)[65:70]
#data diagnozy i data remisji

legend.labs <- list( 
  "PTEN.ABN" = c( 'PTEN.WT', 'PTEN.ABN'),
  "PTEN.MUT" = c('PTEN.WT', 'PTEN.MUT'),
  "PTEN.DEL" = c( 'PTEN.WT', 'PTEN.DEL'),
  "Risk.group..MRD.d15" = c('standard risk group\n (<0,1% blast cells)\n', 'intermediate risk group\n (0.1 to <10%  blast cells)\n', 'high risk group\n (>10% blast cells)'),
  "Risk.group.MRD.d33." = c('standard risk group\n (<0,1% blast cells)\n', 'intermediate risk group\n (0.1 to <10%  blast cells)\n', 'high risk group\n (>10% blast cells)'),
  "Risk.group.MRD.d78" = c('standard risk group\n (<0,1% blast cells)\n', 'intermediate risk group\n (0.1 to <10%  blast cells)'),
  "PTEN.mutated.only" = c('PTEN.WT', 'PTEN.MUT'),
  "PTEN.deleted.only" = c('PTEN.WT', 'PTEN.DEL'),
  "gLoR_gHoR"= c('gLoR', 'gHoR'),
  "Notch1Fbxw7PTEN_2gr" = c('NOTCH1/FBXW7mut&PTEN.WT', 'NOTCH1/FBXW7/PTEN.WT'),
  "NOTCH1.Sanger" = c('NOTCH1.WT', 'NOTCH1.MUT'),
  "Notch1Fbxw7_2gr" = c('NOTCH1/FBXW7.WT', 'NOTCH1/FBXW7.MUT'),
  "PTEN.mono.vs.biallelic.inactivation" = c("PTEN.WT", "PTEN.biallelic", "PTEN.monoallelic")
)

colnames(data1)[89] <-  "PTEN.mutated.only"
colnames(data1)[90] <-  "PTEN.deleted.only"


# czyli rozpatrujemy zdarzenia od wznowy/rozpoczêcia obserwacji do zdarzenia (uwzglêdniaj¹c zdarzenia konkurencyjne)
# jeœli nast¹pi³a œmieræ, to oznaczmy jako 2
# DFS - czas od uzyskania remisji do wyst¹pienia zdarzenia (tylko wznowa) lub dnia ostatniej obserwacji 
#if data œmierci -> 2
# je¿eli wiemy, ¿e œmieræ nast¹pi³a z innej przyczyny ni¿ nawrót choroby, to mo¿emy to zdarzenie traktowaæ jako ryzyko konkuruj¹ce.
# w przeciwnym razie powinniœmy œmieræ traktowaæ wy³¹cznie jako cenzurowanie - nie mamy informacji o tym, czy pacjent zmar³ z powodu nawrotu, którego nie wykryliœmy,
# czy powód by³ inny
#(czas do nawrotu 1 czas do œmierci 2)

######################################################################################################################
library("cr17")

data1$min_date <- ifelse(data1$date.of.relapse %in% NA, data1$Last.follow.up..date., data1$date.of.relapse)  
indCIR <- ifelse(data1$date.of.relapse %in% NA, 0, 1)
timeCIR <- data1$min_date - data1$Date.of.remission # liczymy od daty remisji
data1$date.of.death
indCIR[(data1$date.of.relapse %in% NA) & data1$date.of.death] <- 2
data1$CIR <- timeCIR
data1$indCIR <- indCIR
ind <- "indCIR"
czas <- "CIR"
live_frame <- data.frame(indCIR, as.Date(as.numeric(data1$date.of.death), origin = "1970-01-01"), as.Date(as.numeric(data1$date.of.relapse), origin = "1970-01-01"), as.Date(as.numeric(data1$Date.of.remission), origin="1970-01-01"), as.Date(as.numeric(data1$Date.of.Dg), origin="1970-01-01"))
result <- data.frame(indCIR, timeCIR)


#i <- "PTEN.ABN"
#i <- "PTEN.deleted.only"
#i <- "PTEN.mutated.only"
#i <- "PTEN.MUT"
i <- "PTEN.mono.vs.biallelic.inactivation"

plik <- paste0('pCIR_', i, ' ')
plik_tiff <- paste0('pCIR_', i, '.tiff')

if(i=="PTEN.mono.vs.biallelic.inactivation"){
  group <- factor(data1[,i],0:2,legend.labs[[i]])
  num <- 4
}else{
  group <- factor(data1[,i],0:1,legend.labs[[i]])
  num <- 3}

event <- factor(data1[,ind],0:2,c('no','relapse', 'death'))

cuminc <- cr17::fitCuminc(time = data1[,czas]/365, risk = event, group = group, cens = "no")

ci <- cuminc
plot(cuminc)
cens = "no"
low <- NULL
up <- NULL
est <- NULL
time <- NULL
group <- NULL
risk <- NULL
if (is.null(cens)) 
  cens <- as.character(risk[1])
timePoints <- 0:10
nrTests <- which(names(ci) == "Tests")
ci <- ci[-nrTests]
aggNames <- names(ci)[num:(nrTests-1)]
toPlot <- lapply(aggNames, function(i) data.frame(time = ci[[i]]$time, est = ci[[i]]$est, var = ci[[i]]$var, 
                                                  group = strsplit(i," ")[[1]][1], risk = strsplit(i, " ")[[1]][2]))
toPlot <- do.call(rbind, toPlot)
risks <- sort(unique(toPlot$risk))
groups <- sort(unique(toPlot$group))

toPlot$lowerBound <- sapply(1:nrow(toPlot), function(x) {
  est <- toPlot[x, "est"]
  var <- toPlot[x, "var"]
  exp(log(est) - 1.96 * sqrt(var)/est)
})
toPlot$upperBound <- sapply(1:nrow(toPlot), function(x) {
  est <- toPlot[x, "est"]
  var <- toPlot[x, "var"]
  exp(log(est) + 1.96 * sqrt(var)/est)
})

plot1 <- ggplot(data = toPlot, aes(time, est, color = group)) + 
                                          geom_step(size = 1) + theme_survminer() + 
                                          theme(legend.text = element_text(size = 14, face = "bold"),
                                           axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=.5,face="bold"),
                                           axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0,face="bold"),
                                           axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
                                           axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold")) +
                                          ggtitle('') + theme(plot.title = element_text(size = 13, 
                                                                   face = "bold", hjust = 0.5), legend.position = c(0.8,0.2)) + 
                                          scale_y_continuous('pCIR') + scale_x_continuous('Years', breaks = timePoints) + 
                                          coord_cartesian(xlim = range(timePoints)) + 
                                          theme(legend.title = element_text(size = 10, face = "bold")) + 
                                          scale_color_manual(values=c("#00BFC4", "#F8766D"))



plot1
group <- factor(data1[,i],0:2,legend.labs[[i]])
p_data <- plot1$data
gr1 = levels(p_data$group)[1]
gr2 = levels(p_data$group)[2]
gr3 = levels(p_data$group)[3]

text_data1 <- tail(p_data[(p_data$group == gr1 & p_data$time <= 5),],1)
text_data2 <- tail(p_data[(p_data$group == gr2 & p_data$time <= 5),],1)
if(!is.na(gr3)){text_data3 <- tail(p_data[(p_data$group == gr3 & p_data$time <= 5),],1)}

p_data$est
p_data$time


text_data <- tail(p_data[(p_data$group == gr1 & p_data$time <= 5),],1)
e <- sum(na.omit(event[group == gr1] == 'relapse'))
n <- sum(na.omit(event[group == gr1] == 'no')) + e
text1 <- paste0(round(text_data$est,2),', SE=', round(text_data$var,3),' (n=', n, ', ', e,' events)')

plot2 <- plot1 + annotate("text", x=5, y=round(text_data$est,2) + 0.02, label= text1, size=5)

text_data <- tail(p_data[(p_data$group == gr2 & p_data$time <= 5),],1)
e <- sum(na.omit(event[group == gr2] == 'relapse'))
n <- sum(na.omit(event[group == gr2] == 'no')) + e
text2 <- paste0(round(text_data$est,2),', SE=', round(text_data$var,3),' (n=', n, ', ', e,' events)')


test_v = round(cuminc$Tests['relapse',]['pv'], 4)
plot3 <- plot2 + annotate("text", x=5, y=round(text_data$est,2) + 0.015, label= text2, size =5) +
 annotate("text", x = 4, y = 0.04, label = paste0("p = ", test_v), size = 5)
plot3

text_data <- tail(p_data[(p_data$group == gr3 & p_data$time <= 5),],1)
e <- sum(na.omit(event[group == gr3] == 'relapse'))
n <- sum(na.omit(event[group == gr3] == 'no')) + e
text3 <- paste0(round(text_data$est,2),', SE=', round(text_data$var,3),' (n=', n, ', ', e,' events)')

plot4 <- plot3 + annotate("text", x=5, y=round(text_data$est,2) + 0.013, label= text3, size=5)
plot4


ppi <- 1200
tiff(plik_tiff, width = 14.7*ppi, height = 8.6*ppi, res = ppi,  compression = "lzw")
plot4
dev.off()


png(plik, width = 820, height = 480)
plot3
dev.off()

################################################ tu!

plot <- plotCuminc(ci = cuminc,
           cens = "no",
           ggtheme = theme_survminer(),
           legendtitle = "") + theme(legend.position = c(0.2, 0.8),
                                     legend.text = element_text(size = 14, face = "bold"))
##########################################################



