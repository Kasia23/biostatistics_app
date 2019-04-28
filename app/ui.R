#data1 <<- readRDS("C:\\Users\\User\\Desktop\\aplikacja_MRD\\data1.rds")
#zmienne_uv <- ifelse((data1[,"NOTCH1.Sanger"]==1 | data1[,"FBXW7.Sanger"]==1), 1, 0 )
#data1$PTEN.ABN.NOTCH1.FBXW7 <- ifelse(zmienne_uv==1, data1[,"PTEN.ABN"], NA)
#data1$Notch1Fbxw7_2gr <- zmienne_uv
#pten_plus <- data1[,"PTEN.ABN"] + 1
#pten_bez_mut <- ifelse(pten_plus==2, 0, pten_plus)
#data1$Notch1Fbxw7PTEN_2gr <- ifelse( (zmienne_uv==1) & (pten_bez_mut==1) , 1, 0)
#data1$PTEN.ABNSIL_TAL1 <- ifelse((data1[,"PTEN.ABN"]==1) | (data1[,"SIL.TAL"]==1), 1, 0)
#data1$gLoR_gHoR <- ifelse(zmienne_uv==1, pten_bez_mut, pten_bez_mut-pten_bez_mut)
#write.table(data1, file="C:\\Users\\User\\Desktop\\aplikacja_MRD\\data1.txt")
#data1 <- read.table("D:\\DEV\\aplikacja_mrd\\data1.txt")


#data1[which(data1["X.ID"]=='T1633'),]["risk.group..SR.IR.HR."] = 2
#data1[which(data1["X.ID"]=='T1701'),]["risk.group..SR.IR.HR."] = 2
#data1[which(data1["X.ID"]=='T0280'),]["indOS"] = 0
#data1[which(data1["X.ID"]=='T0280'),]["indEFS"] = 0

#data1 <<- read.table("C:\\Users\\katar\\OneDrive\\Desktop\\aplikacja_20180305\\data1.txt")
#kolumnaBW <- read.table("C:\\Users\\katar\\OneDrive\\Desktop\\aplikacja_20180305\\kolumnaBW.txt", header = TRUE)
#col_NA <- data.frame(rep(NA,24))
#colnames(col_NA) = "PTEN.mono.vs.biallelic.inactivation"
#kolumnaBW <- rbind(kolumnaBW, col_NA)

#data1$PTEN.WT <- ifelse((data1$PTEN.Sanger==0) & (data1$PTEN.MLPA==0) , 1, 0)
#data1$PTEN.zMUTbezDEL <- ifelse((data1$PTEN.Sanger==1) & (data1$PTEN.MLPA==0) , 1, NA)
#data1$PTEN.bezMUTzDEL <- ifelse((data1$PTEN.Sanger==0) & (data1$PTEN.MLPA==1) , 1, NA)
#data1$PTEN.zMUTbezDEL <- ifelse( data1$PTEN.WT == 1 , 0, data1$PTEN.zMUTbezDEL)
#data1$PTEN.bezMUTzDEL <- ifelse( data1$PTEN.WT == 1 , 0, data1$PTEN.bezMUTzDEL)

#nazwy <- colnames(data1)
#nazwy[23] <- "PTEN.MUT"
#nazwy[25] <- "PTEN.DEL"
#colnames(data1) <- nazwy

#write.table(data1,"C:\\Users\\katar\\OneDrive\\Desktop\\aplikacja_20180305\\data1.txt")
#write.table(data1,"D:\\DEV\\aplikacja_mrd\\data1.txt")

#data1 <<- cbind(data1, kolumnaBW)
#write.table(data1,"C:\\Users\\katar\\OneDrive\\Desktop\\aplikacja_20180305\\data1.txt")
data1 <<- read.table("data1.txt")
library("ggpubr")

#data1 <<- readRDS("data1.rds")
#Sys.setlocale("LC_ALL","Polish")
#pakiety---------#

library("data.table")
library("xlsx")
library("vcd")
library("gplots")
library("survival")
library("survminer")
library("ggplot2")
library("cowplot")
library("foreign")
library("rms")
library("mstate")
library("stringi")

colnames(data1)[89] <-  "PTEN.mutated.only"
colnames(data1)[90] <-  "PTEN.deleted.only"


navbarPage("Wstepna analiza - zaleznosci",
           
           
           
           tabPanel("Tablice kontyngencji",
                  
                    sidebarLayout(
                       sidebarPanel(
                          helpText("Wybierz gen i ceche, ktore chcesz porownac:"),
                          selectInput('gen1', 'mutacja', colnames(data1[ , c(20:53, 73:74, 83:87, 89:91)]), selectize = TRUE),
                          selectInput('cecha1', 'cecha', colnames(data1[ , c(54:60, 62:64, 5:7, 9:11, 13:19, 71:72)]), selectize = TRUE)
                       ),
                       
                       mainPanel(
                          plotOutput("kont"),
                          verbatimTextOutput('tabela'),
                          helpText("Procenty:"),
                          verbatimTextOutput("prop")
                          
                       )
                       
                    )
           ),
           
           tabPanel("Krzywe przezycia",
                    
                    sidebarLayout(
                       sidebarPanel(
                          helpText("Wybierz gen i czas, ktore chcesz analizowac:"),
                          selectInput('gen2', 'mutacja', colnames(data1[ , c(20:53, 73:74, 5:7, 9:11, 13:19, 83:87, 89:91)]), selectize = TRUE),
                          selectInput('czas1', 'czas', colnames(data1[ , c(75,77,79,81)]), selectize = TRUE),
                          selectInput('val', 'Risk.group.protocol', c(1,2), selectize = TRUE)
                       ),
                       
                       mainPanel(
                          
                          tabsetPanel(
                             tabPanel("istotnosc podzialu na protokoly", 
                                      plotOutput("protokoly", height = "750px"),
                                      verbatimTextOutput("pv_protokoly"),
                                      verbatimTextOutput('KM_protokoly')),
                             
                             tabPanel("krzywe przezycia bez stratyfikacji", 
                                      plotOutput("bez_strata", height = "750px"),
                                      verbatimTextOutput("pv_bez_strata"),
                                      verbatimTextOutput("KM_bez_strata")),
                             
                             tabPanel("krzywe przezycia z podzialem na protokoly", 
                                      plotOutput("protokoly1", height = "750px"),
                                      verbatimTextOutput("pv_protokoly1"),
                                      verbatimTextOutput("KM_protokoly1"),
                                      plotOutput("protokoly2", height = "750px"),
                                      verbatimTextOutput("pv_protokoly2"),
                                      verbatimTextOutput("KM_protokoly2")),
                             
                             tabPanel("krzywe przezycia z podzialem na protokoly i stratyfikacja wzgledem grupy ryzyka", 
                                      plotOutput("grupa_ryzyka1", height = "750px"),
                                      verbatimTextOutput("pv_grupa_ryzyka1"), 
                                      verbatimTextOutput("KM_grupa_ryzyka1"), 
                                      plotOutput("grupa_ryzyka2", height = "750px"),
                                      verbatimTextOutput("pv_grupa_ryzyka2"),
                                      verbatimTextOutput("KM_grupa_ryzyka2")),
                             
                             tabPanel("krzywe przezycia ze stratyfikacja wzgledem protokolow i grup ryzyka", 
                                      plotOutput("wszystko", height = "750px"),
                                      verbatimTextOutput("pv_wszystko"),
                                      verbatimTextOutput("KM_wszystko")),
                             
                             tabPanel("krzywe przezycia IR i HR", 
                                      plotOutput("risk_group1", height = "750px"),
                                      verbatimTextOutput("pv_risk_group1"),
                                      verbatimTextOutput("KM_risk_group1"),
                                      plotOutput("risk_group2", height = "750px"),
                                      verbatimTextOutput("pv_risk_group2"),
                                      verbatimTextOutput("KM_risk_group2")),
                             
                             tabPanel("krzywe przezycia dla mrd", 
                                      plotOutput("MRD_15", height = "750px"),
                                      verbatimTextOutput("pv_MRD_15"),
                                      verbatimTextOutput("KM_MRD_15"),
                                      plotOutput("MRD_33", height = "750px"),
                                      verbatimTextOutput("pv_MRD_33"),
                                      verbatimTextOutput("KM_MRD_33"),
                                      plotOutput("MRD_78", height = "750px"),
                                      verbatimTextOutput("pv_MRD_78"),
                                      verbatimTextOutput("KM_MRD_78"))
                             
                             
                          )
                          
                       )
             
                              )
           ),
           
           tabPanel("model PH",
                    
                    sidebarLayout(
                       sidebarPanel(
                        
                           helpText("Wybierz ceche i czas, ktore chcesz analizowac:"),
                           selectInput('gen3', 'cecha', colnames(data1[ , c(20:53, 73:74,54:60, 62:64, 5:7, 9:11, 13:19, 71:72, 83:87, 89:91)]), selectize = TRUE),
                           selectInput('czas2', 'czas', colnames(data1[ , c(75,77,79,81)]), selectize = TRUE)
                           

                       ),
                       
                       mainPanel(
                          
                          tabsetPanel(
                            
                            tabPanel("Model Coxa dla wszystkich zmiennych", tableOutput("coxPH_all")),
                            tabPanel("Model Coxa dla wybranych zmiennych", tableOutput("coxPH")),
                            tabPanel("Model Coxa dla pojedynczych zmiennych", tableOutput("coxPH_stop"))
                             
                              )
                          
                       )
                    )
           ),
           
           
           tabPanel("Dodatkowe statystyki",
                    
                    sidebarLayout(
                      sidebarPanel(
                        helpText("Wybierz gen i czas, ktore chcesz analizowac:"),
                        selectInput('gen4', 'mutacja', colnames(data1[ , c(20:53, 73:74, 5:7, 9:11, 13:19, 83:87, 89:91)]), selectize = TRUE),
                        selectInput('czas3', 'czas', colnames(data1[ , c(75,77,79,81)]), selectize = TRUE)
                        #selectInput('val2', 'Risk.group.protocol', c(1,2), selectize = TRUE)
                      ),
                       
                       mainPanel(
                          
                          tabsetPanel(
                            
                            tabPanel("pCIR", 
                                     plotOutput("pCIR"),
                                     tableOutput("pv_pCIR"),
                                     verbatimTextOutput("pv_pCIR_table"),
                                     plotOutput("pCIR1"),
                                     tableOutput("pv_pCIR1"),
                                     verbatimTextOutput("pv_pCIR1_table"),
                                     plotOutput("pCIR2"),
                                     tableOutput("pv_pCIR2"),
                            verbatimTextOutput("pv_pCIR2_table")),
                            
                            #tabPanel("pCIR - uwzglednienie konkurujacych ryzyk", 
                            #         plotOutput("pCIR_multi"),
                            #         tableOutput("gray_pCIR"),
                            #         plotOutput("pCIR1_multi"),
                            #         tableOutput("gray_pCIR1"),
                            #         plotOutput("pCIR2_multi"),
                            #         tableOutput("gray_pCIR2"))
                            
                            tabPanel("statystyki - median follow-up (dni)", tableOutput("mediana")) #,
                        #    tabPanel("Mediana dla pacjentow leczonych zgodnie z protokolem 2009", tableOutput("mediana2009"))
                               #           helpText("model coxa z wyborem zmiennych"),
                  #            helpText("Dodatkowo w aplikacji nalezy zmienic:"),
                  #           helpText("- szerokosc paneli wyswietlajacych krzywe przezycia,"),
                  #           helpText("- legendy do krzywych przezycia,"),
                  #           helpText("- w przypadku LMO2.RAG2.1 przy podziale na protokoly tylko jedna krzywa - obsluga wyjatkow dla testu logrank")
                          )
                          
                       )
                    )
           )
)