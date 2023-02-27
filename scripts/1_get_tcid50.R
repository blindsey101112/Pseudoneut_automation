library(readxl)
library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)

###Read in .xlsx file
input.files <- list.files("./input/")

luminescence.files <- lapply(input.files,
                             function(lum.file){
                               run.num <- str_extract(lum.file, "\\d+.xlsx")
                               run.num <- sub(".xlsx", "", run.num)
                               
                               if (is.na(run.num)){print(paste0("Not analysing (Wrong format): ", lum.file))}
                               
                               tibble(file=lum.file, run=run.num)
                             }) %>% bind_rows()

luminescence.files <- luminescence.files[!is.na(luminescence.files$run),]

##get sample names
suppressMessages(sample.name.df <- readxl::read_excel(paste0("./input/ADOO_v1_plate_layouts.xlsx"), col_names = F))
plate.ref.index <- lapply(seq(1, ncol(sample.name.df)), 
                          function(col.n){
                            ref.row <- grep("Plate Ref:",unlist(sample.name.df[,col.n]))
                            ref.index <-c()
                            if (length(ref.row)>0){ref.index <- paste0(col.n, "-", ref.row)}
                            return(tibble(ref=ref.index))
                          }) %>% bind_rows()

sample.ids <- lapply(plate.ref.index$ref , 
                     function(plate.id.ref){
                       index.col <- as.numeric(str_split(plate.id.ref, pattern = "-")[[1]][1])
                       index.row <- as.numeric(str_split(plate.id.ref, pattern = "-")[[1]][2])
                       plate.ref <- unlist(sample.name.df[index.row, index.col-1])
                       plate.samples <- unname(unlist(sample.name.df[seq(index.row, index.row+6),index.col+1]))
                       tibble(
                         plate=plate.ref,
                         samples= plate.samples
                       )
                     }) %>% bind_rows()
sample.ids$plate <- sub("_", "-", sample.ids$plate)

##loop through files
loop.log <- lapply(seq(1, nrow(luminescence.files)),
                   function(row.n){
                     file.name <- luminescence.files$file[row.n]
                     run.name <-  str_extract(file.name, "\\d+-\\d+.xlsx")
                     run.name <- sub(".xlsx", "", run.name)
                     print(paste0("Analysing ", run.name))
                     run.samples <- sample.ids$samples[sample.ids$plate==run.name]
                     
                     dir.create(paste0("./output/run",  run.name))
                     dir.create(paste0("./output/run",  run.name, "/norm_virus"))
                     dir.create(paste0("./output/run",  run.name, "/norm_neg"))
                     #excel.file <- readxl::read_excel(paste0("./input/test/B1-1_Luminescence Quick Read 2022.11.09 11_54_37 AM1 1-1.xlsx"), sheet=2)
                     suppressMessages(excel.file <- readxl::read_excel(paste0("./input/", file.name), sheet=2))
                     excel.file <- excel.file[c(9:16), c(6:17)]
                     cell.only.mean <- mean(as.numeric(excel.file[8,])[1:4])
                     negcontrol.mean <- mean(c(as.numeric(excel.file[7,])[3:4]))
                     poscontrol.mean <- mean(c(as.numeric(excel.file[7,])[1:2]))
                     pseudovirus.only.mean <- mean(c(as.numeric(excel.file[7,])[11:12], as.numeric(excel.file[8,])[11:12]))
                     ###Get sample values
                     sample1 <- excel.file[c(1:6),c(1,2)] %>% mutate(sample=1) %>% lapply(as.numeric)
                     sample2 <- excel.file[c(1:6),c(3,4)] %>% mutate(sample=2) %>% lapply(as.numeric)
                     sample3 <- excel.file[c(1:6),c(5,6)] %>% mutate(sample=3) %>% lapply(as.numeric)
                     sample4 <- excel.file[c(1:6),c(7,8)] %>% mutate(sample=4) %>% lapply(as.numeric)
                     sample5 <- excel.file[c(1:6),c(9,10)] %>% mutate(sample=5) %>% lapply(as.numeric)
                     sample6 <- excel.file[c(1:6),c(11,12)] %>% mutate(sample=6)%>% lapply(as.numeric)
                     suppressMessages(sample.values <- bind_rows(sample1, sample2, sample3, sample4, sample5, sample6))
                     colnames(sample.values) <- c("rep1", "rep2", "sample")
                     sample7 <- as.data.frame(t(excel.file[c(7:8),c(5:10)])) %>% mutate(sample=7) 
                     colnames(sample7) <- c("rep1", "rep2", "sample")
                     rownames(sample7) <- NULL
                     sample7 <- as_tibble(sample7) %>% lapply(as.numeric)
                     sample.values <- bind_rows(sample.values, sample7) 
                     sample.values$sample.id <- rep(run.samples, each=6)
                     sample.values$conc <- rep(c(1.30, 1.78, 2.26, 2.73, 3.21, 3.69), 7)
                     rm(sample1, sample2, sample3, sample4, sample5, sample6, sample7)
                     ###get mean of replicates
                     sample.values$mean.value <- (sample.values$rep1 + sample.values$rep2)/2 
                     ###Normalise neutralisation
                     sample.values$normalised.values.pseudo <- 100* ((sample.values$mean.value- pseudovirus.only.mean)/(cell.only.mean- pseudovirus.only.mean))
                     sample.values$normalised.values.neg <- 100* ((sample.values$mean.value- negcontrol.mean)/(cell.only.mean- negcontrol.mean))
                     qc.fail <- F
                     #Cell only should be <1x10^4
                     cell.only.qc <- cell.only.mean <1000
                     #Positive controls should be <5x10^4
                     pos.control.qc <- poscontrol.mean <50000
                     #Negative controls and PV only should be >2x10^6
                     pv.only.control.qc <- pseudovirus.only.mean >2000000
                     neg.control.qc <- negcontrol.mean >2000000
                     
                     if(cell.only.qc + pos.control.qc + pv.only.control.qc + neg.control.qc >0){qc.fail<-T}
                      
                     
                     #Duplicates should be highlighted if max/min>2
                     sample.values <- lapply(seq(1,7), function(i){
                       sample.df <- sample.values %>% filter(sample==i)
                       sample.df$replicate.qc.fail <- F 
                       for (rep in seq(1,6)){
                         rep.values <- unlist(as.vector(sample.df[rep,c(1,2)]))
                         max.min <- max(rep.values)/min(rep.values)
                         if(max.min>2){
                           sample.df$replicate.qc.fail[rep]<-T
                           }
                       }
                       sample.df
                     }) %>% bind_rows()
                     
                     qc.df <- tibble(qc.fail = qc.fail,
                                     cell.only.qc = cell.only.qc,
                                     pos.control.qc = pos.control.qc,
                                     pv.only.control.qc = pv.only.control.qc,
                                     neg.control.qc  = neg.control.qc)
                       
                     write.csv(qc.df, paste0("./output/run",  run.name, "/run", run.name, "_qc_report.csv"), quote = F, row.names = F)
                     write.csv(sample.values, paste0("./output/run",  run.name, "/run", run.name, "_input_data_summary.csv"), quote = F, row.names = F)
                     
                     ### Get TCID50. Look at screenshot. Constrained hill slope.
                     ic50.df <- lapply(seq(1,7), function(i){
                       sample.df <- sample.values %>% filter(sample==i)
                       
                       predict.conc <- seq(0.1, 7, by=0.1)
                       
                       ##pseduovirus normalisation
                       pv.nls.try <- try(
                         {  
                           #define function
                           f <- function(conc,LogIC50,HillSlope) {100/(1+10^((LogIC50-conc)*HillSlope))}
                           #linearise to ger staring values
                           fm0 <- nls(log(normalised.values.pseudo) ~ log(f(conc, LogIC50, HillSlope)), data=sample.df, start = c(LogIC50=2, HillSlope=-1))
                           #run model
                           nls.model.pseudo <- nls(normalised.values.pseudo~f(conc, LogIC50, HillSlope), data=sample.df, 
                                                   start=coef(fm0), 
                                                   algorithm = "port",
                                                   lower = c(LogIC50=0.1, HillSlope=-3),
                                                   upper=c(LogIC50=5, HillSlope=0))
                           
                           #nls.model.pseudo <- nls(normalised.values.pseudo~100/(1+10^((LogIC50-conc)*HillSlope)), data=sample.df, start=c(LogIC50=2.5, HillSlope=-1))
                           nls.pseudo.pars <- nls.model.pseudo$m$getPars()
                           nls.pseudo.prediction <- 100/(1 + 10^((nls.pseudo.pars[1] - predict.conc) * nls.pseudo.pars[2]))
                           nls.pseudo.logic50 <- nls.pseudo.pars[1]
                           nls.pseudo.ic50=10^nls.pseudo.pars[1]
                         }, silent=T)
                       
                       
                       if (class(pv.nls.try) != "try-error"){ 
                         output.plot <- ggplot() +
                           geom_line(aes(x=predict.conc, y=nls.pseudo.prediction), colour="red") +
                           geom_point(aes(x=sample.df$conc, y= sample.df$normalised.values.pseudo),size=3) +
                           geom_point(aes(x=nls.pseudo.logic50, y=50), shape=4, size=3, colour="blue") +
                           ylab("% Neutralisation") +
                           xlab("Log10Conc") +
                           ggtitle(paste0("Sample ", sample.df$sample.id[sample.df$sample==i][1])) +
                           labs(colour="Model")
                         ggsave(paste0("./output/run",  run.name, "/norm_virus/", sample.df$sample.id[sample.df$sample==i][1],".jpeg"))
                         }
                       if (class(pv.nls.try) == "try-error"){ 
                         output.plot <- ggplot() +
                           geom_point(aes(x=sample.df$conc, y= sample.df$normalised.values.pseudo),size=3) +
                           ylab("% Neutralisation") +
                           xlab("Log10Conc") +
                           ggtitle(paste0("Sample ", sample.df$sample.id[sample.df$sample==i])) +
                           labs(colour="Model")
                         nls.pseudo.ic50 <- NA 
                         nls.pseudo.logic50 <- NA
                         ggsave(paste0("./output/run",  run.name, "/norm_virus/", sample.df$sample.id[sample.df$sample==i][1],".jpeg"))
                       }
                       pv.ic50 <- tibble(sample=i, sample.id = sample.df$sample.id[sample.df$sample==i][1], ic50=nls.pseudo.ic50, logic50=nls.pseudo.logic50, norm="pv")
                       
                       
                       
                       ##negative control normalisation
                       neg.nls.try <- try(
                         {  
                           
                           #define function
                           f <- function(conc,LogIC50,HillSlope) {100/(1+10^((LogIC50-conc)*HillSlope))}
                           #linearise to ger staring values
                           fm0 <- nls(log(normalised.values.neg) ~ log(f(conc, LogIC50, HillSlope)), data=sample.df, start = c(LogIC50=2, HillSlope=-1))
                           #run model
                           nls.model.neg <- nls(normalised.values.neg~f(conc, LogIC50, HillSlope), data=sample.df, 
                                                   start=coef(fm0), 
                                                   algorithm = "port",
                                                   lower = c(LogIC50=0.1, HillSlope=-3),
                                                   upper=c(LogIC50=5, HillSlope=0))
                           
                           
                           #nls.model.neg <- nls(normalised.values.neg~100/(1+10^((LogIC50-conc)*HillSlope)), data=sample.df, start=c(LogIC50=2.5, HillSlope=-1))
                           nls.neg.pars <- nls.model.neg$m$getPars()
                           nls.neg.prediction <- 100/(1 + 10^((nls.neg.pars[1] - predict.conc) * nls.neg.pars[2]))
                           nls.neg.logic50 <- nls.neg.pars[1]
                           nls.neg.ic50=10^nls.neg.pars[1]
                         }, silent=T)
                       
                       
                       if (class(neg.nls.try) != "try-error"){ 
                         output.plot <- ggplot() +
                           geom_line(aes(x=predict.conc, y=nls.neg.prediction), colour="red") +
                           geom_point(aes(x=sample.df$conc, y= sample.df$normalised.values.neg),size=3) +
                           geom_point(aes(x=nls.neg.logic50, y=50), shape=4, size=3, colour="blue") +
                           ylab("% Neutralisation") +
                           xlab("Log10Conc") +
                           ggtitle(paste0("Sample ", sample.df$sample.id[sample.df$sample==i][1])) +
                           labs(colour="Model")
                         ggsave(paste0("./output/run",  run.name, "/norm_neg/", sample.df$sample.id[sample.df$sample==i][1],".jpeg"))
                       }
                       if (class(neg.nls.try) == "try-error"){ 
                         output.plot <- ggplot() +
                           geom_point(aes(x=sample.df$conc, y= sample.df$normalised.values.neg),size=3) +
                           ylab("% Neutralisation") +
                           xlab("Log10Conc") +
                           ggtitle(paste0("Sample ", sample.df$sample.id[sample.df$sample==i])) +
                           labs(colour="Model")
                         nls.neg.ic50 <- NA 
                         nls.neg.logic50 <- NA
                         ggsave(paste0("./output/run",  run.name, "/norm_neg/", sample.df$sample.id[sample.df$sample==i][1],".jpeg"))
                       }
                       neg.ic50 <- tibble(sample=i,sample.id=sample.df$sample.id[sample.df$sample==i][1], ic50=nls.neg.ic50, logic50=nls.neg.logic50, norm="neg")
                       
                       ic50 <- bind_rows(pv.ic50, neg.ic50)
                       
                       return(ic50)
                     }) %>% bind_rows()
                     
                     pv.ic50 <- ic50.df %>% filter(norm=="pv") %>% select(-norm)
                     write.csv(pv.ic50, paste0("./output/run",  run.name, "/norm_virus/ic50.csv"), quote = F, row.names = F)
                     
                     neg.ic50 <- ic50.df %>% filter(norm=="neg") %>% select(-norm)
                     write.csv(neg.ic50, paste0("./output/run",  run.name, "/norm_neg/ic50.csv"), quote = F, row.names = F)
                     
                     print(paste0("Number ", row.n, " of ", nrow(luminescence.files), " complete"))
                     }) 






















