library(readxl)
library(tidyr)
library(dplyr)

###Read in .xlsx file
input.files <- list.files("./input/")

luminescence.file <- input.files[1]
excel.file <- readxl::read_excel(paste0("./input/", luminescence.file), sheet=2) 
excel.file <- excel.file[c(9:16), c(6:17)]
rm(luminescence.file)

###Average cell only (H1-4) and pseudo virus only (G/H11-12)
cell.only.mean <- mean(as.numeric(excel.file[8,])[1:4])
pseudovirus.only.mean <- mean(c(as.numeric(excel.file[7,])[11:12], as.numeric(excel.file[8,])[11:12]))

###Get sample values
sample1 <- excel.file[c(1:6),c(1,2)] %>% select(rep1 = 1, rep2 = 2) %>% mutate(sample=1, value=(rep1+rep2)/2)
sample2 <- excel.file[c(1:6),c(3,4)] %>% select(rep1 = 1, rep2 = 2) %>% mutate(sample=2, value=(rep1+rep2)/2)
sample3 <- excel.file[c(1:6),c(5,6)] %>% select(rep1 = 1, rep2 = 2) %>% mutate(sample=3, value=(rep1+rep2)/2)
sample4 <- excel.file[c(1:6),c(7,8)] %>% select(rep1 = 1, rep2 = 2) %>% mutate(sample=4, value=(rep1+rep2)/2)
sample5 <- excel.file[c(1:6),c(9,10)] %>% select(rep1 = 1, rep2 = 2) %>% mutate(sample=5, value=(rep1+rep2)/2)
sample6 <- excel.file[c(1:6),c(11,12)] %>% select(rep1 = 1, rep2 = 2) %>% mutate(sample=6, value=(rep1+rep2)/2)
sample7 <- as.data.frame(t(excel.file[c(7:8),c(5:10)])) %>% select(rep1 = 1, rep2 = 2) %>% mutate(sample=7, value=(rep1+rep2)/2)
rownames(sample7) <- NULL

sample.values <- rbind(sample1, sample2, sample3, sample4, sample5, sample6, sample7)
rm(sample1, sample2, sample3, sample4, sample5, sample6, sample7)

###Linear regresion of pseudo virus and cell only & regress out average read out
normalise.df <- tibble(well=c("pseudovirus", "cell"),
                       value=c(pseudovirus.only.mean, cell.only.mean),
                       percent=c(0,100))

model <- lm(percent ~ value, data = normalise.df)
sample.values$normalised.values <- predict(model, newdata = sample.values)


### Get TCID50. Look at screenshot. Constrained hill slope. 


