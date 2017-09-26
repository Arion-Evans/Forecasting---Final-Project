library(forecast)
library(expsmooth)
library(TSA)
library(dynlm)
library(x12)
library(Hmisc)
library(car)
library(AER)
library(dLagM)
library(knitr) 
library(scales)


monthly <- read.csv("C:/Users/Arion/Desktop/Uni/Forecasting/Final Project/Monthly.csv")
quarterly <- read.csv("C:/Users/Arion/Desktop/Uni/Forecasting/Final Project/Quarterly.csv")
yearly <- read.csv("C:/Users/Arion/Desktop/Uni/Forecasting/Final Project/Yearly.csv")   


# monthly loop
freq <- 12
for(i in 1:nrow(monthly)){
  if(0.05*monthly$N[i]<2){
    series95 <- ts(as.vector(t(as.matrix(monthly[i,7:(monthly$N[i]+4)]))),
                   start = c(monthly$Starting.Year[i], monthly$Starting.Month[i]),
                   frequency = freq)
    series05 <- ts(as.vector(t(as.matrix(monthly[i,(monthly$N[i]+5):(monthly$N[i]+6)]))),
                   start=end(series95)+c(0,1), 
                   frequency = freq)
  }else{
    series95 <- ts(as.vector(t(as.matrix(monthly[i,7:(round(0.95*monthly$N[i])+6)]))),
                   start=c(monthly$Starting.Year[i],monthly$Starting.Month[i]), 
                   frequency = freq)
    
    series05 <- ts(as.vector(t(as.matrix(monthly[i,(round(0.95*monthly$N[i])+7):(monthly$N[i]+6)]))),
                   start=end(series95)+c(0,1), 
                   frequency = freq)
  }
}

# quarterly loop
freq <- 4
for(i in 1:nrow(quarterly)){
  if(0.05*quarterly$N[i]<2){
    series95 <- ts(as.vector(t(as.matrix(quarterly[i,7:(quarterly$N[i]+4)]))),
                   start = c(quarterly$Starting.Year[i], quarterly$Starting.Quarter[i]),
                   frequency = freq)
    series05 <- ts(as.vector(t(as.matrix(quarterly[i,(quarterly$N[i]+5):(quarterly$N[i]+6)]))),
                   start=end(series95)+c(0,1), 
                   frequency = freq)
  }else{
    series95 <- ts(as.vector(t(as.matrix(quarterly[i,7:(round(0.95*quarterly$N[i])+6)]))),
                   start=c(quarterly$Starting.Year[i],quarterly$Starting.Quarter[i]), 
                   frequency = freq)
    
    series05 <- ts(as.vector(t(as.matrix(quarterly[i,(round(0.95*quarterly$N[i])+7):(quarterly$N[i]+6)]))),
                   start=end(series95)+c(0,1), 
                   frequency = freq)
  }
}

# yearly loop

for(i in 1:nrow(yearly)){
  if(0.05*yearly$N[i]<2){
    series95 <- ts(as.vector(t(as.matrix(yearly[i,7:(yearly$N[i]+4)]))),
                   start = yearly$Starting.Year[i])
    series05 <- ts(as.vector(t(as.matrix(yearly[i,(yearly$N[i]+5):(yearly$N[i]+6)]))),
                   start=end(series95)+c(1,0))
  }else{
    series95 <- ts(as.vector(t(as.matrix(yearly[i,7:(round(0.95*yearly$N[i])+6)]))),
                   start=yearly$Starting.Year[i])
    
    series05 <- ts(as.vector(t(as.matrix(yearly[i,(round(0.95*yearly$N[i])+7):(yearly$N[i]+6)]))),
                   start=end(series95)+c(1,0))
  }
}




# Testing
series95 <- ts(as.vector(t(as.matrix(quarterly[45,7:(quarterly$N[45]+4)]))),
               start = c(quarterly$Starting.Year[45], quarterly$Starting.Quarter[45]),
               frequency = 4)
series05 <- ts(as.vector(t(as.matrix(quarterly[45,(quarterly$N[45]+5):(quarterly$N[45]+6)]))),
               start=c((quarterly$Starting.Year[45]+trunc(length(series95)/4)),
                       (quarterly$Starting.Quarter[45]+(length(series95)%%4)%%4)), 
               frequency = 4)

series95 <- ts(as.vector(t(as.matrix(monthly[1,7:(round(0.95*monthly$N[1])+6)]))),
               start=c(monthly$Starting.Year[1],monthly$Starting.Month[1]), 
               frequency = 12)

series05 <- ts(as.vector(t(as.matrix(monthly[1,(round(0.95*monthly$N[1])+7):(monthly$N[1]+6)]))),
               start=end(series95)+c(0,8), 
               frequency = 12)

series95 <- ts(as.vector(t(as.matrix(quarterly[1,7:(round(0.95*quarterly$N[1])+6)]))),
               start=c(quarterly$Starting.Year[1],quarterly$Starting.Quarter[1]), 
               frequency = freq)

series05 <- ts(as.vector(t(as.matrix(quarterly[1,(round(0.95*quarterly$N[1])+7):(quarterly$N[1]+6)]))),
               start=end(series95)+c(0,1), 
               frequency = freq)
series04 <- ts(as.vector(t(as.matrix(quarterly[1,(round(0.95*quarterly$N[1])+7):(quarterly$N[1]+6)]))),
               start=end(series05)+c(0,1), 
               frequency = freq)

series95 <- ts(as.vector(t(as.matrix(yearly[1,7:(round(0.95*yearly$N[1])+6)]))),
               start=yearly$Starting.Year[1])

series05 <- ts(as.vector(t(as.matrix(yearly[1,(round(0.95*yearly$N[1])+7):(yearly$N[1]+6)]))),
               start=end(series95)+c(1,0))

