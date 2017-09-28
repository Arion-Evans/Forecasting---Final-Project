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


#  objects to store results (monthly)
MASE.monthly <- array()
shapiro.monthly <- array()
models.monthly <- list()
model.info.monthly <- list()
MASE.frc.monthly <- array()

# monthly loop
freq <- 12
for(i in 1:nrow(monthly)){
  
  # converting to ts objects leaving at least 2 observations in the 5% series
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
  
  # model fitting
  fit <- expSmooth(series95)
  models.monthly[[i]] <- fit[[1]]
  model.info.monthly[[i]] <- fit[[2]]
  MASE.monthly[i] <- model.info.monthly[[i]][,"MASE"]
  shapiro.monthly[i] <- model.info.monthly[[i]][,"Shapiro-Wilks"]
  
  # forecasts
  forecasts <- forecast(models.monthly[[i]], h = length(series05))$mean
  
  # reversing transformations
  if(model.info.monthly[[i]][,"Series Used"] == "Original"){
    forecasts <- forecasts - model.info.monthly[[i]][,"Added Value"]
  }
  if(model.info.monthly[[i]][,"Series Used"] == "Box-Cox"){
    forecasts <- invBoxCox(forecasts, model.info.monthly[[i]][,"Lambda"]) - model.info.monthly[[i]][,"Added Value"]
  }
  if(model.info.monthly[[i]][,"Series Used"] == "1st Difference"){
    comb <- ts.union(diff(series95) , forecasts - model.info.monthly[[i]][,"Added Value"])
    ts.combined.diff  = pmin(comb[,1], comb[,2], na.rm = TRUE)
    back.series = diffinv(ts.combined.diff, xi = series95[1],lag = 1)
    forecasts = window(back.series,start = start(series05))
  }
  if(model.info.monthly[[i]][,"Series Used"] == "Seasonal Difference"){
    comb <- ts.union(diff(series95, lag = freq) , forecasts - model.info.monthly[[i]][,"Added Value"])
    ts.combined.diff  = pmin(comb[,1], comb[,2], na.rm = TRUE)
    back.series = diffinv(ts.combined.diff, xi = series95[1:freq],lag = freq)
    forecasts = window(back.series,start = start(series05))
  }
  
  MASE.frc.monthly[i] <- as.numeric(MASE.custom(series05, forecasts))
  
  
}

#  objects to store results (quarterly)
MASE.quarterly <- array()
shapiro.quarterly <- array()
models.quarterly <- list()
model.info.quarterly <- list()
MASE.frc.quarterly <- array()

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
  
  # model fitting
  fit <- expSmooth(series95)
  models.quarterly[[i]] <- fit[[1]]
  model.info.quarterly[[i]] <- fit[[2]]
  MASE.quarterly[i] <- model.info.quarterly[[i]][,"MASE"]
  shapiro.quarterly[i] <- model.info.quarterly[[i]][,"Shapiro-Wilks"]
  
  # forecasts
  forecasts <- forecast(models.quarterly[[i]], h = length(series05))$mean
  
  # reversing transformations
  if(model.info.quarterly[[i]][,"Series Used"] == "Original"){
    forecasts <- forecasts - model.info.quarterly[[i]][,"Added Value"]
  }
  if(model.info.quarterly[[i]][,"Series Used"] == "Box-Cox"){
    forecasts <- invBoxCox(forecasts, model.info.quarterly[[i]][,"Lambda"]) - model.info.quarterly[[i]][,"Added Value"]
  }
  if(model.info.quarterly[[i]][,"Series Used"] == "1st Difference"){
    comb <- ts.union(diff(series95) , forecasts - model.info.quarterly[[i]][,"Added Value"])
    ts.combined.diff  = pmin(comb[,1], comb[,2], na.rm = TRUE)
    back.series = diffinv(ts.combined.diff, xi = series95[1],lag = 1)
    forecasts = window(back.series,start = start(series05))
  }
  if(model.info.quarterly[[i]][,"Series Used"] == "Seasonal Difference"){
    comb <- ts.union(diff(series95, lag = freq) , forecasts - model.info.quarterly[[i]][,"Added Value"])
    ts.combined.diff  = pmin(comb[,1], comb[,2], na.rm = TRUE)
    back.series = diffinv(ts.combined.diff, xi = series95[1:freq],lag = freq)
    forecasts = window(back.series,start = start(series05))
  }
  
  MASE.frc.quarterly[i] <- as.numeric(MASE.custom(series05, forecasts))
}

#  objects to store results (yearly)
MASE.yearly <- array()
shapiro.yearly <- array()
models.yearly <- list()
model.info.yearly <- list()
MASE.frc.yearly <- array()

# yearly loop
freq = 1
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
  
  # model fitting
  fit <- expSmooth(series95)
  models.yearly[[i]] <- fit[[1]]
  model.info.yearly[[i]] <- fit[[2]]
  MASE.yearly[i] <- model.info.yearly[[i]][,"MASE"]
  shapiro.yearly[i] <- model.info.yearly[[i]][,"Shapiro-Wilks"]
  
  # forecasts
  forecasts <- forecast(models.yearly[[i]], h = length(series05))$mean
  
  # reversing transformations
  if(model.info.yearly[[i]][,"Series Used"] == "Original"){
    forecasts <- forecasts - model.info.yearly[[i]][,"Added Value"]
  }
  if(model.info.yearly[[i]][,"Series Used"] == "Box-Cox"){
    forecasts <- invBoxCox(forecasts, model.info.yearly[[i]][,"Lambda"]) - model.info.yearly[[i]][,"Added Value"]
  }
  if(model.info.yearly[[i]][,"Series Used"] == "1st Difference"){
    comb <- ts.union(diff(series95) , forecasts - model.info.yearly[[i]][,"Added Value"])
    ts.combined.diff  = pmin(comb[,1], comb[,2], na.rm = TRUE)
    back.series = diffinv(ts.combined.diff, xi = series95[1],lag = 1)
    forecasts = window(back.series,start = start(series05))
  }
  if(model.info.yearly[[i]][,"Series Used"] == "Seasonal Difference"){
    comb <- ts.union(diff(series95, lag = freq) , forecasts - model.info.yearly[[i]][,"Added Value"])
    ts.combined.diff  = pmin(comb[,1], comb[,2], na.rm = TRUE)
    back.series = diffinv(ts.combined.diff, xi = series95[1:freq],lag = freq)
    forecasts = window(back.series,start = start(series05))
  }
  
  MASE.frc.yearly[i] <- as.numeric(MASE.custom(series05, forecasts))
}

# Results
results <- data.frame(matrix(NA, nrow = 3, ncol = 3), row.names = c("Monthly","Quarterly","Yearly"))
colnames(results) <- c("Mean Model MASE","Mean Forecast MASE","Mean Shapiro")
results[1,1] <- mean(MASE.monthly)
results[2,1] <- mean(MASE.quarterly)
results[3,1] <- mean(MASE.yearly)
results[1,2] <- mean(MASE.frc.monthly)
results[2,2] <- mean(MASE.frc.quarterly)
results[3,2] <- mean(MASE.frc.yearly)
results[1,3] <- mean(shapiro.monthly)
results[2,3] <- mean(shapiro.quarterly)
results[3,3] <- mean(shapiro.yearly)



