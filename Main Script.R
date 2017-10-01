library(forecast)


# read data
monthly <- read.csv("Data/Monthly.csv")
quarterly <- read.csv("Data/Quarterly.csv")
yearly <- read.csv("Data/Yearly.csv")   


#  objects to store results (monthly)
models.monthly <- list()
model.info.monthly <- list()
forecasts.monthly <- list()
series95.monthly <- list()
series05.monthly <- list()

# array to store overall model fitting info
overall.info <- array(NA,
                      dim = c(303,11), 
                      dimnames = list(NULL,c("Series Category","Season Frequency","Model Fitted","Series Used",
                                             "Model MASE","Forecast MASE","Shapiro-Wilks","Ljung-Box","AIC","AICc","BIC")))
count <- 1

# monthly loop
freq <- 12
for(i in 1:nrow(monthly)){
  
  # converting to ts objects leaving at least 2 observations in the 5% series
  if((0.05*monthly$N[i]<3) && (monthly[i,monthly$N[i]+5] == monthly[i,monthly$N[i]+6])){
    series95.monthly[[i]] <- ts(as.vector(t(as.matrix(monthly[i,7:(monthly$N[i]+3)]))),
                   start = c(monthly$Starting.Year[i], monthly$Starting.Month[i]),
                   frequency = freq)
    series05.monthly[[i]] <- ts(as.vector(t(as.matrix(monthly[i,(monthly$N[i]+4):(monthly$N[i]+6)]))),
                   start=end(series95.monthly[[i]])+c(0,1), 
                   frequency = freq)
  }else if(0.05*monthly$N[i]<3){
    series95.monthly[[i]] <- ts(as.vector(t(as.matrix(monthly[i,7:(monthly$N[i]+4)]))),
                   start = c(monthly$Starting.Year[i], monthly$Starting.Month[i]),
                   frequency = freq)
    series05.monthly[[i]] <- ts(as.vector(t(as.matrix(monthly[i,(monthly$N[i]+5):(monthly$N[i]+6)]))),
                   start=end(series95.monthly[[i]])+c(0,1), 
                   frequency = freq)
  }else{
    series95.monthly[[i]] <- ts(as.vector(t(as.matrix(monthly[i,7:(round(0.95*monthly$N[i])+6)]))),
                   start=c(monthly$Starting.Year[i],monthly$Starting.Month[i]), 
                   frequency = freq)
    
    series05.monthly[[i]] <- ts(as.vector(t(as.matrix(monthly[i,(round(0.95*monthly$N[i])+7):(monthly$N[i]+6)]))),
                   start=end(series95.monthly[[i]])+c(0,1), 
                   frequency = freq)
  }
  
  # model fitting
  fit <- expSmooth(series95.monthly[[i]])
  models.monthly[[i]] <- fit[[1]]
  model.info.monthly[[i]] <- fit[[2]]
  overall.info[count,"Model MASE"] <- model.info.monthly[[i]][,"MASE"]
  overall.info[count,"Shapiro-Wilks"] <- model.info.monthly[[i]][,"Shapiro-Wilks"]
  overall.info[count,"Ljung-Box"] <- model.info.monthly[[i]][,"Ljung-Box"]
  
  # forecasts
  forecasts.monthly[[i]] <- forecast(models.monthly[[i]], h = length(series05.monthly[[i]]))$mean
  
  # reversing transformations
  if(model.info.monthly[[i]][,"Series Used"] == "Original"){
    forecasts.monthly[[i]] <- forecasts.monthly[[i]] - model.info.monthly[[i]][,"Added Value"]
  }
  if(model.info.monthly[[i]][,"Series Used"] == "Box-Cox"){
    forecasts.monthly[[i]] <- invBoxCox(forecasts.monthly[[i]], model.info.monthly[[i]][,"Lambda"]) - model.info.monthly[[i]][,"Added Value"]
  }
  if(model.info.monthly[[i]][,"Series Used"] == "1st Difference"){
    comb <- ts.union(diff(series95.monthly[[i]]) , forecasts.monthly[[i]] - model.info.monthly[[i]][,"Added Value"])
    ts.combined.diff  = pmin(comb[,1], comb[,2], na.rm = TRUE)
    back.series = diffinv(ts.combined.diff, xi = series95.monthly[[i]][1],lag = 1)
    forecasts.monthly[[i]] = window(back.series,start = start(series05.monthly[[i]]))
  }
  if(model.info.monthly[[i]][,"Series Used"] == "Seasonal Difference"){
    comb <- ts.union(diff(series95.monthly[[i]], lag = freq) , forecasts.monthly[[i]] - model.info.monthly[[i]][,"Added Value"])
    ts.combined.diff  = pmin(comb[,1], comb[,2], na.rm = TRUE)
    back.series = diffinv(ts.combined.diff, xi = series95.monthly[[i]][1:freq],lag = freq)
    forecasts.monthly[[i]] = window(back.series,start = start(series05.monthly[[i]]))
  }
  
  overall.info[count,"Forecast MASE"] <- as.numeric(MASE.custom(series05.monthly[[i]], forecasts.monthly[[i]]))
  
  # storing overall info 
  overall.info[count,"Series Category"] <- as.character(monthly[i,"Category"])
  overall.info[count,"Season Frequency"] <- freq
  overall.info[count,"Model Fitted"] <- models.monthly[[i]]$method
  overall.info[count,"Series Used"] <- as.character(model.info.monthly[[i]][,"Series Used"])
  overall.info[count,"AIC"] <- if(grepl("ETS",models.monthly[[i]]$method)){models.monthly[[i]]$aic}else{models.monthly[[i]]$model$aic}
  overall.info[count,"AICc"] <- if(grepl("ETS",models.monthly[[i]]$method)){models.monthly[[i]]$aicc}else{models.monthly[[i]]$model$aicc}
  overall.info[count,"BIC"] <- if(grepl("ETS",models.monthly[[i]]$method)){models.monthly[[i]]$bic}else{models.monthly[[i]]$model$bic}
  
  count <- count + 1
}


#  objects to store results (quarterly)
models.quarterly <- list()
model.info.quarterly <- list()
forecasts.quarterly <- list()
series95.quarterly <- list()
series05.quarterly <- list()

# quarterly loop
freq <- 4
for(i in 1:nrow(quarterly)){
  if((0.05*quarterly$N[i]<3) && (quarterly[i,quarterly$N[i]+5] == quarterly[i,quarterly$N[i]+6])){
    series95.quarterly[[i]] <- ts(as.vector(t(as.matrix(quarterly[i,7:(quarterly$N[i]+3)]))),
                                  start = c(quarterly$Starting.Year[i], quarterly$Starting.Month[i]),
                                  frequency = freq)
    series05.quarterly[[i]] <- ts(as.vector(t(as.matrix(quarterly[i,(quarterly$N[i]+4):(quarterly$N[i]+6)]))),
                                  start=end(series95.quarterly[[i]])+c(0,1), 
                                  frequency = freq)
  }else if(0.05*quarterly$N[i]<3){
    series95.quarterly[[i]] <- ts(as.vector(t(as.matrix(quarterly[i,7:(quarterly$N[i]+4)]))),
                                  start = c(quarterly$Starting.Year[i], quarterly$Starting.Quarter[i]),
                                  frequency = freq)
    series05.quarterly[[i]] <- ts(as.vector(t(as.matrix(quarterly[i,(quarterly$N[i]+5):(quarterly$N[i]+6)]))),
                                  start=end(series95.quarterly[[i]])+c(0,1), 
                                  frequency = freq)
  }else{
    series95.quarterly[[i]] <- ts(as.vector(t(as.matrix(quarterly[i,7:(round(0.95*quarterly$N[i])+6)]))),
                                  start=c(quarterly$Starting.Year[i],quarterly$Starting.Quarter[i]), 
                                  frequency = freq)
    
    series05.quarterly[[i]] <- ts(as.vector(t(as.matrix(quarterly[i,(round(0.95*quarterly$N[i])+7):(quarterly$N[i]+6)]))),
                                  start=end(series95.quarterly[[i]])+c(0,1), 
                                  frequency = freq)
  }
  
  # model fitting
  fit <- expSmooth(series95.quarterly[[i]])
  models.quarterly[[i]] <- fit[[1]]
  model.info.quarterly[[i]] <- fit[[2]]
  overall.info[count,"Model MASE"] <- model.info.quarterly[[i]][,"MASE"]
  overall.info[count,"Shapiro-Wilks"] <- model.info.quarterly[[i]][,"Shapiro-Wilks"]
  overall.info[count,"Ljung-Box"] <- model.info.quarterly[[i]][,"Ljung-Box"]
  
  # forecasts
  forecasts.quarterly[[i]] <- forecast(models.quarterly[[i]], h = length(series05.quarterly[[i]]))$mean
  
  # reversing transformations
  if(model.info.quarterly[[i]][,"Series Used"] == "Original"){
    forecasts.quarterly[[i]] <- forecasts.quarterly[[i]] - model.info.quarterly[[i]][,"Added Value"]
  }
  if(model.info.quarterly[[i]][,"Series Used"] == "Box-Cox"){
    forecasts.quarterly[[i]] <- invBoxCox(forecasts.quarterly[[i]], model.info.quarterly[[i]][,"Lambda"]) - model.info.quarterly[[i]][,"Added Value"]
  }
  if(model.info.quarterly[[i]][,"Series Used"] == "1st Difference"){
    comb <- ts.union(diff(series95.quarterly[[i]]) , forecasts.quarterly[[i]] - model.info.quarterly[[i]][,"Added Value"])
    ts.combined.diff  = pmin(comb[,1], comb[,2], na.rm = TRUE)
    back.series = diffinv(ts.combined.diff, xi = series95.quarterly[[i]][1],lag = 1)
    forecasts.quarterly[[i]] = window(back.series,start = start(series05.quarterly[[i]]))
  }
  if(model.info.quarterly[[i]][,"Series Used"] == "Seasonal Difference"){
    comb <- ts.union(diff(series95.quarterly[[i]], lag = freq) , forecasts.quarterly[[i]] - model.info.quarterly[[i]][,"Added Value"])
    ts.combined.diff  = pmin(comb[,1], comb[,2], na.rm = TRUE)
    back.series = diffinv(ts.combined.diff, xi = series95.quarterly[[i]][1:freq],lag = freq)
    forecasts.quarterly[[i]] = window(back.series,start = start(series05.quarterly[[i]]))
  }
  
  overall.info[count,"Forecast MASE"] <- as.numeric(MASE.custom(series05.quarterly[[i]], forecasts.quarterly[[i]]))
  
  # storing overall info 
  overall.info[count,"Series Category"] <- as.character(quarterly[i,"Category"])
  overall.info[count,"Season Frequency"] <- freq
  overall.info[count,"Model Fitted"] <- models.quarterly[[i]]$method
  overall.info[count,"Series Used"] <- as.character(model.info.quarterly[[i]][,"Series Used"])
  overall.info[count,"AIC"] <- if(grepl("ETS",models.quarterly[[i]]$method)){models.quarterly[[i]]$aic}else{models.quarterly[[i]]$model$aic}
  overall.info[count,"AICc"] <- if(grepl("ETS",models.quarterly[[i]]$method)){models.quarterly[[i]]$aicc}else{models.quarterly[[i]]$model$aicc}
  overall.info[count,"BIC"] <- if(grepl("ETS",models.quarterly[[i]]$method)){models.quarterly[[i]]$bic}else{models.quarterly[[i]]$model$bic}
  
  count <- count + 1
}


#  objects to store results (yearly)
models.yearly <- list()
model.info.yearly <- list()
forecasts.yearly <- list()
series95.yearly <- list()
series05.yearly <- list()

# yearly loop
freq = 1
for(i in 1:nrow(yearly)){
  if((0.05*yearly$N[i]<3) && (yearly[i,yearly$N[i]+5] == yearly[i,yearly$N[i]+6])){
    series95.yearly[[i]] <- ts(as.vector(t(as.matrix(yearly[i,7:(yearly$N[i]+3)]))),
                   start = c(yearly$Starting.Year[i], yearly$Starting.Month[i]),
                   frequency = freq)
    series05.yearly[[i]] <- ts(as.vector(t(as.matrix(yearly[i,(yearly$N[i]+4):(yearly$N[i]+6)]))),
                   start=end(series95.yearly[[i]])+c(0,1), 
                   frequency = freq)
  }else if(0.05*yearly$N[i]<3){
    series95.yearly[[i]] <- ts(as.vector(t(as.matrix(yearly[i,7:(yearly$N[i]+4)]))),
                   start = yearly$Starting.Year[i])
    series05.yearly[[i]] <- ts(as.vector(t(as.matrix(yearly[i,(yearly$N[i]+5):(yearly$N[i]+6)]))),
                   start=end(series95.yearly[[i]])+c(1,0))
  }else{
    series95.yearly[[i]] <- ts(as.vector(t(as.matrix(yearly[i,7:(round(0.95*yearly$N[i])+6)]))),
                   start=yearly$Starting.Year[i])
    
    series05.yearly[[i]] <- ts(as.vector(t(as.matrix(yearly[i,(round(0.95*yearly$N[i])+7):(yearly$N[i]+6)]))),
                   start=end(series95.yearly[[i]])+c(1,0))
  }
  
  # model fitting
  fit <- expSmooth(series95.yearly[[i]])
  models.yearly[[i]] <- fit[[1]]
  model.info.yearly[[i]] <- fit[[2]]
  overall.info[count,"Model MASE"] <- model.info.yearly[[i]][,"MASE"]
  overall.info[count,"Shapiro-Wilks"] <- model.info.yearly[[i]][,"Shapiro-Wilks"]
  overall.info[count,"Ljung-Box"] <- model.info.yearly[[i]][,"Ljung-Box"]
  
  # forecasts
  forecasts.yearly[[i]] <- tryCatch(forecast(models.yearly[[i]], h = length(series05.yearly[[i]]))$mean,
                                    error = function(cond){
                                      return(forecast(ets(series95.yearly[[i]], "ZZZ"), h = length(series05.yearly[[i]]))$mean)
                                      })
  
  # reversing transformations
  if(model.info.yearly[[i]][,"Series Used"] == "Original"){
    forecasts.yearly[[i]] <- forecasts.yearly[[i]] - model.info.yearly[[i]][,"Added Value"]
  }
  if(model.info.yearly[[i]][,"Series Used"] == "Box-Cox"){
    forecasts.yearly[[i]] <- invBoxCox(forecasts.yearly[[i]], model.info.yearly[[i]][,"Lambda"]) - model.info.yearly[[i]][,"Added Value"]
  }
  if(model.info.yearly[[i]][,"Series Used"] == "1st Difference"){
    comb <- ts.union(diff(series95.yearly[[i]]) , forecasts.yearly[[i]] - model.info.yearly[[i]][,"Added Value"])
    ts.combined.diff  = pmin(comb[,1], comb[,2], na.rm = TRUE)
    back.series = diffinv(ts.combined.diff, xi = series95.yearly[[i]][1],lag = 1)
    forecasts.yearly[[i]] = window(back.series,start = start(series05.yearly[[i]]))
  }
  if(model.info.yearly[[i]][,"Series Used"] == "Seasonal Difference"){
    comb <- ts.union(diff(series95.yearly[[i]], lag = freq) , forecasts.yearly[[i]] - model.info.yearly[[i]][,"Added Value"])
    ts.combined.diff  = pmin(comb[,1], comb[,2], na.rm = TRUE)
    back.series = diffinv(ts.combined.diff, xi = series95.yearly[[i]][1:freq],lag = freq)
    forecasts.yearly[[i]] = window(back.series,start = start(series05.yearly[[i]]))
  }
  
  overall.info[count,"Forecast MASE"] <- as.numeric(MASE.custom(series05.yearly[[i]], forecasts.yearly[[i]]))
  
  # storing overall info 
  overall.info[count,"Series Category"] <- as.character(yearly[i,"Category"])
  overall.info[count,"Season Frequency"] <- freq
  overall.info[count,"Model Fitted"] <- models.yearly[[i]]$method
  overall.info[count,"Series Used"] <- as.character(model.info.yearly[[i]][,"Series Used"])
  overall.info[count,"AIC"] <- if(grepl("ETS",models.yearly[[i]]$method)){models.yearly[[i]]$aic}else{models.yearly[[i]]$model$aic}
  overall.info[count,"AICc"] <- if(grepl("ETS",models.yearly[[i]]$method)){models.yearly[[i]]$aicc}else{models.yearly[[i]]$model$aicc}
  overall.info[count,"BIC"] <- if(grepl("ETS",models.yearly[[i]]$method)){models.yearly[[i]]$bic}else{models.yearly[[i]]$model$bic}
  
  count <- count + 1
}

# Results
overall.info <- as.data.frame(overall.info)
overall.info[, c(2,5:11)] <- sapply(overall.info[, c(2,5:11)], as.character)
overall.info[, c(2,5:11)] <- sapply(overall.info[, c(2,5:11)], as.numeric)
colnames(overall.info) <- c("Series Category","Season Frequency","Model Fitted","Series Used","Model MASE",
                            "Forecast MASE","Shapiro-Wilks","Ljung-Box","AIC","AICc","BIC")

# counting non-normal and correlated residuals
overall.info$`Non-Normal Std Residuals` <- ifelse(overall.info$`Shapiro-Wilks` < 0.05, 1, 0)
overall.info$`Correlated Std Residuals` <- ifelse(overall.info$`Ljung-Box` < 0.05, 1, 0)

# filtering by series season frequency
overall.info.monthly <- overall.info[overall.info[,"Season Frequency"] == 12,]
overall.info.quarterly <- overall.info[overall.info[,"Season Frequency"] == 4,]
overall.info.yearly <- overall.info[overall.info[,"Season Frequency"] == 1,]


results <- data.frame(matrix(NA, nrow = 3, ncol = 6), row.names = c("Monthly","Quarterly","Yearly"))
colnames(results) <- c("Mean Model MASE","Mean Forecast MASE","Mean Shapiro","Mean Ljung",
                       "Non-Normal Std Residuals","Correlated Std Residuals")
results[1,1] <- mean(overall.info.monthly[,"Model MASE"])
results[2,1] <- mean(overall.info.quarterly[,"Model MASE"])
results[3,1] <- mean(overall.info.yearly[,"Model MASE"])
results[1,2] <- mean(overall.info.monthly[,"Forecast MASE"])
results[2,2] <- mean(overall.info.quarterly[,"Forecast MASE"])
results[3,2] <- mean(overall.info.yearly[,"Forecast MASE"])
results[1,3] <- mean(overall.info.monthly[,"Shapiro-Wilks"])
results[2,3] <- mean(overall.info.quarterly[,"Shapiro-Wilks"])
results[3,3] <- mean(overall.info.yearly[,"Shapiro-Wilks"])
results[1,4] <- mean(overall.info.monthly[,"Ljung-Box"])
results[2,4] <- mean(overall.info.quarterly[,"Ljung-Box"])
results[3,4] <- mean(overall.info.yearly[,"Ljung-Box"])
results[1,5] <- mean(sum(overall.info.monthly[,"Non-Normal Std Residuals"]))
results[2,5] <- mean(sum(overall.info.quarterly[,"Non-Normal Std Residuals"]))
results[3,5] <- mean(sum(overall.info.yearly[,"Non-Normal Std Residuals"]))
results[1,6] <- mean(sum(overall.info.monthly[,"Correlated Std Residuals"]))
results[2,6] <- mean(sum(overall.info.quarterly[,"Correlated Std Residuals"]))
results[3,6] <- mean(sum(overall.info.quarterly[,"Correlated Std Residuals"]))


