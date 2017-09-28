library(forecast)


# read data
monthly <- read.csv("Data/Monthly.csv")
quarterly <- read.csv("Data/Quarterly.csv")
yearly <- read.csv("Data/Yearly.csv")   


#  objects to store results (monthly)
MASE.monthly <- array()
shapiro.monthly <- array()
models.monthly <- list()
model.info.monthly <- list()
MASE.frc.monthly <- array()
forecasts.monthly <- list()
series95.monthly <- list()
series05.monthly <- list()

# array to store overall model fitting info
overall.info <- array(NA,
                      dim = c(303,9), 
                      dimnames = list(NULL,c("Series Category","Season Frequency","Model Fitted","Series Used","MASE",
                                             "Shapiro-Wilks","AIC","AICc","BIC")))
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
  MASE.monthly[i] <- model.info.monthly[[i]][,"MASE"]
  shapiro.monthly[i] <- model.info.monthly[[i]][,"Shapiro-Wilks"]
  
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
  
  MASE.frc.monthly[i] <- as.numeric(MASE.custom(series05.monthly[[i]], forecasts.monthly[[i]]))
  
  # storing overall info 
  overall.info[count,"Series Category"] <- as.character(monthly[i,"Category"])
  overall.info[count,"Season Frequency"] <- freq
  overall.info[count,"Model Fitted"] <- models.monthly[[i]]$method
  overall.info[count,"Series Used"] <- as.character(model.info.monthly[[i]][,"Series Used"])
  overall.info[count,"MASE"] <- MASE.monthly[i]
  overall.info[count,"Shapiro-Wilks"] <- shapiro.monthly[i]
  overall.info[count,"AIC"] <- if(grepl("ETS",models.monthly[[i]]$method)){models.monthly[[i]]$aic}else{models.monthly[[i]]$model$aic}
  overall.info[count,"AICc"] <- if(grepl("ETS",models.monthly[[i]]$method)){models.monthly[[i]]$aicc}else{models.monthly[[i]]$model$aicc}
  overall.info[count,"BIC"] <- if(grepl("ETS",models.monthly[[i]]$method)){models.monthly[[i]]$bic}else{models.monthly[[i]]$model$bic}
  
  count <- count + 1
}

#  objects to store results (quarterly)
MASE.quarterly <- array()
shapiro.quarterly <- array()
models.quarterly <- list()
model.info.quarterly <- list()
MASE.frc.quarterly <- array()
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
  MASE.quarterly[i] <- model.info.quarterly[[i]][,"MASE"]
  shapiro.quarterly[i] <- model.info.quarterly[[i]][,"Shapiro-Wilks"]
  
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
  
  MASE.frc.quarterly[i] <- as.numeric(MASE.custom(series05.quarterly[[i]], forecasts.quarterly[[i]]))
  
  # storing overall info 
  overall.info[count,"Series Category"] <- as.character(quarterly[i,"Category"])
  overall.info[count,"Season Frequency"] <- freq
  overall.info[count,"Model Fitted"] <- models.quarterly[[i]]$method
  overall.info[count,"Series Used"] <- as.character(model.info.quarterly[[i]][,"Series Used"])
  overall.info[count,"MASE"] <- MASE.quarterly[i]
  overall.info[count,"Shapiro-Wilks"] <- shapiro.quarterly[i]
  overall.info[count,"AIC"] <- if(grepl("ETS",models.quarterly[[i]]$method)){models.quarterly[[i]]$aic}else{models.quarterly[[i]]$model$aic}
  overall.info[count,"AICc"] <- if(grepl("ETS",models.quarterly[[i]]$method)){models.quarterly[[i]]$aicc}else{models.quarterly[[i]]$model$aicc}
  overall.info[count,"BIC"] <- if(grepl("ETS",models.quarterly[[i]]$method)){models.quarterly[[i]]$bic}else{models.quarterly[[i]]$model$bic}
  
  count <- count + 1
}

#  objects to store results (yearly)
MASE.yearly <- array()
shapiro.yearly <- array()
models.yearly <- list()
model.info.yearly <- list()
MASE.frc.yearly <- array()
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
  MASE.yearly[i] <- model.info.yearly[[i]][,"MASE"]
  shapiro.yearly[i] <- model.info.yearly[[i]][,"Shapiro-Wilks"]
  
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
  
  MASE.frc.yearly[i] <- as.numeric(MASE.custom(series05.yearly[[i]], forecasts.yearly[[i]]))
  
  # storing overall info 
  overall.info[count,"Series Category"] <- as.character(yearly[i,"Category"])
  overall.info[count,"Season Frequency"] <- freq
  overall.info[count,"Model Fitted"] <- models.yearly[[i]]$method
  overall.info[count,"Series Used"] <- as.character(model.info.yearly[[i]][,"Series Used"])
  overall.info[count,"MASE"] <- MASE.yearly[i]
  overall.info[count,"Shapiro-Wilks"] <- shapiro.yearly[i]
  overall.info[count,"AIC"] <- if(grepl("ETS",models.yearly[[i]]$method)){models.yearly[[i]]$aic}else{models.yearly[[i]]$model$aic}
  overall.info[count,"AICc"] <- if(grepl("ETS",models.yearly[[i]]$method)){models.yearly[[i]]$aicc}else{models.yearly[[i]]$model$aicc}
  overall.info[count,"BIC"] <- if(grepl("ETS",models.yearly[[i]]$method)){models.yearly[[i]]$bic}else{models.yearly[[i]]$model$bic}
  
  count <- count + 1
}

# Results
overall.info <- as.data.frame(overall.info)
overall.info[, c(2,5:9)] <- sapply(overall.info[, c(2,5:9)], as.character)
overall.info[, c(2,5:9)] <- sapply(overall.info[, c(2,5:9)], as.numeric)
colnames(overall.info) <- c("Series Category","Season Frequency","Model Fitted","Series Used","MASE",
                          "Shapiro-Wilks","AIC","AICc","BIC")

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


