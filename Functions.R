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


expSmooth <- function(ts) {
  # array to store models
  models <- list()
  model.info <- array(NA,
                      dim = c(700,6), 
                      dimnames = list(NULL,c("ID","MASE","Shapiro-Wilks","Transformed","Lambda","Added Value")))
  
  # model options
  exponential <- c(TRUE,FALSE)
  damped <- c(TRUE,FALSE)
  method <- c("ANN","ANA","MNN","MNA","MNM")
  methodDamp <- c("AAN","AAA","MAN","MMN","MAA","MAM","MMM")  
  opt.crit <- c("lik","amse","mse","sigma","mae")
  bounds <- c("both","usual","admissible")
  
  # time series list
  ts.list <- list()
  ts.list[[1]] <- ts
  
  # adjusting negative series & box-cox transformations
  min_value <- min(ts)
  add_value <- 0
  if(min_value <= 0){
    add_value <- abs(min_value)+2
    ts_add <- ts + add_value
    lambda <- BoxCox.lambda(ts_add, method = "loglik")
    bc <- BoxCox(ts_add,lambda)
    ts.list[[2]] <- bc
    ts.list[[3]] <- ts_add
  }else{
    lambda <- BoxCox.lambda(ts, method = "loglik")
    bc <- BoxCox(ts,lambda)
    ts.list[[2]] <- bc
  }
  
  
  # array row counter
  count <- 1
  
    for(j in 1:2){
      model <- holt(ts.list[[j]], damped = TRUE, exponential = FALSE)
      models[[count]] <- model
      model.info[count,"ID"] <- count
      model.info[count,"MASE"] <- accuracy(model$model)[1,"MASE"]
      model.info[count,"Shapiro-Wilks"] <- shapiro.test(model$model$residuals)$p.value
      model.info[count,"Added Value"] <- 0
      if(j == 1){
        model.info[count,"Transformed"] <- FALSE
        model.info[count,"Lambda"] <- NA
      }else{
        model.info[count,"Transformed"] <- TRUE
        model.info[count,"Lambda"] <- lambda
      }
      
      count <- count + 1
      
      model <- ses(ts.list[[j]])
      models[[count]] <- model
      model.info[count,"ID"] <- count
      model.info[count,"MASE"] <- accuracy(model$model)[1,"MASE"]
      model.info[count,"Shapiro-Wilks"] <- shapiro.test(model$model$residuals)$p.value
      model.info[count,"Added Value"] <- 0
      if(j == 1){
        model.info[count,"Transformed"] <- FALSE
        model.info[count,"Lambda"] <- NA
      }else{
        model.info[count,"Transformed"] <- TRUE
        model.info[count,"Lambda"] <- lambda
      }
      
      count <- count + 1
    }
  
  
    for(j in 1:2){
      if(exponential[j] == TRUE && min_value <= 0){
        model <- holt(ts.list[[j+1]], damped = FALSE, exponential = exponential[j])
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- accuracy(model$model)[1,"MASE"]
        model.info[count,"Shapiro-Wilks"] <- shapiro.test(model$model$residuals)$p.value
        model.info[count,"Added Value"] <- 0
        if(j == 2){
          model.info[count,"Transformed"] <- FALSE
          model.info[count,"Lambda"] <- NA
        }else{
          model.info[count,"Transformed"] <- TRUE
          model.info[count,"Lambda"] <- lambda
        }
        
        count <- count + 1
      }else{
        model <- holt(ts.list[[j]], damped = FALSE, exponential = exponential[j])
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- accuracy(model$model)[1,"MASE"]
        model.info[count,"Shapiro-Wilks"] <- shapiro.test(model$model$residuals)$p.value
        model.info[count,"Added Value"] <- 0
        if(j == 1){
          model.info[count,"Transformed"] <- FALSE
          model.info[count,"Lambda"] <- NA
        }else{
          model.info[count,"Transformed"] <- TRUE
          model.info[count,"Lambda"] <- lambda
        }
        
        count <- count + 1
      }
      
      
      model <- hw(ts.list[[j]], damped = damped[j], exponential = FALSE, 
                  seasonal = "additive")
      models[[count]] <- model
      model.info[count,"ID"] <- count
      model.info[count,"MASE"] <- accuracy(model$model)[1,"MASE"]
      model.info[count,"Shapiro-Wilks"] <- shapiro.test(model$model$residuals)$p.value
      model.info[count,"Added Value"] <- 0
      if(j == 1){
        model.info[count,"Transformed"] <- FALSE
        model.info[count,"Lambda"] <- NA
      }else{
        model.info[count,"Transformed"] <- TRUE
        model.info[count,"Lambda"] <- lambda
      }
      
      count <- count + 1
    }
  
  
  # option grid
  expand.de <- expand.grid(damped, exponential, stringsAsFactors = FALSE)
  
  for(i in 1:4){
    for(j in 1:2){
      if(min_value <= 0){
        model <- hw(ts.list[[j+1]], damped = expand.de[i,1], exponential = expand.de[i,2], 
                    seasonal = "multiplicative")
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- accuracy(model$model)[1,"MASE"]
        model.info[count,"Shapiro-Wilks"] <- shapiro.test(model$model$residuals)$p.value
        model.info[count,"Added Value"] <- 0
        if(j == 2){
          model.info[count,"Transformed"] <- FALSE
          model.info[count,"Lambda"] <- NA
        }else{
          model.info[count,"Transformed"] <- TRUE
          model.info[count,"Lambda"] <- lambda
        }
        
        count <- count + 1
      }else{
        model <- hw(ts.list[[j]], damped = expand.de[i,1], exponential = expand.de[i,2], 
                    seasonal = "multiplicative")
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- accuracy(model$model)[1,"MASE"]
        model.info[count,"Shapiro-Wilks"] <- shapiro.test(model$model$residuals)$p.value
        model.info[count,"Added Value"] <- 0
        if(j == 1){
          model.info[count,"Transformed"] <- FALSE
          model.info[count,"Lambda"] <- NA
        }else{
          model.info[count,"Transformed"] <- TRUE
          model.info[count,"Lambda"] <- lambda
        }
        
        count <- count + 1
      }
    }
  }
  
  # options grid
  expand.mob <- expand.grid(method, opt.crit, bounds, stringsAsFactors = FALSE)
  
  for(i in 1:nrow(expand.mob)){
    for(j in 1:2){
      if(grepl("M",expand.mob[i,1]) && min_value <= 0){
        model <- ets(ts.list[[j+1]],model = expand.mob[i,1], opt.crit = expand.mob[i,2],
                     bounds = expand.mob[i,3])
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- accuracy(model$model)[1,"MASE"]
        model.info[count,"Shapiro-Wilks"] <- shapiro.test(model$model$residuals)$p.value
        model.info[count,"Added Value"] <- 0
        if(j == 2){
          model.info[count,"Transformed"] <- FALSE
          model.info[count,"Lambda"] <- NA
        }else{
          model.info[count,"Transformed"] <- TRUE
          model.info[count,"Lambda"] <- lambda
        }
        
        count <- count + 1
      }else{
        model <- ets(ts.list[[j]],model = expand.mob[i,1], opt.crit = expand.mob[i,2],
                     bounds = expand.mob[i,3])
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- accuracy(model$model)[1,"MASE"]
        model.info[count,"Shapiro-Wilks"] <- shapiro.test(model$model$residuals)$p.value
        model.info[count,"Added Value"] <- 0
        if(j == 1){
          model.info[count,"Transformed"] <- FALSE
          model.info[count,"Lambda"] <- NA
        }else{
          model.info[count,"Transformed"] <- TRUE
          model.info[count,"Lambda"] <- lambda
        }
        
        count <- count + 1
      }
    }
  }
  
  # options grid
  expand.mdob <- expand.grid(methodDamp, damped, opt.crit, bounds, stringsAsFactors = FALSE)
  
  
  for(i in 1:nrow(expand.mdob)){
    for(j in 1:2){
      if(grepl("M",expand.mdob[i,1]) && min_value <= 0){
        model <- ets(ts.list[[j+1]],model = expand.mdob[i,1], damped = expand.mdob[i,2], 
                     bounds = expand.mdob[i,4])
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- accuracy(model$model)[1,"MASE"]
        model.info[count,"Shapiro-Wilks"] <- shapiro.test(model$model$residuals)$p.value
        model.info[count,"Added Value"] <- 0
        if(j == 2){
          model.info[count,"Transformed"] <- FALSE
          model.info[count,"Lambda"] <- NA
        }else{
          model.info[count,"Transformed"] <- TRUE
          model.info[count,"Lambda"] <- lambda
        }
        
        count <- count + 1
      }else{
        model <- ets(ts.list[[j]],model = expand.mdob[i,1], damped = expand.mdob[i,2], 
                     bounds = expand.mdob[i,4])
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- accuracy(model$model)[1,"MASE"]
        model.info[count,"Shapiro-Wilks"] <- shapiro.test(model$model$residuals)$p.value
        model.info[count,"Added Value"] <- 0
        if(j == 1){
          model.info[count,"Transformed"] <- FALSE
          model.info[count,"Lambda"] <- NA
        }else{
          model.info[count,"Transformed"] <- TRUE
          model.info[count,"Lambda"] <- lambda
        }
        
        count <- count + 1
      }
      
    }
  }
  
  # extracting models with best MASE and residuals
  best.MASE.info <- model.info[order(model.info[,"MASE"]),][1,]
  best.shapiro.info <- model.info[order(-model.info[,"Shapiro-Wilks"]),][1,]
  
  if(best.MASE.info["Shapiro-Wilks"] > 0.05){
    best.model <- models[[best.MASE.info["ID"]]]
    best.model.info <- best.MASE.info
  }else if(best.shapiro.info["Shapiro-Wilks"] < 0.05){
    best.model <- models[[best.MASE.info["ID"]]]
    best.model.info <- best.MASE.info
  }else if((best.shapiro.info["MASE"] - best.MASE.info["MASE"]) > 0.05){
    best.model <- models[[best.MASE.info["ID"]]]
    best.model.info <- best.MASE.info
  }else{
    best.model <- models[[best.shapiro.info["ID"]]]
    best.model.info <- best.shapiro.info
  }
  
  result <- list(best.model, best.model.info)
  return(result)
  
}