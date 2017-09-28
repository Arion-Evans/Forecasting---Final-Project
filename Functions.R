library(forecast)


expSmooth <- function(ts) {
  # array to store models
  models <- list()
  model.info <- array(NA,
                      dim = c(1500,6), 
                      dimnames = list(NULL,c("ID","MASE","Shapiro-Wilks","Series Used","Lambda","Added Value")))
  
  # model options
  exponential <- c(TRUE,FALSE)
  damped <- c(TRUE,FALSE)
  method <- c("ANN","ANA","MNN","MNA","MNM")
  methodDamp <- c("AAN","AAA","MAN","MMN","MAA","MAM","MMM")  
  opt.crit <- c("lik","amse","mse","sigma","mae")
  bounds <- c("both","usual","admissible")
  
  # time series lists
  ts.pos.list <- list()
  ts.neg.list <- list()
  ts.neg.list[[1]] <- ts
  
  
  # adjusting negative series & box-cox transformations
  min_value <- min(ts)
  add_value <- 0
  if(min_value <= 0){
    add_value <- abs(min_value) + 2
    ts_add <- ts + add_value
    lambda <- BoxCox.lambda(ts_add, method = "loglik")
    bc <- BoxCox(ts_add,lambda)
    ts.pos.list[[1]] <- ts_add
    ts.pos.list[[2]] <- bc
    ts.neg.list[[2]] <- bc
  }else{
    lambda <- BoxCox.lambda(ts, method = "loglik")
    bc <- BoxCox(ts,lambda)
    ts.pos.list[[1]] <- ts
    ts.pos.list[[2]] <- bc
    ts.neg.list[[2]] <- bc
  }
  
  # differencing
  ts.diff <- diff(ts)
  ts.neg.list[[3]] <- ts.diff
  
  ts.seas.diff <- diff(ts, lag = frequency(ts))
  ts.neg.list[[4]] <- ts.seas.diff
  
  add_value_seas_diff <- abs(min(ts.seas.diff)) + 2
  ts.seas.diff.add <- ts.seas.diff + add_value_seas_diff
  ts.pos.list[[4]] <- ts.seas.diff.add
  
  add_value_diff <- abs(min(ts.diff)) + 2
  ts.diff.add <- ts.diff + add_value_diff
  ts.pos.list[[3]] <- ts.diff.add
    
  # ts info
  added.value.neg <- c(0,add_value,0,0)
  added.value.pos <- c(add_value,add_value,add_value_diff,add_value_seas_diff)
  lambda.array <- c(0,lambda,0,0)
  ts.type <- c("Original","Box-Cox","1st Difference","Seasonal Difference")
  
  # array row counter
  count <- 1
  
    for(j in 1:length(ts.neg.list)){
      model <- try(holt(ts.neg.list[[j]], damped = TRUE, exponential = FALSE), silent = TRUE)
      models[[count]] <- model
      model.info[count,"ID"] <- count
      model.info[count,"MASE"] <- try(accuracy(model$model)[1,"MASE"], silent = TRUE)
      model.info[count,"Shapiro-Wilks"] <- try(shapiro.test(model$model$residuals)$p.value, silent = TRUE)
      model.info[count,"Series Used"] <- ts.type[j]
      model.info[count,"Lambda"] <- lambda.array[j]
      model.info[count,"Added Value"] <- added.value.neg[j]
      
      count <- count + 1
      
      model <- try(ses(ts.neg.list[[j]]), silent = TRUE)
      models[[count]] <- model
      model.info[count,"ID"] <- count
      model.info[count,"MASE"] <- try(accuracy(model$model)[1,"MASE"], silent = TRUE)
      model.info[count,"Shapiro-Wilks"] <- try(shapiro.test(model$model$residuals)$p.value, silent = TRUE)
      model.info[count,"Series Used"] <- ts.type[j]
      model.info[count,"Lambda"] <- lambda.array[j]
      model.info[count,"Added Value"] <- added.value.neg[j]
      
      count <- count + 1
    }
  
  
    for(i in 1:2){
      for(j in 1:length(ts.pos.list)){
      if(exponential[i] == TRUE){
        model <- try(holt(ts.pos.list[[j]], damped = FALSE, exponential = exponential[i]), silent = TRUE)
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- try(accuracy(model$model)[1,"MASE"], silent = TRUE)
        model.info[count,"Shapiro-Wilks"] <- try(shapiro.test(model$model$residuals)$p.value, silent = TRUE)
        model.info[count,"Series Used"] <- ts.type[j]
        model.info[count,"Lambda"] <- lambda.array[j]
        model.info[count,"Added Value"] <- added.value.pos[j]
        
        count <- count + 1
      }else{
        model <- try(holt(ts.neg.list[[j]], damped = FALSE, exponential = exponential[i]), silent = TRUE)
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- try(accuracy(model$model)[1,"MASE"], silent = TRUE)
        model.info[count,"Shapiro-Wilks"] <- try(shapiro.test(model$model$residuals)$p.value, silent = TRUE)
        model.info[count,"Series Used"] <- ts.type[j]
        model.info[count,"Lambda"] <- lambda.array[j]
        model.info[count,"Added Value"] <- added.value.neg[j]
        
        count <- count + 1
      }
      
      
      model <- try(hw(ts.neg.list[[j]], damped = damped[i], exponential = FALSE, 
                  seasonal = "additive"), silent = TRUE)
      models[[count]] <- model
      model.info[count,"ID"] <- count
      model.info[count,"MASE"] <- try(accuracy(model$model)[1,"MASE"], silent = TRUE)
      model.info[count,"Shapiro-Wilks"] <- try(shapiro.test(model$model$residuals)$p.value, silent = TRUE)
      model.info[count,"Series Used"] <- ts.type[j]
      model.info[count,"Lambda"] <- lambda.array[j]
      model.info[count,"Added Value"] <- added.value.neg[j]
      
      count <- count + 1
      }
    }
  
  
  # option grid
  expand.de <- expand.grid(damped, exponential, stringsAsFactors = FALSE)
  
  for(i in 1:nrow(expand.de)){
    for(j in 1:length(ts.pos.list)){
      
        model <- try(hw(ts.pos.list[[j]], damped = expand.de[i,1], exponential = expand.de[i,2], 
                    seasonal = "multiplicative"), silent = TRUE)
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- try(accuracy(model$model)[1,"MASE"], silent = TRUE)
        model.info[count,"Shapiro-Wilks"] <- try(shapiro.test(model$model$residuals)$p.value, silent = TRUE)
        model.info[count,"Series Used"] <- ts.type[j]
        model.info[count,"Lambda"] <- lambda.array[j]
        model.info[count,"Added Value"] <- added.value.pos[j]
        
        count <- count + 1
    }
  }
  
  # options grid
  expand.mob <- expand.grid(method, opt.crit, bounds, stringsAsFactors = FALSE)
  
  for(i in 1:nrow(expand.mob)){
    for(j in 1:length(ts.pos.list)){
      if(grepl("M",expand.mob[i,1])){
        model <- try(ets(ts.pos.list[[j]],model = expand.mob[i,1], opt.crit = expand.mob[i,2],
                     bounds = expand.mob[i,3]), silent = TRUE)
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- try(accuracy(model)[1,"MASE"], silent = TRUE)
        model.info[count,"Shapiro-Wilks"] <- try(shapiro.test(model$residuals)$p.value, silent = TRUE)
        model.info[count,"Series Used"] <- ts.type[j]
        model.info[count,"Lambda"] <- lambda.array[j]
        model.info[count,"Added Value"] <- added.value.pos[j]
        
        count <- count + 1
      }else{
        model <- try(ets(ts.neg.list[[j]],model = expand.mob[i,1], opt.crit = expand.mob[i,2],
                     bounds = expand.mob[i,3]), silent = TRUE)
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- try(accuracy(model)[1,"MASE"], silent = TRUE)
        model.info[count,"Shapiro-Wilks"] <- try(shapiro.test(model$residuals)$p.value, silent = TRUE)
        model.info[count,"Series Used"] <- ts.type[j]
        model.info[count,"Lambda"] <- lambda.array[j]
        model.info[count,"Added Value"] <- added.value.neg[j]
        
        count <- count + 1
      }
    }
  }
  
  # options grid
  expand.mdob <- expand.grid(methodDamp, damped, opt.crit, bounds, stringsAsFactors = FALSE)
  
  
  for(i in 1:nrow(expand.mdob)){
    for(j in 1:length(ts.pos.list)){
      if(grepl("M",expand.mdob[i,1])){
        model <- try(ets(ts.pos.list[[j]],model = expand.mdob[i,1], damped = expand.mdob[i,2], 
                     bounds = expand.mdob[i,4]), silent = TRUE)
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- try(accuracy(model)[1,"MASE"], silent = TRUE)
        model.info[count,"Shapiro-Wilks"] <- try(shapiro.test(model$residuals)$p.value, silent = TRUE)
        model.info[count,"Series Used"] <- ts.type[j]
        model.info[count,"Lambda"] <- lambda.array[j]
        model.info[count,"Added Value"] <- added.value.pos[j]
        
        count <- count + 1
      }else{
        model <- try(ets(ts.neg.list[[j]],model = expand.mdob[i,1], damped = expand.mdob[i,2], 
                     bounds = expand.mdob[i,4]), silent = TRUE)
        models[[count]] <- model
        model.info[count,"ID"] <- count
        model.info[count,"MASE"] <- try(accuracy(model)[1,"MASE"], silent = TRUE)
        model.info[count,"Shapiro-Wilks"] <- try(shapiro.test(model$residuals)$p.value, silent = TRUE)
        model.info[count,"Series Used"] <- ts.type[j]
        model.info[count,"Lambda"] <- lambda.array[j]
        model.info[count,"Added Value"] <- added.value.neg[j]
        
        count <- count + 1
      }
      
    }
  }
  
  # formatting array
  model.info <- as.data.frame(model.info)
  model.info <- model.info[complete.cases(model.info),]
  model.info[, c(1:3,5:6)] <- sapply(model.info[, c(1:3,5:6)], as.character)
  model.info[, c(1:3,5:6)] <- sapply(model.info[, c(1:3,5:6)], as.numeric)
  colnames(model.info) <- c("ID","MASE","Shapiro-Wilks","Series Used","Lambda","Added Value")
  
  # best MASE model
  best.MASE.info <- model.info[order(model.info[,"MASE"]),][1,]
  
  # models with more than 0.5 shapiro p-value
  best.residual.models <- model.info[model.info[,"Shapiro-Wilks"] > 0.4,]
  best.residual.models <- best.residual.models[complete.cases(best.residual.models),]
  
  # select best shapiro/MASE combo
  if(nrow(best.residual.models) > 0){
    best.shapiro.info <- best.residual.models[order(best.residual.models[,"MASE"]),][1,]
  }else{
    best.shapiro.info <- model.info[order(-model.info[,"Shapiro-Wilks"]),][1,]
  }
  
  MASE.diff <- best.shapiro.info[,"MASE"] - best.MASE.info[,"MASE"]
  shapiro.diff <- best.shapiro.info[,"Shapiro-Wilks"] - best.MASE.info[,"Shapiro-Wilks"]
  
  if(best.MASE.info[,"Shapiro-Wilks"] > 0.4){
    best.model <- models[[best.MASE.info[,"ID"]]]
    best.model.info <- best.MASE.info
  }else if(best.shapiro.info[,"Shapiro-Wilks"] < 0.05){
    best.model <- models[[best.MASE.info[,"ID"]]]
    best.model.info <- best.MASE.info
  }else if(MASE.diff > 0.05){
    best.model <- models[[best.MASE.info[,"ID"]]]
    best.model.info <- best.MASE.info
  }else if((MASE.diff < 0.02) && (shapiro.diff > 0.1)){
    best.model <- models[[best.shapiro.info[,"ID"]]]
    best.model.info <- best.shapiro.info
  }else{
    best.model <- models[[best.MASE.info[,"ID"]]]
    best.model.info <- best.MASE.info
  }
  
  result <- list(best.model, best.model.info)
  return(result)
  
}


MASE.custom <- function(observed , fitted ){ 
  
  Y.t = observed 
  n = length(fitted) 
  e.t = Y.t - fitted 
  sum = 0  
  
  for (i in 2:n){   
    sum = sum + abs(Y.t[i] - Y.t[i-1] ) 
  } 
  
  q.t = e.t / (sum/(n-1)) 
  MASE = mean(abs(q.t))  
  
  return(list(MASE = MASE))
}

invBoxCox <- function(x, lambda){
  if (lambda == 0){
    exp(x)
  }else{
    (lambda*x + 1)^(1/lambda)
  } 
}
  