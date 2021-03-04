chkres <- function(model, variable1 = NA, variable2 = NA) {
  require(RVAideMemoire)
  sresid <- resid(model, type = "deviance")
  hist(sresid)
  fitted.glmm <- fitted(model, level=1)        # Extract the fitted (predicted) values
  plot(sresid ~ fitted.glmm)                   # Check for homoscedasticity
  par(mfrow=c(2,2))
  plot(model)                                  # gives a heteroscedasticity plot #fans out, not good
  require(arm)
  par(mfrow=c(1,1))
  binnedplot(fitted(model),resid(model))
  
  # check if a plotting variable is specified; if not, chooses the first explanatory variable
  if(is.na(variable1)){
    if(isS4(model)){
      variable1 <- model@frame[[2]]
    } else {
      variable1 <- model$data[,which(colnames(model$data) == model$formula[[3]])]
    }
  }
  
  plot(sresid ~ variable1)
  
  # if a second plotting variable is specified, plots that too
  suppressWarnings(  if(!is.na(variable2)){ plot(sresid ~ variable2) } )
  
}





chkres.PQL <- function(model,variable1 = FALSE,variable2 = FALSE) {
  require(RVAideMemoire)
  plot(model)
  presid <- resid(model, type = "pearson")
  hist(presid)
  fitted.glmm <- fitted(model, level=1)        # Extract the fitted (predicted) values
  plot(presid ~ fitted.glmm)                   # Check for homoscedasticity
  if (variable1[1] != FALSE){ plot(presid ~ variable1) }
  if (variable2[1] != FALSE){ plot(presid ~ variable2) }
}


chkres.zi <- function(model,variable1 = FALSE,variable2 = FALSE) {
  require(RVAideMemoire)
  presid <- residuals(model)
  hist(presid)
  fitted.glmm <- fitted(model, level=1)        # Extract the fitted (predicted) values
  plot(presid ~ fitted.glmm)                   # Check for homoscedasticity
  require(arm)
  binnedplot(fitted.glmm,presid)
  if (variable1[1] != FALSE){ plot(presid ~ variable1) }
  if (variable2[1] != FALSE){ plot(presid ~ variable2) }
  
}

