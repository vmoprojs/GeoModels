####################################################
### File name: GeoTest.r
####################################################

### Procedures are in alphabetical order.

### Statistical hypothesis testing for nested models
GeoTests <- function(object1, object2, ..., statistic)
  {
      # check the type of test:
      Istest <- function(statistic)
      {
        Istest=NULL
          Istest <- switch(statistic,
                           Wald=1,
                           Wilks=2,
                           WilksCB=3,
                           WilksPSS=4,
                           WilksRJ=5,
                           WilksS=6
                           #Rao=7
                           )
          return(Istest)
      }
      # compute gradient, sensitivity and variability under the null:
      UnderNull <- function(model1, model2)
      {
    
          return(list(score=model1$score, sensmat=model1$sensmat, varimat=model1$varimat))
      }
      # compute the statistic:
StatiTest <- function(df, model1, model2, statistic)
      {

          model1$param=unlist(model1$param)
          model2$param=unlist(model2$param)
          model1$fixed=unlist(model1$fixed)
          model2$fixed=unlist(model2$fixed)
          namesparam <- names(model1$param)[!names(model1$param)%in%names(model2$param)]
         
  ### composite likelihood case 
  if(model1$likelihood %in% c("Marginal","Conditional","Difference"))
  {
          # compute the Wald-type statistic:
          if(statistic=='Wald'){
              mm2=model2$fixed[namesparam]
              if(any(is.na(mm2))) mm2=rep(0,df)
              theta <- model1$param[namesparam]-mm2
              varcov <- model1$varcov[namesparam,namesparam]# Restricted variance-covariance matrix
              W <- t(theta)%*%solve(varcov)%*%theta
              nu <- df}
          if(statistic=="WilksCB"){
              W <- 2*(model1$logCompLik-model2$logCompLik)
              mm2=model2$fixed[namesparam]
              if(any(is.na(mm2))) mm2=rep(0,df)
              theta <- model1$param[namesparam]-mm2
              G <- solve(model1$varcov)[namesparam,namesparam]
              H <- model1$sensmat[namesparam,namesparam]
              W <- W*((t(theta)%*%G%*%theta)/(t(theta)%*%H%*%theta))
              nu <- df}
          if(any(statistic==c("WilksPSS","WilksRJ","WilksS"))){
              null <- UnderNull(model1,model2)

              H <- null$sensmat# H - sensitivity matrix
              Hi <- solve(H)# inverse of H
              J <- null$varimat# J - variability matrix
              Ji <- solve(J)# inverse of J
              G <- H%*%Ji%*%H# compute the Godambe matrix
              varcov <- solve(G)# variance-covariance matrix

              #if(statistic=="WilksPSS"){
              #    W <- 2*(model1$logCompLik-model2$logCompLik)
              #    varcov <- varcov[namesparam,namesparam]
              #    G <- try(solve(varcov),silent=TRUE)
              #    score <- null$score[namesparam]# select the score related to the null
              #    Hi <- Hi[namesparam,namesparam]# select the inversed of sens related to the null
              #    Rao <- t(score)%*%Hi%*%G%*%Hi%*%score# compute the Rao statistic
              #    W <- W*Rao/(t(score)%*%Hi%*%score)
              #    nu <- df}
              if(statistic=="WilksRJ"){
                  W <- 2*(model1$logCompLik-model2$logCompLik)
                  lambda <- eigen(solve(Hi[namesparam, namesparam])%*%varcov[namesparam, namesparam])$values
                  W <- W/mean(lambda)
                  nu <- df}
              if(statistic=="WilksS"){
                  W <- 2*(model1$logCompLik-model2$logCompLik)
                  lambda <- eigen(solve(Hi[namesparam, namesparam])%*%varcov[namesparam, namesparam])$values
                  slambda <- sum(lambda)
                  nu <- slambda^2/sum(lambda^2)
                  W <- nu*W/slambda}
          #if(statistic=="Rao"){
          #        varcov <- varcov[namesparam,namesparam]
          #        G <- try(solve(varcov),silent=TRUE)
          #        score <- null$score[namesparam]# select the score related to the null
          #        Hi <- Hi[namesparam,namesparam]# select the inversed of sens related to the null
          #        W <- t(score)%*%Hi%*%G%*%Hi%*%score# compute the Rao statistic
          #        nu <- df}
      
                }
        }
  ### full likelihood case
        if(model1$likelihood %in% c("Full"))
  {
      if(statistic=='Wald'){

                 mm2=model2$fixed[namesparam]
              if(any(is.na(mm2))) mm2=rep(0,df)

              theta <- model1$param[namesparam]-mm2
              
              varcov <- as.matrix(model1$varcov[namesparam,namesparam])
           
    
              W <- t(theta)%*%solve(varcov)%*%theta
              nu <- df}
       if(statistic=="Wilks"){
              W <- 2*(model1$logCompLik-model2$logCompLik)
              nu <- df}

  }
          return(list(W=W,nu=nu))
      }
###########################################################
################# Start ###################################
###########################################################
      # check fitted models:
      if(any(missing(object1), missing(object2)))
         stop("Models one and two must be specified\n")
      # check statistics:
      if(missing(statistic))
         stop('Insert the type of statistic use in the hypothesis test\n')
      if(is.null(Istest(statistic)))
          stop("The name of test does not match with one those available\n")



      # check if there are multipl fitted models:
      objects <- as.list(substitute(list(...)))[-1]
      objects <- sapply(objects,function(x) deparse(x))
      if(!length(objects)) objects <- NULL
      # build a sequence of fitted models:
      models <- c(deparse(substitute(object1)),
                  deparse(substitute(object2)),
                  objects)
      nummod <- length(models)# number of models
      numparam <- NULL# parameters size
      numstat <- length(statistic)# number of statistics
      lmodels <- vector("list", nummod)# define a list of models
      W <- double(nummod - 1)# statistics
      pvalue <- double(nummod - 1)# p-values
      df <- double(nummod - 1)# degrees of freedoms (of the tests)
      nu <- double(nummod - 1)# adjusted df
      # Compute hypothesis tests:
      for(i in 1:nummod)
      {   # consider the ith model:
          model <- get(models[i], envir=parent.frame())
          if(!inherits(model, "GeoFit"))
              stop("use HypoTest only with 'GeoFit' objects\n")
          numparam <- c(numparam, length(unlist(model$param)))
          lmodels[[i]] <- model
          if(!is.matrix(lmodels[[i]]$varcov))
              stop("one of the fitted models does not have a valid variance-covariance matrix\n")
          if(i>1){
            j <- i-1
            if((!all(names(lmodels[[j]]) %in% names(lmodels[[i]]))) &&
               (!identical(lmodels[[j]]$model,lmodels[[i]]$model))) 
              stop("models are not nested\n")
            # Define the degrees of freedom:
            df[j] <- length(unlist(lmodels[[j]]$param))-length(unlist(lmodels[[i]]$param))
            if(df[j] <= 0) stop("model are not nested\n")
            stat <- StatiTest(df[j],lmodels[[j]],lmodels[[i]],statistic)
            nu[j] <- stat$nu
            W[j] <- stat$W
            pvalue[j] <- pchisq(W[j], df = nu[j], lower.tail = FALSE)
          }
      }
###########################################################
###########################################################
###########################################################
      #print a table with the hypothesis testing:
      table <- data.frame(numparam, c(NA, df), c(NA, nu),c(NA, W), c(NA, pvalue))
      dimnames(table) <- list(models, c("Num.Par", "Diff.Par", "Df","Chisq", "Pr(>chisq)"))
      structure(table, heading = c("Statistical Hypothesis Test Table\n"),
                class = c("data.frame"))
  }
