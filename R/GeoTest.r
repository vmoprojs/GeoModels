######################################################
### Authors: Moreno Bevilacqua, Víctor Morales Oñate.
### Emails:  moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Institutions: 
### Universidad de Valparaiso
### File name: GeoTests.r
### Description:
### This file contains a set of procedures
### for the computation of composite likelihood-based
### statistics and tests.
### Last change: 28/03/2017.
######################################################

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
         # fixed <- model1$fixed
         # if(!is.null(fixed)){
         #     start <- model2$fixed[!names(model2$fixed)%in%names(fixed)]
         #     namesparam <- names(start)
         #     fixed <- as.list(model1$fixed)
         #     start <- as.list(c(start,model2$param))}
         # else{
         #     fixed <- NULL
         #     start <- as.list(c(model2$param,model2$fixed))
         #     namesparam <- names(model2$fixed)}
          # derive the initial parameters: 

         # initparam <- StartParam(model2$coordx,model2$coordy,model2$coordt,
         #                        model2$coordx_dyn, model2$corrmodel,
         #                        model2$data,model2$distance,"Fitting",fixed,model2$grid,
         #                        model2$likelihood,model2$maxdist,model2$maxtime,
         #                        model2$model,model2$n,NULL,NULL,NULL,model2$radius,
         #                        start,NULL,NULL,
         #                        model2$type,model2$type,TRUE,model2$vartype,
         #                        model2$weigthed,model2$winconst,model2$winstp,
         #                        model2$winconst_t,model2$winstp_t,model2$X)            
         
          #numparam <- initparam$numparam
         # dimmat <- numparam^2
         # dmat <- numparam*(numparam+1)/2
         # score <- double(numparam)
         # sensmat <- double(dmat)# H - sensitivity matrix
         # varimat <- double(dmat)# J - variability matrix
         # eps <- (.Machine$double.eps)^(1/3)
         # spacetime <- model2$ spacetime
         # bivariate <- model2$bivariate
         # param <- c(initparam$param,initparam$fixed)
         # paramcorr <- param[initparam$namescorr]
         # nuisance <- param[initparam$namesnuis]
          
          #sel=substr(names(nuisance),1,4)=="mean"
          #mm=as.numeric(nuisance[sel])   ## mean paramteres
          #other_nuis=as.numeric(nuisance[!sel]) 
          #nuisance=c(mm,other_nuis)
          #num_betas=ncol(model2$X)

          
          #GD=.C('GodambeMat',as.double(mm),as.integer(bivariate),as.double(model2$coordx),as.double(model2$coordy),
           # as.double(model2$coordt),as.integer(initparam$corrmodel),as.double(model2$data),as.integer(initparam$distance),as.double(eps),
           # as.integer(initparam$flagcorr),as.integer(initparam$flagnuis), as.integer(model2$grid),as.integer(initparam$likelihood),
           #  as.double(c((model2$X)%*%mm)),as.integer(initparam$model),as.double(model2$n),as.integer(num_betas),
           #  as.integer(numparam),as.integer(initparam$numparamcorr),as.double(paramcorr),as.double(nuisance),
           #  score=score,sensmat=sensmat,as.integer(spacetime),as.integer(initparam$type),
           #  varimat=varimat,as.integer(initparam$vartype),as.double(model2$winconst),as.double(model2$winstp),
           #  as.double(model2$winconst_t),as.double(model2$winstp_t),
           #  as.integer(model2$weigthed),c(model2$X),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)

         # score=GD$score;
         # varimat=GD$varimat 

        #  H <- matrix(rep(0,dimmat),ncol=numparam)# H - sensitivity matrix
        #  J <- matrix(rep(0,dimmat),ncol=numparam)# J - variability matrix
        #  namesmat <- c(initparam$namesnuis[as.logical(initparam$flagnuis)],
        #                initparam$namescorr[as.logical(initparam$flagcorr)])
        #  names(score) <- namesmat
        #  score <- score[initparam$namesparam]
        #  dimnames(H) <- list(namesmat,namesmat)
        #  dimnames(J) <- list(namesmat,namesmat)
        #  if(numparam>1){
        #      H[lower.tri(H,diag=TRUE)] <- sensmat
        #      H <- t(H)
        #      H[lower.tri(H,diag=TRUE)] <- sensmat
        #      J[lower.tri(J,diag=TRUE)] <- varimat
        #      J <- t(J)
        #      J[lower.tri(J,diag=TRUE)] <- varimat
        #      H <- H[initparam$namesparam,initparam$namesparam]
        #      J <- J[initparam$namesparam,initparam$namesparam]
        #    }
        #  else{
        #      H[1,1] <- sensmat
        #      J[1,1] <- varimat}
        #  .C('DeleteGlobalVar', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
          #return(list(score=score, sensmat=model1$sensmat, varimat=J))
          return(list(score=model1$score, sensmat=model1$sensmat, varimat=model1$varimat))
      }
      # compute the statistic:
      StatiTest <- function(df, model1, model2, statistic)
      {
          namesparam <- names(model1$param)[!names(model1$param)%in%names(model2$param)]
  ### composite likelihood case 
  if(model1$likelihood %in% c("Marginal","Conditional","Difference"))
  {
          # compute the Wald-type statistic:
          if(statistic=='Wald'){
              theta <- model1$param[namesparam]-model2$fixed[namesparam]
              varcov <- model1$varcov[namesparam,namesparam]# Restricted variance-covariance matrix
              W <- t(theta)%*%solve(varcov)%*%theta
              print(model1$param[namesparam])
              print(model2$fixed[namesparam])
              print(theta)
              print(df)
              print(varcov)
              nu <- df}
          if(statistic=="WilksCB"){
              W <- 2*(model1$logCompLik-model2$logCompLik)
              theta <- model1$param[namesparam]-model2$fixed[namesparam]
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

              if(statistic=="WilksPSS"){
                  W <- 2*(model1$logCompLik-model2$logCompLik)
                  varcov <- varcov[namesparam,namesparam]
                  G <- try(solve(varcov),silent=TRUE)
                  score <- null$score[namesparam]# select the score related to the null
                  Hi <- Hi[namesparam,namesparam]# select the inversed of sens related to the null
                  Rao <- t(score)%*%Hi%*%G%*%Hi%*%score# compute the Rao statistic
                  W <- W*Rao/(t(score)%*%Hi%*%score)
                  nu <- df}
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
              theta <- model1$param[namesparam]-model2$fixed[namesparam]
              varcov <- model1$varcov[namesparam,namesparam]# Restricted variance-covariance matrix
              W <- t(theta)%*%solve(varcov)%*%theta
              nu <- df}
       if(statistic=="Wilks"){
              W <- 2*(model1$logCompLik-model2$logCompLik)
              nu <- df}

  }
          return(list(W=W,nu=nu))
      }
      ### START THE MAIN BODY OF THE PROCEDURE
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
          numparam <- c(numparam, length(model$param))
          lmodels[[i]] <- model
          if(!is.matrix(lmodels[[i]]$varcov))
              stop("one of the fitted models does not have a valid variance-covariance matrix\n")
          if(i>1){
            j <- i-1
            if((!all(names(lmodels[[j]]) %in% names(lmodels[[i]]))) &&
               (!identical(lmodels[[j]]$model,lmodels[[i]]$model))) 
              stop("models are not nested\n")
            # Define the degrees of freedom:
            df[j] <- length(lmodels[[j]]$param)-length(lmodels[[i]]$param)
            if(df[j] <= 0) stop("model are not nested\n")
            stat <- StatiTest(df[j],lmodels[[j]],lmodels[[i]],statistic)
            nu[j] <- stat$nu
            W[j] <- stat$W
            pvalue[j] <- pchisq(W[j], df = nu[j], lower.tail = FALSE)
          }
      }
      ### END THE MAIN BODY OF THE PROCEDURE
      #print a table with the hypothesis testing:
      table <- data.frame(numparam, c(NA, df), c(NA, nu),c(NA, W), c(NA, pvalue))
      dimnames(table) <- list(models, c("Num.Par", "Diff.Par", "Df","Chisq", "Pr(>chisq)"))
      structure(table, heading = c("Statistical Hypothesis Test Table\n"),
                class = c("data.frame"))
  }
