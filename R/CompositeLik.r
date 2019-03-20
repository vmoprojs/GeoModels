####################################################
### Authors:  Moreno Bevilacqua, Víctor Morales Oñate.
### Email: moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: CompLik.r
### Description:
### This file contains a set of procedures
### for maximum composite-likelihood fitting of
### random fields.
### Last change: 28/03/2017.
####################################################




### Optim call for Composite log-likelihood maximization

CompLik <- function(bivariate, coordx, coordy ,coordt,coordx_dyn,corrmodel, data, distance, flagcorr, flagnuis, fixed, GPU,grid,
                           likelihood, local,lower, model, n, namescorr, namesnuis, namesparam,
                           numparam, numparamcorr, optimizer, onlyvar, param, spacetime, type,
                           upper, varest, vartype, weigthed, winconst, winstp,winconst_t, winstp_t, ns,X,sensitivity)
  {
    ### Define the object function:
    comploglik <- function(param,coordx, coordy ,coordt, corrmodel, data, fixed, fan, n, namescorr, namesnuis,namesparam,weigthed,X,ns,NS,GPU,local)
      {
        names(param) <- namesparam
        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
        mm=as.numeric(nuisance[sel])   ## mean paramteres
        other_nuis=as.numeric(nuisance[!sel])   ## or nuis parameters (nugget sill skew df)
        result <- .C(as.character(fan),as.integer(corrmodel),as.double(coordx),as.double(coordy),as.double(coordt), 
                    as.double(data), 
                   as.integer(n),as.double(paramcorr), as.integer(weigthed), 
                   res=double(1),as.double(c(X%*%mm)),as.double(0),as.double(other_nuis),
                    as.integer(ns),as.integer(NS),as.integer(local),as.integer(GPU),
                    PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)$res   

        #result <- dotCall64::.C64(as.character(fan), 
         #            SIGNATURE = c("integer",rep("double",4),"integer","double","integer",
          #                          rep("double",4),rep("integer",4)), 
           #          INTENT = c("r","r","r","r","r","r","r","r","w","r","r","r","r","r","r","r"), 
            #         NAOK = TRUE,PACKAGE='GeoModels',
             #        corrmodel,coordx,coordy,coordt,data, n,paramcorr, weigthed, 
              #                    res = vector_dc("numeric",1),c(X%*%mm),0,other_nuis,ns,NS,local,GPU)$res
              print(result)
         return(-result)
      }
     comploglik_biv <- function(param,coordx, coordy ,coordt, corrmodel, data, fixed, fan, n, namescorr, namesnuis,namesparam,weigthed,X,ns,NS,GPU,local)
      {

        names(param) <- namesparam
        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
        mm=as.numeric(nuisance[sel])
        #mm1=c(rep(mm[1],dimat/2),rep(mm[2],dimat/2))
        mm1=c(rep(mm[1],ns[1]),rep(mm[2],ns[2]))
        other_nuis=as.numeric(nuisance[!sel]) 
        result <- .C(fan,as.integer(corrmodel),as.double(coordx),as.double(coordy),as.double(coordt), as.double(data),as.integer(n), 
                     as.double(paramcorr), as.integer(weigthed), res=double(1),as.double(c(X*mm1)),
                    as.double(0),as.double(other_nuis),as.integer(ns),as.integer(NS),as.integer(local),as.integer(GPU),
                     PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)$res

        return(-result)
      }
   ############################################################################################################################################
   ############################################################################################################################################
   ############################################################################################################################################
    numcoord=length(coordx);numtime=1;spacetime_dyn=FALSE
    if(spacetime) numtime=length(coordt)
    if(bivariate) numtime=2
    if(!is.null(coordx_dyn)) spacetime_dyn=TRUE
    dimat <- numcoord*numtime#
    NS=cumsum(ns)
    if(is.null(dim(X)))
    {
       X=as.matrix(rep(1,dimat))
       if((spacetime||bivariate)&& spacetime_dyn)  X=as.matrix(rep(1,NS[numtime]))
    }
    num_betas=ncol(X)   
    fname <- NULL; hessian <- FALSE

    if(all(model==1,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss'
    if(all(model==1,likelihood==3,type==1)) fname <- 'Comp_Diff_Gauss'
    if(all(model==1,likelihood==3,type==2)) {fname <- 'Comp_Pair_Gauss'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==2,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomGauss'
                                             if(varest & vartype==2) hessian <- TRUE}
    if(all(model==11,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomGauss'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==19,likelihood==3,type==2)){ namesnuis=c(namesnuis,"z")
                                              fixed<- c(fixed, list(z=min(n)))
                                              fname <- 'Comp_Pair_Binom2Gauss'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==14,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomnegGauss'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==16,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomnegGauss'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==15,likelihood==3,type==2)){ fname <- 'Comp_Pair_PoisbinGauss'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==17,likelihood==3,type==2)){ fname <- 'Comp_Pair_PoisbinnegGauss'
                                              if(varest & vartype==2) hessian <- TRUE}    
    if(all(model==13,likelihood==3,type==2)){ fname <- 'Comp_Pair_WrapGauss'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==10,likelihood==3,type==2)){ fname <- 'Comp_Pair_SkewGauss'
                                              if(varest & vartype==2) hessian <- TRUE}
     if(all(model==21,likelihood==3,type==2)){ fname <- 'Comp_Pair_Gamma'
                                              if(varest & vartype==2) hessian <- TRUE}
      if(all(model==21,likelihood==1,type==2)){ fname <- 'Comp_Cond_Gamma'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==33,likelihood==3,type==2)){ fname <- 'Comp_Pair_Kumaraswamy'
                                              if(varest & vartype==2) hessian <- TRUE}
   if(all(model==26,likelihood==3,type==2)){ fname <- 'Comp_Pair_Weibull'
                                              if(varest & vartype==2) hessian <- TRUE}    
   if(all(model==26,likelihood==1,type==2)){ fname <- 'Comp_Cond_Weibull'
                                              if(varest & vartype==2) hessian <- TRUE}                                       
    if(all(model==24,likelihood==3,type==2)){ fname <- 'Comp_Pair_LogLogistic'
                                              if(varest & vartype==2) hessian <- TRUE}     
    if(all(model==25,likelihood==3,type==2)){ fname <- 'Comp_Pair_Logistic'
                                              if(varest & vartype==2) hessian <- TRUE}                                                                              
    if(all(model==23,likelihood==3,type==2)){ fname <- 'Comp_Pair_2Gamma'
                                              if(varest & vartype==2) hessian <- TRUE}                                        
    if(all(model==22,likelihood==3,type==2)){ fname <- 'Comp_Pair_LogGauss';
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==18,likelihood==3,type==2)){ fname <- 'Comp_Pair_SkewTGauss'
                                              if(varest & vartype==2) hessian <- TRUE}  
    if(all(model==27,likelihood==3,type==2)){ fname <- 'Comp_Pair_TWOPIECET'
                                              if(varest & vartype==2) hessian <- TRUE} 
    if(all(model==29,likelihood==3,type==2)){ fname <- 'Comp_Pair_TWOPIECEGauss'
                                              if(varest & vartype==2) hessian <- TRUE} 
    if(all(model==31,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomTWOPIECEGauss'
                                              if(varest & vartype==2) hessian <- TRUE} 
    if(all(model==32,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomnegTWOPIECEGauss'
                                              if(varest & vartype==2) hessian <- TRUE}
     if(all(model==12,likelihood==3,type==2)){ fname <- 'Comp_Pair_T'
                                              if(varest & vartype==2) hessian <- TRUE}   
     if(all(model==20,likelihood==3,type==2)){ fname <- 'Comp_Pair_SinhGauss'
                                              if(varest & vartype==2) hessian <- TRUE}
                                              #print(fname)
    if(sensitivity)hessian=TRUE
    if(spacetime) fname <- paste(fname,"_st",sep="")
    if(bivariate) fname <- paste(fname,"_biv",sep="")
    fname <- paste(fname,"2",sep="")
    path.parent <- getwd()
    # if(!is.null(GPU)) 
    # {
    #   path <- system.file("CL", "Kernel.cl", package = "GeoModels")
    #   path <- gsub("/Kernel.cl","/",path);setwd(path)
    #   fname <- paste(fname,"_OCL",sep="")
    #   .C("create_binary_kernel",  as.integer(GPU),as.character(fname),  PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)
    # }
    if(!is.null(GPU))
    {
      fname <- paste(fname,"_OCL",sep="")
      # cat("fname de Composit.r: ",fname,"\n")
      
      path <- system.file("CL", paste(fname,".cl",sep = ""), package = "GeoModels")
      path <- gsub(paste("/",paste(fname,".cl",sep = ""),sep = ""),"/",path)
      # .C("create_binary_kernel",  as.integer(GPU),as.character(fname),  PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)
      setwd(path)
      .C("create_binary_kernel",  as.integer(GPU),as.character(fname),  PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)
    }
   
    if(grid)    {a=expand.grid(coordx,coordy);coordx=a[,1];coordy=a[,2]; }
    else {      if((spacetime||bivariate)&&(!spacetime_dyn)){
                                  data=c(t(data))
                                  coordx=rep(coordx,numtime);coordy=rep(coordy,numtime)
                }
                if((spacetime||bivariate)&&(spacetime_dyn)) data=unlist(data)          
    }
   if(spacetime||bivariate)   NS=c(0,NS)[-(length(ns)+1)]



   if(!onlyvar){
  ##############################.  spatial or space time ############################################
   if(!bivariate)           {
    if(length(param)==1) {
         optimizer="optimize"  
     CompLikelihood <- optimize(f=comploglik, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed, fan=fname,  lower=lower, n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, maximum = FALSE,
                              upper=upper,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)}
   if(length(param)>1) {
    if(optimizer=='L-BFGS-B')
      CompLikelihood <- optim(par=param,fn=comploglik, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, control=list(fnscale=1,
                              factr=1, pgtol=1e-14, maxit=100000), data=data, fixed=fixed,
                              fan=fname, hessian=hessian, lower=lower, method=optimizer,n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                              upper=upper,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
    if(optimizer=='BFGS')
        CompLikelihood <- optim(par=param, fn=comploglik,  coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, control=list(fnscale=1,
                             reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=hessian, method=optimizer,n=n,namescorr=namescorr,
                                  namesnuis=namesnuis,namesparam=namesparam,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
      if(optimizer=='Nelder-Mead')
        CompLikelihood <- optim(par=param, fn=comploglik,  coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, control=list(fnscale=1,
                             reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=hessian, method=optimizer,n=n,namescorr=namescorr,
                                  namesnuis=namesnuis,namesparam=namesparam,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
    if(optimizer=='nlm')
    CompLikelihood <- nlm(f=comploglik,p=param,steptol = 1e-4, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                               fan=fname,hessian=hessian,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                               weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
                               }}
     ############################## bivariate  ############################################                           
    if(bivariate)           {
     if(length(param)==1)
        {
         optimizer="optimize" 
       CompLikelihood <- optimize(f=comploglik_biv, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed,fan=fname,  lower=lower,n=n,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,maximum = FALSE,
                              upper=upper,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)}
      if(length(param)>1) {   
    if(optimizer=='L-BFGS-B')
      CompLikelihood <- optim(param,comploglik_biv, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, control=list(fnscale=1,
                              factr=1, pgtol=1e-14, maxit=100000), data=data, fixed=fixed,
                              fan=fname, hessian=hessian, lower=lower, method=optimizer,n=n,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,
                              upper=upper,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
       if(optimizer=='BFGS')
      CompLikelihood <- optim(param,comploglik_biv, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, control=list(fnscale=1,
                              reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=hessian, method=optimizer,n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam ,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
   if(optimizer=='Nelder-Mead')
      CompLikelihood <- optim(param,comploglik_biv, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, control=list(fnscale=1,
                              reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=hessian, method=optimizer,n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam ,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
   if(optimizer=='nlm') 
        CompLikelihood <- nlm( f=comploglik_biv,p=param, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                               fan=fname,hessian=hessian,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                               weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
                               }
    }                    
      ########################################################################################   
      ########################################################################################
    # check the optimisation outcome
      if(optimizer=='Nelder-Mead'){
        CompLikelihood$value = -CompLikelihood$value
        names(CompLikelihood$par)<- namesparam
        if(CompLikelihood$convergence == 0)
        CompLikelihood$convergence <- 'Successful'
        else
        if(CompLikelihood$convergence == 1)
        CompLikelihood$convergence <- 'Iteration limit reached'
        else
        CompLikelihood$convergence <- "Optimization may have failed"
        if(CompLikelihood$value==-1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
    }
    if(optimizer=='L-BFGS-B'||  optimizer=='BFGS'){
        CompLikelihood$value = -CompLikelihood$value
        names(CompLikelihood$par)<- namesparam
        if(CompLikelihood$convergence == 0)
        CompLikelihood$convergence <- 'Successful'
        else
        if(CompLikelihood$convergence == 1)
        CompLikelihood$convergence <- 'Iteration limit reached'
        else
        CompLikelihood$convergence <- "Optimization may have failed"
        if(CompLikelihood$value==-1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
    }
    if(optimizer=='nlm'){
        CompLikelihood$par <- CompLikelihood$estimate
        names(CompLikelihood$par)<- namesparam
        CompLikelihood$value <- -CompLikelihood$minimum
        if(CompLikelihood$code == 1)
        CompLikelihood$convergence <- 'Successful'
        else
        if(CompLikelihood$code == 4)
        CompLikelihood$convergence <- 'Iteration limit reached'
        else
        CompLikelihood$convergence <- "Optimization may have failed"
        if(CompLikelihood$value==-1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
    }
    if(optimizer=='optimize'){
    param<-CompLikelihood$minimum
    CompLikelihood$par<-param  
    names(CompLikelihood$par)<- namesparam
    maxfun <- -CompLikelihood$objective
    CompLikelihood$value <- maxfun
    CompLikelihood$convergence <- 'Successful'
    }
  }
    else {
          CompLikelihood=as.list(0)
          names(CompLikelihood)="value"
          CompLikelihood$par <- param
          CompLikelihood$claic <- NULL;CompLikelihood$clbic <- NULL;
          CompLikelihood$convergence <- 'Successful'
          CompLikelihood$value=0
          if(hessian) 
          {
               if(!bivariate)  
                CompLikelihood$hessian=numDeriv::hessian(func=comploglik,x=param,method="Richardson",  coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed,fan=fname,n=n,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,
                              weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
               if(bivariate)  
               CompLikelihood$hessian=numDeriv::hessian(func=comploglik_biv,x=param,method="Richardson",coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                             data=data, fixed=fixed,fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                             weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
               rownames(CompLikelihood$hessian)=namesparam
               colnames(CompLikelihood$hessian)=namesparam
               print(CompLikelihood$hessian)
          }
  }



#####################################
if(!is.matrix(CompLikelihood$hessian)||!is.numeric(sum(CompLikelihood$hessian)))
  {
if(!bivariate)  

CompLikelihood$hessian=numDeriv::hessian(func=comploglik,x=CompLikelihood$par,method="Richardson",  coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed,fan=fname,n=n,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,
                              weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
if(bivariate)  
CompLikelihood$hessian=numDeriv::hessian(func=comploglik_biv,x=CompLikelihood$par,method="Richardson",coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                             data=data, fixed=fixed,fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                               weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
rownames(CompLikelihood$hessian)=namesparam
colnames(CompLikelihood$hessian)=namesparam
  }

####################################
        if( (CompLikelihood$convergence!='Successful')||CompLikelihood$value==-1e+15) { print("Optimization failed: try with other starting values ")}
          else{
    if(varest)
          {
            # The sensitivity (H) and the variability (J) matrices
            # are estimated by the sample estimators contro-parts:
            dimmat <- numparam^2
            dmat <- numparam*(numparam+1)/2
            eps <- (.Machine$double.eps)^(1/3)
            param <- c(CompLikelihood$par, fixed)
            score <- double(numparam)
            paramcorr <- param[namescorr]
            nuisance <- param[namesnuis]
            sel=substr(names(nuisance),1,4)=="mean"
            mm=as.numeric(nuisance[sel])   ## mean paramteres
            other_nuis=as.numeric(nuisance[!sel]) 
            nuisance=c(mm,other_nuis)
            sensmat <- double(dmat);varimat <- double(dmat)
            # Set the window parameter:

            GD=.C('GodambeMat',as.double(mm),as.integer(bivariate),as.double(coordx),as.double(coordy),
              as.double(coordt),as.integer(corrmodel), as.double(data),as.integer(distance),as.double(eps),
              as.integer(flagcorr), as.integer(flagnuis),as.integer(grid),as.integer(likelihood),
              as.double(c(X%*%mm)),as.integer(model),as.double(n),as.integer(num_betas),
              as.integer(numparam),as.integer(numparamcorr),as.integer(length(paramcorr)),as.double(paramcorr),as.double(nuisance),
              score=score,sensmat=sensmat,as.integer(spacetime),as.integer(type),
              varimat=varimat,as.integer(vartype),as.double(winconst),as.double(winstp),as.double(winconst_t),as.double(winstp_t),
              as.integer(weigthed),c(t(X)),as.integer(ns),as.integer(NS),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)

            #if(!sum(GD$varimat)) print("Std error estimation failed")
            # Set score vectore:

            CompLikelihood$winconst<-winconst
            CompLikelihood$winstp<-winstp
            CompLikelihood$score <- GD$score
            # Set sensitivity matrix:
            CompLikelihood$sensmat <- matrix(rep(0,dimmat),ncol=numparam)
            # Set variability matrix:
            CompLikelihood$varimat <- matrix(rep(0,dimmat),ncol=numparam)
            namesgod <- c(namesnuis[as.logical(flagnuis)], namescorr[as.logical(flagcorr)])
            names(CompLikelihood$score) <- namesgod
            CompLikelihood$score <- CompLikelihood$score[namesparam]
            dimnames(CompLikelihood$sensmat) <- list(namesgod, namesgod)
            dimnames(CompLikelihood$varimat) <- list(namesgod, namesgod)
            if(numparam>1){
             # CompLikelihood$sensmat[lower.tri(CompLikelihood$sensmat, diag=TRUE)] <- GD$sensmat
             # CompLikelihood$sensmat <- t(CompLikelihood$sensmat)
             # CompLikelihood$sensmat[lower.tri(CompLikelihood$sensmat, diag=TRUE)] <- GD$sensmat
              CompLikelihood$varimat[lower.tri(CompLikelihood$varimat, diag=TRUE)] <- GD$varimat
              CompLikelihood$varimat <- t(CompLikelihood$varimat)
              CompLikelihood$varimat[lower.tri(CompLikelihood$varimat, diag=TRUE)] <- GD$varimat
              CompLikelihood$sensmat <- CompLikelihood$sensmat[namesparam, namesparam]
              CompLikelihood$varimat <- CompLikelihood$varimat[namesparam, namesparam]}
            else {CompLikelihood$sensmat[1,1] <- sensmat
                  CompLikelihood$varimat[1,1] <- varimat}
            if(hessian) CompLikelihood$sensmat=CompLikelihood$hessian

            print(CompLikelihood$sensmat)
            print(CompLikelihood$varimat)
      
            icholsensmat <- try(chol(CompLikelihood$sensmat), silent = TRUE)
            isensmat <- try(chol2inv(icholsensmat), silent = TRUE)

             if(!is.matrix(isensmat) || !is.matrix(CompLikelihood$varimat))
              {
                warning("observed information matrix is singular")
                CompLikelihood$varcov <- 'none'
                CompLikelihood$stderr <- 'none'
              }
            else
              { 
                penalty <- crossprod(CompLikelihood$varimat,isensmat)
                CompLikelihood$claic <- -2 * CompLikelihood$value + 2*sum(diag(penalty))
                CompLikelihood$clbic <- -2 * CompLikelihood$value + log(dimat)*sum(diag(penalty))
                CompLikelihood$varcov <- crossprod(isensmat,penalty)
                dimnames(CompLikelihood$varcov) <- list(namesparam, namesparam)
                CompLikelihood$stderr <- diag(CompLikelihood$varcov)
                if(any(CompLikelihood$stderr < 0))
                  CompLikelihood$stderr <- 'none'
                else
                  CompLikelihood$stderr <- sqrt(CompLikelihood$stderr)
              }
        }
    setwd(path.parent)
      }
      if(hessian) CompLikelihood$sensmat=CompLikelihood$hessian
    if(!is.null(GPU)) gc()
    return(CompLikelihood)
  }

