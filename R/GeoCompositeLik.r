####################################################
### File name: CompLik.r
####################################################



CompLik <- function(copula,bivariate, coordx, coordy ,coordt,coordx_dyn,corrmodel, data, distance, flagcorr, flagnuis, fixed, GPU,grid,
                           likelihood, local,lower, model, n, namescorr, namesnuis, namesparam,
                           numparam, numparamcorr, optimizer, onlyvar, parallel, param, spacetime, type,
                           upper, varest, vartype, weigthed, winconst, winstp,winconst_t, winstp_t, ns, X,sensitivity,MM)
  {
    ### Define the object function:
    comploglik <- function(param,coordx, coordy ,coordt, corrmodel, data, fixed, fan, n, namescorr, 
                              namesnuis,namesparam,weigthed,X,ns,NS,GPU,local,MM)
      {
        names(param) <- namesparam
        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
        mm=as.numeric(nuisance[sel])
        Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM

       other_nuis=as.numeric(nuisance[!sel])   ## or nuis parameters (nugget sill skew df)
        res=double(1)
       # result <- .C(as.character(fan),as.integer(corrmodel),as.double(coordx),as.double(coordy),as.double(coordt), as.double(data), 
        #           as.integer(n),as.double(paramcorr), as.integer(weigthed), 
         #          res=res,as.double(c(X%*%mm)),as.double(0),as.double(other_nuis),
          #          as.integer(ns),as.integer(NS),as.integer(local),as.integer(GPU),
           #         PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)$res   

        result=dotCall64::.C64(as.character(fan),
        SIGNATURE = c("integer","double","double","double","double",
                         "integer","double","integer","double","double","double","double",
                          "integer","integer","integer","integer"),  
                         corrmodel ,coordx,coordy , coordt ,  data , n , paramcorr ,  weigthed , 
                   res=res, Mean,0,other_nuis ,
                     ns , NS , local ,GPU,
         INTENT =    c("r","r","r","r","r","r","r","r","rw","r", "r","r","r","r","r","r"),
         PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$res
 
         return(-result)
      }
     comploglik_biv <- function(param,coordx, coordy ,coordt, corrmodel, data, fixed, fan, n, namescorr, namesnuis,namesparam,weigthed,X,ns,NS,GPU,local,MM)
      {

        names(param) <- namesparam
        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]
        sel1=substr(names(nuisance),1,6)=="mean_1"
        mm1=as.numeric(nuisance[sel1])
        sel2=substr(names(nuisance),1,6)=="mean_2"
        mm2=as.numeric(nuisance[sel2])
        sel=substr(names(nuisance),1,4)=="mean"
        X1=as.matrix(X[1:ns[1],]);X2=as.matrix(X[(ns[1]+1):(ns[2]+ns[1]),]); 
        other_nuis=as.numeric(nuisance[!sel]) 

        res=double(1)
        #result <- .C(fan,as.integer(corrmodel),as.double(coordx),as.double(coordy),as.double(coordt), as.double(data),as.integer(n), 
        #            as.double(paramcorr), as.integer(weigthed), res=res,as.double(c(X1%*%mm1,X2%*%mm2)),
        #           as.double(0),as.double(other_nuis),as.integer(ns),as.integer(NS),as.integer(local),as.integer(GPU),
        #           PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)$res

        result=dotCall64::.C64(as.character(fan),
        SIGNATURE = c("integer","double","double","double","double",
                         "integer","double","integer","double","double","double","double",
                          "integer","integer","integer","integer"),  
                         corrmodel ,coordx,coordy , coordt ,  data , n , paramcorr ,  weigthed , 
                   res=res, c(X1%*%mm1,X2%*%mm2),0,other_nuis ,
                     ns , NS , local ,GPU,
         INTENT =    c("r","r","r","r","r","r","r","r","rw","r", "r","r","r","r","r","r"),
         PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$res
        return(-result)
      }
   ############################################################################################################################################
   ############################################################################################################################################
   ############################################################################################################################################
 
    numcoord=length(coordx);numtime=1;spacetime_dyn=FALSE
    if(spacetime) numtime=length(coordt)
    if(bivariate) numtime=2
    if(!is.null(coordx_dyn)) spacetime_dyn=TRUE
    if(!spacetime_dyn)  dimat <- numcoord*numtime#
    if(spacetime_dyn)  dimat <- sum(ns)

    NS=cumsum(ns)
    if(is.null(dim(X)))
    {
       X=as.matrix(rep(1,dimat))
       if((spacetime||bivariate)&& spacetime_dyn)  X=as.matrix(rep(1,NS[numtime]))
    }
    else(if(bivariate) X=rbind(X,X))

  
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
      if(all(model==42,likelihood==3,type==2)){ fname <- 'Comp_Pair_Kumaraswamy2'
                                              if(varest & vartype==2) hessian <- TRUE}
   if(all(model==26,likelihood==3,type==2)){ fname <- 'Comp_Pair_Weibull'
                                              if(varest & vartype==2) hessian <- TRUE}    
   if(all(model==26,likelihood==1,type==2)){ fname <- 'Comp_Cond_Weibull'
                                              if(varest & vartype==2) hessian <- TRUE}    
    if(all(model==28,likelihood==3,type==2)){ fname <- 'Comp_Pair_Beta'
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
    if(all(model==39,likelihood==3,type==2)){ fname <- 'Comp_Pair_TWOPIECEBIMODAL'
                                              if(varest & vartype==2) hessian <- TRUE} 
    if(all(model==29,likelihood==3,type==2)){ fname <- 'Comp_Pair_TWOPIECEGauss'
                                              if(varest & vartype==2) hessian <- TRUE} 
    if(all(model==31,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomTWOPIECEGauss'
                                              if(varest & vartype==2) hessian <- TRUE} 
    if(all(model==32,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomnegTWOPIECEGauss'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==12,likelihood==3,type==2)){ fname <- 'Comp_Pair_T'
                                              if(varest & vartype==2) hessian <- TRUE}  
    if(all(model==34,likelihood==3,type==2)){ fname <- 'Comp_Pair_Tukeyh' 
                                              if(varest & vartype==2) hessian <- TRUE} 
    if(all(model==41,likelihood==3,type==2)){ fname <- 'Comp_Pair_Gauss_misp_Tukeygh' 
                                              if(varest & vartype==2) hessian <- TRUE} 
    if(all(model==40,likelihood==3,type==2)){ fname <- 'Comp_Pair_Tukeyhh' 
                                              if(varest & vartype==2) hessian <- TRUE} 
    if(all(model==36,likelihood==3,type==2)){ fname <- 'Comp_Pair_Gauss_misp_Pois'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==35,likelihood==3,type==2)){ fname <- 'Comp_Pair_Gauss_misp_T'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==37,likelihood==3,type==2)){ fname <- 'Comp_Pair_Gauss_misp_SkewT'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==20,likelihood==3,type==2)){ fname <- 'Comp_Pair_SinhGauss'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==38,likelihood==3,type==2)){ fname <- 'Comp_Pair_TWOPIECETukeyh'
                                              if(varest & vartype==2) hessian <- TRUE} 
    if(all(model==30,likelihood==3,type==2)){ fname <- 'Comp_Pair_Pois'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==46,likelihood==3,type==2)){ fname <- 'Comp_Pair_PoisGamma'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==46,likelihood==3,type==2)){ fname <- 'Comp_Pair_PoisGamma'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==43,likelihood==3,type==2)){ fname <- 'Comp_Pair_PoisZIP'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==44,likelihood==3,type==2)){ fname <- 'Comp_Pair_Gauss_misp_PoisZIP'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(all(model==45,likelihood==3,type==2)){ fname <- 'Comp_Pair_BinomnegGaussZINB'
                                              if(varest & vartype==2) hessian <- TRUE}
    if(sensitivity) hessian=TRUE
    if(spacetime) fname <- paste(fname,"_st",sep="")
    if(bivariate) fname <- paste(fname,"_biv",sep="")
    fname <- paste(fname,"2",sep="")

    if(!is.null(copula))
    {
        if(copula=="Gaussian") fname <- paste(fname,"GCop",sep="")
        if(copula=="Clayton")     fname <- paste(fname,"BCop",sep="")
    }
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
      cat("fname de GeoComposite.r: ",fname,"\n")
      
      path <- system.file("CL", paste(fname,".cl",sep = ""), package = "GeoModels")
      path <- gsub(paste("/",paste(fname,".cl",sep = ""),sep = ""),"/",path)
      # .C("create_binary_kernel",  as.integer(GPU),as.character(fname),  PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)
      setwd(path)
      .C("create_binary_kernel",  as.integer(GPU),as.character(fname),  PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)
    }
    
    if((spacetime||bivariate)&&(!spacetime_dyn)){
                                  data=c(t(data))
                                  coordx=rep(coordx,numtime);coordy=rep(coordy,numtime)
                }
                if((spacetime||bivariate)&&(spacetime_dyn)) data=unlist(data)          
    
   if(spacetime||bivariate)   NS=c(0,NS)[-(length(ns)+1)]
   if(is.null(GPU)) GPU=0




   if(!onlyvar){
  ##############################.  spatial or space time ############################################
   if(!bivariate)           {
    if(length(param)==1) {
         optimizer="optimize"  
     CompLikelihood <- optimize(f=comploglik, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed, fan=fname,  lower=lower, n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, maximum = FALSE,
                              upper=upper,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)
          }
   if(length(param)>1) {
    if(optimizer=='L-BFGS-B'&&!parallel)
      CompLikelihood <- optim(par=param,fn=comploglik, 
                              control=list(factr=1e-10,pgtol=1e-14, maxit=100000), 
                              coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                              fan=fname, lower=lower, method='L-BFGS-B',n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                              upper=upper,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU, hessian=FALSE,MM=MM)
      if(optimizer=='L-BFGS-B'&&parallel){
        ncores=max(1, parallel::detectCores() - 1)
        if(Sys.info()[['sysname']]=="Windows") cl <- parallel::makeCluster(ncores,type = "PSOCK")
        else                                   cl <- parallel::makeCluster(ncores,type = "FORK")
        parallel::setDefaultCluster(cl = cl)
        CompLikelihood <- optimParallel::optimParallel(par=param,fn=comploglik, 
                              control=list(pgtol=1e-14, maxit=100000,factr=1e-10),
                              coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                               data=data, fixed=fixed,fan=fname, lower=lower, method='L-BFGS-B',n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                              upper=upper,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU, hessian=FALSE,MM=MM
                               )
         parallel::setDefaultCluster(cl=NULL)
         parallel::stopCluster(cl)
         }
    if(optimizer=='BFGS') 
        CompLikelihood <- optim(par=param, fn=comploglik,  coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                           control=list(factr=1e-10,
                             reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='BFGS',n=n,namescorr=namescorr,
                                  namesnuis=namesnuis,namesparam=namesparam,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)

      if(optimizer=='Nelder-Mead')
        CompLikelihood <- optim(par=param, fn=comploglik,  coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
          control=list( reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='Nelder-Mead',n=n,namescorr=namescorr,
                                  namesnuis=namesnuis,namesparam=namesparam,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)
    
 #    if(optimizer=='lbfgsb3')
 #      CompLikelihood <-lbfgsb3::lbfgsb3(prm=param,fn=comploglik,gr=NULL,lower=lower,upper=upper,
 #                             #control=list(factr=1e-10,pgtol=1e-14, maxit=100000), 
 #                             coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
 #                             fan=fname,n=n,
 #                             namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
 #                             weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)

    if(optimizer=='nmk')
      CompLikelihood <-dfoptim::nmk(par=param, fn=comploglik, control = list(maxfeval=100000,tol=1e-10),
                           coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed, fan=fname,
                           n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)
    if(optimizer=='nmkb')
    {
      CompLikelihood <-dfoptim::nmkb(par=param, fn=comploglik, control = list(maxfeval=100000,tol=1e-10),
                         lower=lower,upper=upper,
                         coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed, fan=fname,
                         n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)
    }
    if(optimizer=='nlm')
    CompLikelihood <- nlm(f=comploglik,p=param,steptol = 1e-4, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                               fan=fname,hessian=FALSE,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                               iterlim=100000, weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)
  
    if(optimizer=='nlminb')
     CompLikelihood <-nlminb(objective=comploglik,start=param,coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                                control = list( iter.max=100000),
                              lower=lower,upper=upper,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                               weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)

         if(optimizer=='multinlminb'){
       CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn=comploglik,
              coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                               weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM,
                                    lower=lower,upper=upper,method = "nlminb", nbtrials = 500, 
                              control = list( iter.max=100000),
                           typerunif = "sobol"#,nbclusters=2,
                     )
                               }
     if(optimizer=='multiNelder-Mead'){
       CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn=comploglik,
         coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam,MM=MM, 
                               weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,lower=lower,upper=upper,
          method = "Nelder-Mead", nbtrials = 500, 
                              control=list( reltol=1e-14, maxit=100000),
                           typerunif = "sobol"#,nbclusters=4,
                     )
  }

    if(optimizer=='ucminf')   
      CompLikelihood <-ucminf::ucminf(par=param, fn=comploglik, hessian=as.numeric(hessian),   
                        control=list( maxeval=100000),MM=MM,
                            coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed, fan=fname,
                            n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)
                               
    }}
     ############################## bivariate  ############################################                           
    if(bivariate)           {
      
    
     if(length(param)==1)
        {
         optimizer="optimize" 
       CompLikelihood <- optimize(f=comploglik_biv, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed,fan=fname,  lower=lower,n=n,MM=MM,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,maximum = FALSE,
                              upper=upper,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)}
      if(length(param)>1) {   
    if(optimizer=='L-BFGS-B'&&!parallel){
      CompLikelihood <- optim(param,comploglik_biv, control=list(pgtol=1e-14, maxit=100000),
                              method='L-BFGS-B',hessian=FALSE,lower=lower, upper=upper,MM=MM,
                              coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed,fan=fname,n=n,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,
                             weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU )}
     #  if(optimizer=='lbfgsb3')
     # CompLikelihood <- lbfgsb3c::lbfgsb3c(param,comploglik_biv, 
     #                         #control=list(pgtol=1e-14, maxit=100000),
     #                         lower=lower, upper=upper,
     #                         coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
     #                         data=data, fixed=fixed,fan=fname,n=n,
     #                         namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,
     #                        weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU )
    if(optimizer=='L-BFGS-B'&&parallel){
           ncores=max(1, parallel::detectCores() - 1)
           if(Sys.info()[['sysname']]=="Windows") cl <- parallel::makeCluster(ncores,type = "PSOCK")
           else                                   cl <- parallel::makeCluster(ncores,type = "FORK")
           parallel::setDefaultCluster(cl = cl)
           CompLikelihood <- optimParallel::optimParallel(param,comploglik_biv, 
                              control=list(pgtol=1e-14, maxit=100000),
                              coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel,,MM=MM,
                              data=data, fixed=fixed,
                              fan=fname,  n=n, namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,
                              weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU, 
                              lower=lower,upper=upper,
                              hessian=FALSE)
              parallel::setDefaultCluster(cl=NULL)
              parallel::stopCluster(cl)
         }
      if(optimizer=='BFGS')
      CompLikelihood <- optim(param,comploglik_biv, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, control=list(
                              reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='BFGS',n=n,MM=MM,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam ,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
   if(optimizer=='Nelder-Mead')
      CompLikelihood <- optim(param,comploglik_biv, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, control=list(
                              reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='Nelder-Mead',n=n,MM=MM,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam ,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
    if(optimizer=='nmk')
        CompLikelihood <- dfoptim::nmk(par=param, fn=comploglik_biv, control = list(maxfeval=100000,tol=1e-10),
                             coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed, fan=fname,MM=MM,
                            n=n, namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam ,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
    if(optimizer=='nmkb')
        CompLikelihood <- dfoptim::nmkb(par=param, fn=comploglik_biv, control = list(maxfeval=100000,tol=1e-10),
                             lower=lower,upper=upper,
                             coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed, fan=fname,MM=MM,
                            n=n, namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam ,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
     if(optimizer=='nlm') 
        CompLikelihood <- nlm( f=comploglik_biv,p=param, coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                               fan=fname,hessian=FALSE,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam,MM=MM, 
                               weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
    if(optimizer=='ucminf') 
         CompLikelihood <-ucminf::ucminf(par=param, fn=comploglik_biv, hessian=as.numeric(hessian),   
                        control=list( maxeval=100000), 
                        coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                        fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam,MM=MM,
                        weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
    if(optimizer=='nlminb') 
        CompLikelihood <- nlminb( objective=comploglik_biv,start=param, 
                                     control = list( iter.max=100000),
                              lower=lower,upper=upper,
                               coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam,MM=MM,
                               weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU)
          if(optimizer=='multinlminb')
       CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn=comploglik_biv,
         coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                               weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM,
                                    lower=lower,upper=upper,method = "nlminb", nbtrials = 500, 
                              control = list( iter.max=100000),
                           typerunif = "sobol")
          if(optimizer=='multiNelder-Mead')
       CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn=comploglik_biv,
         coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                               weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM,
                                    lower=lower,upper=upper,method = "Nelder-Mead", nbtrials = 500, 
                              control = list( iter.max=100000),
                           typerunif = "sobol")
    
   }}                    
      ########################################################################################   
      ########################################################################################
    # check the optimisation outcome
      if(optimizer=='Nelder-Mead'||optimizer=='multiNelder-Mead'){
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
        if(optimizer=='nmk'||optimizer=='nmkb'){
        CompLikelihood$value = -CompLikelihood$value
        names(CompLikelihood$par)<- namesparam
        if(CompLikelihood$convergence == 0)
        CompLikelihood$convergence <- 'Successful'
        else CompLikelihood$convergence <- "Optimization may have failed"
        if(CompLikelihood$value==-1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
    }
      if(optimizer=='ucminf'){
        CompLikelihood$value = -CompLikelihood$value
        names(CompLikelihood$par)<- namesparam
        if(CompLikelihood$convergence == 1||CompLikelihood$convergence == 2||CompLikelihood$convergence == 4)
        CompLikelihood$convergence <- 'Successful'
        else
        if(CompLikelihood$convergence == 3)
        CompLikelihood$convergence <- 'Iteration limit reached'
        else
        CompLikelihood$convergence <- "Optimization may have failed"
        if(CompLikelihood$value==-1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
    }
    if(optimizer=='L-BFGS-B'||optimizer=='BFGS'||optimizer=='lbfgsb3c'){
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
        if(CompLikelihood$code == 1|| CompLikelihood$code == 2)
        CompLikelihood$convergence <- 'Successful'
        else
        if(CompLikelihood$code == 4)
        CompLikelihood$convergence <- 'Iteration limit reached'
        else
        CompLikelihood$convergence <- "Optimization may have failed"
        if(CompLikelihood$value==-1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
    }

    if(optimizer=='nlminb'||optimizer=='multinlminb' ){
        CompLikelihood$par <- CompLikelihood$par
        names(CompLikelihood$par)<- namesparam
        CompLikelihood$value <- -CompLikelihood$objective
        if(CompLikelihood$convergence == 0) { CompLikelihood$convergence <- 'Successful' }
        else {CompLikelihood$convergence <- "Optimization may have failed" }
        if(CompLikelihood$objective==-1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
    }
    if(optimizer=='optimize'){
    param<-CompLikelihood$minimum
    CompLikelihood$par<-param  
    names(CompLikelihood$par)<- namesparam
    maxfun <- -CompLikelihood$objective
    CompLikelihood$value <- maxfun
    CompLikelihood$convergence <- 'Successful'
    }
  } ##### end if !onlyvar
    else {
          CompLikelihood=as.list(0)
          names(CompLikelihood)="value"
          CompLikelihood$par <- param
          CompLikelihood$claic <- NULL;CompLikelihood$clbic <- NULL;
          CompLikelihood$convergence <- 'Successful'
          if(!bivariate)CompLikelihood$value = - comploglik(param=CompLikelihood$par ,  coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed, fan=fname,
                             n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)
          else CompLikelihood$value = -comploglik_biv(param=CompLikelihood$par ,  coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, data=data, fixed=fixed, fan=fname,
                             n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)
          

          if(hessian) 
          {
               if(!bivariate)  
                CompLikelihood$hessian=numDeriv::hessian(func=comploglik,x=param,method="Richardson",  coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed,fan=fname,n=n,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,
                              weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)
               if(bivariate)  
               CompLikelihood$hessian=numDeriv::hessian(func=comploglik_biv,x=param,method="Richardson",coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                             data=data, fixed=fixed,fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                             weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)
               rownames(CompLikelihood$hessian)=namesparam
               colnames(CompLikelihood$hessian)=namesparam
          }
  }

#####################################
#if((sensitivity||varest)&&(is.null(CompLikelihood$hessian)||min(eigen(CompLikelihood$hessian)$values)<0))
if((sensitivity||varest))
  {
if(!bivariate)  

CompLikelihood$hessian=numDeriv::hessian(func=comploglik,x=CompLikelihood$par,method="Richardson",  coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                              data=data, fixed=fixed,fan=fname,n=n,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,
                              weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)
if(bivariate)  
CompLikelihood$hessian=numDeriv::hessian(func=comploglik_biv,x=CompLikelihood$par,method="Richardson",coordx=coordx, coordy=coordy, coordt=coordt,corrmodel=corrmodel, 
                             data=data, fixed=fixed,fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, 
                               weigthed=weigthed,X=X,ns=ns,NS=NS,local=local,GPU=GPU,MM=MM)
rownames(CompLikelihood$hessian)=namesparam
colnames(CompLikelihood$hessian)=namesparam
  }


####################################
        if( (CompLikelihood$convergence!='Successful')||CompLikelihood$value==-1e+15)  print("Optimization failed: try with other starting values ")
          else{
    if(varest)
          {
        
            dimmat <- numparam^2
            dmat <- numparam*(numparam+1)/2
            eps <- (.Machine$double.eps)^(1/3)
            param <- c(CompLikelihood$par, fixed)
            score <- double(numparam)
            paramcorr <- param[namescorr]
            nuisance <- param[namesnuis]
            sel=substr(names(nuisance),1,4)=="mean"
            mm=as.numeric(nuisance[sel])   ## mean paramteres

            if(bivariate){
            sel1=substr(names(nuisance),1,6)=="mean_1"
            mm1=as.numeric(nuisance[sel1])
            sel2=substr(names(nuisance),1,6)=="mean_2"
            mm2=as.numeric(nuisance[sel2])
            mm=c(mm1,mm2)}

            num_betas=length(mm)
            other_nuis=as.numeric(nuisance[!sel]) 
            nuisance=c(mm,other_nuis)
            sensmat <- double(dmat);varimat <- double(dmat)


            # Set the window parameter:
           if(length(winconst)==1) winconst=c(winconst,0)

 
            GD=.C('GodambeMat',as.double(mm),as.integer(bivariate),as.double(coordx),as.double(coordy),
              as.double(coordt),as.integer(corrmodel), as.double(data),as.integer(distance),as.double(eps),
              as.integer(flagcorr), as.integer(flagnuis),as.integer(grid),as.integer(likelihood),
              as.double(c(X%*%mm)),as.integer(model),as.double(n),as.integer(num_betas),
              as.integer(numparam),as.integer(numparamcorr),as.integer(length(paramcorr)),as.double(paramcorr),as.double(nuisance),
              score=score,sensmat=sensmat,as.integer(spacetime),as.integer(type),
              varimat=varimat,as.integer(vartype),as.double(winconst),as.double(winstp),as.double(winconst_t),as.double(winstp_t),
              as.integer(weigthed),c(t(X)),as.integer(ns),as.integer(NS),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)
            
  
            if(!sum(GD$varimat)) print("Std error estimation failed")
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
            
            #namesgod=namesparam
            #names(CompLikelihood$score )=namesparam
#print(namesgod)
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

