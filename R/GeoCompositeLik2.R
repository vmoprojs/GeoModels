####################################################
### File name: CompLik2.r
####################################################

### Optim call for Composite log-likelihood maximization

CompLik2 <- function(copula,bivariate, coordx, coordy ,coordz,coordt,coordx_dyn,corrmodel, data, distance, flagcorr, flagnuis, fixed, GPU,grid,
                           likelihood, local,lower, model, n, namescorr, namesnuis, namesparam,
                           numparam, numparamcorr, optimizer, onlyvar, parallel, param, spacetime, type,
                           upper, varest,  weigthed, ns, X,sensitivity,
                           colidx,rowidx,neighb,MM,aniso)
  {


comploglik2 <- function(param,colidx,rowidx, corrmodel, coords,data1,data2,fixed, fan, n, namescorr, 
                              namesnuis,namesparam,namesaniso,weigthed,X,GPU,local,MM,aniso,type_cop,cond_pair)
      {
        names(param) <- namesparam
        param <- c(param, fixed)
        paramcorr <- as.numeric(param[namescorr])
        nuisance <- param[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"

        if(is.null(MM)){ mm=as.numeric(nuisance[sel]) ### linear mean
                         Mean=c(X%*%mm)
                       }
        else           Mean=MM                     ### fixed mean

        other_nuis=as.numeric(nuisance[!sel])   ## or nuis parameters (nugget sill skew df)         

############################################
#if(!type_cop) { # not copula models 

        if(aniso){
            anisopar<-param[namesaniso]
            coords1=GeoAniso(coords, anisopars=anisopar)
          c1=c(t(coords1[colidx,]));c2=c(t(coords1[rowidx,]))
          result=dotCall64::.C64(as.character(fan),
         SIGNATURE = c("integer","double","double","double","double", "integer","integer","double","integer","double","double","double","double","integer","integer","integer","integer"),  
                        corrmodel,c1,c2,data1, data2, n[colidx],n[rowidx],paramcorr,weigthed, res= dotCall64::numeric_dc(1),Mean[colidx], Mean[rowidx], other_nuis,local,GPU,
                        type_cop,cond_pair,
          INTENT =    c(rep("r",9),"w",rep("r",7)),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$res


        }
         else
         {
            #print(fan)
                  
         result=dotCall64::.C64(as.character(fan),
         SIGNATURE = c("integer","double","double", "integer","integer","double","integer","double","double","double","double","integer","integer","integer","integer"),  
                        corrmodel,data1, data2, n[colidx],n[rowidx],paramcorr,weigthed, res= dotCall64::numeric_dc(1),Mean[colidx], Mean[rowidx], other_nuis,local,GPU,type_cop,cond_pair,
          INTENT =    c(rep("r",7),"w",rep("r",7)),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$res
       }
return(-result)
      }
##################################################################
##################################################################
comploglik_biv2 <- function(param,colidx,rowidx, corrmodel, coords,data1,data2,fixed, fan, n, 
                          namescorr, namesnuis,namesparam,namesaniso,weigthed,X,GPU,local,MM,aniso,type_cop,cond_pair)
      {

        names(param) <- namesparam
        param <- c(param, fixed)
        paramcorr <- param[namescorr]
        nuisance <- param[namesnuis]
        anisopar<-param[namesaniso]

        sel1=substr(names(nuisance),1,6)=="mean_1"
        mm1=as.numeric(nuisance[sel1])
        sel2=substr(names(nuisance),1,6)=="mean_2"
        mm2=as.numeric(nuisance[sel2])
        sel=substr(names(nuisance),1,4)=="mean"
        X1=as.matrix(X[1:ns[1],]);X2=as.matrix(X[(ns[1]+1):(ns[2]+ns[1]),]); 
        other_nuis=as.numeric(nuisance[!sel]) 
        Mean=c(X1%*%mm1,X2%*%mm2)
        res=double(1)
       
        result=dotCall64::.C64(as.character(fan),
          SIGNATURE = c("integer","double","double", "integer","double","integer","double","double","double","double","integer","integer"),  
                        corrmodel,data1, data2, n,paramcorr,weigthed, res=res,Mean[colidx],Mean[rowidx],other_nuis,local,GPU,
          INTENT =    c(rep("r",6),"w", rep("r", 5)),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$res
        return(-result)
      }

   ##################################################################################################
   ############### starting functtion ###############################################################
   ##################################################################################################

    numcoord=length(coordx);numtime=1;spacetime_dyn=FALSE
    if(spacetime) numtime=length(coordt)
    if(bivariate) numtime=2
    if(!is.null(coordx_dyn)) spacetime_dyn=TRUE

    dimat=numcoord*numtime#
    if((spacetime||bivariate)) dimat <- sum(ns)
    NS=cumsum(ns)
    if(is.null(dim(X))){X=as.matrix(rep(1,dimat))}
    
    fname <- NULL; hessian <- FALSE
    if(all(model==1,likelihood==4,type==2)) fname <- 'Comp_Diff_Gauss'
    
    namesaniso=c("angle","ratio")

   
   if(length(n)==1) n=rep(n,dimat)
####################### conditional ##############################################
    if(all(model==1, likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss'                          
    if(all(model==21,likelihood==1,type==2)) fname <- 'Comp_Cond_Gamma'                                 
    if(all(model==26,likelihood==1,type==2)) fname <- 'Comp_Cond_Weibull'
    if(all(model==22,likelihood==1,type==2)) fname <- 'Comp_Cond_LogGauss'
    if(all(model==28,likelihood==1,type==2)) fname <- 'Comp_Cond_Beta'
    if(all(model==50,likelihood==1,type==2)) fname <- 'Comp_Cond_Beta2' 
    if(all(model==13,likelihood==1,type==2)) fname <- 'Comp_Cond_WrapGauss'   
    if(all(model==33,likelihood==1,type==2)) fname <- 'Comp_Cond_Kumaraswamy'
    if(all(model==42,likelihood==1,type==2)) fname <- 'Comp_Cond_Kumaraswamy2'                                          
    if(all(model==20,likelihood==1,type==2)) fname <- 'Comp_Cond_SinhGauss'                                           
    if(all(model==38,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECETukeyh'                                          
    if(all(model==27,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECET'                                          
    if(all(model==29,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECEGauss'                                         
    if(all(model==10,likelihood==1,type==2)) fname <- 'Comp_Cond_SkewGauss'                                        
    if(all(model==12,likelihood==1,type==2)) fname <- 'Comp_Cond_T'                                          
    if(all(model==34,likelihood==1,type==2)) fname <- 'Comp_Cond_Tukeyh'                                          
    if(all(model==40,likelihood==1,type==2)) fname <- 'Comp_Cond_Tukeyhh'                                         
    if(all(model==35,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss_misp_T'                                         
    if(all(model==27,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECET'                                         
    if(all(model==39,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECEBIMODAL'                                        
    if(all(model==29,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECEGauss'                                        
    if(all(model==31,likelihood==1,type==2)) fname <- 'Comp_Cond_BinomTWOPIECEGauss'                                        
    if(all(model==32,likelihood==1,type==2)) fname <- 'Comp_Cond_BinomnegTWOPIECEGauss'                                       
    if(all(model==38,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECETukeyh'                                        
    if(all(model==30,likelihood==1,type==2)) fname <- 'Comp_Cond_Pois'                                      
    if(all(model==46,likelihood==1,type==2)) fname <- 'Comp_Cond_PoisGamma'                                       
    if(all(model==36,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss_misp_Pois'
    if(all(model==2,likelihood==1,type==2))  fname <- 'Comp_Cond_BinomGauss'                                 
    if(all(model==11,likelihood==1,type==2)&length(n)==1) fname <- 'Comp_Cond_BinomGauss'                                 
    if(all(model==51,likelihood==1,type==2)&length(n)>=1) fname <- 'Comp_Cond_BinomNNGauss_misp'                                
    if(all(model==11,likelihood==1,type==2)&length(n)>1) fname <- 'Comp_Cond_BinomNNGauss'
    if(all(model==49,likelihood==1,type==2)&length(n)==1) fname <- 'Comp_Cond_BinomLogi'                        
    if(all(model==49,likelihood==1,type==2)&length(n)>1) fname <- 'Comp_Cond_BinomNNLogi'                                                              
    if(all(model==14,likelihood==1,type==2)) fname <- 'Comp_Cond_BinomnegGauss'                                     
    if(all(model==16,likelihood==1,type==2)) fname <- 'Comp_Cond_BinomnegGauss'                                      
    if(all(model==54,likelihood==1,type==2)) fname <- 'Comp_Cond_BinomnegBinary'                                      
    if(all(model==27,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECET'                                        
    if(all(model==39,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECEBIMODAL'                                       
    if(all(model==29,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECEGauss'                                       
    if(all(model==38,likelihood==1,type==2)) fname <- 'Comp_Cond_TWOPIECETukeyh'                                      
    if(all(model==43,likelihood==1,type==2)) fname <- 'Comp_Cond_PoisZIP'                                     
    if(all(model==57,likelihood==1,type==2)) fname <- 'Comp_Cond_PoisGammaZIP'                                                                                 
    if(all(model==44,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss_misp_PoisZIP'                                       
    if(all(model==45,likelihood==1,type==2)) fname <- 'Comp_Cond_BinomnegGaussZINB'                                       
    if(all(model==24,likelihood==1,type==2)) fname <- 'Comp_Cond_LogLogistic'                                        
    if(all(model==25,likelihood==1,type==2)) fname <- 'Comp_Cond_Logistic'                                       
    if(all(model==41,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss_misp_Tukeygh'                                         
    if(all(model==37,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss_misp_SkewT'                                        
    if(all(model==47,likelihood==1,type==2)) fname <- 'Comp_Cond_Gauss_misp_PoisGamma'                          
###################### pairwise ###############################################
    if(all(model==1,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss'                                        
    if(all(model==2,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomGauss'                                       
    if(all(model==11,likelihood==3,type==2)&length(n)==1) fname <- 'Comp_Pair_BinomGauss'                                       
    if(all(model==51,likelihood==3,type==2)&length(n)>=1) fname <- 'Comp_Pair_BinomNNGauss_misp'                                       
    if(all(model==11,likelihood==3,type==2)&length(n)>1) fname <- 'Comp_Pair_BinomNNGauss'                                       
    if(all(model==49,likelihood==3,type==2)&length(n)==1) fname <- 'Comp_Pair_BinomLogi'                                       
    if(all(model==49,likelihood==3,type==2)&length(n)>1) fname <- 'Comp_Pair_BinomNNLogi'                                      
    if(all(model==19,likelihood==3,type==2)) { namesnuis=c(namesnuis,"z")
                                              fixed<- c(fixed, list(z=min(n)))
                                              fname <- 'Comp_Pair_Binom2Gauss'}                                      
    if(all(model==14,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomnegGauss'                                      
    if(all(model==16,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomnegGauss'                                      
    if(all(model==54,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomnegBinary'                                      
    if(all(model==15,likelihood==3,type==2)) fname <- 'Comp_Pair_PoisbinGauss'                                    
    if(all(model==17,likelihood==3,type==2)) fname <- 'Comp_Pair_PoisbinnegGauss'                                         
    if(all(model==13,likelihood==3,type==2)) fname <- 'Comp_Pair_WrapGauss'                                     
    if(all(model==10,likelihood==3,type==2)) fname <- 'Comp_Pair_SkewGauss'
    if(all(model==21,likelihood==3,type==2)) fname <- 'Comp_Pair_Gamma'                                         
    if(all(model==33,likelihood==3,type==2)) fname <- 'Comp_Pair_Kumaraswamy'                                         
    if(all(model==42,likelihood==3,type==2)) fname <- 'Comp_Pair_Kumaraswamy2'                                         
    if(all(model==28,likelihood==3,type==2)) fname <- 'Comp_Pair_Beta'                                        
    if(all(model==50,likelihood==3,type==2)) fname <- 'Comp_Pair_Beta2'                                        
    if(all(model==26,likelihood==3,type==2)) fname <- 'Comp_Pair_Weibull'                                                                             
    if(all(model==24,likelihood==3,type==2)) fname <- 'Comp_Pair_LogLogistic'                                            
    if(all(model==25,likelihood==3,type==2)) fname <- 'Comp_Pair_Logistic'                                                                                                                      
    if(all(model==23,likelihood==3,type==2)) fname <- 'Comp_Pair_2Gamma'                                                                               
    if(all(model==22,likelihood==3,type==2)) fname <- 'Comp_Pair_LogGauss';                                       
    if(all(model==18,likelihood==3,type==2)) fname <- 'Comp_Pair_SkewTGauss'                                         
    if(all(model==27,likelihood==3,type==2)) fname <- 'Comp_Pair_TWOPIECET'                                          
    if(all(model==39,likelihood==3,type==2)) fname <- 'Comp_Pair_TWOPIECEBIMODAL'                                         
    if(all(model==29,likelihood==3,type==2)) fname <- 'Comp_Pair_TWOPIECEGauss'                                         
    if(all(model==31,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomTWOPIECEGauss'                                        
    if(all(model==32,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomnegTWOPIECEGauss'                                       
    if(all(model==12,likelihood==3,type==2)) fname <- 'Comp_Pair_T'                                         
    if(all(model==34,likelihood==3,type==2)) fname <- 'Comp_Pair_Tukeyh'                                         
    if(all(model==40,likelihood==3,type==2)) fname <- 'Comp_Pair_Tukeyhh'                                        
    if(all(model==41,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss_misp_Tukeygh'                                        
    if(all(model==36,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss_misp_Pois'                                      
    if(all(model==35,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss_misp_T'                                       
    if(all(model==37,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss_misp_SkewT'                                         
    if(all(model==20,likelihood==3,type==2)) fname <- 'Comp_Pair_SinhGauss'                                        
    if(all(model==38,likelihood==3,type==2)) fname <- 'Comp_Pair_TWOPIECETukeyh'                                         
    if(all(model==30,likelihood==3,type==2)) fname <- 'Comp_Pair_Pois'                                        
    if(all(model==46,likelihood==3,type==2)) fname <- 'Comp_Pair_PoisGamma'                                        
    if(all(model==43,likelihood==3,type==2)) fname <- 'Comp_Pair_PoisZIP'                                         
    if(all(model==57,likelihood==3,type==2)) fname <- 'Comp_Pair_PoisGammaZIP'                                         
    if(all(model==44,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss_misp_PoisZIP'                                        
    if(all(model==45,likelihood==3,type==2)) fname <- 'Comp_Pair_BinomnegGaussZINB'                                        
    if(all(model==47,likelihood==3,type==2)) fname <- 'Comp_Pair_Gauss_misp_PoisGamma'
 #############################################################################
    if(sensitivity) hessian=TRUE
    if(spacetime) fname <- paste(fname,"_st",sep="")
    if(bivariate) fname <- paste(fname,"_biv",sep="")

    type_cop=0
    cond_pair=0

   
  if(!is.null(copula))
    {
        fname <- paste(fname,"Cop",sep="");
        if(copula=="Gaussian") {type_cop=1; }
        if(copula=="Clayton")  {type_cop=2; }
        if(all(likelihood==1,type==2)) {cond_pair=1;  fname=gsub("Cond", "Pair", fname)}
        if(all(likelihood==3,type==2)) {cond_pair=0}
    }

    fname <- paste(fname,"2mem",sep="")

    if(aniso) fname <- paste(fname,"_aniso",sep="")
     
   
    #path.parent <- getwd()

  
    # if(!is.null(GPU)) 
    # {
    #   path <- system.file("CL", "Kernel.cl", package = "GeoModels")
    #   path <- gsub("/Kernel.cl","/",path);setwd(path)
    #   fname <- paste(fname,"_OCL",sep="")
    #   .C("create_binary_kernel",  as.integer(GPU),as.character(fname),  PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)
    # }
    if(!is.null(GPU))
    {
      # fname <- paste(fname,"_OCL",sep="")
      # #cat("fname de Composit.r: ",fname,"\n")
      # 
      # path <- system.file("CL", paste(fname,".cl",sep = ""), package = "GeoModels")
      # path <- gsub(paste("/",paste(fname,".cl",sep = ""),sep = ""),"/",path)
      # # .C("create_binary_kernel",  as.integer(GPU),as.character(fname),  PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)
      # setwd(path)
      # .C("create_binary_kernel",  as.integer(GPU),as.character(fname),  PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)
    }

  
     if((spacetime||bivariate)&&(!spacetime_dyn))    data=c(t(data))
     if((spacetime||bivariate)&&(spacetime_dyn))     data=unlist(data)          
     if(spacetime||bivariate)   NS=c(0,NS)[-(length(ns)+1)]

###### selectin data with indexes from composite likelihood
 
   if(is.null(neighb)) {colidx=colidx+1; rowidx=rowidx+1}  #updating if #using "my distances from C" 
   data1=data[colidx]; data2=data[rowidx]                  ##using "RANN distances" 
 

   if((model==11||model==49||model==51)&&length(n)>1) {n1=n[colidx];n2=n[rowidx];n=c(n1,n2)} ## for binomials type models
   if(is.null(GPU)) GPU=0



##################
#if(!bivariate){
#if(is.null(MM)) lname="comploglik2"
#else            lname="comploglik2MM"
#}

lname="comploglik2"


 
coords=cbind(coordx,coordy,coordz)

if(!onlyvar){
  
#ptm=proc.time()
  ##############################.  spatial or space time ############################################
   if(!bivariate)           {
    if(length(param)==1) {
         optimizer="optimize"  
     CompLikelihood <- optimize(f=eval(as.name(lname)), colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords,
                              data1=data1,data2=data2, fixed=fixed, fan=fname,  lower=lower, n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,
                               maximum = FALSE,
                              upper=4,weigthed=weigthed,X=X, local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)}
   if(length(param)>1) {

    if(optimizer=='L-BFGS-B'&&!parallel)
      CompLikelihood <- optim(par=param,fn=eval(as.name(lname)), 
                              control=list(factr=1e-10,pgtol=1e-14, maxit=100000), 
                              colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords,data1=data1,data2=data2, fixed=fixed,
                              fan=fname, lower=lower, method='L-BFGS-B',n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
                              upper=upper,weigthed=weigthed,X=X, local=local,GPU=GPU, hessian=FALSE,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
    
    if(optimizer=='L-BFGS-B'&&parallel){
        #ncores=max(1, parallel::detectCores() - 1)
        ncores=length(param) * 2 + 1
        if(Sys.info()[['sysname']]=="Windows") cl <- parallel::makeCluster(ncores,type = "PSOCK")
        else                                   cl <- parallel::makeCluster(ncores,type = "FORK")
        parallel::setDefaultCluster(cl = cl)
        CompLikelihood <- optimParallel::optimParallel(par=param,fn=eval(as.name(lname)), 
                              #control=list(pgtol=1e-14, maxit=100000,factr = 1e8), # factr = 1e-10
                                control=list(factr=1e-10,pgtol=1e-14, maxit=100000), 
                              colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                               data1=data1,data2=data2,fixed=fixed,fan=fname, lower=lower,n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
                         parallel = list(forward = FALSE,loginfo=FALSE),
                          upper=upper,weigthed=weigthed,X=X, local=local,GPU=GPU, hessian=FALSE,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
         parallel::setDefaultCluster(cl=NULL)
         parallel::stopCluster(cl)
         }
    if(optimizer=='BFGS') 
        CompLikelihood <- optim(par=param, fn=eval(as.name(lname)),  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                           control=list(factr=1e-10,
                             reltol=1e-14, maxit=100000),data1=data1,data2=data2, fixed=fixed, fan=fname,
                              hessian=FALSE, method='BFGS',n=n,namescorr=namescorr,
                                  namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,weigthed=weigthed,X=X, local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
    if(optimizer=='SANN') 
        CompLikelihood <- optim(par=param, fn=eval(as.name(lname)),  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                           control=list(factr=1e-10,
                             reltol=1e-14, maxit=100000),data1=data1,data2=data2, fixed=fixed, fan=fname,
                              hessian=FALSE, method='SANN',n=n,namescorr=namescorr,
                                  namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,weigthed=weigthed,X=X, local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)

      if(optimizer=='Nelder-Mead')
        CompLikelihood <- optim(par=param, fn=eval(as.name(lname)),  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
          control=list( reltol=1e-14, maxit=100000), data1=data1,data2=data2, fixed=fixed, fan=fname,
                              hessian=FALSE, method='Nelder-Mead',n=n,namescorr=namescorr,
                                  namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,weigthed=weigthed,X=X, local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)

 #if(optimizer=='multinlminb'){
  #     CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn=eval(as.name(lname)),
   #     colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords, data1=data1,data2=data2, fixed=fixed,
    #                           fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
     #                          weigthed=weigthed,X=X, local=local,GPU=GPU,MM=MM,aniso=aniso,
      #    lower=lower,upper=upper,method = "nlminb", nbtrials = 400, typerunif = "sobol",
       #                       control = list( iter.max=100000))
  #}
 #if(optimizer=='multiNelder-Mead'){
  #     CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn=eval(as.name(lname)),
   #     colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,data1=data1,data2=data2, fixed=fixed,
    #                           fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
     #                          weigthed=weigthed,X=X, local=local,GPU=GPU,lower=lower,upper=upper,MM=MM,aniso=aniso,
      #    method = "Nelder-Mead", nbtrials = 500, 
       #                       control=list( reltol=1e-14, maxit=100000),
        #                   typerunif = "sobol"#,nbclusters=4,
         #            )
  #}


    #if(optimizer=='nmk')
     # CompLikelihood <-dfoptim::nmk(par=param, fn=eval(as.name(lname)), control = list(maxfeval=100000,tol=1e-10),
      #                    colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords,data1=data1,data2=data2,fixed=fixed, fan=fname,
       #                    n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,weigthed=weigthed,X=X, local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
    #if(optimizer=='nmkb')
    #{
     # CompLikelihood <-dfoptim::nmkb(par=param, fn=eval(as.name(lname)), control = list(maxfeval=100000,tol=1e-10),
      #                   lower=lower,upper=upper,
       #                  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,data1=data1,data2=data2, fixed=fixed, fan=fname,
        #                 n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,weigthed=weigthed,X=X, local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
    #}
    if(optimizer=='nlm')
    CompLikelihood <- nlm(f=eval(as.name(lname)),p=param,steptol = 1e-4, colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,data1=data1,data2=data2, fixed=fixed,
                               fan=fname,hessian=FALSE,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso, 
                               iterlim=100000, weigthed=weigthed,X=X, local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
  

     if(optimizer=='bobyqa')   
  {
     CompLikelihood <-minqa::bobyqa(par=param, fn=eval(as.name(lname)),lower=lower,upper=upper,  
                        control = list( maxfun=100000),
                        colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords, data1=data1,data2=data2, fixed=fixed,fan=fname,n=n,namescorr=namescorr, 
                       namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso, weigthed=weigthed,X=X, local=local,GPU=GPU,MM=MM,
                        aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
    }

    if(optimizer=='nlminb'){

    # tt1 <- proc.time() 


     CompLikelihood <-nlminb(objective=eval(as.name(lname)),start=param,colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords, data1=data1,data2=data2, fixed=fixed,
                                control = list( iter.max=100000),
                              lower=lower,upper=upper,hessian=FALSE,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
                               weigthed=weigthed,X=X, local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
     # tt1 <- proc.time()-tt1;
    }
 

                           

    #  if(optimizer=='sa')
    #    CompLikelihood <- optimization::optim_sa(start=param, fun=eval(as.name(lname)), maximization = FALSE, 
     #     lower=lower,upper=upper,
    #       colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, 
     #     #control=list( reltol=1e-14, maxit=100000), 
      #    data1=data1,data2=data2, fixed=fixed, fan=fname,
       #                     n=n,namescorr=namescorr,
        #                          namesnuis=namesnuis,namesparam=namesparam,weigthed=weigthed,X=X, local=local,GPU=GPU,MM=MM,aniso=aniso)  

       # bb=data1*data2/length(data1)                     
    }}
######################################################################################
############################## bivariate  ############################################ 
######################################################################################                          
    if(bivariate)           { 
     if(length(param)==1)
        {
         optimizer="optimize" 
       CompLikelihood <- optimize(f=comploglik_biv2,  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                              data1=data1,data2=data2, fixed=fixed,fan=fname,  lower=lower,n=n,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,namesaniso=namesaniso,maximum = FALSE,
                              upper=upper,weigthed=weigthed,X=X,local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)}
      if(length(param)>1) {   
    if(optimizer=='L-BFGS-B'&&!parallel){
      CompLikelihood <- optim(param,comploglik_biv2, control=list(pgtol=1e-14, maxit=100000),
                              method='L-BFGS-B',hessian=FALSE,lower=lower, upper=upper,
                               colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                              data1=data1,data2=data2, fixed=fixed,fan=fname,n=n,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,namesaniso=namesaniso,
                                   parallel = list(forward = FALSE,loginfo=FALSE),
                             weigthed=weigthed,X=X,local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)}
  
    if(optimizer=='L-BFGS-B'&&parallel){
           ncores=max(1, parallel::detectCores() - 1)
           if(Sys.info()[['sysname']]=="Windows") cl <- parallel::makeCluster(ncores,type = "PSOCK")
           else                                   cl <- parallel::makeCluster(ncores,type = "FORK")
           parallel::setDefaultCluster(cl = cl)
           CompLikelihood <- optimParallel::optimParallel(param,comploglik_biv2, 
                              control=list(pgtol=1e-14, maxit=100000),
                               colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                              data1=data1,data2=data2, fixed=fixed,
                              fan=fname,  n=n, namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,namesaniso=namesaniso,
                              weigthed=weigthed,X=X,local=local,GPU=GPU, 
                              lower=lower,upper=upper,
                              hessian=FALSE,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
              parallel::setDefaultCluster(cl=NULL)
              parallel::stopCluster(cl)
         }
      if(optimizer=='BFGS')
      CompLikelihood <- optim(param,comploglik_biv2,  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,control=list(
                              reltol=1e-14, maxit=100000), data1=data1,data2=data2, fixed=fixed, fan=fname,
                              hessian=FALSE, method='BFGS',n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam ,namesaniso=namesaniso,weigthed=weigthed,X=X,local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
   if(optimizer=='Nelder-Mead')
      CompLikelihood <- optim(param,comploglik_biv2,  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,control=list(
                              reltol=1e-14, maxit=100000), data1=data1,data2=data2, fixed=fixed, fan=fname,
                              hessian=FALSE, method='Nelder-Mead',n=n,
                              namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam ,namesaniso=namesaniso,weigthed=weigthed,X=X,local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
    if(optimizer=='nlm') 
        CompLikelihood <- nlm( f=comploglik_biv2,p=param,  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords, data1=data1,data2=data2, fixed=fixed,
                               fan=fname,hessian=FALSE,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
                               weigthed=weigthed,X=X,local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
    
   if(optimizer=='bobyqa')   {
        CompLikelihood <- minqa::bobyqa(fn=comploglik_biv2,par=param, 
                                     control = list(iter.max=100000),
                              lower=lower,upper=upper,
                                colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,data1=data1,data2=data2, fixed=fixed,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
                               weigthed=weigthed,X=X,local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
                      }
    if(optimizer=='nlminb') 
        CompLikelihood <- nlminb( objective=comploglik_biv2,start=param, 
                                     control = list(iter.max=100000),
                              lower=lower,upper=upper, hessian=FALSE,
                                colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,data1=data1,data2=data2, fixed=fixed,
                               fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
                               weigthed=weigthed,X=X,local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
     #if(optimizer=='multinlminb'){
      # CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn=comploglik_biv2,
       # colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords, data1=data1,data2=data2, fixed=fixed,
        #                       fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
         #                      weigthed=weigthed,X=X,local=local,GPU=GPU,,MM=MM,aniso=aniso,
          #                          lower=lower,upper=upper,method = "nlminb", nbtrials = 500, 
           #                   control = list( iter.max=100000),
            #               typerunif = "sobol")
             #                  }
    # if(optimizer=='multiNelder-Mead'){
     #  CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn=comploglik_biv2,
      #  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel, coords=coords, data1=data1,data2=data2, fixed=fixed,
       #                        fan=fname,n=n,namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
        #                       weigthed=weigthed,X=X, local=local,GPU=GPU,lower=lower,upper=upper,,MM=MM,aniso=aniso,
         # method = "Nelder-Mead", nbtrials = 500, 
          #                    control=list( reltol=1e-14, maxit=100000),
           #                typerunif = "sobol")
  #}

   }
 }  
 

      ########################################################################################   
      ########################################################################################
    # check the optimisation outcome
      if(optimizer=='Nelder-Mead'||optimizer=='multiNelder-Mead'||optimizer=='SANN'){
        CompLikelihood$value = -CompLikelihood$value
        names(CompLikelihood$par)<- namesparam
        if(CompLikelihood$convergence == 0)
        CompLikelihood$convergence <- 'Successful'
        else
        if(CompLikelihood$convergence == 1)
        CompLikelihood$convergence <- 'Iteration limit reached'
        else
        CompLikelihood$convergence <- "Optimization may have failed"
        if(CompLikelihood$value>=1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
        CompLikelihood$counts=as.numeric(CompLikelihood$counts[1])
    }
        if(optimizer=='nmk'||optimizer=='nmkb'){
        CompLikelihood$value = -CompLikelihood$value
        names(CompLikelihood$par)<- namesparam
        if(CompLikelihood$convergence == 0)
        CompLikelihood$convergence <- 'Successful'
        else CompLikelihood$convergence <- "Optimization may have failed"
        if(CompLikelihood$value>=1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
        CompLikelihood$counts=as.numeric(CompLikelihood$feval)
    }
    #  if(optimizer=='ucminf'){
     #   CompLikelihood$value = -CompLikelihood$value
      #  names(CompLikelihood$par)<- namesparam
      #  if(CompLikelihood$convergence == 1||CompLikelihood$convergence == 2||CompLikelihood$convergence == 4)
      #  CompLikelihood$convergence <- 'Successful'
      #  else
      #  if(CompLikelihood$convergence == 3)
      #  CompLikelihood$convergence <- 'Iteration limit reached'
      #  else
      #  CompLikelihood$convergence <- "Optimization may have failed"
      #  if(CompLikelihood$value>=1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
   # }
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
        if(CompLikelihood$value>=1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
        CompLikelihood$counts=as.numeric(CompLikelihood$counts[1])
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
        if(CompLikelihood$value>= 1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
        CompLikelihood$counts=as.numeric(CompLikelihood$iterations)
    }


     if(optimizer=='bobyqa'){

     
        CompLikelihood$par <- CompLikelihood$par
        names(CompLikelihood$par)<- namesparam
        CompLikelihood$value <- -CompLikelihood$fval
        if(CompLikelihood$ierr == 0) { CompLikelihood$convergence <- 'Successful' }
        else {CompLikelihood$convergence <- "Optimization may have failed" }
        if(CompLikelihood$fval >= 1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
        CompLikelihood$counts=as.numeric(CompLikelihood$feval)
    }



    if(optimizer=='nlminb'||optimizer=='multinlminb'){

     
        CompLikelihood$par <- CompLikelihood$par
        names(CompLikelihood$par)<- namesparam
        CompLikelihood$value <- -CompLikelihood$objective
        if(CompLikelihood$convergence == 0) { CompLikelihood$convergence <- 'Successful' }
        else {CompLikelihood$convergence <- "Optimization may have failed" }
        if(CompLikelihood$objective>= 1.0e8) CompLikelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
        CompLikelihood$counts=as.numeric(CompLikelihood$iterations)
    }
    if(optimizer=='optimize'){
    param<-CompLikelihood$minimum
    CompLikelihood$par<-param  
    names(CompLikelihood$par)<- namesparam
    maxfun <- -CompLikelihood$objective
    CompLikelihood$value <- maxfun
    CompLikelihood$convergence <- 'Successful'
       CompLikelihood$counts=NULL
    }
  } ##### end if !onlyvar
    else {
          CompLikelihood=as.list(0)
          names(CompLikelihood)="value"
          CompLikelihood$par <- param
          CompLikelihood$claic <- NULL;CompLikelihood$clbic <- NULL;
          CompLikelihood$convergence <- 'Successful'
          if(!bivariate) CompLikelihood$value = - comploglik2(param=CompLikelihood$par ,  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                              data1=data1,data2=data2, fixed=fixed, fan=fname,type_cop=type_cop,
                             n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,weigthed=weigthed,X=X, local=local,GPU=GPU)
          else CompLikelihood$value = -comploglik_biv2(param=CompLikelihood$par ,  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                data1=data1,data2=data2, fixed=fixed, fan=fname,type_cop=type_cop,
                             n=n,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,namesaniso=namesaniso,weigthed=weigthed,X=X, local=local,GPU=GPU,MM=MM,aniso=aniso)

          if(hessian) 
          {
               if(!bivariate)  
                CompLikelihood$hessian=numDeriv::hessian(func=comploglik2,x=as.numeric(param),method="Richardson",  colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                              data1=data1,data2=data2,fixed=fixed,fan=fname,n=n,type_cop=type_cop,
                              namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,namesaniso=namesaniso,
                              weigthed=weigthed,X=X,local=local,GPU=GPU,MM=MM,aniso=aniso)
               if(bivariate)  
               CompLikelihood$hessian=numDeriv::hessian(func=comploglik_biv2,x=as.numeric(param),method="Richardson",colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                             data1=data1,data2=data2, fixed=fixed,fan=fname,n=n,type_cop=type_cop,
                             namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
                             weigthed=weigthed,X=X,local=local,GPU=GPU,MM=MM,aniso=aniso)
               rownames(CompLikelihood$hessian)=namesparam
               colnames(CompLikelihood$hessian)=namesparam
          }
  }

#####################################
#if((sensitivity||varest)&&(is.null(CompLikelihood$hessian)||min(eigen(CompLikelihood$hessian)$values)<0))
if((sensitivity||varest))
  {

 
if(!bivariate)  { if(is.null(CompLikelihood$hessian)) {CompLikelihood$hessian=matrix(c(1,2,3,2),2,2) }
                  if( min(eigen( CompLikelihood$hessian )$values) <0)
                                  {
                                 CompLikelihood$hessian=numDeriv::hessian(func=comploglik2,x=as.numeric(CompLikelihood$par),method="Richardson",   colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                                                        data1=data1,data2=data2, fixed=fixed,fan=fname,n=n,
                                                        namescorr=namescorr, namesnuis=namesnuis, namesparam=namesparam,namesaniso=namesaniso,
                                                        weigthed=weigthed,X=X,local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)
                                   }
                 }

if(bivariate)   CompLikelihood$hessian=numDeriv::hessian(func=comploglik_biv2,x=as.numeric(CompLikelihood$par),method="Richardson", colidx=colidx,rowidx=rowidx,corrmodel=corrmodel,  coords=coords,
                          data1=data1,data2=data2, fixed=fixed,fan=fname,n=n,
                          namescorr=namescorr, namesnuis=namesnuis,namesparam=namesparam, namesaniso=namesaniso,
                          weigthed=weigthed,X=X,local=local,GPU=GPU,MM=MM,aniso=aniso,type_cop=type_cop,cond_pair=cond_pair)

rownames(CompLikelihood$hessian)=namesparam
colnames(CompLikelihood$hessian)=namesparam
  }


####################################
       if( (CompLikelihood$convergence!='Successful')||CompLikelihood$value==-1e+15)  print("Optimization may have failed: try with other starting values ")
          else{
    if(varest){
           stop("subsambling is not working")
            }
      }
      if(hessian) CompLikelihood$sensmat=CompLikelihood$hessian
    if(!is.null(GPU)) gc()
    return(CompLikelihood)
  }

