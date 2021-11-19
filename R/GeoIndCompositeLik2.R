####################################################
### File name: CompIndLik2.r
####################################################

### Optim call for Indipendence Composite log-likelihood maximization

CompIndLik2 <- function(bivariate, coordx, coordy ,coordt,coordx_dyn, data, flagcorr, flagnuis, fixed,grid,
                           lower, model, n, namescorr, namesnuis, namesparam,
                           numparam,  optimizer, onlyvar, parallel, param, spacetime, type,
                           upper,namesupper, varest, ns, X,sensitivity,copula,MM)
  {


##################################
one_log_tukeyh=function(data,mm,sill,tail) {
          q = (data - mm)/sqrt(sill);
          aa=tail*q*q
          x = sign(q)*sqrt(lamW::lambertW0(aa)/tail);
          extra = 1/((1+lamW::lambertW0(aa)));
          return(log(dnorm(x)* x  * extra/(q*sqrt(sill))))
 }
#################################
indloglik<- function(fan,data,mm,nuis){


## gaussian and misspecified gaussian
    if(fan=="Ind_Pair_Gauss")                  { sill=nuis[1];
                                                 res=sum(dnorm(data,mean =mm,sd=sqrt(sill), log = TRUE))
                                               }

    if(fan=="Ind_Pair_Gauss_misp_T")           {    df=1/nuis[1];sill=nuis[2]
                                                    vv=sill*df/(df-2)
                                                  res=sum(dnorm(data,mean =mm,sd=sqrt(vv), log = TRUE))
                                               }

    if(fan=="Ind_Pair_Gauss_misp_SkewT")        {    df=1/nuis[1];sill=nuis[2];skew=nuis[3]; 
                                                     kk=gamma(0.5*(df-1))/gamma(0.5*(df)) 
                                                     me=  sqrt(sill)*sqrt(df/pi)*kk*skew  ;    
                                                     vv=  sill*(df/(df-2) - me^2)
                                                    res=sum(dnorm(data,mean=(mm+me),sd=sqrt(vv), log = TRUE)) 
                                                }

    if(fan=="Ind_Pair_Gauss_misp_Tukeygh")      {   sill =nuis[1];eta  = nuis[2];tail = nuis[3]; 
                                                    eta2=eta*eta; u=1-tail;
                                                    me=sqrt(sill)*(exp(eta2/(2*u))-1)/(eta*sqrt(u));
                                                    vv=sill*((exp(2*eta2/(1-2*tail))-2*exp(eta2/(2*(1-2*tail)))+1)/(eta2*sqrt(1-2*tail)))-me*me;
                                                  res=sum(dnorm(data, mean=mm+me,sd=sqrt(vv),  log = TRUE))
                                                }

## non gaussian over R
    if(fan=="Ind_Pair_T")                     {  sill=nuis[2]; df=1/nuis[1]
                                                res=sum(dt((data-mm)/sqrt(sill), df, log = TRUE)-0.5*log(sill))
                                              }  
    if(fan=="Ind_Pair_Logistic")              {  sill=nuis[1]; 
                                                res=sum(dlogis(data,mm, sqrt(sill) ,log = TRUE))
                                              }  
    if(fan=="Ind_Pair_SkewGauss")   {
                                                sk=nuis[2];sill=nuis[1]
                                                omega=sk*sk + sill
                                                #alpha=sk/(sill^0.5)
                                               #res=sum( sn::dsn((data-mm)/sqrt(omega),xi=0,omega= 1,alpha= alpha,log=TRUE)-0.5*log(omega)) 
                                                q=data-mm
                                               res=sum(log(2)-0.5*log(omega)+dnorm(q/(sqrt(omega)),log=TRUE)+pnorm( (sk*q)/(sqrt(sill)*sqrt(omega)),log.p=TRUE))
                                   }


   if(fan=="Ind_Pair_SinhGauss")   {           
                                                skew=nuis[2];sill=nuis[1];tail=nuis[3]
                                                q=(data-mm)/(sqrt(sill));
                                                b1=tail*asinh(q)-skew;Z1=sinh(b1);
                                                res=sum(-0.5*log(q^2+1)-0.5*log(2*pi*sill)+log(cosh(b1))+log(tail)-Z1*Z1/2);
                                       }

  if(fan=="Ind_Pair_Tukeyh")        { 
                                           sill=nuis[1];tail=nuis[2]
                                           res=sum(one_log_tukeyh(data,mm,sill,tail)) 
                                    }
   if(fan=="Ind_Pair_Tukeyhh")        { sill=nuis[1];tail1=nuis[3];tail2=nuis[2];
    res=sum(  log(  exp(one_log_tukeyh(data,mm,sill,tail2))*I(data>=mm) +
                    exp(one_log_tukeyh(data,mm,sill,tail1))*I(data<mm) )) 
                                      }
   if(fan=="Ind_Pair_TWOPIECEGauss") { sill=nuis[1];eta=nuis[2]
                                        y=(data-mm)/sqrt(sill)
                    res=sum( log( dnorm(y/(1-eta))*I(y>=0) + dnorm(y/(1+eta))*I(y<0)) -0.5*log(sill))
                                      }  
   if(fan=="Ind_Pair_TWOPIECETukeyh") { sill=nuis[1];eta=nuis[2];tail=nuis[3]
                                        y=(data-mm)/sqrt(sill);
                    res= sum(  log( exp(one_log_tukeyh(y/(1-eta),0,1,tail))*I(y>=0) + 
                                    exp(one_log_tukeyh(y/(1+eta),0,1,tail))*I(y<0) )     -0.5*log(sill) )
                                       }
   if(fan=="Ind_Pair_TWOPIECET") { sill=nuis[2];eta=nuis[3];tail=1/nuis[1]
                                  y=(data-mm)/sqrt(sill);
                    res=sum( log(  dt(y/(1-eta),df=tail)*I(y>=0) + dt(y/(1+eta),df=tail)*I(y<0)) -0.5*log(sill))
                                       }
## non gaussian over R^+

if(fan== "Ind_Pair_Gamma")                 {
                                              shape=nuis[2]
                                              res=sum(dgamma(data, shape=shape/2, scale = 1/(shape/(2*exp(mm))), log = TRUE))
                                            }
if(fan== "Ind_Pair_Weibull")              {

                                              shape=nuis[2] 
                                              res=sum(dweibull(data, shape=shape,scale=exp(mm)/(gamma(1+1/shape)) , log = TRUE))
                                            }
 if(fan== "Ind_Pair_LogGauss")              {
                                              sill=nuis[1]
                                              c1=mm-sill/2
                                              res=sum(dlnorm(data, meanlog =c1, sdlog = sqrt(sill), log = TRUE))
                                            }
if(fan== "Ind_Pair_LogLogistic")            {  shape=nuis[2]
                                              ci=gamma(1+1/shape)*gamma(1-1/shape)
                                              res=sum(actuar::dllogis(data, shape, scale = exp(mm)/ci, log = TRUE))
                                            }
## non Gassian bounded support
  if(fan== "Ind_Pair_Beta2")                {    mmax=nuis[4];mmin=nuis[3]
                                                 shape=nuis[2]
                                    
                                                 me=1/(1+exp(-mm))
                                                 res=sum(dbeta((data-mmin)/(mmax-mmin), me*shape, (1-me)*shape,log=TRUE)-log(mmax-mmin))
                                            }
 if(fan== "Ind_Pair_Kumaraswamy2")                {  
                                             mmax=nuis[4];mmin=nuis[3]
                                             shape=nuis[2]
                                             q=(data-mmin)/(mmax-mmin);k=1-q^shape
                                             m1=1/(1+exp(-mm));
                                             shapei=log(0.5)/log1p(-m1^shape);
                                             res=sum(log(shapei)+log(shape)+(shape-1)*log(q)+(shapei-1)*log(k)-log(mmax-mmin))            
                                             }

#### discrete
    if(fan=="Ind_Pair_Pois")              {mm=exp(mm);res=sum(dpois(data, mm, log = TRUE)) }
    if(fan=="Ind_Pair_Gauss_misp_Pois")   {mm=exp(mm);res=sum(dnorm(data, mean = mm, sd =sqrt(mm), log = TRUE))}
    if(fan=="Ind_Pair_BinomGauss")        res=sum(dbinom(data, n, pnorm(mm), log = TRUE))
    if(fan=="Ind_Pair_BinomGauss_misp")   { pp=pnorm(mm); mm=n*pp;vv=mm*(1-pp)
                                           res=sum(dnorm(data, mean =mm , sd =sqrt(vv), log = TRUE))}
    if(fan=="Ind_Pair_BinomnegGauss")     res=sum(dnbinom(data, n, pnorm(mm), log = TRUE))
    if(fan=="Ind_Pair_PoisGamma")         {mm=exp(mm);res=sum(dnbinom(data, nuis[2], mu=mm, log = TRUE))}
    if(fan=="Ind_Pair_Gauss_misp_PoisGamma") {mm=exp(mm);res=sum(dnorm(data, mean = mm, sd =sqrt(mm*(1+mm/nuis[2])), log = TRUE))}
    if(fan=="Ind_Pair_PoisZIP")           {mm=exp(mm);pp=pnorm(nuis[3])
                                            res1=sum(log(pp+(1-pp)*dpois(data[data==0],mm[data==0],0)));
                                            res2=sum(log(1-pp)+dpois(data[data!=0],mm[data!=0],1)); res=res1+res2}
    if(fan=="Ind_Pair_Gauss_misp_PoisZIP"){mm=exp(mm);pp=pnorm(nuis[3])
                                            res1=sum(log(pp+(1-pp)*dnorm(data[data==0],mean=mm[data==0],sd =sqrt(mm[data==0]), log = FALSE)));
                                            res2=sum(log(1-pp)+dnorm(data[data!=0],mean=mm[data!=0], sd =sqrt(mm[data!=0]), log = TRUE)); res=res1+res2}
    if(fan=="Ind_Pair_BinomnegGaussZINB") {pm=pnorm(mm);pp=pnorm(nuis[3])
                                            res1=sum(log(pp+(1-pp)*dnbinom(data[data==0],n,pm[data==0],log = FALSE)));
                                            res2=sum(log(1-pp)+dnbinom(data[data!=0],n, pm[data!=0], log = TRUE)); res=res1+res2}
### .......
return(-res)
}

 compindloglik2 <- function(param,  data,fixed, fan, n, 
                              namesnuis,namesparam,X,MM)
      {


        names(param) <- namesparam
        param <- c(param, fixed)
        nuisance <- param[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
        mm=as.numeric(nuisance[sel])   ## mean paramteres
   
        other_nuis=as.numeric(nuisance[!sel])   ## or nuis parameters (nugget sill skew df)
        if((is.null(MM))) Mean=c(X%*%mm)
        else Mean=c(MM)
        Mean=c(X%*%mm)
    
        result=indloglik(fan,data,Mean,other_nuis)
        return(result)
      }

 compindloglik_biv2 <- function(param, data1,data2,fixed, fan, n,   ## to do...
                           namesnuis,namesparam,X,MM)
      {

        names(param) <- namesparam
        param <- c(param, fixed)
        nuisance <- param[namesnuis]
        sel1=substr(names(nuisance),1,6)=="mean_1"
        mm1=as.numeric(nuisance[sel1])
        sel2=substr(names(nuisance),1,6)=="mean_2"
        mm2=as.numeric(nuisance[sel2])
        sel=substr(names(nuisance),1,4)=="mean"
        X1=as.matrix(X[1:ns[1],]);X2=as.matrix(X[(ns[1]+1):(ns[2]+ns[1]),]); 
        other_nuis=as.numeric(nuisance[!sel]) 
        Mean=c(X1%*%mm1,X2%*%mm2)
        res=double(1)
     
        result=4
        return(-result)
      }

   ##################################################################################################
   ############### starting function ###############################################################
   ##################################################################################################


    numcoord=length(coordx);numtime=1;spacetime_dyn=FALSE
    if(spacetime) numtime=length(coordt)
    if(bivariate) numtime=2
    if(!is.null(coordx_dyn)) spacetime_dyn=TRUE

    dimat=numcoord*numtime#
    if((spacetime||bivariate)) dimat <- sum(ns)
    NS=cumsum(ns)
    if(is.null(dim(X))){X=as.matrix(rep(1,dimat))}
    #else(if(bivariate) X=rbind(X,X))


    fname <- NULL; hessian <- FALSE
     


###################### pairwise ###############################################
    if( model==1 ) fname <- 'Ind_Pair_Gauss'                #
    if( model==2 ) fname <- 'Ind_Pair_BinomGauss'           #
    if( model==14 ) fname <- 'Ind_Pair_BinomnegGauss'       #                                 
    if( model==16 ) fname <- 'Ind_Pair_BinomnegGauss'       #                                   
    #if( model==15 ) fname <- 'Ind_Pair_PoisbinGauss'   
    #if( model==17 ) fname <- 'Ind_Pair_PoisbinnegGauss'
    if( model==13) fname <- 'Ind_Pair_WrapGauss'            #
    if( model==10) fname <- 'Ind_Pair_SkewGauss'            #
    if( model==21 ) fname <- 'Ind_Pair_Gamma'                #
   # if( model==33 ) fname <- 'Ind_Pair_Kumaraswamy'                                           
    if( model==42 ) fname <- 'Ind_Pair_Kumaraswamy2'         #                            
    #if( model==28 ) fname <- 'Ind_Pair_Beta'                                     
    if( model==50 ) fname <- 'Ind_Pair_Beta2'                #      
    if( model==26 ) fname <- 'Ind_Pair_Weibull'              #                                                                  
    if( model==24 ) fname <- 'Ind_Pair_LogLogistic'          #
    if( model==25 ) fname <- 'Ind_Pair_Logistic'             #                                                                                                                                                                                
    if( model==22 ) fname <- 'Ind_Pair_LogGauss';            #                                  
    if( model==27 ) fname <- 'Ind_Pair_TWOPIECET'            #                             
    #if( model==39 ) fname <- 'Ind_Pair_TWOPIECEBIMODAL'
    if( model==29 ) fname <- 'Ind_Pair_TWOPIECEGauss'                #
    if( model==12 ) fname <- 'Ind_Pair_T'                            #                            
    if( model==34 ) fname <- 'Ind_Pair_Tukeyh'                       #                                     
    if( model==40 ) fname <- 'Ind_Pair_Tukeyhh'                      #                      
    if( model==41 ) fname <- 'Ind_Pair_Gauss_misp_Tukeygh'           #                                 
    if( model==36 ) fname <- 'Ind_Pair_Gauss_misp_Pois'              #                                  
    if( model==35 ) fname <- 'Ind_Pair_Gauss_misp_T'                 #                                  
    if( model==37 ) fname <- 'Ind_Pair_Gauss_misp_SkewT'             #                                 
    if( model==20 ) fname <- 'Ind_Pair_SinhGauss'                    #                                  
    if( model==38 ) fname <- 'Ind_Pair_TWOPIECETukeyh'               #
    if( model==30 ) fname <- 'Ind_Pair_Pois'                         #
    if( model==46 ) fname <- 'Ind_Pair_PoisGamma'
    if( model==43 ) fname <- 'Ind_Pair_PoisZIP'
    if( model==44 ) fname <- 'Ind_Pair_Gauss_misp_PoisZIP'
    if( model==45 ) fname <- 'Ind_Pair_BinomnegGaussZINB'
    if( model==47 ) fname <- 'Ind_Pair_Gauss_misp_PoisGamma'
    if( model==11 ) fname <- 'Ind_Pair_BinomGauss'                   #                    
    if( model==51 ) fname <- 'Ind_Pair_BinomGauss_misp'              #
    #if( model==49) fname <- 'Ind_Pair_BinomLogi'
   
########################################################################                                            
    if(sensitivity) hessian=TRUE
    
     if((spacetime||bivariate)&&(!spacetime_dyn))    data=c(t(data))
     if((spacetime||bivariate)&&(spacetime_dyn))     data=unlist(data)          
     if(spacetime||bivariate)   NS=c(0,NS)[-(length(ns)+1)]


####     

tot=c(param,fixed) ## all the parameters
## deleting corr para
tot=tot[is.na(pmatch(names(tot),namescorr))]
## deleting nugget
tot=tot[names(tot)!='nugget'];namesnuis=namesnuis[namesnuis!="nugget"]
param=tot[pmatch(namesparam,names(tot))];param=param[!is.na(names(param))];param=param[!is.na(param)];


if(model  %in%  c(2,14,16,21,42,50,26,24,25,30,46,43,11))  ## model where sill must be fixed = to 1
 {param=param[names(param)!='sill'];a=1; names(a)="sill";
  if(is.null(unlist(fixed['sill']))) fixed=c(fixed,a)}


if(!is.null(copula))
    {if(copula=="Clayton")
           {param=param[names(param)!='nu'];a=2; names(a)="nu";
           if(is.null(unlist(fixed['nu']))) fixed=c(fixed,a)}}

namesparam=names(param)
  

###updating upper and lower bound if necessary
sel=pmatch(namesparam,namesupper)
lower=lower[sel]
upper=upper[sel]

 ###
   if(!onlyvar){
  ##############################.  spatial or space time ############################################
   if(!bivariate)           {
    if(length(param)==1) {
      
         optimizer="optimize"  
         if(is.na(lower)||is.na(upper))  {
            if(model %in% c(2,14,16,45,11)) {lower=-5;upper=5}
            else                            {lower=-1e+10;upper=1e+10}
           }

     CompLikelihood <- optimize(f= compindloglik2,    
                              data=data, fixed=fixed, fan=fname,  lower=lower, n=n,
                               namesnuis=namesnuis,namesparam=namesparam, maximum = FALSE,
                              upper= upper,  X=X,MM=MM)}
   if(length(param)>1) {
    if(optimizer=='L-BFGS-B'&&!parallel)
      CompLikelihood <- optim(par=param,fn= compindloglik2, 
                              control=list(factr=1e-10,pgtol=1e-14, maxit=100000), 
                                data=data, fixed=fixed,
                              fan=fname, lower=lower, method='L-BFGS-B',n=n,
                               namesnuis=namesnuis,namesparam=namesparam, 
                              upper=upper,  X=X,MM=MM,   hessian=FALSE)
      if(optimizer=='L-BFGS-B'&&parallel){
        ncores=max(1, parallel::detectCores() - 1)
        if(Sys.info()[['sysname']]=="Windows") cl <- parallel::makeCluster(ncores,type = "PSOCK")
        else                                   cl <- parallel::makeCluster(ncores,type = "FORK")
        parallel::setDefaultCluster(cl = cl)
        CompLikelihood <- optimParallel::optimParallel(par=param,fn= compindloglik2, 
                              control=list(pgtol=1e-14, maxit=100000,factr=1e-10),
                                 
                               data=data,fixed=fixed,fan=fname, lower=lower, method='L-BFGS-B',n=n,
                               namesnuis=namesnuis,namesparam=namesparam, 
                              upper=upper,  X=X,MM=MM,  hessian=FALSE)
         parallel::setDefaultCluster(cl=NULL)
         parallel::stopCluster(cl)
         }
    if(optimizer=='BFGS') 
        CompLikelihood <- optim(par=param, fn= compindloglik2,     
                           control=list(factr=1e-10,
                             reltol=1e-14, maxit=100000),data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='BFGS',n=n,
                              namesnuis=namesnuis,namesparam=namesparam,  X=X,MM=MM)
   if(optimizer=='SANN'){ 
   print(param)
    CompLikelihood <- optim(par=param, fn= compindloglik2,     
                           control=list(factr=1e-10,
                             reltol=1e-14, maxit=100000),data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='SANN',n=n,
                              namesnuis=namesnuis,namesparam=namesparam,  X=X,MM=MM)
  }
      if(optimizer=='Nelder-Mead')
        CompLikelihood <- optim(par=param, fn= compindloglik2,     
          control=list( reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='Nelder-Mead',n=n,
                                  namesnuis=namesnuis,namesparam=namesparam,  X=X,MM=MM)
 if(optimizer=='multinlminb'){

       CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn= compindloglik2,
           data=data, fixed=fixed,
                               fan=fname,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                                 X=X,MM=MM,  
          lower=lower,upper=upper,method = "nlminb", nbtrials = 400, typerunif = "sobol",
                              control = list( iter.max=100000),
                           )#,nbclusters=4,
                     
      
  }

 if(optimizer=='multiNelder-Mead'){
       CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn= compindloglik2,
           data=data, fixed=fixed,
                               fan=fname,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                                 X=X,MM=MM,  lower=lower,upper=upper,
          method = "Nelder-Mead", nbtrials = 500, 
                              control=list( reltol=1e-14, maxit=100000),
                           typerunif = "sobol"#,nbclusters=4,
                     )
    
  }


    if(optimizer=='nmk')
      CompLikelihood <-dfoptim::nmk(par=param, fn= compindloglik2, control = list(maxfeval=100000,tol=1e-10),
                            data=data,fixed=fixed, fan=fname,
                           n=n,namesnuis=namesnuis,namesparam=namesparam,  X=X,MM=MM)
    if(optimizer=='nmkb')
    {
      CompLikelihood <-dfoptim::nmkb(par=param, fn= compindloglik2, control = list(maxfeval=100000,tol=1e-10),
                         lower=lower,upper=upper,
                            data=data, fixed=fixed, fan=fname,
                         n=n,namesnuis=namesnuis,namesparam=namesparam,  X=X,MM=MM)
    }
    if(optimizer=='nlm')
    CompLikelihood <- nlm(f= compindloglik2,p=param,steptol = 1e-4,    data=data, fixed=fixed,
                               fan=fname,hessian=FALSE,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                               iterlim=100000,   X=X)
  
    if(optimizer=='nlminb')
     CompLikelihood <-nlminb(objective= compindloglik2,start=param,   data=data, fixed=fixed,
                                control = list( iter.max=100000),
                              lower=lower,upper=upper,
                               fan=fname,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                                 X=X,MM=MM)
    if(optimizer=='ucminf')   
      CompLikelihood <-ucminf::ucminf(par=param, fn= compindloglik2, hessian=as.numeric(hessian),   
                        control=list( maxeval=100000),
                               data=data,fixed=fixed, fan=fname,
                            n=n,namesnuis=namesnuis,namesparam=namesparam,  X=X,MM=MM)
                               
    }}
######################################################################################
############################## bivariate  ############################################ 
######################################################################################                          
    if(bivariate)           { 
     if(length(param)==1)
        {
         optimizer="optimize" 

       CompLikelihood <- optimize(f= compindloglik_biv2,     
                              data=data, fixed=fixed,fan=fname,  lower=lower,n=n,
                               namesnuis=namesnuis, namesparam=namesparam,maximum = FALSE,
                              upper=upper,  X=X,MM=MM )}
      if(length(param)>1) {   
    if(optimizer=='L-BFGS-B'&&!parallel){
      CompLikelihood <- optim(param, compindloglik_biv2, control=list(pgtol=1e-14, maxit=100000),
                              method='L-BFGS-B',hessian=FALSE,lower=lower, upper=upper,
                                  
                              data=data, fixed=fixed,fan=fname,n=n,
                               namesnuis=namesnuis, namesparam=namesparam,
                               X=X ,MM=MM )}
     #  if(optimizer=='lbfgsb3')
     # CompLikelihood <- lbfgsb3c::lbfgsb3c(param, compindloglik_biv, 
     #                         #control=list(pgtol=1e-14, maxit=100000),
     #                         lower=lower, upper=upper,
     #                             
     #                         data=data, fixed=fixed,fan=fname,n=n,
     #                          namesnuis=namesnuis, namesparam=namesparam,
     #                          X=X  )
    if(optimizer=='L-BFGS-B'&&parallel){
           ncores=max(1, parallel::detectCores() - 1)
           if(Sys.info()[['sysname']]=="Windows") cl <- parallel::makeCluster(ncores,type = "PSOCK")
           else                                   cl <- parallel::makeCluster(ncores,type = "FORK")
           parallel::setDefaultCluster(cl = cl)
           CompLikelihood <- optimParallel::optimParallel(param, compindloglik_biv2, 
                              control=list(pgtol=1e-14, maxit=100000),
                                  
                              data=data, fixed=fixed,
                              fan=fname,  n=n,  namesnuis=namesnuis, namesparam=namesparam,
                                X=X,MM=MM, 
                              lower=lower,upper=upper,
                              hessian=FALSE)
              parallel::setDefaultCluster(cl=NULL)
              parallel::stopCluster(cl)
         }
      if(optimizer=='BFGS')
      CompLikelihood <- optim(param, compindloglik_biv2,     control=list(
                              reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='BFGS',n=n,
                               namesnuis=namesnuis,namesparam=namesparam ,  X=X,MM=MM )
   if(optimizer=='Nelder-Mead')
      CompLikelihood <- optim(param, compindloglik_biv2,     control=list(
                              reltol=1e-14, maxit=100000), data=data, fixed=fixed, fan=fname,
                              hessian=FALSE, method='Nelder-Mead',n=n,
                               namesnuis=namesnuis,namesparam=namesparam ,  X=X,MM=MM )
    if(optimizer=='nmk')
        CompLikelihood <- dfoptim::nmk(par=param, fn= compindloglik_biv2, control = list(maxfeval=100000,tol=1e-10),
                                 data=data, fixed=fixed, fan=fname,
                            n=n,  namesnuis=namesnuis,namesparam=namesparam ,  X=X ,MM=MM)
    if(optimizer=='nmkb')
        CompLikelihood <- dfoptim::nmkb(par=param, fn= compindloglik_biv2, control = list(maxfeval=100000,tol=1e-10),
                             lower=lower,upper=upper,
                                 data=data, fixed=fixed, fan=fname,
                            n=n,  namesnuis=namesnuis,namesparam=namesparam ,  X=X ,MM=MM)
     if(optimizer=='nlm') 
        CompLikelihood <- nlm( f= compindloglik_biv2,p=param,     data=data, fixed=fixed,
                               fan=fname,hessian=FALSE,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                                 X=X,MM=MM )
    if(optimizer=='ucminf') 
         CompLikelihood <-ucminf::ucminf(par=param, fn= compindloglik_biv2, hessian=as.numeric(hessian),   
                        control=list( maxeval=100000), 
                            data=data, fixed=fixed,
                        fan=fname,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                          X=X )
    if(optimizer=='nlminb') 
        CompLikelihood <- nlminb( objective= compindloglik_biv2,start=param, 
                                     control = list( iter.max=100000),
                              lower=lower,upper=upper,
                                   data=data, fixed=fixed,
                               fan=fname,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                                 X=X ,MM=MM)
     if(optimizer=='multinlminb'){
       CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn= compindloglik_biv2,
           data=data, fixed=fixed,
                               fan=fname,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                                 X=X,MM=MM,
                                    lower=lower,upper=upper,method = "nlminb", nbtrials = 500, 
                              control = list( iter.max=100000),
                           typerunif = "sobol")
                               }
     if(optimizer=='multiNelder-Mead'){
       CompLikelihood <- mcGlobaloptim::multiStartoptim(objectivefn= compindloglik_biv2,
           data=data, fixed=fixed,
                               fan=fname,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                                 X=X,MM=MM,  lower=lower,upper=upper,
          method = "Nelder-Mead", nbtrials = 500, 
                              control=list( reltol=1e-14, maxit=100000),
                           typerunif = "sobol"#,nbclusters=4,
                     )
  }

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

    if(optimizer=='nlminb'||optimizer=='multinlminb'){
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
          if(!bivariate) CompLikelihood$value = -  compindloglik2(param=CompLikelihood$par ,     
                              data=data, fixed=fixed, fan=fname,
                             n=n,namesnuis=namesnuis,namesparam=namesparam,  X=X,MM=MM)
          else CompLikelihood$value = - compindloglik_biv2(param=CompLikelihood$par ,     
                data=data, fixed=fixed, fan=fname,
                             n=n,namesnuis=namesnuis,namesparam=namesparam,  X=X,MM=MM)

          if(hessian) 
          {
               if(!bivariate)  
                CompLikelihood$hessian=numDeriv::hessian(func= compindloglik2,x=param,method="Richardson",     
                              data=data,fixed=fixed,fan=fname,n=n,
                               namesnuis=namesnuis, namesparam=namesparam,
                                X=X,MM=MM )
               if(bivariate)  
               CompLikelihood$hessian=numDeriv::hessian(func= compindloglik_biv2,x=param,method="Richardson",   
                             data=data, fixed=fixed,fan=fname,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                               X=X ,MM=MM)
               rownames(CompLikelihood$hessian)=namesparam
               colnames(CompLikelihood$hessian)=namesparam
          }
  }

#####################################
if((sensitivity||varest))
  {
if(!bivariate)  

CompLikelihood$hessian=numDeriv::hessian(func= compindloglik2,x=CompLikelihood$par,method="Richardson",      
                              data=data, fixed=fixed,fan=fname,n=n,
                               namesnuis=namesnuis, namesparam=namesparam,
                                X=X ,MM=MM)
if(bivariate)  
CompLikelihood$hessian=numDeriv::hessian(func= compindloglik_biv2,x=CompLikelihood$par,method="Richardson",    
                          data=data, fixed=fixed,fan=fname,n=n, namesnuis=namesnuis,namesparam=namesparam, 
                                 X=X ,MM=MM)
rownames(CompLikelihood$hessian)=namesparam
colnames(CompLikelihood$hessian)=namesparam
  }


if(hessian) CompLikelihood$sensmat=CompLikelihood$hessian

    return(CompLikelihood)
  }

