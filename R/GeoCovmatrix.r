####################################################
### Authors:  Moreno Bevilacqua, Víctor Morales Oñate.
### Email: moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: GeoCovmatrix.r
### Description:
### This file contains a set of procedures
### for computing a covariance (tapered) matrix for a given
### space(time) covariance model.
### Last change: 28/03/2019.
####################################################



### decomposition of a square  matrix
MatDecomp<-function(mtx,method)    {
        if(method=="cholesky")  {
            mat.decomp <- try(chol(mtx), silent=TRUE)
            if (inherits(mat.decomp , "try-error")) return (FALSE)
        }
        if(method=="svd")      {
            mat.decomp <- svd(mtx)
            cov.logdeth <- try(sum(log(sqrt(mat.decomp$d))), silent=TRUE)
            if (inherits(cov.logdeth, "try-error"))  return (FALSE)
        }
        return(mat.decomp)
    }
### square root of a square matrix
MatSqrt<-function(mat.decomp,method)    {  
        if(method=="cholesky")  varcov.sqrt <- mat.decomp
        if(method=="svd")       varcov.sqrt <- sqrt(diag(mat.decomp$d))%*%t(mat.decomp$u)
        return(varcov.sqrt)
    }  
### inverse a square matrix given a decomposition
MatInv<-function(mat.decomp,method)    {

        if(method=="cholesky")  varcov.inv <- chol2inv(mat.decomp)
        if(method=="svd")       
        { 
              tol = sqrt(.Machine$double.eps)
              e <- mat.decomp$d;e[e > tol] <- 1/e[e > tol] 
              varcov.inv<-mat.decomp$v %*% diag(e,nrow=length(e)) %*% t(mat.decomp$u) 
        } 
        return(varcov.inv)
    } 
### determinant a square matrix given a decomposition    
MatLogDet<-function(mat.decomp,method)    {
        if(method=="cholesky")  det.mat <- 2*sum(log(diag(mat.decomp)))
        if(method=="svd")       det.mat <- sum(log(mat.decomp$d))
        return(det.mat)
    }      


######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

GeoCovmatrix <- function(coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL,corrmodel, distance="Eucl", grid=FALSE,
                       maxdist=NULL, maxtime=NULL, model="Gaussian", n=1, param, radius=6371, 
                       sparse=FALSE,taper=NULL, tapsep=NULL, type="Standard",X=NULL)

{
  ########################################################################################################
  ##########  Internal function: computing covariance matrix #############################################
  ########################################################################################################

    Cmatrix <- function(bivariate, coordx, coordy, coordt,corrmodel, dime, n, ns, NS, nuisance, numpairs,
                           numpairstot, model, paramcorr, setup, radius, spacetime, spacetime_dyn,type,X)
    {
   
#####################################################f
    if(model %in% c(1,20))   ## gaussian case  or sinhgaussian case (sinh no implementerd)
    {
      
        if(type=="Standard")  {

        fname <-"CorrelationMat2"
      if(spacetime) fname <- "CorrelationMat_st_dyn2"
        if(bivariate) fname <- "CorrelationMat_biv_dyn2"
        cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),
          as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),as.double(radius), 
          as.integer(ns),as.integer(NS),
          PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)

     #   cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",
     #  "integer","double","double","double","integer","integer"),         
     #  corr=double(numpairstot),
     #  coordx,coordy,coordt,corrmodel, nuisance,paramcorr,radius,ns,NS,
     #        INTENT = c("w", "r", "r", "r", "r","r","r","r", "r", "r"),
     #        NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0)

        corr=cr$corr

        if(!bivariate)                  {
        # Builds the covariance matrix:
          varcov <- (nuisance['nugget'] + nuisance['sill']) * diag(dime)
          corr <- corr * nuisance['sill']
          varcov[lower.tri(varcov)] <- corr
          varcov <- t(varcov)
          varcov[lower.tri(varcov)] <- corr }
        if(bivariate)      {
          varcov<-diag(dime)
          varcov[lower.tri(varcov,diag=TRUE)] <- corr
          varcov <- t(varcov)
          varcov[lower.tri(varcov,diag=TRUE)] <- corr
        }}
       if(type=="Tapering")  {
        fname <- "CorrelationMat_tap"
        if(spacetime) fname <- "CorrelationMat_st_tap"
        #if(spacetime) fname <- "CorrelationMat_st_tap_dyn2"
       if(bivariate) fname <- "CorrelationMat_biv_tap"
        #if(bivariate) fname <- "CorrelationMat_biv_tap_dyn2"

        cr=.C(fname,  corr=double(numpairs), as.double(coordx),as.double(coordy),as.double(coordt),
          as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),as.double(radius),as.integer(ns),
           as.integer(NS),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)

    #      cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",
    #   "integer","double","double","double","integer","integer"),        
    #   #corr=vector_dc("numeric", numpairstot), 
    #   corr=double(numpairs),
    #   coordx,coordy,coordt,corrmodel, nuisance,paramcorr,radius,ns,NS,
    #         INTENT = c("w", "r", "r", "r", "r","r","r","r", "r", "r"),
    #         NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0) 


        corr=cr$corr
        if(!bivariate){
             vcov <- corr*nuisance['sill']; vcov[vcov==(nuisance['sill'])] <- nuisance['sill']+nuisance['nugget']
             varcov <- new("spam",entries=vcov*setup$taps,colindices=setup$ja,
                             rowpointers=setup$ia,dimension=as.integer(rep(dime,2)))}
          else  {varcov <-new("spam",entries=corr*setup$taps,colindices=setup$ja,
                         rowpointers=setup$ia,dimension=as.integer(rep(dime,2))) }              
         }
       }  


###############################################################
   if(model %in% c(2,11,19)){ #  binomial Gaussian type 1  and 2
            sel=substr(names(nuisance),1,4)=="mean"
            mm=as.numeric(nuisance[sel]) 
            mu = X%*%mm
            other_nuis=as.numeric(nuisance[!sel])   ## other nuis parameters (nugget sill skew df)
            fname <-"CorrelationMat_bin2"
            ##if(spacetime) fname <- "CorrelationMat_st_bin2"
            if(spacetime) fname <- "CorrelationMat_st_dyn_bin2"
            #print(other_nuis)
           # if(bivariate) fname <- "CorrelationMat_biv_bin_dyn2"
            if(bivariate) fname <- "CorrelationMat_biv_bin_dyn2"
            cr=.C(fname, corr=double(numpairstot+dime),  as.double(coordx),as.double(coordy),as.double(coordt),
              as.integer(corrmodel), as.double(c(mu)),as.integer(min(n)), as.double(other_nuis), as.double(paramcorr),as.double(radius),
              as.integer(ns), as.integer(NS),
              PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
            corr=cr$corr
            if(!bivariate)                  {
                # Builds the covariance matrix:
                varcov <-  diag(dime) 
                varcov[lower.tri(varcov,diag=TRUE)] <- corr
                varcov <- t(varcov)
                varcov[lower.tri(varcov,diag=TRUE)] <- corr 
                if(model %in% c(2,11)) {
                 pg=pnorm(mu)
                 diag(varcov)=pg*(1-pg)*n}
                }  

                 ## covariance matrix  min(n)=k for the type 2
            #   if(bivariate)      {
            #    varcov<-diag(dime)
            #    varcov[lower.tri(varcov,diag=T)] <- corr
            #    varcov <- t(varcov)
            #    varcov[lower.tri(varcov,diag=T)] <- corr
            #}
        }
###############################################################
         if(model==14)  {  ##  geometric Gaussian
          ##nuisance <- param[namesnuis]
            sel=substr(names(nuisance),1,4)=="mean"
            mm=as.numeric(nuisance[sel]) 
            mu = X%*%mm
            other_nuis=as.numeric(nuisance[!sel])   ## other nuis parameters (nugget sill skew df)
           fname <-"CorrelationMat_geom2"
            #print(other_nuis)
            #if(spacetime) fname <- "CorrelationMat_st_geom2"
            if(spacetime) fname <- "CorrelationMat_st_dyn_geom2"
           # if(bivariate) fname <- "CorrelationMat_biv_geom_dyn2"
            if(bivariate) fname <- "CorrelationMat_biv_geom_dyn2"
      
            cr=.C(fname, corr=double(numpairstot+dime),  as.double(coordx),as.double(coordy),as.double(coordt),
              as.integer(corrmodel), as.double(c(mu)), as.double(other_nuis), as.double(paramcorr),as.double(radius),
              as.integer(ns), as.integer(NS),
              PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
            corr=cr$corr
            if(!bivariate)                  {
                # Builds the covariance matrix:
                varcov <-  diag(dime) 
                varcov[lower.tri(varcov,diag=TRUE)] <- corr
                varcov <- t(varcov)
                varcov[lower.tri(varcov,diag=TRUE)] <- corr 
                pg=pnorm(mu)
                diag(varcov)=(1-pg)/pg^2
                } ## cov matrix
                ##pp=pnorm(X%*%mm) # pr of success
                ##varcov=varcov*((1-pp)/pp^2) ## covariance matrix
            #   if(bivariate)      {
            #    varcov<-diag(dime)
            #    varcov[lower.tri(varcov,diag=T)] <- corr
            #    varcov <- t(varcov)
            #    varcov[lower.tri(varcov,diag=T)] <- corr
            #}
        } 
      if(model==16)  {  ##  binomial negative
            sel=substr(names(nuisance),1,4)=="mean"
            mm=as.numeric(nuisance[sel]) 
            mu = X%*%mm
            other_nuis=as.numeric(nuisance[!sel])   ## other nuis parameters (nugget sill skew df)
            fname <-"CorrelationMat_binneg2"
            #if(spacetime) fname <- "CorrelationMat_st_binneg2"
            if(spacetime) fname <- "CorrelationMat_st_dyn_binneg2"
           # if(bivariate) fname <- "CorrelationMat_biv_binneg_dyn2"
           # if(bivariate&&spacetime_dyn) fname <- "CorrelationMat_biv_binneg_dyn2"
            cr=.C(fname, corr=double(numpairstot+dime),  as.double(coordx),as.double(coordy),as.double(coordt),
              as.integer(corrmodel), as.double(c(mu)),as.integer(min(n)), as.double(other_nuis), as.double(paramcorr),as.double(radius),
              as.integer(ns), as.integer(NS), PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
            corr=cr$corr
            if(!bivariate)                  {
                # Builds the covariance matrix:
                varcov <-  diag(dime) 
                varcov[lower.tri(varcov,diag=TRUE)] <- corr
                varcov <- t(varcov)
                varcov[lower.tri(varcov,diag=TRUE)] <- corr 
                pg=pnorm(mu)
                diag(varcov)=n*(1-pg)/pg^2
              }  
                 ## covariance matrix  min(n)=k for the type 2
            #   if(bivariate)      {
            #    varcov<-diag(dime)
            #    varcov[lower.tri(varcov,diag=T)] <- corr
            #    varcov <- t(varcov)
            #    varcov[lower.tri(varcov,diag=T)] <- corr
            #}
        }
###############################################################
            if(model==9)  {  ## TukeyGH
 h=nuisance['tail']
 g=nuisance['skew']
              fname <-"CorrelationMat2"
      if(spacetime) fname <- "CorrelationMat_st_dyn2"
        if(bivariate) fname <- "CorrelationMat_biv_dyn2"
         cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),
          as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),as.double(radius), 
          as.integer(ns),as.integer(NS),
          PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)

    #     cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",
    #   "integer","double","double","double","integer","integer"),        
    #   #corr=vector_dc("numeric", numpairstot), 
    #   corr=double(numpairstot),
    #   coordx,coordy,coordt,corrmodel, nuisance,paramcorr,radius,ns,NS,
    #         INTENT = c("w", "r", "r", "r", "r","r","r","r", "r", "r"),
    #         NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0)

            corr=cr$corr
            if(!bivariate)                

  if(!g&&!h){     #ok ## gaussian case
              vv=1 # variance
              varcov <-  diag(dime)
               varcov[lower.tri(varcov)] <- corr
              varcov <- t(varcov)
              varcov[lower.tri(varcov)] <- corr 
              varcov=varcov*vv*as.numeric(nuisance['sill']) } 

  if(g&&!h){  #  
              vv=( -exp(g^2)+exp(g^2*2))*g^(-2)
              varcov <-  diag(dime)
              corr <- (( -exp(g^2)+exp(g^2*(1+corr)))*g^(-2))/vv
              varcov[lower.tri(varcov)] <- corr
              varcov <- t(varcov)
              varcov[lower.tri(varcov)] <- corr 
              varcov=varcov*vv*as.numeric(nuisance['sill']) } ## cov matrix
   if(!g&&h){ ##
              vv=(1-2*h)^(-1.5) # variance
              varcov <-  diag(dime)
              corr <- (-corr/((1+h*(corr-1))*(-1+h+h*corr)*(1+h*(-2+h-h*corr^2))^0.5))/vv
              varcov[lower.tri(varcov)] <- corr
              varcov <- t(varcov)
              varcov[lower.tri(varcov)] <- corr 
              varcov=varcov*vv*as.numeric(nuisance['sill']) } ## cov matrix

  if(h&&g){ # ok
              varcov <-  diag(dime)
              vv=(exp(g^2*(2)/(1-h*2))-2*exp(1/((1-h)^2-h^2)*(g^2/2))+1)/(g^2*((1-h)^2-h^2)^(0.5))-((exp(g^2/(2*(1-h)))-1)/(g*(1-h)^0.5))^2
              corr <-((exp(g^2*(1+corr)/(1-h*(1+corr)))-2*exp((1-h*(1-corr^2))/((1-h)^2-h^2*corr^2)*(g^2/2))+1)/(g^2*((1-h)^2-corr^2*h^2)^(0.5))-((exp(g^2/(2*(1-h)))-1)/(g*(1-h)^0.5))^2) /vv
              varcov[lower.tri(varcov)] <- corr
              varcov <- t(varcov)
              varcov[lower.tri(varcov)] <- corr  
              varcov=varcov*vv*as.numeric(nuisance['sill']) 
            } ## 
                

            #   if(bivariate)      {
            #    varcov<-diag(dime)
            #    varcov[lower.tri(varcov,diag=T)] <- corr
            #    varcov <- t(varcov)
            #    varcov[lower.tri(varcov,diag=T)] <- corr
            #}
        
        } 

###############################################################
if(model==34)  {  ## TukeyH  
 h=nuisance['tail']
              fname <-"CorrelationMat2"
      if(spacetime) fname <- "CorrelationMat_st_dyn2"
        if(bivariate) fname <- "CorrelationMat_biv_dyn2"
           cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),
          as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),as.double(radius), 
          as.integer(ns),as.integer(NS),
          PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)

   #         cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",
   #    "integer","double","double","double","integer","integer"),        
   #    #corr=vector_dc("numeric", numpairstot), 
   #    corr=double(numpairstot),
   #    coordx,coordy,coordt,corrmodel, nuisance,paramcorr,radius,ns,NS,
   #          INTENT = c("w", "r", "r", "r", "r","r","r","r", "r", "r"),
   #          NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0)
   
            corr=cr$corr
            if(!bivariate)                

  if(!h){     #ok ## gaussian case
              vv=1 # variance
              varcov <-  diag(dime)
               varcov[lower.tri(varcov)] <- corr
              varcov <- t(varcov)
              varcov[lower.tri(varcov)] <- corr 
              varcov=varcov*vv*as.numeric(nuisance['sill']) } 
   if(h){ ##ok
              vv=(1-2*h)^(-1.5) # variance
              varcov <-  diag(dime)
              corr <- (-corr/((1+h*(corr-1))*(-1+h+h*corr)*(1+h*(-2+h-h*corr^2))^0.5))/vv
              varcov[lower.tri(varcov)] <- corr
              varcov <- t(varcov)
              varcov[lower.tri(varcov)] <- corr 
              varcov=varcov*vv*as.numeric(nuisance['sill']) 
        } 
                
            #   if(bivariate)      {
            #    varcov<-diag(dime)
            #    varcov[lower.tri(varcov,diag=T)] <- corr
            #    varcov <- t(varcov)
            #    varcov[lower.tri(varcov,diag=T)] <- corr
            #}
        
        } 
######################################################################        
if(model==12)   ##  student case 
    {
        fname <-"CorrelationMat2"
        #if(spacetime) fname <- "CorrelationMat_st2"
        if(spacetime) fname <- "CorrelationMat_st_dyn2"
        #if(bivariate) fname <- "CorrelationMat_biv2"
        if(bivariate) fname <- "CorrelationMat_biv_dyn2"
            cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),
          as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),as.double(radius),
          as.integer(ns), as.integer(NS),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)

    # cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",
    #   "integer","double","double","double","integer","integer"),        
    #   #corr=vector_dc("numeric", numpairstot), 
    #   corr=double(numpairstot),
    #   coordx,coordy,coordt,corrmodel, nuisance,paramcorr,radius,ns,NS,
    #         INTENT = c("w", "r", "r", "r", "r","r","r","r", "r", "r"),
    #         NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0) 

        nu=as.numeric(1/nuisance['df'])    
  if(nu<170) corr=((nu-2)*gamma((nu-1)/2)^2*Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,cr$corr^2))*cr$corr)/(2*gamma(nu/2)^2)
  else corr=exp(log(nu-2)+2*lgamma(0.5*(nu-1))-log(2)-2*lgamma(nu/2)+log(Re(hypergeo::hypergeo(0.5,0.5, nu/2,cr$corr^2)))+log(cr$corr))
  if(!bivariate) {
        # Builds the covariance matrix:
        varcov <-  diag(dime)
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr   
        vv=(nu)/(nu-2)
        varcov=varcov*vv*as.numeric(nuisance['sill'])
        }
        }
############################################################### 
  if(model==18)   ##  skew student case 
    {
        fname <-"CorrelationMat2"
        #if(spacetime) fname <- "CorrelationMat_st2"
        if(spacetime) fname <- "CorrelationMat_st_dyn2"
        #if(bivariate) fname <- "CorrelationMat_biv2"
        if(bivariate) fname <- "CorrelationMat_biv_dyn2"
         cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),
          as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),as.double(radius),
          as.integer(ns), as.integer(NS),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)  
        cc=cr$corr

        nu=as.numeric(1/nuisance['df']); sk=as.numeric(nuisance['skew'])

        sk2=sk^2; KK=2*sk2/pi; D1=(nu-1)/2;D2=nu/2;
        CC=(pi*(nu-2)*gamma(D1)^2) /(2*( pi*gamma(D2)^2 *(1+sk2) - sk2*(nu-2)*gamma(D1)^2) );
        corr2= (1/(-1+1/KK))*(  sqrt(1-cc^2) + cc*asinh(cc) - 1 )+(1-sk2)*cc/(1-KK);
        
        corr=CC*( Re(hypergeo::hypergeo(0.5,0.5 ,nu/2 ,cc^2)) * ((1+sk2*(1-2/pi))*corr2 + KK)-KK )


  if(!bivariate) {
        # Builds the covariance matrix:
        varcov <-  diag(dime)
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr   
        vv=(nu)/(nu-2)  -    (nu*sk2/pi)*(gamma(D1)/gamma(D2))^2
        varcov=varcov*vv*as.numeric(nuisance['sill'])
        }
        }
###############################################################           
 
if(model==27)   ##  two piece student case case
    {
        fname <-"CorrelationMat2"
        #if(spacetime) fname <- "CorrelationMat_st2"
        if(spacetime) fname <- "CorrelationMat_st_dyn2"
        #if(bivariate) fname <- "CorrelationMat_biv2"
        if(bivariate) fname <- "CorrelationMat_biv_dyn2"

             cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),
          as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),as.double(radius),
          as.integer(ns), as.integer(NS),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE) 

    # cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",
    #   "integer","double","double","double","integer","integer"),        
    #   #corr=vector_dc("numeric", numpairstot), 
    #   corr=double(numpairstot),
    #   coordx,coordy,coordt,corrmodel, nuisance,paramcorr,radius,ns,NS,
    #         INTENT = c("w", "r", "r", "r", "r","r","r","r", "r", "r"),
    #         NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0)    
          nu=as.numeric(1/nuisance['df']); sk=as.numeric(nuisance['skew'])
          corr2=cr$corr^2;sk2=sk^2
          a1=Re(hypergeo::hypergeo(0.5,0.5,nu/2,corr2))
          a2=cr$corr*asin(cr$corr) + (1-corr2)^(0.5)
          ll=qnorm((1-sk)/2)
          p11=pbivnorm::pbivnorm(ll,ll, rho = cr$corr, recycle = TRUE)
          a3=3*sk2 + 2*sk + 4*p11 - 1
          KK=( nu*(nu-2)*gamma((nu-1)/2)^2) / (nu*pi*gamma(nu/2)^2*(3*sk2+1)-4*sk2*nu*(nu-2)*gamma((nu-1)/2)^2 )
          corr= KK*(a1*a2*a3-4*sk2);
         # print(corr)
  if(!bivariate) {
        # Builds the covariance matrix:
        varcov <-  diag(dime)
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr   
        vv=nu*(3*sk^2+1)/(nu-2)-(4*sk^2*nu*gamma((nu-1)/2)^2)/(pi*gamma(nu/2)^2)
        varcov=varcov*vv*nuisance['sill']
        }
    ## todo
    #if(bivariate)      {
     #     varcov<-diag(dime)
      #    varcov[lower.tri(varcov,diag=T)] <- corr
       #   varcov <- t(varcov)
        #  varcov[lower.tri(varcov,diag=T)] <- corr
        #}
      ###  
}
###############################################################  
if(model==29)   ##  two piece gaussian case
    {
        fname <-"CorrelationMat2"
        #if(spacetime) fname <- "CorrelationMat_st2"
        if(spacetime) fname <- "CorrelationMat_st_dyn2"
        #if(bivariate) fname <- "CorrelationMat_biv2"
        if(bivariate) fname <- "CorrelationMat_biv_dyn2"

 cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),
          as.integer(corrmodel), as.double(nuisance), as.double(paramcorr),as.double(radius),
          as.integer(ns), as.integer(NS),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)     
     

     #cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",
     #  "integer","double","double","double","integer","integer"),        
     #  #corr=vector_dc("numeric", numpairstot), 
     #  corr=double(numpairstot),
     #  coordx,coordy,coordt,corrmodel, nuisance,paramcorr,radius,ns,NS,
     #        INTENT = c("w", "r", "r", "r", "r","r","r","r", "r", "r"),
     #        NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0)    

          corr2=sqrt(1-cr$corr^2)
          xx=as.numeric(nuisance['skew']); xx2=xx^2
          ll=qnorm((1-xx)/2)
          p11=pbivnorm::pbivnorm(ll,ll, rho = cr$corr, recycle = TRUE)
          KK=3*xx2+2*xx+ 4*p11 - 1
          corr=(2*((corr2 + cr$corr*atan(cr$corr/corr2))*KK)- 8*xx2)/(3*pi*xx2  -  8*xx2   +pi   )
         # print(corr)
  if(!bivariate) {
        # Builds the covariance matrix:
        varcov <-  diag(dime)
        varcov[lower.tri(varcov)] <- corr
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr   
        vv=(1+3*xx2)-8*xx2/pi
        varcov=varcov*vv*nuisance['sill']
        }
    ## todo
    #if(bivariate)      {
     #     varcov<-diag(dime)
      #    varcov[lower.tri(varcov,diag=T)] <- corr
       #   varcov <- t(varcov)
        #  varcov[lower.tri(varcov,diag=T)] <- corr
        #}
      ###  
}


###############################################################           
if(model==21)   ##  gamma case
    {
         sel=substr(names(nuisance),1,4)=="mean"
         mm=as.numeric(nuisance[sel]) 
         mu = X%*%mm
         other_nuis=0# not necessary as.numeric(nuisance[!sel])   ## other nuis parameters (nugget sill skew df)
        fname <-"CorrelationMat2"
       # if(spacetime) fname <- "CorrelationMat_st2"
        if(spacetime) fname <- "CorrelationMat_st_dyn2"
       #if(bivariate) fname <- "CorrelationMat_biv2"
        if(bivariate) fname <- "CorrelationMat_biv_dyn2"

          cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),
          as.integer(corrmodel), as.double(other_nuis), as.double(paramcorr),as.double(radius),
          as.integer(ns), as.integer(NS),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)      
  
    # cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",
    #   "integer","double","double","double","integer","integer"),        
    #   #corr=vector_dc("numeric", numpairstot), 
    #   corr=double(numpairstot),
    #   coordx,coordy,coordt,corrmodel, nuisance,paramcorr,radius,ns,NS,
    #         INTENT = c("w", "r", "r", "r", "r","r","r","r", "r", "r"),
    #         NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0)     
  
  corr=cr$corr^2   ### gamma correlation
  if(!bivariate) {
        # Builds the covariance matrix:
        varcov <-  diag(1,dime)
        varcov[lower.tri(varcov)] <- corr*(1-nuisance['nugget'])
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <-corr*(1-nuisance['nugget'])  
        vv=exp(mu)^2 * 2/nuisance['shape']
        V=vv%*%t(vv)
        varcov=varcov*sqrt(V)
        }
    ## todo
    #if(bivariate)      {
     #     varcov<-diag(dime)
      #    varcov[lower.tri(varcov,diag=T)] <- corr
       #   varcov <- t(varcov)
        #  varcov[lower.tri(varcov,diag=T)] <- corr
        #}
      ###  
}
###############################################################
if(model==30)   ##  poisson case
    {

         sel=substr(names(nuisance),1,4)=="mean"
         mm=as.numeric(nuisance[sel]) 
         mu = X%*%mm
         other_nuis=0# not necessary as.numeric(nuisance[!sel])   ## other nuis parameters (nugget sill skew df)
        fname <-"CorrelationMat2"
       # if(spacetime) fname <- "CorrelationMat_st2"
        if(spacetime) fname <- "CorrelationMat_st_dyn2"
       #if(bivariate) fname <- "CorrelationMat_biv2"
        if(bivariate) fname <- "CorrelationMat_biv_dyn2"

         cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),
          as.integer(corrmodel), as.double(other_nuis), as.double(paramcorr),as.double(radius),
          as.integer(ns), as.integer(NS),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE) 

 # cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",
 #      "integer","double","double","double","integer","integer"),        
  #     #corr=vector_dc("numeric", numpairstot), 
  #     corr=double(numpairstot),
  #     coordx,coordy,coordt,corrmodel, nuisance,paramcorr,radius,ns,NS,
  #           INTENT = c("w", "r", "r", "r", "r","r","r","r", "r", "r"),
  #           NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0)      

corr=cr$corr*(1-nuisance['nugget'])         
corr2=corr^2   ### gamma correlation
corriinv=1/(1-corr2)
  if(!bivariate) {
        # Builds the covariance matrix:
        varcov <-  diag(1,dime);  varcov1 <-  diag(1,dime);
        varcov[lower.tri(varcov)] <- corr2
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <-corr2 
        varcov1[lower.tri(varcov1)] <- corriinv
        varcov1 <- t(varcov1)
        varcov1[lower.tri(varcov1)] <-corriinv
        vv=exp(mu);V=sqrt(vv%*%t(vv)); z=2*V*varcov1
        varcov=varcov*(1-exp(-z)*(besselI(z,0)+besselI(z,1))) ## correlation matrix 
        varcov=varcov*V     ## covarianvce matrix 
        diag(varcov)=diag(V) ## fixing diagonal
        }
    ## todo
    #if(bivariate)      {
     #     varcov<-diag(dime)
      #    varcov[lower.tri(varcov,diag=T)] <- corr
       #   varcov <- t(varcov)
        #  varcov[lower.tri(varcov,diag=T)] <- corr
        #}
      ###  
}
###############################################################           
if(model==26)   ##  weibull case
    {
         sel=substr(names(nuisance),1,4)=="mean"
         mm=as.numeric(nuisance[sel]) 
         mu = X%*%mm
         other_nuis=0# not necessary as.numeric(nuisance[!sel])   ## other nuis parameters (nugget sill skew df)
        fname <-"CorrelationMat2"
        #if(spacetime) fname <- "CorrelationMat_st2"
        if(spacetime) fname <- "CorrelationMat_st_dyn2"
        #if(bivariate) fname <- "CorrelationMat_biv2"
        if(bivariate) fname <- "CorrelationMat_biv_dyn2"

                 cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),
          as.integer(corrmodel), as.double(other_nuis), as.double(paramcorr),as.double(radius),
          as.integer(ns), as.integer(NS),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)  

   # cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",
   #    "integer","double","double","double","integer","integer"),        
     #  #corr=vector_dc("numeric", numpairstot), 
    #   corr=double(numpairstot),
      # coordx,coordy,coordt,corrmodel, nuisance,paramcorr,radius,ns,NS,
      #       INTENT = c("w", "r", "r", "r", "r","r","r","r", "r", "r"),
      #       NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0) 
           bcorr=    gamma(1+1/nuisance['shape'])^2/
           (gamma(1+2/nuisance['shape'])-gamma(1+1/nuisance['shape'])^2)    
         # weibull correlations                

      corr=bcorr*((1-cr$corr^2)^(1+2/nuisance['shape'])*Re(hypergeo::hypergeo(1+1/nuisance['shape'],1+1/nuisance['shape'] ,1 ,cr$corr^2))-1) 
  if(!bivariate) {
        # Builds the covariance matrix:
        varcov <-  diag(1,dime)
        varcov[lower.tri(varcov)] <- corr*(1-nuisance['nugget'])
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr*(1-nuisance['nugget'])
        vv=exp(mu)^2 * (gamma(1+2/nuisance['shape'])/gamma(1+1/nuisance['shape'])^2-1)
        V=vv%*%t(vv)
        varcov=varcov*sqrt(V)
        }
    ## todo
    #if(bivariate)      {
     #     varcov<-diag(dime)
      #    varcov[lower.tri(varcov,diag=T)] <- corr
       #   varcov <- t(varcov)
        #  varcov[lower.tri(varcov,diag=T)] <- corr
        #}
      ###  
}
############################################################### 
    if(model==24)   ## log-logistic case
    {
         sel=substr(names(nuisance),1,4)=="mean"
         mm=as.numeric(nuisance[sel]) 
         mu = X%*%mm
         other_nuis=0# not necessary as.numeric(nuisance[!sel])   ## other nuis parameters (nugget sill skew df)
        fname <-"CorrelationMat2"
        #if(spacetime) fname <- "CorrelationMat_st2"
        if(spacetime) fname <- "CorrelationMat_st_dyn2"
        #if(bivariate) fname <- "CorrelationMat_biv2"
        if(bivariate) fname <- "CorrelationMat_biv_dyn2"
           cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),
          as.integer(corrmodel), as.double(other_nuis), as.double(paramcorr),as.double(radius),
          as.integer(ns), as.integer(NS),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE) 

    # cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",
    #   "integer","double","double","double","integer","integer"),        
    #   #corr=vector_dc("numeric", numpairstot), 
    #   corr=double(numpairstot),
    #   coordx,coordy,coordt,corrmodel, nuisance,paramcorr,radius,ns,NS,
    #         INTENT = c("w", "r", "r", "r", "r","r","r","r", "r", "r"),
    #         NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0)        
    corr= ((pi*sin(2*pi/nuisance['shape']))/(2*nuisance['shape']*
           (sin(pi/nuisance['shape']))^2-pi*sin(2*pi/nuisance['shape'])))*
             (Re(hypergeo::hypergeo(-1/nuisance['shape'],-1/nuisance['shape'] ,1 ,cr$corr^2))*
                 Re(hypergeo::hypergeo(1/nuisance['shape'],1/nuisance['shape'] ,1 ,cr$corr^2)) -1)
  if(!bivariate) {
        varcov <-  diag(1,dime)
        varcov[lower.tri(varcov)] <- corr*(1-nuisance['nugget'])
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- corr*(1-nuisance['nugget'])  
        vv=(exp(mu))^2*
        (2*nuisance['shape']*sin(pi/nuisance['shape'])^2/(pi*sin(2*pi/nuisance['shape']))-1)
        V=vv%*%t(vv)
        varcov=varcov*sqrt(V)
        }
      }
############################################################### 
if(model==22)  {  ## Log Gaussian
            sel=substr(names(nuisance),1,4)=="mean"
            mm=as.numeric(nuisance[sel]) 
            mu = X%*%mm
            other_nuis=as.numeric(nuisance[!sel])   ## other nuis parameters (nugget sill skew df)
         fname <-"CorrelationMat2"
        #if(spacetime) fname <- "CorrelationMat_st2"
        if(spacetime) fname <- "CorrelationMat_st_dyn2"
        #if(bivariate) fname <- "CorrelationMat_biv2"
        if(bivariate) fname <- "CorrelationMat_biv_dyn2"
        vvar=nuisance['sill']
        nuisance['sill']=1

            cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),
          as.integer(corrmodel), as.double(other_nuis), as.double(paramcorr),as.double(radius),
          as.integer(ns), as.integer(NS),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)  

  #cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",
  #     "integer","double","double","double","integer","integer"),        
  #     #corr=vector_dc("numeric", numpairstot), 
  #     corr=double(numpairstot),
  #     coordx,coordy,coordt,corrmodel, nuisance,paramcorr,radius,ns,NS,
  #           INTENT = c("w", "r", "r", "r", "r","r","r","r", "r", "r"),
  #           NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0)      
            corr=cr$corr*(1-nuisance['nugget'])
            corr=exp(vvar*corr)
            #print(corr[1:10])
       if(!bivariate)                  {
           ### cov matrix of the GRF
          varcov <- 1 * diag(dime)
          varcov[lower.tri(varcov)] <- corr
          varcov <- t(varcov)
          varcov[lower.tri(varcov)] <- corr
          diag(varcov)=rep(exp(vvar),dime)
          ones=rep(1,dime)
          varcov=varcov- ones%*%t(ones)
          M=exp(mu-vvar*0.5)%*%t(exp(mu-vvar*0.5))   #exp(2*mu)*(exp(sill*corr)-1) 
          varcov=M*varcov        #ok
          #print(M[1:5,1:5])
         # print( (exp(varcov)-diag(dime))[1:3,1:3])
         # print(varcov[1:5,1:5])
             }
            #   if(bivariate)      {
            #    varcov<-diag(dime)
            #    varcov[lower.tri(varcov,diag=T)] <- corr
            #    varcov <- t(varcov)
            #    varcov[lower.tri(varcov,diag=T)] <- corr
            #}
        } 
###############################################################
         if(model==10)  {  ##  skew Gaussian case

          if(!bivariate){ 
        fname <-"CorrelationMat2"
       # if(spacetime) fname <- "CorrelationMat_st2"
        if(spacetime) fname <- "CorrelationMat_st_dyn2"
      cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),as.integer(corrmodel),
             as.double(nuisance), as.double(paramcorr),as.double(radius),as.integer(ns), as.integer(NS),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
        
 # cr=dotCall64::.C64(fname,SIGNATURE = c("double","double","double","double",
 #      "integer","double","double","double","integer","integer"),        
 #      #corr=vector_dc("numeric", numpairstot), 
 #      corr=double(numpairstot),
 #      coordx,coordy,coordt,corrmodel, nuisance,paramcorr,radius,ns,NS,
 #            INTENT = c("w", "r", "r", "r", "r","r","r","r", "r", "r"),
 #            NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0) 
              corr2=cr$corr^2; sk=as.numeric(nuisance['skew']); sk2=sk^2; vv=as.numeric(nuisance['sill'])
              corr=((2*sk2/pi)*(sqrt(1-corr2) + cr$corr*asin(cr$corr)-1) + cr$corr*vv)/(vv+sk2*(1-2/pi));
                # Builds the covariance matrix:
             varcov <-  diag(1,dime)
             varcov[lower.tri(varcov)] <- corr*(1-nuisance['nugget'])
             varcov <- t(varcov)
             varcov[lower.tri(varcov)] <- corr*(1-nuisance['nugget'])
             vv=nuisance['sill']+nuisance['skew']^2*(1-2/pi)
             varcov=varcov*vv
            }
            if(bivariate)      {
                fname <- "CorrelationMat_biv_skew_dyn2"
                cr=.C(fname, corr=double(numpairstot),  as.double(coordx),as.double(coordy),as.double(coordt),as.integer(corrmodel),
                  as.double(nuisance), as.double(paramcorr),as.double(radius),as.integer(ns),as.integer(NS),PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)
                corr=cr$corr
                varcov <- diag(dime)
                varcov[lower.tri(varcov,diag=T)] <- corr
                varcov <- t(varcov)
                varcov[lower.tri(varcov,diag=T)] <- corr
            }
      } 
      return(varcov)
    }
  #############################################################################################
  #################### end internal function ##################################################
  #############################################################################################
    # Check the user input
    spacetime<-CheckST(CkCorrModel(corrmodel))
    bivariate<-CheckBiv(CkCorrModel(corrmodel))
    ## setting zero mean and nugget if no mean or nugget is fixed
    if(!bivariate){
    if(is.null(param$mean)) param$mean<-0
    if(is.null(param$nugget)) param$nugget<-0  }
    else{
    if(is.null(param$mean_1)) param$mean_1<-0
    if(is.null(param$mean_2)) param$mean_2<-0
    if(is.null(param$nugget_1)) param$nugget_1<-0 
    if(is.null(param$nugget_2)) param$nugget_2<-0 }
    unname(coordt)
    if(is.null(coordx_dyn)){
    unname(coordx);unname(coordy)}
    #if the covariance is compact supported  and option sparse is used
    #then set the code as a tapering and an object spam is returned
    if(sparse) {
    covmod=CkCorrModel(corrmodel)
    if(covmod %in% c(10,11,13,15,19,
                     63,64,65,66,67,68,
                     69,70,71,72,73,74,75,76,77,
                     111,112,129,113,114,131,
                     115,116,120))
    {
      type="Tapering"
      if(bivariate){
      #taper="unit_matrix_biv"
      if(covmod %in% c(111,113,115)) maxdist=c(param$scale,param$scale,param$scale) 
      if(covmod %in% c(112,114,116)) maxdist=c(param$scale_1,param$scale_12,param$scale_2)
      if(covmod %in% c(120,129,131)) maxdist=c(param$scale_1,0.5*(param$scale_1+param$scale_2),param$scale_2)  
      }
      if(spacetime)
      {#taper="unit_matrix_st"
      maxdist=param$scale_s;maxtime=param$scale_t
       if(covmod==63||covmod==65||covmod==67) {  tapsep=c(param$power2_s,param$power_t,param$scale_s,param$scale_t,param$sep) }
       if(covmod==64||covmod==66||covmod==68) {  tapsep=c(param$power_s,param$power2_t,param$scale_s,param$scale_t,param$sep) }
    }
      if(!(spacetime||bivariate)){
        #taper="unit_matrix";
        maxdist=param$scale}
  }
  taper=corrmodel
}
    checkinput <- CkInput(coordx, coordy, coordt, coordx_dyn, corrmodel, NULL, distance, "Simulation",
                             NULL, grid, NULL, maxdist, maxtime,  model=model, n,  NULL,
                              param, radius, NULL, taper, tapsep,  "Standard", NULL, NULL, NULL,X)

    if(!is.null(checkinput$error)) stop(checkinput$error)
    spacetime_dyn=FALSE
    if(!is.null(coordx_dyn))  spacetime_dyn=TRUE
    # Initialising the parameters:
    initparam <- StartParam(coordx, coordy, coordt,coordx_dyn, corrmodel, NULL, distance, "Simulation",
                           NULL, grid, NULL, maxdist, maxtime, model, n, 
                           param, NULL, NULL, radius, NULL, taper, tapsep,  type, type,
                           NULL, NULL, FALSE, NULL, NULL,NULL,NULL,X)

    if(grid) cc=expand.grid(initparam$coordx,initparam$coordy)
    else     cc=cbind(initparam$coordx,initparam$coordy)  


    if(!spacetime_dyn) dime=initparam$numcoord*initparam$numtime 
    else               dime=sum(initparam$ns)

    if(!initparam$bivariate) numpairstot=dime*(dime-1)*0.5
    if(initparam$bivariate)  numpairstot=dime*(dime-1)*0.5+dime
    if(!is.null(initparam$error)) stop(initparam$error)
    setup<-initparam$setup
    if(initparam$type=="Tapering")
    {
      if(initparam$spacetime) fname= "CorrelationMat_st_tap"
      if(!initparam$spacetime) fname= "CorrelationMat_tap"
      if(initparam$bivariate) fname= "CorrelationMat_biv_tap"
      corr <- double(initparam$numpairs)
      #tapmod <- setup$tapmodel
      ### unit taperssss ####
      if(sparse){
           if(spacetime) tapmod=230
           if(bivariate) tapmod=147
           if(!(spacetime||bivariate)) tapmod=36
           }
      else(tapmod=CkCorrModel(taper))
    #######################
    tp=.C(fname,tapcorr=double(initparam$numpairs),as.double(cc[,1]),as.double(cc[,2]),as.double(initparam$coordt),as.integer(tapmod),
      as.double(1),as.double(tapsep),as.double(1),
      PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)
        setup$taps<-tp$tapcorr
    }
    if(is.null(X))  initparam$X=as.matrix(rep(1,dime))
    if(spacetime||bivariate){
          initparam$NS=cumsum(initparam$ns);
            if(spacetime_dyn){  initparam$NS=c(0,initparam$NS)[-(length(initparam$ns)+1)]}
            else{               initparam$NS=rep(0,initparam$numtime)}
    }
    if(is.null(initparam$NS)) initparam$NS=0
    covmatrix<- Cmatrix(initparam$bivariate,cc[,1],cc[,2],initparam$coordt,initparam$corrmodel,dime,n,initparam$ns,
                        initparam$NS,
                        initparam$param[initparam$namesnuis],
                        initparam$numpairs,numpairstot,initparam$model,
                        initparam$param[initparam$namescorr],setup,initparam$radius,initparam$spacetime,spacetime_dyn,initparam$type,initparam$X)

    initparam$param=initparam$param[names(initparam$param)!='mean']
    if(sparse==TRUE) if(!spam::is.spam(covmatrix)) covmatrix=spam::as.spam(covmatrix)
    if(type=="Tapering") sparse=TRUE

    # Delete the global variables:
    .C('DeleteGlobalVar', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
    if(initparam$bivariate)   initparam$numtime=2
    # Return the objects list:
    CovMat <- list(bivariate =  initparam$bivariate,
                   coordx = initparam$coordx,
                   coordy = initparam$coordy,
                   coordt = initparam$coordt,
                   coordx_dyn = coordx_dyn,
                   covmatrix=covmatrix,
                   corrmodel = corrmodel,
                   distance = distance,
                   grid=   grid,
                   nozero=initparam$setup$nozero,
                   maxdist = maxdist,
                   maxtime = maxtime,
                   n=n,
                   ns=initparam$ns,
                   NS=initparam$NS,
                   model=initparam$model,
                   namescorr = initparam$namescorr,
                   namesnuis = initparam$namesnuis,
                   namessim = initparam$namessim,
                   numblock = initparam$numblock,
                   numcoord = initparam$numcoord,
                   numcoordx = initparam$numcoordx,
                   numcoordy = initparam$numcoordy,
                   numtime = initparam$numtime,
                   param = initparam$param,
                   radius = initparam$radius,
                   setup=setup,
                   spacetime = initparam$spacetime,
                   sparse=sparse,
                   tapmod=taper,
                   tapsep=tapsep,
                   X=initparam$X)
    structure(c(CovMat, call = call), class = c("CovMat"))
}

