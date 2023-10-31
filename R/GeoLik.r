####################################################
### File name: GeoLik.r
####################################################


### Optim call for log-likelihood maximization 
Lik <- function(copula,bivariate,coordx,coordy,coordt,coordx_dyn,corrmodel,data,fixed,flagcor,flagnuis,grid,lower,
                       mdecomp,model,namescorr,namesnuis,namesparam,numcoord,numpairs,numparamcor,numtime,
                       optimizer,onlyvar,parallel,param,radius,setup,spacetime,sparse,varest,taper,type,upper,ns,X,neighb,MM,aniso)
{
 ######### computing upper trinagular of covariance matrix   
    matr <- function(corrmat,corr,coordx,coordy,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
    {

        #cc <- .C(corrmat,cr=corr,as.double(coordx),as.double(coordy),as.double(coordt),as.integer(corrmodel),as.double(nuisance),
        #as.double(paramcorr),as.double(radius),as.integer(ns),as.integer(NS),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)$cr

           cc=dotCall64::.C64(as.character(corrmat),
         SIGNATURE = c("double","double","double","double", "integer","double","double","double","integer","integer"),  
                          cr=corr, coordx, coordy, coordt, corrmodel, nuisance,paramcorr,radius, ns,NS,
         INTENT =    c("rw","r","r","r","r","r","r","r","r", "r"),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$cr
        return(cc)
    }
   ######### computing upper trinagular of covariance matrix   it is for poisson
    matr2 <- function(corrmat,corr,coordx,coordy,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius,model,mu)
    {
       # cc <- .C(corrmat,cr=corr,as.double(coordx),as.double(coordy),as.double(coordt),as.integer(corrmodel),
       #  as.double(c(mu)), as.integer(1), as.double(nuisance['nugget']),
       # as.double(paramcorr),as.double(radius),as.integer(ns),as.integer(NS),as.integer(model),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)$cr
        
        hh=1;nn=nuisance['nugget'];mm=c(mu)
  cc=dotCall64::.C64(as.character(corrmat),
         SIGNATURE = c("double","double","double","double", "integer","double", "integer","double","double","double","integer","integer","integer"),  
                          cr=corr, coordx, coordy, coordt, corrmodel, mm,hh,nn,paramcorr,radius, ns,NS,model,
         INTENT =    c("rw","r","r","r","r","r","r","r","r","r","r","r","r"),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$cr
        return(cc)
    }
    ### START Defining the objective functions
######### Restricted log-likelihood for multivariate normal density:
    LogNormDenRestr <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        varcov <- (nuisance['sill'])*ident
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova
        #  decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(varcov,mdecomp)
        if(is.logical(decompvarcov)) return(llik)       
        logdetvarcov <- MatLogDet(decompvarcov,mdecomp) 
        ivarcov <- MatInv(decompvarcov,mdecomp)
        sumvarcov <- sum(ivarcov)
        p <- ivarcov-array(rowSums(ivarcov),c(dimat,1))%*%colSums(ivarcov)/sumvarcov
        llik <- 0.5*(const+logdetvarcov+log(sumvarcov)+crossprod(t(crossprod(stdata,p)),stdata))
        return(llik)
    }
######### Restricted log-likelihood for bivariate multivariate normal density:
         LogNormDenRestr_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        ident[lower.tri(ident,diag=T)] <- cova
        ident <- t(ident)
        ident[lower.tri(ident,diag=T)] <- cova
        #  decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(ident,mdecomp)
        if(is.logical(decompvarcov)) return(llik)       
        logdetvarcov <- MatLogDet(decompvarcov,mdecomp) 
        ivarcov <- MatInv(decompvarcov,mdecomp)
        sumvarcov <- sum(ivarcov)
        p <- ivarcov-array(rowSums(ivarcov),c(dimat,1))%*%colSums(ivarcov)/sumvarcov
        llik <- 0.5*(const+ logdetvarcov +log(sumvarcov)+crossprod(t(crossprod(stdata,p)),stdata))
        return(llik)
    }
    
    
######### Tapering 2 log-likelihood for multivariate normal density:
    LogNormDenTap <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        lliktap <- 1.0e8
        # Computes the vector of the correlations:
        varcovtap <- new("spam",entries=cova*setup$taps,colindices=setup$ja,
        rowpointers=setup$ia,dimension=as.integer(rep(dimat,2)))
        cholvctap <- spam::update.spam.chol.NgPeyton(setup$struct, varcovtap)
        #if(class(cholvctap)=="try-error") return(lliktap)
        logdet <- c(spam::determinant(cholvctap)$modulus)
        inv <- spam::solve.spam(cholvctap)
        slot(varcovtap,"entries") <- inv[setup$idx]*setup$taps
        lliktap= 0.5*(const+2*logdet+drop(t(stdata)%*%varcovtap%*%stdata))
        return(lliktap)
    }

######### Tapering 1 sas log-likelihood 
LogshDenTap1<- function(const,cova,ident,dimat,mdecomp,nuisance,setup,sill,stdata)
{

    varcovtap <- new("spam",entries=setup$taps*cova,colindices=setup$ja,
    rowpointers=setup$ia,dimension=as.integer(rep(dimat,2)))
   
          #cholvctap <-spam::chol.spam(varcovtap)
           cholvctap <- spam::update.spam.chol.NgPeyton(setup$struct, varcovtap)
          logdet <- 2*c(spam::determinant(cholvctap)$modulus)

    ################################################
        skew=as.numeric(nuisance["skew"])
        delta=as.numeric(nuisance["tail"])
        Z=sinh(delta * asinh(stdata)-skew)
        C=delta*sqrt((1+Z^2)/(stdata^2+1))
    llik <- 0.5*( const + const*log(sill)/log(2*pi)+logdet - 2*sum(log(C)) +sum(Z* spam::solve.spam(cholvctap, Z)))
                    
    return(llik)
}


######### Tapering 1 normal log-likelihood 
LogNormDenTap1 <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
{

    varcovtap <- new("spam",entries=cova*setup$taps,colindices=setup$ja,
    rowpointers=setup$ia,dimension=as.integer(rep(dimat,2)))
    cholvctap <- spam::update.spam.chol.NgPeyton(setup$struct, varcovtap)
    #cholvctap <-spam::chol.spam(varcovtap)
    logdet <- c(spam::determinant(cholvctap)$modulus)
    lliktap= 0.5*(const+2*logdet+sum(stdata* spam::solve.spam(cholvctap, stdata)))
    return(lliktap)
}

    
######### Tapering 2 log-likelihood for bivariate multivariate normal density:
    LogNormDenTap_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        lliktap <- 1.0e8
        # Computes the vector of the correlations:
        varcovtap <- try(new("spam",entries=cova*setup$taps,colindices=setup$ja,
        rowpointers=setup$ia,dimension=as.integer(rep(dimat,2))),silent=TRUE)
        cholvctap <- try(spam::chol.spam(varcovtap),silent=TRUE)
        if(inherits(cholvctap,"try-error")) {return(lliktap)}
        logdet <- c(spam::determinant(cholvctap)$modulus)
        inv <- spam::solve.spam(cholvctap)
        slot(varcovtap,"entries") <- inv[setup$idx]*setup$taps
        lliktap= 0.5*(const+2*logdet+drop(t(stdata)%*%varcovtap%*%stdata))
        return(lliktap)
    }
######### Tapering 1 log-likelihood for bivariate multivariate normal density:
    LogNormDenTap1_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        lliktap <- 1.0e8
        # Computes the vector of the correlations:
        cova[cova==(nuisance['sill'])] <- nuisance['sill']+nuisance['nugget']
        varcovtap <- new("spam",entries=cova*setup$taps,colindices=setup$ja,
        rowpointers=setup$ia,dimension=as.integer(rep(dimat,2)))
        cholvctap <- spam::update.spam.chol.NgPeyton(setup$struct, varcovtap)
        if(inherits(cholvctap, "try-error"))  {return(lliktap)}
        logdet <- c(spam::determinant.spam.chol.NgPeyton(cholvctap)$modulus)
        lliktap <- 0.5*(const+2*logdet+sum(stdata* spam::solve.spam(cholvctap, stdata)))
        return(lliktap)
    }

######### Standard log-likelihood function for multivariate normal density with sparse alg matrices
    LogNormDenStand_spam <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        varcov <- (nuisance['sill'])*ident
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova
        mcov=spam::as.spam(varcov)#,silent=TRUE);if(class(mcov)=="try-error")   return(llik)
        cholS <- spam::chol.spam(mcov)# ,silent=TRUE);if(class(cholS)=="try-error")   return(llik)         
        llik=0.5*( const+2*c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus) 
                                    + sum(stdata* spam::solve.spam(cholS,stdata)))
        return(llik)
    }    
######### Standard log-likelihood function for multivariate normal density:
    LogNormDenStand <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        varcov <- (nuisance['sill'])*ident
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova    
        # decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(varcov,mdecomp)
        if(is.logical(decompvarcov)) return(llik)  
        logdetvarcov <- MatLogDet(decompvarcov,mdecomp) 
        llik <- 0.5*(const+logdetvarcov+  sum((backsolve(decompvarcov, stdata, transpose = TRUE))^2))
        #llik <- 0.5*(const+logdetvarcov +  
         #    sum(stdata * backsolve(decompvarcov, forwardsolve(decompvarcov, stdata, transpose = TRUE, upper.tri = TRUE))))
        return(llik)
    }
    LogNormDenStand22 <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        varcov <- cova      
        # decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(varcov,mdecomp)
        if(is.logical(decompvarcov)) return(llik)  
        logdetvarcov <- MatLogDet(decompvarcov,mdecomp) 
       # invarcov <- MatInv(decompvarcov,mdecomp)
       # llik <- 0.5*(const+logdetvarcov+crossprod(t(crossprod(stdata,invarcov)),stdata))
         llik <- 0.5*(const+logdetvarcov+  sum((backsolve(decompvarcov, stdata, transpose = TRUE))^2))
        return(llik)
    }
    ######### CVV mdecomp:
    CVV <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        varcov <- (nuisance['sill'])*ident
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova      
        # decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(varcov,mdecomp)
        if(is.logical(decompvarcov)) return(llik)  
        invarcov <- MatInv(decompvarcov,mdecomp)
        D=diag(1/diag(invarcov))
        M=crossprod(invarcov,D);C=tcrossprod(M,M)
        llik <- mean(crossprod(t(crossprod(stdata,C)),stdata))
        return(llik)
    }

    ######### CVV mdecomp in the bivariate case
CVV_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        ident[lower.tri(ident,diag=T)] <- cova
        ident <- t(ident)
        ident[lower.tri(ident,diag=T)] <- cova
        # decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(ident,mdecomp)
        if(is.logical(decompvarcov)) return(llik)
        invarcov <- MatInv(decompvarcov,mdecomp)
        D=diag(1/diag(invarcov))
        M=crossprod(invarcov,D);C=tcrossprod(M,M)
        llik <- mean(crossprod(t(crossprod(stdata,C)),stdata))
        return(llik)
    }

 ######## Standard log-likelihood function for log gaussian random fields       
  LogNormDenStand_LG <- function(const,cova,ident,dimat,mdecomp,nuisance,det,sill,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        varcov <- (sill)*ident
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova      
        # decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(varcov,mdecomp)
        if(is.logical(decompvarcov)) return(llik)  
        logdetvarcov <- MatLogDet(decompvarcov,mdecomp) 
       # invarcov <- MatInv(decompvarcov,mdecomp)
       # llik <- 0.5*(const+logdetvarcov+crossprod(t(crossprod(stdata,invarcov)),stdata))
         llik <- 0.5*(const+logdetvarcov+2*det+  sum((backsolve(decompvarcov, stdata, transpose = TRUE))^2))
        return(llik)
    }

######## Standard log-likelihood function for tukeyH random fields

    LogNormDenStand_TukeyH <- function(const,cova,ident,dimat,mdecomp,nuisance,sill,setup,stdata)
    {
        llik <- 1.0e8
   # Computes the covariance matrix:
        varcov <- ident
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova      
        # decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(varcov,mdecomp)
        if(is.logical(decompvarcov)) return(llik)  
        logdetvarcov <- MatLogDet(decompvarcov,mdecomp) 
################################################
        delta=nuisance["tail"]
        vv=sqrt(VGAM::lambertW(delta*stdata^2)/delta)
        IL=sign(stdata)*vv
        IW=1/(stdata*(1+VGAM::lambertW(delta*stdata^2)))
        llik <- 0.5*( const*log(sill)/log(2*pi) + 
                      const + logdetvarcov + sum((backsolve(decompvarcov, IL, transpose = TRUE))^2)
                      - 2*sum(log(IL*IW)))
        return(llik)
    }


 LogNormDenStand_Tukey2H <- function(const,cova,ident,dimat,mdecomp,nuisance,sill,setup,stdata)
    {
        llik <- 1.0e8
   # Computes the covariance matrix:
        varcov <- ident
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova      
        # decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(varcov,mdecomp)
        if(is.logical(decompvarcov)) return(llik)  
        logdetvarcov <- MatLogDet(decompvarcov,mdecomp) 
################################################
        delta1=nuisance["tail1"]
        delta2=nuisance["tail2"]
a=which(stdata>=0); b=which(stdata<0)
stmas=stdata[stdata>=0]; stmenos=stdata[stdata<0]

g1<-(VGAM::lambertW(delta1*(stmas)^2)/(delta1))^(1/2)
jac1<- g1/(stmas*(1+VGAM::lambertW(delta1*(stmas)^2)))

g2<-(VGAM::lambertW(delta2*(stmenos)^2)/(delta2))^(1/2)
jac2<- -g2/(stmenos*(1+VGAM::lambertW(delta2*(stmenos)^2)))

pp=data.frame(rbind(cbind(a,g1),cbind(b,g2)))
qq=data.frame(rbind(cbind(a,jac1),cbind(b,jac2)))

g=pp[with(pp, order(pp$a)), ] 
jac=qq[with(qq, order(qq$a)), ] 
tau_inv=sign(stdata)*c(g$g1)

llik <- 0.5*( const*log(sill)/log(2*pi) + 
                      const + logdetvarcov + sum((backsolve(decompvarcov, c(tau_inv), transpose = TRUE))^2)- 2*sum(log(jac$jac1)))
        return(llik)
    }

######## Standard log-likelihood function for SH random fields
    LogNormDenStand_SH <- function(const,cova,ident,dimat,mdecomp,nuisance,sill,setup,stdata)
    {
        llik <- 1.0e8
   # Computes the covariance matrix:
        varcov <- ident
        varcov[lower.tri(varcov)] <- cova
        varcov <- t(varcov)
        varcov[lower.tri(varcov)] <- cova      
        # decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(varcov,mdecomp)
        if(is.logical(decompvarcov)) return(llik)  
        logdetvarcov <- MatLogDet(decompvarcov,mdecomp) 
################################################
        skew=as.numeric(nuisance["skew"])
        delta=as.numeric(nuisance["tail"])
        Z=sinh(delta * asinh(stdata)-skew)
        C=delta*sqrt((1+Z^2)/(stdata^2+1))
        llik <- 0.5*( const + const*log(sill)/log(2*pi)
                      +logdetvarcov - 2*sum(log(C))
                      +sum((backsolve(decompvarcov, Z, transpose = TRUE))^2))
        return(llik)
    }
       
       
######### Standard log-likelihood function for multivariate bivariate normal density:
 LogNormDenStand_biv <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        ident[lower.tri(ident,diag=T)] <- cova
        ident <- t(ident)
        ident[lower.tri(ident,diag=T)] <- cova
        # decomposition of the covariance matrix:
        decompvarcov <- MatDecomp(ident,mdecomp)
        if(is.logical(decompvarcov)) return(llik)
        logdetvarcov <- MatLogDet(decompvarcov,mdecomp) 
        #invarcov <- MatInv(decompvarcov,mdecomp)
        #llik <- 0.5*(const+logdetvarcov+crossprod(t(crossprod(stdata,invarcov)),stdata))
        llik <- 0.5*(const+logdetvarcov+  sum((backsolve(decompvarcov, stdata, transpose = TRUE))^2))
        return(llik)

        
    }
    ######### Standard log-likelihood function for multivariate bivariate normal density with sparse alg matrices
 LogNormDenStand_biv_spam <- function(const,cova,ident,dimat,mdecomp,nuisance,setup,stdata)
    {
        llik <- 1.0e8
        # Computes the covariance matrix:
        ident[lower.tri(ident,diag=T)] <- cova
        ident <- t(ident)
        ident[lower.tri(ident,diag=T)] <- cova
        # decomposition of the covariance matrix:
        mcov=try(spam::as.spam(ident),silent=TRUE)
        if(inherits(mcov,"try-error")){return(llik)}
        cholS <- try(chol(mcov) ,silent=T);if(inherits(cholS,"try-error")){ return(llik)}
        llik=0.5*( const+2*c(spam::determinant.spam.chol.NgPeyton(cholS)$modulus) 
                                    + sum(stdata* spam::solve.spam(cholS,stdata)))
         return(llik)
    }
    ### END Defining the objective functions
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
        # Call to the objective functions:
         loglik_loggauss <- function(param,const,coordx,coordy,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM)
    {

        llik <- 1.0e8
        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"

        mm=as.numeric(nuisance[sel])
        Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM
      
        # Computes the vector of the correlations:
        if(nuisance['nugget']<0||nuisance['nugget']>=1) return(llik)
        sill=nuisance['sill']
        if(sill<0) return(llik)
        nuisance['sill']=1
        corr=matr(corrmat,corr,coordx,coordy,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
       ## if(corr[1]==-2||is.nan(corr[1])) return(llik)
        # Computes the correlation matrix:
        cova <- corr*(1-nuisance['nugget'])
        #KK=exp(sill)/2
        # Computes the log-likelihood
        loglik_u <- do.call(what="LogNormDenStand_LG",
            args=list(stdata=(log(data)-(c(Mean-sill*0.5))),const=const,cova=cova,dimat=dimat,ident=ident,
            mdecomp=mdecomp,nuisance=nuisance,det=sum(1/(data)),sill=sill,setup=setup))
        return(loglik_u)
      }
################################################################################################
   loglik_tukey2h <- function(param,const,coordx,coordy,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM)
    {

        llik <- 1.0e8

        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
         mm=as.numeric(nuisance[sel])
                 Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM
   
        #if(nuisance['tail']<0||nuisance['tail1']>0.5||nuisance['tail2']>0.5||nuisance['nugget']<0||nuisance['nugget']>=1) return(llik)
        # Computes the vector of the correlations:
        sill=nuisance['sill']
        nuisance['sill']=1
         corr=matr(corrmat,corr,coordx,coordy,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
         corr= corr*(1-nuisance['nugget'])
        loglik_u <- do.call(what="LogNormDenStand_Tukey2H",
            args=list(stdata=((data-c(Mean))/(sqrt(sill))),const=const,cova=corr,dimat=dimat,ident=ident,
            mdecomp=mdecomp,nuisance=nuisance,sill=sill,setup=setup))
        return(loglik_u)
      }
################################################################################################
    loglik_tukeyh <- function(param,const,coordx,coordy,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM)
    {

        llik <- 1.0e8

        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
          mm=as.numeric(nuisance[sel])
               Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM
      
     if(nuisance['tail']<0||nuisance['tail']>0.5||nuisance['nugget']<0||nuisance['nugget']>=1||nuisance['sill']<0) return(llik)
        # Computes the vector of the correlations:
        sill=nuisance['sill']
        nuisance['sill']=1
         corr=matr(corrmat,corr,coordx,coordy,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
         corr= corr*(1-nuisance['nugget'])
        loglik_u <- do.call(what="LogNormDenStand_TukeyH",
            args=list(stdata=((data-c(Mean))/(sqrt(sill))),const=const,cova=corr,dimat=dimat,ident=ident,
            mdecomp=mdecomp,nuisance=nuisance,sill=(sill),setup=setup))

        return(loglik_u)
      }


################################################################################################
   # Call to the objective functions:
        loglik_miss_skewT<- function(param,const,coordx,coordy,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM)
    {
        llik <- 1.0e8
        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
                  mm=as.numeric(nuisance[sel])
          Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM
   
        # Computes the vector of the correlations:
        corr=matr(corrmat,corr,coordx,coordy,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
        nu=1/nuisance['df']; eta2=nuisance['skew']^2
        #if(nu<2||abs(nuisance['skew'])>1)  return(llik)
        w=sqrt(1-eta2);
        KK=2*eta2/pi
        D1=(nu-1)/2; D2=nu/2
        CorSkew<-(2*eta2/(pi*w^2+eta2*(pi-2)))*(sqrt(1-corr^2)+corr*asin(corr)-1)+w^2*corr/(w^2+eta2*(1-2/pi))   
        corr3<-(pi*(nu-2)*gamma(D1)^2/(2*(pi*gamma(D2)^2-eta2*(nu-2)*gamma(D1)^2)))*(Re(hypergeo::hypergeo(0.5,0.5,D2,corr^2))*((1-KK)*CorSkew+KK)-KK)
    #   if(is.nan(corr[1])||nuisance['sill']<0||nuisance['nugget']<0||nuisance['nugget']>1) return(llik)
        cova <- corr3*nuisance['sill'] *(1-nuisance['nugget'])
       #nuisance['nugget']=0
      loglik_u <- do.call(what="LogNormDenStand",args=list(stdata=(data-c(Mean)),const=const,cova=cova,dimat=dimat,ident=ident,
            mdecomp=mdecomp,nuisance=nuisance,setup=setup))
        return(loglik_u)
      }
################################################################################################
    # Call to the objective functions:
    loglik_miss_T <- function(param,const,coordx,coordy,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM)
    {
        llik <- 1.0e8
        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
           mm=as.numeric(nuisance[sel])
                 Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM
      
        # Computes the vector of the correlations:
        corr=matr(corrmat,corr,coordx,coordy,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
        df=1/nuisance['df']
        #if(df<2)  return(llik)
         #if(df<170) corr=(df-2)*gamma((df-1)/2)^2/(2*gamma(df/2)^2)* corr *Re(hypergeo::hypergeo(0.5,0.5,df/2,corr^2)) 
         #else      
    corr=exp(log(df-2)+2*lgamma(0.5*(df-1))-(log(2)+2*lgamma(df/2))+log(Re(hypergeo::hypergeo(0.5,0.5, df/2,corr^2)))+log(corr))
        #if(is.nan(corr[1])||nuisance['sill']<0||nuisance['nugget']<0||nuisance['nugget']>1) return(llik)
    cova <- corr*(nuisance['sill'])*(1-nuisance['nugget'])
        #nuisance['nugget']=0
      loglik_u <- do.call(what="LogNormDenStand",args=list(stdata=(data-c(Mean)),const=const,cova=cova,dimat=dimat,ident=ident,
            mdecomp=mdecomp,nuisance=nuisance,setup=setup))
        return(loglik_u)
      }
################################################################################################
    # Call to the objective functions:
    loglik_miss_Pois <- function(param,const,coordx,coordy,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM)
    {
        llik <- 1.0e8
        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
           mm=as.numeric(nuisance[sel])
           Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM
      
        # Computes the vector of the correlations:
        mu=Mean
        model=30
        corr=matr2(corrmat,corr,coordx,coordy,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius,model,mu)
        cova <-  ident
        cova[lower.tri(cova)] <- corr   
        cova <- t(cova)
        cova[lower.tri(cova)] <- corr 
        diag(cova)=exp(mu)
    if(nuisance['nugget']<0||nuisance['nugget']>1) return(llik)
        #nuisance['nugget']=0
      loglik_u <- do.call(what="LogNormDenStand22",args=list(stdata=data-c(exp(mu)),
           const=const,cova=cova,dimat=dimat,ident=ident,
            mdecomp=mdecomp,nuisance=nuisance,setup=setup))
        return(loglik_u)
    
      }
      ################################################################################################ 
loglik_sh <- function(param,const,coordx,coordy,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM)
    {

        llik <- 1.0e8

        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
          mm=as.numeric(nuisance[sel])
           Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM
         # Computes the vector of the correlations:
        sill=as.numeric(nuisance['sill'])
        nuisance['sill']=1
        if(as.numeric(nuisance['tail'])<=0||as.numeric(nuisance['sill'])<=0||as.numeric(nuisance['nugget'])<0||as.numeric(nuisance['nugget'])>=1) return(llik)
        corr=matr(corrmat,corr,coordx,coordy,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
        #if(is.nan(corr[1])) return(llik)
        corr= corr*(1-nuisance['nugget'])
        loglik_u <- do.call(what=fname,#"LogNormDenStand_SH",
            args=list(stdata=((data-c(Mean))/(sqrt(sill))),const=const,cova=corr,dimat=dimat,ident=ident,
            mdecomp=mdecomp,nuisance=nuisance,sill=sill,setup=setup))
        return(loglik_u)
      }
################################################################################################
    # Call to the objective functions:
    loglik <- function(param,const,coordx,coordy,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM)
    {

        llik <- 1.0e8
        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
        mm=as.numeric(nuisance[sel])
        Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM
              # Computes the vector of the correlations:
        corr=matr(corrmat,corr,coordx,coordy,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
        if(is.nan(corr[1])||nuisance['sill']<0||nuisance['nugget']<0||nuisance['nugget']>1) return(llik)
        cova <- corr*nuisance['sill']*(1-nuisance['nugget'])
        
      loglik_u <- do.call(what=fname,args=list(stdata=data-c(Mean),const=const,cova=cova,dimat=dimat,ident=ident,
            mdecomp=mdecomp,nuisance=nuisance,setup=setup))

        return(loglik_u)
      }
  #####################################################    

 loglikvecchia <- function(param,vecchia.approx,data,fixed,dimat,
                    model,namescorr,namesnuis,namesparam,X,MM)
    {
        llik <- 1.0e8
        names(param) <- namesparam
        # Set the parameter vector:
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]
        sel=substr(names(nuisance),1,4)=="mean"
        mm=as.numeric(nuisance[sel])
        Mean=c(X%*%mm)
        if(!is.null(MM)) Mean=MM
   
          nuggets=as.numeric(nuisance['nugget'])+ 0.00001
        data=c(data-X%*%mm)
        ppar=as.numeric(c(nuisance['sill'], paramcorr[1], paramcorr[2]))
    if(ppar[2]<0|| ppar[3]<0||nuisance['sill']<0||nuisance['nugget']<0||nuisance['nugget']>1){return(llik)}
        loglik_u=GPvecchia::vecchia_likelihood(data,vecchia.approx,
                         covparms=ppar,nuggets=nuggets,covmodel ="matern")
        return(-loglik_u)
      }

      
################################################################################################                     
     loglik_biv <- function(param,const,coordx,coordy,coordt,corr,corrmat,corrmodel,data,dimat,fixed,fname,
                       grid,ident,mdecomp,model,namescorr,namesnuis,namesparam,radius,setup,X,ns,NS,MM)
      {

        # Set the parameter vector:
        names(param) <- namesparam
        pram <- c(param, fixed)
        paramcorr <- pram[namescorr]
        nuisance <- pram[namesnuis]

        sel1=substr(names(nuisance),1,6)=="mean_1"
        mm1=as.numeric(nuisance[sel1])
        sel2=substr(names(nuisance),1,6)=="mean_2"
        mm2=as.numeric(nuisance[sel2])
       
        X1=as.matrix(X[1:ns[1],]);X2=as.matrix(X[(ns[1]+1):(ns[2]+ns[1]),]); 
        mm=as.double(c(X1%*%mm1,X2%*%mm2))
        # Standardizes the data:
         stdata <- data-mm 
      # Computes the vector of the correlations
         corr=matr(corrmat,corr,coordx,coordy,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
      # Computes the log-likelihood
       loglik_b <- do.call(what=fname,args=list(stdata=stdata,const=const,cova=corr,ident=ident,dimat=dimat,
            mdecomp=mdecomp,nuisance=nuisance,setup=setup))

        return(loglik_b)
      }


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    ### START the main code of the function:

    spacetime_dyn=FALSE; NS=0;fname=NULL
    if(!is.null(coordx_dyn)) spacetime_dyn=TRUE
    if(grid)     {a=expand.grid(coordx,coordy);coordx=a[,1];coordy=a[,2]; }

    ####################################
    if(!spacetime_dyn) dimat <- numcoord*numtime# length of data
    if(spacetime_dyn)  dimat =sum(ns)

    if(is.null(dim(X))) {
    if(!bivariate) X=as.matrix(rep(1,dimat))  # matrix of covariates
    if( bivariate) {X=as.matrix(rep(1,ns[1]+ns[2]));# X=rbind(X,X)
        }}
    else{
    if(!bivariate) num_betas=ncol(X)  
    if( bivariate) num_betas=c(ncol(X),ncol(X)) }
     
    corrmat<-"CorrelationMat"# set the type of correlation matrix     
    if(model==36) corrmat<-"CorrelationMat_dis"# set the type of correlation matrix 
    if(spacetime)  { corrmat<-"CorrelationMat_st_dyn"
                     if(model==36) corrmat="CorrelationMat_st_dyn_dis"
                    } 
    if(bivariate)  corrmat<-"CorrelationMat_biv_dyn"  
     if(spacetime||bivariate){
          NS=cumsum(ns);
          if(spacetime_dyn){ data=t(unlist(data));NS=c(0,NS)[-(length(ns)+1)]}
          else {data=matrix(t(data),nrow=1);NS=rep(0,numtime)}  
    }       
    ####################################
    dd=0
     
    if(bivariate)  dd=dimat
    numpairstot <- dimat*(dimat-1)/2+dd
    const<-dimat*log(2*pi)# set the likelihood constant

     if(is.null(neighb))
     {
     corr<-double(numpairstot)# initialize the correlation
     ident <- diag(dimat)# set the identity matrix}
     }

 #########################################################################  
 ######################################################################### 
if(model==1||model==20){  ## gaussian case
    lname <- 'loglik'
    if(!is.null(neighb))  lname <-'loglikvecchia'
    if(bivariate)  {lname <- 'loglik_biv'}
  
    # detects the type of likelihood:
    if(type==3){
        if(bivariate) fname<-"LogNormDenRestr_biv"
        else          fname<-"LogNormDenRestr"
        const <- const-log(2*pi)}
    if(type==4) { if(bivariate)   fname <- 'LogNormDenStand_biv'
                  else            fname <- 'LogNormDenStand'
                }
    if(type==5||type==6){  #tapering
        corrmat<-"CorrelationMat_tap"
        if(spacetime) corrmat<-"CorrelationMat_st_tap"
        if(bivariate) {if(type==5) fname <- 'LogNormDenTap_biv'
                       if(type==6) fname <- 'LogNormDenTap1_biv'
                      corrmat<-"CorrelationMat_biv_tap" }
        else          {if(type==5) fname <- 'LogNormDenTap'
                       if(type==6)  fname <- 'LogNormDenTap1'
                                   
                      }             
        corr <- double(numpairs)
        tapcorr <- double(numpairs)

           tcor=matr(corrmat,tapcorr,coordx,coordy,coordt,setup$tapmodel,c(0,0,1),1,ns,NS,radius)

           tape <- new("spam",entries=tcor,colindices=setup$ja,rowpointers=setup$ia,dimension=as.integer(rep(dimat,2)))
           setup$struct <- try(spam::chol.spam(tape,silent=TRUE))
           setup$taps<-tcor

        }
        if(type==8)
        { if(bivariate) fname<-"CVV_biv"
          else          fname<-"CVV"
       }
     if(sparse) fname <- paste(fname,"_spam",sep="")
 }
 ############################################################################
############################################################################
hessian=FALSE
 if(model==20){   ## SAS case
     lname <- 'loglik_sh'
    if(bivariate)  {lname <- 'loglik_biv_sh'}
    fname='LogNormDenStand_SH'
    if(type==6) fname <- 'LogshDenTap1'

}
 
 if(model==34){   ## Tukeyh case
     lname <- 'loglik_tukeyh'
    if(bivariate)  {lname <- 'loglik_biv_tukeyh'}

}

 if(model==40){   ## Tukey2h case
     lname <- 'loglik_tukey2h'
    if(bivariate)  {lname <- 'loglik_biv_tukey2h'}
   
}


 if(model==35){   ## gaussian misspecified t
     lname <- 'loglik_miss_T'
    if(bivariate)  {lname <- 'loglik_biv_miss_T'}

}

 if(model==36){   ## Poisson misspecified t
     lname <- 'loglik_miss_Pois'
    if(bivariate)  {lname <- 'loglik_biv_miss_Pois'}

}


 if(model==37){   ## gaussian misspecified skewt
     lname <- 'loglik_miss_skewT'
    if(bivariate)  {lname <- 'loglik_biv_miss_skewT'}

}


 if(model==22){   ## loggaussian  case
     lname <- 'loglik_loggauss'
    if(bivariate)  {lname <- 'loglik_biv_loggauss'}

}


 if(type!=5&&type!=6){ corrmat <- paste(corrmat,"2",sep="") }



##################

if(!onlyvar){   # performing optimization
    maxit=10000
    # Optimize the log-likelihood:
   if(length(param)==1)
        {
         optimizer="optimize"         
  Likelihood <- optimize(f=eval(as.name(lname)),const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,
                          fname=fname,grid=grid,ident=ident,lower=lower,maximum = FALSE,mdecomp=mdecomp,
                          model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,
                          upper=upper,radius=radius,setup=setup,X=X,ns=ns,NS=NS,MM=MM)
        }
  if(length(param)>1)
        {
  #### vecchia  approxxx
    if(!is.null(neighb)) { 
            locs=cbind(coordx,coordy)

            #tt2 <- proc.time()
            vecchia.approx=GPvecchia::vecchia_specify(locs,m=neighb,ordering="maxmin")#,conditioning="NN",cond.yz="SGV")
          #tt2<- proc.time()-tt2;print(tt2[3])

          if(optimizer=="nlminb"){
             tt3 <- proc.time()
             Likelihood <- nlminb(objective=eval(as.name(lname)),start=param,vecchia.approx=vecchia.approx,
                             control = list( iter.max=100000),dimat=dimat,
                         lower=lower,upper=upper, hessian=hessian,
                           data=t(data),fixed=fixed,
                          model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,X=X,MM=MM)
             #tt3<- proc.time()-tt3;print(tt3[3]/as.numeric(Likelihood$iterations)) 
         }
       
              if(optimizer=="nmkb")
                  Likelihood <- dfoptim::nmkb(par=param, fn=eval(as.name(lname)), control = list(maxfeval=100000,tol=1e-10),
                        lower=lower,upper=upper,  vecchia.approx=vecchia.approx,
                        dimat=dimat,data=t(data),fixed=fixed,
                          model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,X=X,MM=MM)
                        
            if(optimizer=="Nelder-Mead")
                  Likelihood <- optim(param,eval(as.name(lname)),vecchia.approx=vecchia.approx,
                             control=list(reltol=1e-14, maxit=maxit),dimat=dimat,hessian=hessian, data=t(data),fixed=fixed,
                             model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,X=X,MM=MM)

             #Likelihood <- optim(param,eval(as.name(lname)),const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,
              #            corrmodel=corrmodel,control=list(
               #              reltol=1e-14, maxit=maxit),data=t(data),dimat=dimat,
                #         fixed=fixed,method="Nelder-Mead",
                 #         model=model,namescorr=namescorr,hessian=hessian,
                  #        namesnuis=namesnuis,namesparam=namesparam,X=X)

                       }
        else{  ### no vecchia
if(optimizer=='L-BFGS-B'&&!parallel)
                        Likelihood <- optim(param,fn=eval(as.name(lname)),const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,control=list(
                          pgtol=1e-14,maxit=maxit),data=t(data),dimat=dimat,fixed=fixed,
                          fname=fname,grid=grid,ident=ident,lower=lower,mdecomp=mdecomp,method=optimizer,
                          model=model,namescorr=namescorr,hessian=hessian,
                          namesnuis=namesnuis,upper=upper,namesparam=namesparam,radius=radius,setup=setup,X=X,ns=ns,NS=NS,MM=MM) 

   if(optimizer=='L-BFGS-B'&&parallel){

            #ncores=max(1, parallel::detectCores() - 1)
        ncores=length(param) * 2 + 1
        if(Sys.info()[['sysname']]=="Windows") cl <- parallel::makeCluster(ncores,type = "PSOCK")
        else                                   cl <- parallel::makeCluster(ncores,type = "FORK")
        parallel::setDefaultCluster(cl = cl)
                          Likelihood <- optimParallel::optimParallel(param,fn=eval(as.name(lname)),const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,control=list(pgtol=1e-14, maxit=100000,factr = 1e8),
                          data=t(data),dimat=dimat,fixed=fixed,
                          fname=fname,grid=grid,ident=ident,lower=lower,mdecomp=mdecomp,method=optimizer,
                          model=model,namescorr=namescorr,hessian=hessian,  parallel = list(forward = FALSE),
                          namesnuis=namesnuis,upper=upper,namesparam=namesparam,radius=radius,setup=setup,X=X,ns=ns,NS=NS,MM=MM)
     parallel::setDefaultCluster(cl=NULL)
     parallel::stopCluster(cl)
  }
  if(optimizer=='BFGS')
     Likelihood <- optim(param,eval(as.name(lname)),const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,control=list(
                        pgtol=1e-14,maxit=maxit),data=t(data),dimat=dimat,
                         fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,method=optimizer,
                          model=model,namescorr=namescorr,hessian=hessian,
                          namesnuis=namesnuis,namesparam=namesparam,radius=radius,setup=setup,X=X,ns=ns,NS=NS,MM=MM)
  if(optimizer=='Nelder-Mead')
                   Likelihood <- optim(param,eval(as.name(lname)),const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,control=list(
                             reltol=1e-14, maxit=maxit),data=t(data),dimat=dimat,
                         fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,method=optimizer,
                          model=model,namescorr=namescorr,hessian=hessian,
                          namesnuis=namesnuis,namesparam=namesparam,radius=radius,setup=setup,X=X,ns=ns,NS=NS,MM=MM)




  if(optimizer=='nlm')
                      Likelihood <- nlm(eval(as.name(lname)),param,const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,hessian=hessian,
                          model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,radius=radius,setup=setup,iterlim = maxit,X=X,ns=ns,NS=NS,MM=MM)
  if(optimizer=='nlminb')
  { #print(lname);print(fname);print(mdecomp)
                       Likelihood <-nlminb(objective=eval(as.name(lname)),start=param,
                             control = list( iter.max=100000),
                         lower=lower,upper=upper, hessian=hessian,
                          const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,
                          corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
                          model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,radius=radius,
                          setup=setup,X=X,ns=ns,NS=NS,MM=MM)
                   }
  # if(optimizer=='multinlminb'){
   #                    Likelihood <-mcGlobaloptim::multiStartoptim(objectivefn=eval(as.name(lname)),
    #                      const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,
     #                     corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
      #                    model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,radius=radius,
       #                   setup=setup,X=X,ns=ns,NS=NS,MM=MM,
        #                     lower=lower,upper=upper,method = "nlminb", nbtrials = 500, 
         #                     control = list( iter.max=100000),#
          #                 typerunif = "sobol")#,nbclusters=2,
           #        }
    # if(optimizer=='multiNelder-Mead'){
     #                  Likelihood <-mcGlobaloptim::multiStartoptim(objectivefn=eval(as.name(lname)),
      #                    const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,
       #                   corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
        #                  model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,radius=radius,
         #                 setup=setup,X=X,ns=ns,NS=NS,MM=MM,
          #                   lower=lower,upper=upper,method = "Nelder-Mead", nbtrials = 500, 
           #                   control = list( iter.max=100000),
            #               typerunif = "sobol")#,nbclusters=2,
             #      }                 
    #if(optimizer=='ucminf')    
     #                   Likelihood <-ucminf::ucminf(par=param, fn=eval(as.name(lname)), hessian=as.numeric(hessian),  
      #                  control=list( maxeval=100000),
       #                 const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,MM=MM,
        #                  corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
         #                 model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,
          #                radius=radius,setup=setup,X=X,ns=ns,NS=NS)
    if(optimizer=='nmk')    
                        Likelihood <- dfoptim::nmk(par=param, fn=eval(as.name(lname)), control = list(maxfeval=100000,tol=1e-10),
                        const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,MM=MM,
                          corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
                          model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,
                          radius=radius,setup=setup,X=X,ns=ns,NS=NS)
   if(optimizer=='nmkb')    
                        Likelihood <- dfoptim::nmkb(par=param, fn=eval(as.name(lname)), control = list(maxfeval=100000,tol=1e-10),
                        lower=lower,upper=upper,
                        const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,MM=MM,
                          corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
                          model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,
                          radius=radius,setup=setup,X=X,ns=ns,NS=NS)
    }
}



 if(optimizer %in% c('Nelder-Mead','L-BFGS-B','BFGS','nmk','nmkb','multiNelder-Mead'))
                   {names(Likelihood$par)=namesparam
                    param <- Likelihood$par
                   maxfun <- -Likelihood$value
                   Likelihood$value <- maxfun
    if(optimizer %in% c('Nelder-Mead','L-BFGS-B','BFGS'))  Likelihood$counts=as.numeric(Likelihood$counts[1])
     if(optimizer %in% c('nmk','nmkb'))                    Likelihood$counts=as.numeric(Likelihood$feval)
               }

    #if(optimizer=='ucminf'){
    #               names(Likelihood$par)=namesparam
    #               param <- Likelihood$par
    #               maxfun <- -Likelihood$value
     #              Likelihood$value <- maxfun
     # }
       if(optimizer %in% c('nlminb','multinlminb')){
                   names(Likelihood$par)=namesparam
                   param <- Likelihood$par
                   maxfun <- -Likelihood$objective
                   Likelihood$value <- maxfun
                   Likelihood$counts=as.numeric(Likelihood$iterations)
      }
    if(optimizer=='nlm')
           { names(Likelihood$estimate)=namesparam
             param <- Likelihood$estimate
             #names(param)<-namesparam
             maxfun <- -as.numeric(Likelihood$minimum)
             Likelihood$value <- maxfun
             Likelihood$param <- param
             Likelihood$counts=as.numeric(Likelihood$iterations)
          }
     if(optimizer=='optimize')  
          {param<-Likelihood$minimum
           names(param)<-namesparam
           maxfun <- -Likelihood$objective
           Likelihood$value <- maxfun
           Likelihood$param<-param}
    numparam<-length(param)# parameter size
    Likelihood$claic <- NULL; Likelihood$clbic <- NULL;
    if(type==3||type==4) 
             {Likelihood$claic <- -2*(maxfun-numparam) 
              Likelihood$clbic <- -2*maxfun+numparam*log(dimat)    }
    ### Some checks of the output from the optimization procedure:
    if(optimizer=='Nelder-Mead' ||  optimizer=='L-BFGS-B'||  optimizer=='BFGS'){
    if(Likelihood$convergence == 0)
      Likelihood$convergence <- 'Successful'
    else
      if(Likelihood$convergence == 1)
        Likelihood$convergence <- 'Iteration limit reached'
      else
        Likelihood$convergence <- 'Optimization may have failed'}
    if(optimizer=='optimize'){  Likelihood$convergence <- 'Successful'}
    if(optimizer=='nmk' || optimizer=='nmkb'){
                   if(Likelihood$convergence == 0) Likelihood$convergence <- 'Successful'
                   else Likelihood$convergence <- 'Optimization may have failed'}
    if(optimizer=='nlm'){
               if(Likelihood$code == 1||Likelihood$code == 2)
               Likelihood$convergence <- 'Successful'
               else
               if(Likelihood$code == 4)
               Likelihood$convergence <- 'Iteration limit reached'
               else
               Likelihood$convergence <- 'Optimization may have failed'}
        if(optimizer=='nlminb'||optimizer=='multinlminb'){
               if(Likelihood$convergence == 0)
               Likelihood$convergence <- 'Successful'
               else
               Likelihood$convergence <- 'Optimization may have failed'}
      # if(optimizer=='ucminf'){
       #        if(Likelihood$convergence== 1||Likelihood$convergence== 2||Likelihood$convergence == 4)
        #       Likelihood$convergence <- 'Successful'
         #      else
          #     Likelihood$convergence <- 'Optimization may have failed'}
    if(maxfun==-1.0e8) Likelihood$convergence <- 'Optimization may have failed: Try with other starting parameters'
     }
     else {Likelihood=as.list(0)   # in the case of computing just the asym variance
           names(Likelihood)="value"
           maxfun=Likelihood$value
           Likelihood$param <- param
           Likelihood$claic <- NULL;Likelihood$clbic <- NULL;
           Likelihood$convergence <- 'None'
           numparam<-length(param)
            }
   


if(varest) 
  {
 #    aa=try(abs(det(Likelihood$hessian)),silent=T)
 #  if(aa<1e-08||is.null(Likelihood$hessian)||min(eigen(Likelihood$hessian)$values)<0)
 # {  

Likelihood$hessian=numDeriv::hessian(func=eval(as.name(lname)),x=Likelihood$par,method="Richardson",  const=const,coordx=coordx,coordy=coordy,
            coordt=coordt,corr=corr,corrmat=corrmat,
            corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
            model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,radius=radius,setup=setup,X=X,ns=ns,NS=NS,MM=MM)

Likelihood$score=numDeriv::grad(func=eval(as.name(lname)),x=Likelihood$par,method="Richardson",  const=const,coordx=coordx,coordy=coordy,
            coordt=coordt,corr=corr,corrmat=corrmat,
            corrmodel=corrmodel,data=t(data),dimat=dimat,fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
            model=model,namescorr=namescorr,namesnuis=namesnuis,namesparam=namesparam,radius=radius,setup=setup,X=X,ns=ns,NS=NS,MM=MM)
rownames(Likelihood$hessian)=namesparam
colnames(Likelihood$hessian)=namesparam
names(Likelihood$score)=namesparam
}
#}


   if(Likelihood$convergence == 'Successful' || Likelihood$convergence =='None')
   {   # if optimization has failed it does not compute stderr

    ### START Computing the asymptotic variance-covariance matrices:  
    if(varest){
    if((model==20||model==22||model==1||model==34||model==35||model==37)&& !(type==5||type==6))      {
       # if(is.null(Likelihood$hessian)) {Likelihood$hessian=numDeriv::hessian(func=eval(as.name(lname)),x=param,method="Richardson",
       #     const=const,coordx=coordx,coordy=coordy,coordt=coordt,corr=corr,corrmat=corrmat,
       #                   corrmodel=corrmodel,data=t(data),dimat=dimat,
       #                   fixed=fixed,fname=fname,grid=grid,ident=ident,mdecomp=mdecomp,
       #                   model=model,namescorr=namescorr,
       #                   namesnuis=namesnuis,namesparam=namesparam,radius=radius,setup=setup,X=X,ns=ns,NS=NS)
       # rownames(Likelihood$hessian)=namesparam; colnames(Likelihood$hessian)=namesparam
        #             }
         aa=try(abs(det(Likelihood$hessian)),silent=T)
         if(aa<1e-08) { 
                        warning("Asymptotic information matrix is singular",immediate.=TRUE)
                        Likelihood$varcov <- NULL  } 
        else 
        {
        
         ii=solve(Likelihood$hessian)
         mm=min(eigen(ii)$values)
            if(mm<=0)   { 
                        warning("Asymptotic information matrix is not positive-definite")
                        Likelihood$varcov <- NULL  } 
            else       {Likelihood$varcov <-  ii
                        Likelihood$stderr <- sqrt(diag(( Likelihood$varcov))) } 
         }
     }
     else
     {
        gnames <- namesparam
        if(!bivariate){
        if(flagnuis[1]) {# cheks if the mean is included
            numparam<-numparam-num_betas
            gnames <- gnames[gnames!="mean"]
         if(num_betas>1) {for(i in 1:(num_betas-1)) {if(flagnuis[i+1]) gnames <- gnames[gnames!=paste("mean",i,sep="")]}}
            }}
        else  
        {
        if(flagnuis[1]) {# cheks if the mean1 is included
            numparam<-numparam-1
            gnames <- gnames[gnames!="mean_1"]}
         if(flagnuis[2]) {# cheks if the mean2 is included
            numparam<-numparam-1
            gnames <- gnames[gnames!="mean_2"]}
         }
        param<-c(param,fixed)# update with the fixed
        paramcorr<-param[namescorr]# correlation components
        numparamcorr<-length(paramcorr)# correlation size
        nuisance<-param[namesnuis]# nuisance components
        eps<-(.Machine$double.eps)^(1/3)
        numfish<-numparam*(numparam-1)/2+numparam# set variance-covariance matrix size
        
        # computing the off diagonal elements of the correlation matrix
        corr<-matr(corrmat,corr,coordx,coordy,coordt,corrmodel,nuisance,paramcorr,ns,NS,radius)
        if(!bivariate) varian<-corr*nuisance['sill']# computes covariance components
        else           varian<-corr
        dname<-"DCorrelationMat"
        if(spacetime) dname<-"DCorrelationMat_st"
        if(bivariate) dname<-"DCorrelationMat_biv"
        numparamcorr<-length(paramcorr[flagcor==1])# set the effective number of the corr param
        namescorr<-namescorr[flagcor==1]# set the effective names of the corr param
        # ML and REML cases:
        if(type==3||type==4){
            # Computing variance-covariance matrix of the random field:
            if(!bivariate)   {varcov<-(nuisance['sill'])*ident;yesdiag=FALSE}  # computes variance components
            else             {varcov=ident;yesdiag=TRUE}
            varcov[lower.tri(varcov,diag=yesdiag)]<-varian
            varcov<-t(varcov)
            varcov[lower.tri(varcov,diag=yesdiag)]<-varian
            decompvarcov <- MatDecomp(varcov,mdecomp)
            if(is.logical(decompvarcov)) {stop("Covariance Matrix is not positive definite")}
            invar <- MatInv(decompvarcov,mdecomp)
            fish<-double(numfish)
            # Restricted likelihood case
            if(type==3) P<-invar-array(rowSums(invar),c(dimat,1))%*%colSums(invar)/sum(invar)
            # set array of gradient matrices
            gradient<-array(0,dim=c(dimat,numparam,dimat))# vector derivatives
            colnames(gradient) <- gnames
            dcorr <- double(numpairstot*numparamcorr)
             dname <- paste(dname,"2",sep="")  
            # correlation gradient vector
            #dc=.C(dname,as.integer(corrmodel),as.double(coordx),as.double(coordy),as#.double(coordt),dr=dcorr,as.double(eps),
               #as.integer(flagcor),as.integer(numparamcorr),as#.double(paramcorr),cr=corr,PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)
               dc <- list()
               dcorr<-dc$dr; corr<-dc$cr
            dim(dcorr) <- c(numpairstot,numparamcorr)
            if(!bivariate){
            # Computing the gradient matrices:
            if(flagnuis[num_betas+1])  # nugget parameter
                 gradient[,namesnuis[num_betas+1],] <- ident
            # gradient matrix, derivatives with respect to the sill (correlation matrix)
            if(flagnuis[num_betas+2]){   # variance parameter
                R<-ident
                R[lower.tri(R)]<-corr
                R<-t(R)
                R[lower.tri(R)]<-corr
                gradient[,namesnuis[num_betas+2],] <- R}}
            # gradient matrix, derivatives with respect to the correlation parameters
            if(numparamcorr) {
            for(i in 1:numparamcorr){
                grad<-matrix(0,dimat,dimat)# set gradient matrix
                grad[lower.tri(grad,diag=yesdiag)]<-dcorr[,i]
                grad<-t(grad)
                grad[lower.tri(grad,diag=yesdiag)]<-dcorr[,i]
                if(!bivariate)  grad<-diag(nuisance['nugget'],dimat)+grad*nuisance['sill']
                if(flagcor[namescorr][i])  gradient[,namescorr[i],] <- grad  }
            }
            i<-1
            # Computing the gradient matrices:
            k<-1
            # Computing off diagonal elements of the asymptotic Fisher/Godambe information matrix:
            for(i in 1:numparam)
                for(j in i:numparam){
                    if(type==3) fish[k]<-0.5*sum(diag(crossprod(P,gradient[,i,])%*%crossprod(P,gradient[,j,])))# REML case
                    if(type==4) fish[k]<-0.5*sum(diag(crossprod(invar,gradient[,i,])%*%crossprod(invar,gradient[,j,])))# ML case
                    k<-k+1}
            # Building Fisher/Godambe information matrix:
            fisher<-diag(0,numparam)# full and restricted likelihood cases
            fisher[lower.tri(fisher,diag=TRUE)] <- fish
            fisher<- t(fisher)
            fisher[lower.tri(fisher,diag=TRUE)] <- fish
            # Adding mean to the asymptotic Fisher/Godambe information matrix:
            if(!bivariate){
            if(flagnuis[1]){          #mean parameter qua!!!!!!!!
                    if(num_betas==1){
                         zeros<-rep(0,numparam)
                         if(type==4) fishmean<-sum(invar)
                         fisher<-rbind(c(fishmean,zeros),cbind(zeros,fisher))
                         }
                    if(num_betas>1){  zeros=matrix(0,ncol=numparam,nrow=num_betas)
                         if(type==4) fishmean<-t(X)%*%invar%*%X
                         fisher<-rbind(cbind(fishmean,zeros),cbind(t(zeros),fisher))
                          }
             }
         }
            else {
                  zeros<-rep(0,numparam)
                  if(flagnuis[1]){  #mean_1 parameter
                                  #unozeros<-rep(c(1,0),numcoord*numtime);
                                  unozeros<-c(rep(1,ns[1]),rep(0,ns[2]));
                                  fishmean1=t(unozeros)%*%invar%*%unozeros;
                                  fisher<-rbind(c(fishmean1,zeros),cbind(zeros,fisher));zeros=c(0,zeros)
                                 }
                  if(flagnuis[2]){   #mean_2 parameter
                                  #zerounos<-rep(c(0,1),numcoord*numtime);
                                  zerounos<-c(rep(0,ns[1]),rep(1,ns[2]));
                                  fishmean2=t(zerounos)%*%invar%*%zerounos;
                                  fisher<-rbind(c(fishmean2,zeros),cbind(zeros,fisher)) }
                  }
            cholfisher<- try(chol(fisher),silent=T)
            invfisher <- try(chol2inv(cholfisher),silent=TRUE)
            if(!is.matrix(invfisher)) invfisher<-NULL
            Likelihood$sensmat <- NULL
            Likelihood$varimat <- NULL
        }
     
        # Computing of the asymptotic variance covariance matrix: case TAPERING
        if(type==5||type==6){
            # define the variance-covariance vector
            if(!bivariate) {sel<- varian==(nuisance['sill'])
                            varian[sel] <- nuisance['sill']}
            # define the sparse variance-covariance matrix
            spamvar <- new("spam",entries=varian,colindices=setup$ja,rowpointers=setup$ia,
                           dimension=as.integer(rep(dimat,2)))
            # define tha variance-covariance  matrix
            varcov <- as.matrix(spamvar)  #  variance-covariance tapered matrix as matrix
            covtap<-spamvar
            slot(covtap,"entries") <- varian*setup$taps
            cholcovtap  <- try(spam::update.spam.chol.NgPeyton(setup$struct,covtap),silent=TRUE)
            invtap <- as.matrix(spam::solve.spam(cholcovtap))
            covtap <- as.matrix(covtap)   #  variance-covariance tapered matrix as matrix
            HH <- double(numfish)# upper triang matrix (with diag) of sensitivity matrix
            JJ <- double(numfish)# upper triang matrix (with diag) of variability matrix
            # computing the matrix of derivatives
            gradient <- array(0,c(numpairs,numparam))# vector derivatives
            colnames(gradient) <- gnames
            dname<-"DCorrelationMat_tap"
            if(spacetime) dname<-"DCorrelationMat_st_tap"
            if(bivariate) dname<-"DCorrelationMat_biv_tap"
            dcorr <- double(numpairs*numparamcorr)
            # correlation gradient vector
            #gr <- #.C(dname,as.integer(corrmodel),as.double(coordx),as.double(coordy),as#.double(coordt),dr=dcorr,as.double(eps),
               #as.integer(flagcor),as.integer(numparamcorr),as.double(paramcorr),cr=corr,PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)
               gr <- list()
            dcorr <- gr$dr
            corr <- gr$cr
            dim(dcorr) <- c(numpairs,numparamcorr)
            if(!bivariate) {if(flagnuis[num_betas+1]) gradient[,namesnuis[num_betas+1]] <- ident[setup$idx]  # nugget parameter
                            if(flagnuis[num_betas+2]) gradient[,namesnuis[num_betas+2]] <- corr              # variance parameter
                            }
            if(!bivariate){ dcorr[sel]<-dcorr[sel]*(nuisance['sill'])
                            dcorr[-sel]<-dcorr[-sel]*nuisance['sill']
                          }
             gradient[,namescorr] <- dcorr
            # computing matrix derivatives
            k <- 1
            # computing uppper triangular of the Fisher matrix
            H<-diag(0,numparam)# sensitivity matrix
            J<-diag(0,numparam)# variability matrix
            gradtapI <- gradtapJ <- bigI <- bigJ <- spamvar
            for(i in 1:numparam){
                for(j in i:numparam){
                    slot(gradtapI,"entries") <- gradient[,i]*setup$taps
                    slot(gradtapJ,"entries") <- gradient[,j]*setup$taps
                    aI=tcrossprod(as.matrix(gradtapI),invtap);aJ=tcrossprod(as.matrix(gradtapJ),invtap)
                    HH[k] <- -0.5*sum(diag(aI%*%aJ))
                    slot(bigI, "entries") <- (invtap%*%aI)[setup$idx]*setup$taps
                    slot(bigJ, "entries") <- (invtap%*%aJ)[setup$idx]*setup$taps
                    JJ[k] <- 0.5*sum(diag(tcrossprod(as.matrix(bigI),varcov)%*%tcrossprod(as.matrix(bigJ),varcov)))
                    k <- k+1}}
            # Building Godambe information  matrix
            H[lower.tri(H,diag=TRUE)] <- HH
            H<- t(H)
            H[lower.tri(H,diag=TRUE)] <- HH
            J[lower.tri(J,diag=TRUE)] <- JJ
            J<- t(J)
            J[lower.tri(J,diag=TRUE)] <- JJ
            # Adding the mean parameter
            if(!bivariate)  {
          if(flagnuis[1]){  ##mean parameter
                slot(spamvar,"entries") <- invtap[setup$idx]*setup$taps
                if(num_betas==1){
                   zeros <- rep(0,numparam)
                   fishmH <- -sum(spamvar)
                   H <- rbind(c(fishmH,zeros),cbind(zeros,H))
                   fishmJ <- sum(spamvar%*%varcov%*%spamvar)
                   J <- rbind(c(fishmJ,zeros),cbind(zeros,J))
                     }
                if(num_betas>1){

                   zeros=matrix(0,ncol=numparam,nrow=num_betas)
                   fishmH=-t(X)%*%spamvar%*%X
                   H<-rbind(cbind(fishmH,zeros),cbind(t(zeros),H))
                   fishmJ <-t(X)%*%(spamvar%*%varcov%*%spamvar)%*%X
                   J<-rbind(cbind(fishmJ,zeros),cbind(t(zeros),J))
             }}}
             else {
              if(flagnuis[1]){   ##mean_1 parameter
               slot(spamvar,"entries") <- invtap[setup$idx]*setup$taps
                  zeros <- rep(0,numparam)
                  unozeros<-rep(c(1,0),numcoord*numtime);
                  fishmH <- -t(unozeros)%*%spamvar%*%unozeros;
                  H <- rbind(c(fishmH,zeros),cbind(zeros,H))
                  fishmJ <- t(unozeros)%*%(spamvar%*%varcov%*%spamvar)%*%unozeros;
                  J <- rbind(c(fishmJ,zeros),cbind(zeros,J))
                  zeros=c(0,zeros)
              }
              if(flagnuis[2]){ ##mean_2 parameter
                  slot(spamvar,"entries") <- invtap[setup$idx]*setup$taps
                  zerosuno<-rep(c(0,1),numcoord*numtime);
                  fishmH <- -t(zerosuno)%*%spamvar%*%zerosuno;
                  H <- rbind(c(fishmH,zeros),cbind(zeros,H))
                  fishmJ <- t(zerosuno)%*%(spamvar%*%varcov%*%spamvar)%*%zerosuno;
                  J <- rbind(c(fishmJ,zeros),cbind(zeros,J))
              }
          }
            cholH <- try(chol(-H),silent = TRUE)
            invH=try(-chol2inv(cholH),silent = TRUE)
              Likelihood$sensmat <- H
              Likelihood$varimat <- J
            if(!is.matrix(invH) || !is.matrix(H)){
                invfisher<-NULL
                Likelihood$claic <- NULL;Likelihood$clbic <- NULL; }
            else{
                prJH=crossprod(J,invH)
                invfisher<-crossprod(invH,prJH)   # godambe matrix
                # invfisher<-invH%*%J%*%invH   # godambe matrix
                Likelihood$claic <- -2*(maxfun-sum(diag(prJH)))# penalty in the TIC case
                Likelihood$clbic <- -2*(maxfun)+log(dimat)*sum(diag(prJH))
                }
        }
        ### END Computing the asymptotic variance-covariance matrices
        Likelihood$varcov <- invfisher
        #Checks if the resulting variance and covariance matrix:
        mm=min(eigen(Likelihood$varcov)$values)
        if(is.null(Likelihood$varcov)||mm<0){
            if(mm<0)                       warning("Asymptotic information matrix is not positive-definite")
            if(is.null(Likelihood$varcov)) warning("Asymptotic information matrix is singular")
            Likelihood$varcov <- 'none'
            Likelihood$stderr <- 'none'}
        else
        {
            dimnames(Likelihood$varcov)<-list(namesparam,namesparam)
            Likelihood$stderr<-diag(Likelihood$varcov)
        if(any(Likelihood$stderr < 0)) Likelihood$stderr <- 'none'
        else{
            Likelihood$stderr<-sqrt(Likelihood$stderr)
            names(Likelihood$stderr)<-namesparam}}
    }}}
    ### END the main code of the function:
if(varest)
    if(is.null(Likelihood$varcov)){
                Likelihood$varcov <- 'none';Likelihood$stderr <- 'none'}

    return(Likelihood)
  }
