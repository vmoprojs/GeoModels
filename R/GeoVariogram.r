####################################################
### File name: GeoVariogram.r
####################################################

### Procedures are in alphabetical order.

GeoVariogram <- function(data, coordx, coordy=NULL, coordt=NULL, coordx_dyn=NULL,cloud=FALSE, distance="Eucl",
                       grid=FALSE, maxdist=NULL,neighb=NULL, maxtime=NULL, memdist=FALSE,numbins=NULL,
                       radius=6371, type='variogram',bivariate=FALSE)
  {
    call <- match.call()
    corrmodel <- 'exponential'
    ### Check the parameters given in input:
    if(is.null(type))
      type <- 'variogram'
    # Set the type of model:
    if(type=='variogram'){
        model <- 'Gaussian'
        fname <- 'Binned_Variogram'}

    
    # Checks if its a spatial or spatial-temporal random field:
    if(bivariate) coordt=c(0,1)
    if(!is.null(coordt))
      if(is.numeric(coordt))
        if(length(coordt)>1) corrmodel <- 'gneiting'
    # Checks the input:
    checkinput <- CkInput(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, "Fitting", NULL, grid,
                             'None', maxdist, maxtime, model,NULL, 'Nelder-Mead', NULL,
                             radius,  NULL, NULL,NULL, 'GeoWLS', FALSE, 'SubSamp', FALSE,NULL,NULL)
                             

    # Checks if there are errors in the input:
    if(!is.null(checkinput$error))
      stop(checkinput$error)
    ### START -- Specific checks of the Empirical Variogram:
    if(!is.null(cloud) & !is.logical(cloud))
      stop('insert a logical value (TRUE/FALSE) for the cloud parameter\n')
    
    if(!is.null(numbins) & !is.integer(numbins))
      if(numbins < 0)
        stop('insert a positive integer value for the number of bins\n')
    if(is.null(numbins))
      numbins <- 13
    ### END -- Specific checks of the Empirical Variogram
    if(bivariate)  corrmodel <- 'Bi_matern_sep'
  

    n=1
    initparam <- StartParam(coordx, coordy, coordt,coordx_dyn, corrmodel, data,distance, "Fitting",
                           NULL, grid, 'None', maxdist,neighb,
                           maxtime, model, n, NULL, NULL, FALSE, radius, 
                           NULL, NULL, NULL, 'GeoWLS', 'GeoWLS', FALSE,
                           'SubSamp', FALSE, 1, 1,1,1,NULL,NULL,FALSE)
    spacetime_dyn=NULL
    coordx=initparam$coordx;coordy=initparam$coordy;coordt=initparam$coordt                 
    # Checks if there are inconsistences:
    if(!is.null(initparam$error))
      stop(initparam$error)
    numvario <- numbins-1
    if(cloud){
        numbins <- numvario <- initparam$numpairs
        fname <- 'Cloud_Variogram'}
    ### Estimation of the empirical spatial or spatial-temporal variogram:
    bins <- double(numbins) # spatial bins
    moments <- double(numvario) # vector of spatial moments
    lenbins <- integer(numvario) # vector of spatial bin sizes
    bint <- NULL
    lenbinst <- NULL
    lenbint <- NULL
    variogramst <- NULL
    variogramt <- NULL  
  #***********************************************************************************************#
  #***********************************************************************************************#
  #***********************************************************************************************#
    if(initparam$bivariate){
               n_var=initparam$numtime
               spacetime_dyn=FALSE
               if(!is.null(coordx_dyn)) spacetime_dyn=TRUE
               ns=initparam$ns
               NS=cumsum(ns)
               if(!spacetime_dyn){
                                  data=c(t(data))
                                  coordx=rep(coordx,n_var)
                                  coordy=rep(coordy,n_var)
                         }
               if(spacetime_dyn) data=unlist(data)
               NS=c(0,NS)[-(length(ns)+1)]
               moments_marg<-double(n_var*numvario)   # vect of square differences for each component (n_var) 11  e 22
               lenbins_marg<-integer(n_var*numvario)  #
               moments_cross<-double(0.5*n_var*(n_var-1)*numvario)  # vect of square differences for cross components (12)
               lenbins_cross<-integer(0.5*n_var*(n_var-1)*numvario) #

               DEV=.C("Binned_Variogram_biv2", bins=bins, as.double(coordx),as.double(coordy),as.double(coordt),as.double(data),
               lenbins_cross=lenbins_cross, moments_cross=moments_cross, as.integer(numbins),lenbins_marg=lenbins_marg,
               moments_marg=moments_marg,as.integer(ns),as.integer(NS),
               PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
               
              #DEV=dotCall64::.C64("Binned_Variogram_biv2", 
               #      SIGNATURE = c("double","double","double","double","double",
                #                "integer","double","integer","integer",
                    #            "double","integer","integer"),
              #  bins=bins, coordx,coordy,coordt,data,
              # lenbins_cross=lenbins_cross, moments_cross=moments_cross, numbins,lenbins_marg=lenbins_marg,
              # moments_marg=moments_marg,ns,NS,
                # INTENT = c("w","r","r","r","r","w","w","r","w","w","r","r"),NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0)

               bins=DEV$bins
               lenbins_cross=DEV$lenbins_cross
               moments_cross=DEV$moments_cross
               lenbins_marg=DEV$lenbins_marg
               moments_marg=DEV$moments_marg
               m_11=moments_marg[1:numvario];m_22=moments_marg[(numvario+1):(2*numvario)];m_12=moments_cross[1:numvario];
               l_11=lenbins_marg[1:numvario];l_22=lenbins_marg[(numvario+1):(2*numvario)];l_12=lenbins_cross[1:numvario];
               indbin_marg <- l_11>0;indbin_cross <- l_12>0  #
               bins<- bins[indbin_marg]
               numbins <-length(bins)
               m_11 <- m_11[indbin_marg];m_22 <- m_22[indbin_marg];l_11 <- l_11[indbin_marg];l_22 <- l_22[indbin_marg]
               m_12 <- m_12[indbin_cross];l_12 <- l_12[indbin_cross]
               variograms_11 <- m_11/l_11;variograms_22 <- m_22/l_22
               variograms_12 <- m_12/l_12   
               centers <-   bins[1:(numbins[1]-1)]+diff(bins)/2
               lenbins=rbind(l_11,l_22);lenbinst=l_12
               variograms=rbind(variograms_11,variograms_22) 
               variogramst=variograms_12
    }
  #***********************************************************************************************#
  #***********************************************************************************************#
  #***********************************************************************************************#
    if(initparam$spacetime){
      numtime=initparam$numtime
      spacetime_dyn=FALSE
      if(!is.null(coordx_dyn)) spacetime_dyn=TRUE
      ns=initparam$ns
      NS=cumsum(ns)
      numbint <- initparam$numtime-1 # number of temporal bins
      bint <- double(numbint)        # temporal bins
      momentt <- double(numbint)     # vector of temporal moments
      lenbint <- integer(numbint)    # vector of temporal bin sizes
      numbinst <- numvario*numbint   # number of spatial-temporal bins
      binst <- double(numbinst)      # spatial-temporal bins
      momentst <- double(numbinst)   # vector of spatial-temporal moments
      lenbinst <- integer(numbinst)  # vector of spatial-temporal bin sizes
      #if(cloud) fname <- 'Cloud_Variogram_st' else 
      fname <- 'Binned_Variogram_st'
      if(grid)     {a=expand.grid(coordx,coordy);coordx=a[,1];coordy=a[,2]; }
      else{
      if(!spacetime_dyn){
                                  data=c(t(data))
                                  coordx=rep(coordx,numtime)
                                  coordy=rep(coordy,numtime)
                         }
      if(spacetime_dyn) data=unlist(data)
      }
         NS=c(0,NS)[-(length(ns)+1)]

      if(spacetime_dyn) fname <- paste(fname,"2_dyn",sep="") 
      if(!spacetime_dyn) fname <- paste(fname,"2",sep="") 
      # Compute the spatial-temporal moments:
if(spacetime_dyn)
      EV=.C("Binned_Variogram_st2_dyn", bins=bins, bint=bint,  as.double(coordx),as.double(coordy),as.double(coordt),as.double((data)),
           lenbins=lenbins,lenbinst=lenbinst,lenbint=lenbint,moments=moments,momentst=momentst,momentt=momentt,
           as.integer(numbins), as.integer(numbint),as.integer(ns),as.integer(NS), PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
if(!spacetime_dyn)
      EV=.C("Binned_Variogram_st2", bins=bins, bint=bint,  as.double(coordx),as.double(coordy),as.double(coordt),as.double((data)),
           lenbins=lenbins,lenbinst=lenbinst,lenbint=lenbint,moments=moments,momentst=momentst,momentt=momentt,
           as.integer(numbins), as.integer(numbint),as.integer(ns),as.integer(NS), PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)

    # EV=dotCall64::.C64(fname, 
             #        SIGNATURE = c("double","double","double","double","double","double",
            #                    "integer","integer","integer","double","double","double",
           #                     "integer","integer","integer","integer"),
          #     bins=bins, bint=bint,coordx,coordy,coordt,data,
         #  lenbins=lenbins,lenbinst=lenbinst,lenbint=lenbint,moments=moments,momentst=momentst,momentt=momentt,
        #   numbins,numbint,ns,NS,
       #          INTENT = c("w","w","r","r","r","r","w","w","w","w","w","w", "r","r","r","r"), 
      #             NAOK = TRUE, PACKAGE = "GeoModels", VERBOSE = 0)
       bins=EV$bins
       bint=EV$bint
       lenbins=EV$lenbins
       lenbint=EV$lenbint
       lenbinst=EV$lenbinst
       moments=EV$moments
       momentt=EV$momentt
       momentst=EV$momentst
       centers <- bins[1:numvario]+diff(bins)/2
     
      indbin <- lenbins>0
      indbint <- lenbint>0
      indbinst <- lenbinst>0
      bins <- bins[indbin]; bint <- bint[indbint]
      centers <- centers[indbin]; 
      moments <- moments[indbin]; lenbins <- lenbins[indbin]
      momentt <- momentt[indbint]; lenbint <- lenbint[indbint]; 
      momentst <- momentst[indbinst]; lenbinst <- lenbinst[indbinst]
      variograms <- moments/lenbins; variogramt <- momentt/lenbint; variogramst <- momentst/lenbinst
    }
  #***********************************************************************************************#
  #***********************************************************************************************#
  #***********************************************************************************************#
    if(!initparam$bivariate&&!initparam$spacetime){  ## spatial univariate case
    if(grid)     {a=expand.grid(coordx,coordy);coordx=a[,1];coordy=a[,2]; }

     if(!memdist)  { 
     fname <- paste(fname,"2",sep="") 
     # Computes the spatial moments
      EV=.C("Binned_Variogram2", bins=bins,  as.double(coordx),as.double(coordy),as.double(coordt),as.double(data), 
        lenbins=lenbins, moments=moments, as.integer(numbins),PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
       }
      else {
       fname="Binned_Variogram2new"
  
         idx=GeoNeighIndex(cbind(coordx,coordy),distance = distance, neighb = neighb, maxdist = maxdist,radius=radius)
         #mm=c(min(idx$lags),max(idx$lags))
         mm=range(idx$lags)
        
         EV=.C("Binned_Variogram2new", bins=bins,  as.integer(length(idx$lags)),as.double(data[idx$colidx]),
                as.double(data[idx$rowidx]), as.double(idx$lags),
        lenbins=lenbins, moments=moments, as.integer(numbins),as.double(mm),PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
       }

   
       bins=EV$bins
       lenbins=EV$lenbins
       moments=EV$moments
      # Computes the spatial variogram:
      indbin <- lenbins>0
      bins <- bins[indbin]
      numbins <- length(bins)
      # check if cloud or binned variogram:
      if(cloud) centers <- bins else centers <- bins[1:(numbins-1)]+diff(bins)/2
      moments <- moments[indbin]
      lenbins <- lenbins[indbin]
      variograms <- moments/lenbins}
    # Start --- compute the extremal coefficient
    if(!memdist).C('DeleteGlobalVar', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
    else .C('DeleteGlobalVar2', PACKAGE='GeoModels', DUP = TRUE, NAOK=TRUE)
 
    GeoVariogram <- list(bins=bins,
                       bint=bint,
                       bivariate=bivariate,
                       cloud=cloud,
                       centers=centers,
                       lenbins=lenbins,
                       lenbinst=lenbinst,
                       lenbint=lenbint,
                       maxdist =maxdist,
                       maxtime = maxtime,
                       spacetime_dyn=spacetime_dyn,
                       variograms=variograms,
                       variogramst=variogramst,
                       variogramt=variogramt,
                       type=type)

    structure(c(GeoVariogram, call = call), class = c("GeoVariogram"))

  }

#########################################################################################################
#########################################################################################################
#########################################################################################################

plot.GeoVariogram <- function(x,...)
  {

if(!inherits(x,"GeoVariogram"))       stop("Enter an object obtained from the function GeoVariogram\n")

opar=par(no.readonly = TRUE)
on.exit(par(opar))

ispatim=bivariate=FALSE 
if(!is.null(x$bint))  ispatim=TRUE
if(x$bivariate)       bivariate=TRUE
  
    # lags associated to empirical variogram estimation
    lags = c(0,x$centers);numlags = length(lags)
    if(ispatim) lagt =c(0,x$bint) else lagt=0
    numlagt = length(lagt)
#########################################################################
    slow=0
    lags_m = seq(slow,max(x$centers),length.out =150)
    if (ispatim) lagt_m =seq(slow,max(x$bint),length.out =150)
    else         lagt_m=0
    numlags_m = length(lags_m)
    numlagt_m = length(lagt_m)
#########################################################################


##########################################
      vario.main = "Spatial semi-variogram"
      vario.ylab = "Semi-Variogram"
        if(ispatim){
            vario.main = "Space-time semi-variogram"
            vario.zlab = "Semi-Variogram"
     }

################################### 
#### bivariate case ###############
###################################
      if(bivariate){
        par(mfrow=c(2,2))
       plot.default(x$centers,x$variograms[1,], main="First semi-variogram",ylim=c(0,max(x$variograms[1,])),
           xlim=c(0,max(x$centers)),
                     xlab="Distance", ylab="Semi-Variogram",...)
       if(min(x$variogramst)>0) {ll1=0;ll=max(x$variogramst)}
       if(min(x$variogramst)<0) {ll1=min(x$variogramst);ll=-min(x$variogramst)}
       plot.default(x$centers,x$variogramst, main="Cross semi-variogram",ylim=c(ll1,ll),
         xlim=c(0,max(x$centers)),
                     xlab="Distance", ylab="Semi-Variogram",...)
       plot.default(x$centers,x$variogramst, main="Cross semivariogram",ylim=c(ll1,ll),
         xlim=c(0,max(x$centers)),
                     xlab="Distance", ylab="Semi-Variogram",...)
       plot.default(x$centers,x$variograms[2,], main="Second semi-variogram",ylim=c(0,max(x$variograms[2,])),
         xlim=c(0,max(x$centers)),
                     xlab="Distance", ylab="Semi-Variogram",...)
      }
################################### 
#### space time case ##############
###################################
    if(ispatim){
        par(mfrow=c(2,2), mai=c(.5,.5,.3,.3), mgp=c(1.4,.5, 0))

plot.default(x$centers, x$variograms, xlab='h', ylab=expression(gamma(h)),
     ylim=c(0, max(x$variograms)), xlim=c(0, max(x$centers)),
     main="Marginal spatial semi-variogram")

plot.default(x$bint, x$variogramt, xlab='t', ylab=expression(gamma(t)),
     ylim=c(0, max(x$variogramt)),xlim=c(0,max(x$bint)),
     main="Marginal temporal semi-variogram")

         evario = matrix(x$variogramst,nrow=length(x$centers),ncol=length(x$bint),byrow=TRUE)
         evario = rbind(c(0,x$variogramt),cbind(x$variograms,evario))
         evario.grid = as.matrix(expand.grid(c(0,x$centers),c(0,x$bint)))
         scatterplot3d::scatterplot3d(evario.grid[,1],evario.grid[,2], c(evario),
                              type="h",highlight.3d=TRUE,cex.axis=.7,cex.lab=.7,
                              main=paste("Empirical",vario.main),xlab="Distance",
                              ylab="Time",zlab=vario.zlab,mar=c(2,2,2,2),mgp=c(0,0,0))
    par(mai=c(.2,.2,.2,.2),mgp=c(1,.3, 0))
     persp(c(0,x$centers), c(0,x$bint), evario,
      xlab="h", ylab="u", zlab=expression(gamma(h,u)),
      ltheta=90, shade=0.75, ticktype="detailed", phi=30,
      theta=30,main="Space-time semi-variogram",cex.axis=.8,
      cex.lab=.8)
            }
############################spatial case#########################################
        if(!ispatim && !bivariate)    plot.default(x$centers, x$variograms,...)

  return(invisible())
  }

