####################################################
### File name: GeoCv.r
####################################################

GeoCV=function(fit, K=100, estimation=TRUE, 
   optimizer=NULL,lower=NULL, upper=NULL,
    n.fold=0.05, local=FALSE,neighb=NULL,maxdist=NULL,maxtime=NULL,
        sparse=FALSE, type_krig="Simple", which=1,parallel=FALSE,ncores=NULL)
{

if(n.fold>0.99||n.fold<0.01) stop("n.fold must be beween 0.01 and 0.99")
print("Cross-validation  kriging can be time consuming ...")
if(ncol(fit$X)==1) {X=Xloc=NULL;tempX=NULL}
mae=rmse=lscore=crps=mad=brie=NULL
space_dyn=FALSE


if(is.null(optimizer)) {optimizer=fit$optimizer;lower=fit$lower;upper=fit$upper}

if(is.list(fit$data)) space_dyn=TRUE
spacetime<-CheckST(CkCorrModel(fit$corrmodel))
bivariate<-CheckBiv(CkCorrModel(fit$corrmodel))
K=round(K)
if(K<2) stop("K must be grater or equal  to 2")
if(K>1000) stop("K is  too large")
if(bivariate)
   {if(!(which==1||which==2))
          stop("which must be 1 or 2")}
if(local) if(is.null(maxdist)&&is.null(neighb)) stop("maxdist or neighb are required
          for local kriging")
i=1
cat("Starting iteration from 1 to",K," ...\n")

space=!spacetime&&!bivariate

dtp=pred=list()

model1=fit$model
if(fit$missp)  ### misspecification
 {if(fit$model=="StudentT")     model1="Gaussian_misp_StudentT"
  if(fit$model=="Poisson")      model1="Gaussian_misp_Poisson"
  if(fit$model=="PoissonZIP")   model1="Gaussian_misp_PoissonZIP"
  if(fit$model=="SkewStudentT") model1="Gaussian_misp_SkewStudenT"
  if(fit$model=="Tukeygh")      model1="Gaussian_misp_Tukeygh"
 }
########################################################################################################################
########### spatial case ###############################################################################################
########################################################################################################################
if(space)
{ 
N=length(fit$data)
coords=cbind(fit$coordx,fit$coordy)
data=fit$data
Mloc=NULL
if(length(fit$fixed$mean)>1)  tempM=fit$fixed$mean
if(!is.null(X)) tempX=fit$X

coremax=parallel::detectCores()
if(is.na(coremax)||coremax==1) parallel=FALSE

###################################
####### not parallel version 
###################################
if(!parallel){

rmse=crps=mae=rmse=double(K)

progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")
pb <- progressr::progressor(along = 1:K)

######################
while(i<=K){
sel_data = sample(1:N,round(N*(1-n.fold)))  
# data and coord used for prediction
if(!is.null(X)) {X=tempX[sel_data,]; Xloc=tempX[-sel_data,]}
else X=Xloc=NULL
if(length(fit$fixed$mean)>1)  { fit$fixed$mean=tempM[sel_data]; Mloc=tempM[-sel_data]}
# data to predict
data_to_pred  = fit$data[-sel_data]
data_to_est=fit$data[sel_data]
coords_est=coords[sel_data,]
coords_to_pred=coords[-sel_data,]
param=append(fit$param,fit$fixed)
########### estimation 
if(estimation) {
          fit_s= GeoFit(data=data_to_est,coordx=coords_est,corrmodel=fit$corrmodel,X=X,
                            likelihood=fit$likelihood,type=fit$type,grid=fit$grid,
                            copula=fit$copula,anisopars=fit$anisopars,est.aniso=fit$est.aniso,
                            model=model1,radius=fit$radius,n=fit$n,
                            local=fit$local,GPU=fit$GPU,
                           maxdist=fit$maxdist, neighb=fit$neighb,distance=fit$distance,
                            optimizer=optimizer, lower=lower,upper=upper,
                           start=fit$param,fixed=fit$fixed)
if(!is.null(fit$anisopars))   {   fit$param$angle=NULL;fit$param$ratio=NULL; fit$fixed$angle=NULL;fit$fixed$ratio=NULL}
            param=append(fit_s$param,fit_s$fixed)
             }
########### prediction
if(!local) { 
             pr=GeoKrig(data=data_to_est, coordx=coords_est,  
	            corrmodel=fit$corrmodel, distance=fit$distance,grid=fit$grid,loc=coords_to_pred, #ok
	            model=fit$model, n=fit$n, mse=TRUE,#ok
                param=param, anisopars=fit$anisopars, type_krig=type_krig,
                radius=fit$radius, sparse=sparse, X=X,Xloc=Xloc,Mloc=Mloc) #ok
             }
if(local) {
              pr=GeoKrigloc(data=data_to_est, coordx=coords_est,  
              corrmodel=fit$corrmodel, distance=fit$distance,grid=fit$grid,loc=coords_to_pred, #ok
              model=fit$model, n=fit$n, mse=TRUE,#ok
              neighb=neighb,maxdist=maxdist,
              param=param, anisopars=fit$anisopars, type_krig=type_krig,
              radius=fit$radius, sparse=sparse, X=X,Xloc=Xloc,Mloc=Mloc) #ok
              }
#pred[i,]=pr$pred
#dtp[i,]=data_to_pred 
pp=GeoScores(data_to_pred,pred=pr$pred,mse=pr$mse,
    score=c("brie","crps","lscore","pe"))

rmse[i]=  pp$rmse
mae[i]=   pp$mae
mad[i]=   pp$mad
lscore[i]=pp$lscore
brie[i]=     pp$brie
crps[i]=  pp$crps
pb(sprintf("i=%g", i))
i=i+1
} 
################end while
} # end no parallel


######################################################################
#######  parallel version 
######################################################################
if(parallel) {

if(is.null(ncores)){ n.cores <- coremax - 1 }
else
{  if(!is.numeric(ncores)) stop("number of cores not valid\n")
   if(ncores>coremax||ncores<1) stop("number of cores not valid\n")
   n.cores=ncores
}

rmse=crps=mae=rmse=double(K)
pp=round(N*(1-n.fold))
data_to_est=matrix(0,nrow=K,ncol=pp)
data_to_pred=matrix(0,nrow=K,ncol=N-pp)
coords_est=coords_to_pred=list()

#######################
#######################
Mloc=Mest=X=Xloc=list()
cat("Selecting",K,"sub-sample...\n")


for(i in 1:K) {
sel_data = sample(1:N,pp)  
if(!is.null(tempX)) {X[[i]]=tempX[sel_data,]; 
                 Xloc[[i]]=tempX[-sel_data,]}
else {X=Xloc=NULL}                 
#### selecting mu and muloc for a fixed constant mean
if(length(fit$fixed$mean)>1)  { Mloc[[i]]=tempM[-sel_data]
                                Mest[[i]]=tempM[sel_data]   
                               }
else {Mloc=NULL} 
# selecting data and coords to pred
coords_to_pred[[i]]=coords[-sel_data,]
coords_est[[i]]=coords[sel_data,]
data_to_pred[i,]  = c(fit$data[-sel_data])
data_to_est[i,]=c(fit$data[sel_data])

}
#### fixed param organized as a matrix (MEST)
if(length(fit$fixed$mean)>1) 
{ 
    FF=fit$fixed
    for(i in 1:K) {FF$mean=Mest[[i]]; MEST=matrix(rep(c(FF),K),ncol=length(FF),byrow=TRUE);colnames(MEST)=names(FF)    }
}
else  {FF=fit$fixed;MEST=matrix(rep(c(FF),K),ncol=length(FF),byrow=TRUE);colnames(MEST)=names(FF)   }
###############################################################################################
###############################################################################################
if(estimation)  
{
cat("Performing",K,"estimations using",n.cores,"cores...\n")
progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")
pb <- progressr::progressor(along = 1:K)
future::plan(multisession, workers = n.cores)
xx=foreach::foreach(i = 1:K,.combine = rbind,
                           .options.future = list(seed = TRUE)) %dofuture% {
         pb(sprintf("i=%g", i))
          a=GeoFit(data=data_to_est[i,],coordx=coords_est[[i]],corrmodel=fit$corrmodel,X=X[[i]],
                            likelihood=fit$likelihood,type=fit$type,grid=fit$grid,
                            copula=fit$copula,anisopars=fit$anisopars,est.aniso=fit$est.aniso,
                            model=model1,radius=fit$radius,n=fit$n,
                           maxdist=fit$maxdist, neighb=fit$neighb,distance=fit$distance,
                            optimizer=optimizer, lower=lower,upper=upper,
                           start=fit$param,fixed=as.list(MEST[i,]))
          a$data=a$coordx=a$coordy=a$coordt=NULL
          c(a$param,a$fixed);
         }   
    future::plan(sequential)
  
}
else
{ ############### not estimation ###################
 ff=append(fit$param,fit$fixed);xx=matrix(rep(c(ff),K),ncol=length(ff),byrow=TRUE) ;colnames(xx)=names(ff)
}   
######################################################################################################

rownames(xx)=NULL

if(!is.null(fit$anisopars))   {   fit$param$angle=NULL;fit$param$ratio=NULL; fit$fixed$angle=NULL;fit$fixed$ratio=NULL}

###############################################################################################
###############################################################################################

cat("Performing",K,"predictions...\n")
  
############ not local ##############  
if(!local) { 
    progressr::handlers(global = TRUE)
    progressr::handlers("txtprogressbar")
    pb <- progressr::progressor(along = 1:K)

    future::plan(multisession, workers = n.cores)
    YY=foreach::foreach(i = 1:K,.combine = rbind,
                           .options.future = list(seed = TRUE)) %dofuture% {
           pb(sprintf("i=%g", i))
             pr=GeoKrig(data=data_to_est[i,], coordx=coords_est[[i]],  
                corrmodel=fit$corrmodel, distance=fit$distance,grid=fit$grid,loc=coords_to_pred[[i]], #ok
                model=fit$model, n=fit$n, mse=TRUE,#ok
                param=as.list(xx[i,]), anisopars=fit$anisopars, type_krig=type_krig,
                radius=fit$radius, sparse=sparse, X=X[[i]],Xloc=Xloc[[i]],Mloc=Mloc[[i]]) #ok
              pr$data=pr$coordx=pr$coordy=pr$coordt=NULL
             c(pr$pred,pr$mse)
             
    }
            res1= YY[,1:(N-pp)]
            res2= YY[,(N-pp+1):(2*(N-pp))]
         
         future::plan(sequential)    
         }
############  local ###################  
if(local) {
    progressr::handlers(global = TRUE)
    progressr::handlers("txtprogressbar")
    pb <- progressr::progressor(along = 1:K)
    future::plan(multisession, workers = n.cores)
    YY=foreach::foreach(i = 1:K,.combine = rbind,
          .options.future = list(seed = TRUE)) %dofuture% {

    pb(sprintf("i=%g", i))
             pr=GeoKrigloc(data=data_to_est[i,], coordx=coords_est[[i]],  neighb=neighb,maxdist=maxdist,
                corrmodel=fit$corrmodel, distance=fit$distance,grid=fit$grid,loc=coords_to_pred[[i]], #ok
                model=fit$model, n=fit$n, mse=TRUE,#ok
                param=as.list(xx[i,]), anisopars=fit$anisopars, type_krig=type_krig,
                radius=fit$radius, sparse=sparse, X=X[[i]],Xloc=Xloc[[i]],Mloc=Mloc[[i]]) #ok
              pr$data=pr$coordx=pr$coordy=pr$coordt=NULL
             c(pr$pred,pr$mse)
             }
            res1= YY[,1:(N-pp)]
            res2= YY[,(N-pp+1):(2*(N-pp))]
         future::plan(sequential)   
          }

##############
for(i in  1:K){
pp=GeoScores(data_to_pred[i,],pred=res1[i,],mse=res2[i,],score=c("brie","crps","lscore","pe"))
rmse[i]=  pp$rmse; mae[i]=   pp$mae; mad[i]=   pp$mad
lscore[i]=pp$lscore; brie[i]=  pp$brie; crps[i]=  pp$crps
}
}  #end parallel


rm(data_to_est,data_to_pred,coords_est,coords_to_pred)
#list(rmse=rmse,mae=mae,lscore=lscore,crps=crps,brie=brie,mad=mad)

} # end spatial

############################################################
########### spatio temporal case ###########################
############################################################

if(spacetime)
{
coords=cbind(fit$coordx,fit$coordy)
ns=fit$ns
coordt=fit$coordt
T=length(coordt)
NT=sum(ns)
if(is.null(X)) X=rep(1,NT)

NS=cumsum(ns)
NS=c(c(0,NS)[-(length(ns)+1)],NT)

if(is.list(fit$data))
{datos=do.call(c,args=c(fit$data))
X <- do.call(rbind,args=c(X))
}else{datos = fit$data}

if(!space_dyn){
data_tot=cc=NULL
for(k in 1:T) {
cc=rbind(cc,coords)
data_tot=rbind(data_tot,cbind(rep(coordt[k],ns[k]),datos[k,]))
}
data_tot=cbind(cc,data_tot,X)
}

if(space_dyn){
  ct=NULL
  for(k in 1:T) {
    ct=c(ct,rep(coordt[k],ns[k]))
  }
  data_tot=cbind(fit$coordx,fit$coordy,ct,datos,X)
}

coremax=parallel::detectCores()
if(is.na(coremax)||coremax==1) parallel=FALSE

if(is.null(ncores)){ n.cores <- coremax - 1 }
else
{  if(!is.numeric(ncores)) stop("number of cores not valid\n")
  if(ncores>coremax||ncores<1) stop("number of cores not valid\n")
  n.cores=ncores
}

#set.seed(round(seed))
if(!parallel){

progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")
pb <- progressr::progressor(along = 1:K)
rmse=crps=mae=rmse=double(K)
while(i<=K){

####################################### 
sel_data = sample(1:NT,round(NT*(1-n.fold))) 
data_sel= data_tot[sel_data,]
data_to_pred=data_tot[-sel_data,]

data_sel_ord=data_sel[order(data_sel[,3]),]
data_to_pred_ord=data_to_pred[order(data_to_pred[,3]),]

DD=ncol(data_sel_ord)
k=1 ; coordx_dynnew=Xnew=Xnew_loc=datanew=coordx_dynnew_loc=list()


utt=unique(data_sel_ord[,3])
utt_1=unique(data_to_pred_ord[,3])


for(k in 1:length(utt_1) ){
  ll=data_to_pred_ord[data_to_pred_ord[,3]==utt_1[k],]

  if(is.matrix(ll)) coordx_dynnew_loc[[k]]=as.matrix(ll)[,1:2]
  else                                  coordx_dynnew_loc[[k]]=matrix(ll[1:2],nrow=1)
  if(!is.null(X))  {if(is.matrix(ll))  Xnew_loc[[k]]=as.matrix(ll)[,5:DD]
                    else               Xnew_loc[[k]]=matrix(ll[5:DD],nrow=1,byrow=T)
                   }
  }

for(k in 1:length(utt) ){
  ss=data_sel_ord[data_sel_ord[,3]==utt[k],]
  datanew[[k]]=as.vector((ss[,4]))
  coordx_dynnew[[k]]=as.matrix(ss[,1:2])
  if(!is.null(X))  Xnew[[k]]=as.matrix(ss[,5:DD])
}
if(ncol(Xnew[[1]])==1) Xnew=NULL

param=append(fit$param,fit$fixed)
if(estimation) {

           fit_s= GeoFit(data=datanew,coordx_dyn=coordx_dynnew,coordt=utt,
                            corrmodel=fit$corrmodel,X=Xnew,
                            likelihood=fit$likelihood,type=fit$type,grid=fit$grid,
                            copula=fit$copula, anisopars=fit$anisopars,est.aniso=fit$est.aniso,
                            model=model1,radius=fit$radius,n=fit$n,
                            local=fit$local,GPU=fit$GPU,
                            maxdist=fit$maxdist, neighb=fit$neighb,maxtime=fit$maxtime,distance=fit$distance,
                            optimizer=optimizer, lower=lower,upper=upper,
                            start=fit$param,fixed=fit$fixed)
                            #start=as.list(fit$param),fixed=as.list(fit$fixed))
            if(!is.null(fit$anisopars))   
                   {   fit_s$param$angle=NULL;fit_s$param$ratio=NULL; fit_s$fixed$angle=NULL;fit_s$fixed$ratio=NULL}
           param=append(fit_s$param,fit_s$fixed)
              }
#####################################
if(!local) {
           pr_st=pr_mse=list()
           for(j in 1:length(utt_1) ){
            if(is.null(Xnew)) {Xnew_loc[j]=list(NULL)}
               pr=GeoKrig(data=datanew,   coordt=utt, coordx_dyn=coordx_dynnew,  #ok
                   corrmodel=fit$corrmodel, distance=fit$distance,grid=fit$grid,loc=coordx_dynnew_loc[[j]], #ok
                   model=fit$model, n=fit$n,  param=param, radius=fit$radius,   time=utt_1[j], mse=TRUE,type_krig=type_krig,
                   X=Xnew,Xloc= Xnew_loc[[j]]) #ok  
               pr_st[[j]]=pr$pred ; pr_mse[[j]]=pr$mse
         }
 }
if(local) {

     
           pr_st=pr_mse=list()
           for(j in 1:length(utt_1) ){
            if(is.null(Xnew)) {Xnew_loc[j]=list(NULL)} 
               pr=GeoKrigloc(data=datanew,   coordt=utt, coordx_dyn=coordx_dynnew,  #ok
                             corrmodel=fit$corrmodel, distance=fit$distance,grid=fit$grid,loc=coordx_dynnew_loc[[j]], #ok
                             model=fit$model, n=fit$n,  param=param, radius=fit$radius,   time=utt_1[j], mse=TRUE,type_krig=type_krig,
                             neighb=neighb,maxdist=maxdist,maxtime=maxtime,
                             X=Xnew,Xloc=Xnew_loc[[j]]) #ok 
                pr_st[[j]]=pr$pred ; pr_mse[[j]]=pr$mse
                   }
}   
#####################################

pp=GeoScores(c(data_to_pred[,4]),pred=as.numeric(unlist(pr_st)),mse=as.numeric(unlist(pr_mse)),
score=c("brie","crps","lscore","pe"))

rmse[i]=  pp$rmse
mae[i]=   pp$mae
mad[i]=   pp$mad
lscore[i]=pp$lscore
brie[i]=     pp$brie
crps[i]=  pp$crps
i=i+1
pb(sprintf("i=%g", i))
}  }### End Not Parallel       
##################################################
##################################################

if(parallel)
{
  # Register parallel backend
  future::plan(multisession)  
  progressr::handlers(global = TRUE)
  progressr::handlers("txtprogressbar")
  # Set up a progressor for the foreach loop
  pb <- progressr::progressor(along = 1:K)
  # Initialize result variables
  rmse <- mae <- mad <- lscore <- brie <- crps <- double(K)
  # Parallelized foreach loop using %dopar%
  results <- foreach::foreach(i = 1:K, .combine = 'rbind', 
             .options.future = list(seed = TRUE)) %dofuture% {
    # Sample and prepare data for each iteration                   
    sel_data <- sample(1:NT, round(NT * (1 - n.fold)))
    data_sel <- data_tot[sel_data, ]
    data_to_pred <- data_tot[-sel_data, ]
    data_sel_ord <- data_sel[order(data_sel[, 3]), ]
    data_to_pred_ord <- data_to_pred[order(data_to_pred[, 3]), ]
    DD <- ncol(data_sel_ord)
    k <- 1
    coordx_dynnew <- Xnew <- Xnew_loc <- datanew <- coordx_dynnew_loc <- list()
    utt <- unique(data_sel_ord[, 3])
    utt_1 <- unique(data_to_pred_ord[, 3])
    # Create dynamic coordinates and data for prediction
    for (k in 1:length(utt_1)) {
      ll <- data_to_pred_ord[data_to_pred_ord[, 3] == utt_1[k], ]
      coordx_dynnew_loc[[k]] <- if (is.matrix(ll)) as.matrix(ll)[, 1:2] else matrix(ll[1:2], nrow = 1)
      if (!is.null(X)) Xnew_loc[[k]] <- if (is.matrix(ll)) as.matrix(ll)[, 5:DD] else matrix(ll[5:DD], nrow = 1, byrow = TRUE)
    }
    
    for (k in 1:length(utt)) {
      ss <- data_sel_ord[data_sel_ord[, 3] == utt[k], ]
      datanew[[k]] <- as.vector((ss[, 4]))
      coordx_dynnew[[k]] <- as.matrix(ss[, 1:2])
      if (!is.null(X)) Xnew[[k]] <- as.matrix(ss[, 5:DD])
    }
    
    if (ncol(Xnew[[1]]) == 1) Xnew <- NULL
    param <- append(fit$param, fit$fixed)
    # Estimation 
    if (estimation) {
      fit_s <- GeoFit(
        data = datanew, coordx_dyn = coordx_dynnew, coordt = utt,
        corrmodel = fit$corrmodel, X = Xnew, likelihood = fit$likelihood,
        type = fit$type, grid = fit$grid, copula = fit$copula,
        anisopars = fit$anisopars, est.aniso = fit$est.aniso,
        model = model1, radius = fit$radius, n = fit$n, local = fit$local,
        GPU = fit$GPU, maxdist = fit$maxdist, neighb = fit$neighb, 
        maxtime = fit$maxtime, distance = fit$distance,
        optimizer = optimizer, lower = lower, upper = upper,
        start = fit$param, fixed = fit$fixed
      )
      if (!is.null(fit$anisopars)) {
        fit_s$param$angle <- NULL
        fit_s$param$ratio <- NULL
        fit_s$fixed$angle <- NULL
        fit_s$fixed$ratio <- NULL
      }
      param <- append(fit_s$param, fit_s$fixed)
    }
    # Kriging 
    pr_st <- pr_mse <- list()
    if (!local) {
      for (j in 1:length(utt_1)) {
        if (is.null(Xnew)) Xnew_loc[j] <- list(NULL)
        pr <- GeoKrig(
          data = datanew, coordt = utt, coordx_dyn = coordx_dynnew,
          corrmodel = fit$corrmodel, distance = fit$distance, grid = fit$grid,
          loc = coordx_dynnew_loc[[j]], model = fit$model, n = fit$n, param = param,
          radius = fit$radius, time = utt_1[j], mse = TRUE, type_krig = type_krig,
          X = Xnew, Xloc = Xnew_loc[[j]]
        )
        pr_st[[j]] <- pr$pred
        pr_mse[[j]] <- pr$mse
      }
    }
    
    # Score calculation
    pp <- GeoScores(
      c(data_to_pred[, 4]), pred = as.numeric(unlist(pr_st)),
      mse = as.numeric(unlist(pr_mse)), score = c("brie", "crps", "lscore", "pe")
    )
    # Store results
    result <- c(pp$rmse, pp$mae, pp$mad, pp$lscore, pp$brie, pp$crps)
    # Update progress bar
    pb(sprintf("i=%g", i))
    return(result)
  }
  
  # Extract the results and assign them to respective variables
  for (i in 1:K) {
    rmse[i] <- results[i, 1]
    mae[i] <- results[i, 2]
    mad[i] <- results[i, 3]
    lscore[i] <- results[i, 4]
    brie[i] <- results[i, 5]
    crps[i] <- results[i, 6]
  }
  # Optional: Reset the plan back to sequential execution
  plan(sequential)
  
} #end  parallel

} ## end spacetime

############################################################
########### spatial bivariate case #########################
############################################################
if(bivariate)
{
ns=fit$ns
if(space_dyn) {data1=fit$data[[1]]
               data2=fit$data[[2]]
               coords1=fit$coordx_dyn[[1]]
               coords2=fit$coordx_dyn[[2]]
               }
if(!space_dyn) {data1=fit$data[1,];data2=fit$data[2,]
                coords=cbind(fit$coordx,fit$coordy)
                coords1=coords; coords2=coords
               }
if(!is.null(X)) {
                if(!is.list(fit$X)){X1=fit$X[1:ns[1],];X2=fit$X[(ns[1]+1):(ns[1]+ns[2]),]}
                if( is.list(fit$X)){X1=fit$X[[1]];X2=fit$X[[2]]}
                }
#set.seed(round(seed))


#param=as.list(c(fit$param,fit$fixed))
param=append(fit$param,fit$fixed)

##########
if(!parallel){

rmse=crps=mae=rmse=double(K)


progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")
pb <- progressr::progressor(along = 1:K)
while(i<=K){

#######################################	
if(which==1) {
	          sel_data = sample(1:ns[1],round(ns[1]*(1-n.fold))) 
              data_to_est=fit$data[sel_data]
              cc1= coords1[sel_data,]
              loc_to_pred   = coords1[-sel_data,]
              d1   = data1[sel_data]
              data_to_pred  = data1[-sel_data]
              if(!is.null(X)) { X=rbind(X1[sel_data,],X2); Xloc=rbind(X1[-sel_data,],X2)}
datanew=list();datanew[[1]]=d1;datanew[[2]]=data2;
coordsnew=list();coordsnew[[1]]=cc1;coordsnew[[2]]=coords2;
          }
if(which==2) {
	          sel_data = sample(1:ns[2],round(ns[2]*(1-n.fold))) 
              data_to_est=fit$data[sel_data]
              cc2= coords2[sel_data,]
              loc_to_pred   = coords2[-sel_data,]
              d2   = data2[sel_data]
              data_to_pred  = data2[-sel_data]
              if(!is.null(X)) {X=rbind(X1,X2[sel_data,]); Xloc=rbind(X1,X2[-sel_data,])}
datanew=list();datanew[[1]]=data1;datanew[[2]]=d2;
coordsnew=list();coordsnew[[1]]=coords1;coordsnew[[2]]=cc2;
            }
#dtp[[i]]=data_to_pred

if(estimation) {
          fit_s= GeoFit(data=data_to_est,coordx=coords[sel_data,],corrmodel=fit$corrmodel,X=X,
                            likelihood=fit$likelihood,c,type=fit$type,grid=fit$grid,
                            model="Gaussian",radius=fit$radius,n=fit$n,
                               copula=fit$copula,
                             local=fit$local,GPU=fit$GPU,
                           maxdist=fit$maxdist, neighb=fit$neighb,distance=fit$distance,
                            optimizer=optimizer, lower=lower,upper=upper,
                           # start=as.list(fit$param),fixed=as.list(fit$fixed))
                               start=fit$param,fixed=fit$fixed)

             if(!is.null(fit$anisopars))   
                   {   fit_s$param$angle=NULL;fit_s$param$ratio=NULL;
                       fit_s$fixed$angle=NULL;fit_s$fixed$ratio=NULL
                        }

            param=append(fit_s$param,fit_s$fixed)
              }
#####################################
if(!local) 
    {
        
        pr=GeoKrig(data=datanew, coordx=NULL,   coordt=NULL, coordx_dyn=coordsnew,  #ok
	       corrmodel=fit$corrmodel, distance=fit$distance,grid=fit$grid,loc=loc_to_pred, #ok
	          model=fit$model, n=fit$n, mse=TRUE,#ok
           param=param, 
           radius=fit$radius, sparse=sparse,   time=NULL,  type_krig=type_krig,
             which=which, X=X,Xloc=Xloc,Mloc=Mloc) 
    }#ok  

if(local) 
   {
          pr=GeoKrigloc(data=datanew, coordx=NULL,   coordt=NULL, coordx_dyn=coordsnew,  #ok
         corrmodel=fit$corrmodel, distance=fit$distance,grid=fit$grid,loc=loc_to_pred, #ok
            model=fit$model, n=fit$n, mse=TRUE,#ok
           neighb=neighb, maxdist=maxdist,
           param=param, 
           radius=fit$radius, sparse=sparse,   time=NULL, type_krig=type_krig,
             which=which, X=X,Xloc=Xloc,Mloc=Mloc) #ok  
   }   


pp=GeoScores(data_to_pred,pred=pr$pred,mse=pr$mse,
    score=c("brie","crps","lscore","pe"))

rmse[i]=  pp$rmse
mae[i]=   pp$mae
mad[i]=   pp$mad

lscore[i]=pp$lscore
brie[i]=     pp$brie
crps[i]=  pp$crps


i=i+1
pb(sprintf("i=%g", i))
}               
} ## end not parallel


if(parallel)
{

    } #end  parallel
} ## end bivariate

return(list(rmse=rmse,mae=mae,mad=mad,brie=brie,crps=crps,lscore=lscore,predicted=pred,data_to_pred=dtp))
}
