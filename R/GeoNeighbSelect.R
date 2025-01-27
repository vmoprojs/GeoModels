####################################################
### File name: GeoNeighbSelect.r
####################################################

GeoNeighbSelect <- function(data, coordx,coordy=NULL, coordz=NULL,coordt=NULL, coordx_dyn=NULL,
  copula=NULL,corrmodel=NULL, distance="Eucl",fixed=NULL,anisopars=NULL,
  est.aniso=c(FALSE,FALSE), grid=FALSE, likelihood='Marginal', 
  lower=NULL,neighb=c(1,2,3,4,5),
  maxtime=Inf, memdist=TRUE,model='Gaussian',n=1, ncores=NULL,
  optimizer='Nelder-Mead', parallel=FALSE, bivariate=FALSE,
  radius=6371, start=NULL,  type='Pairwise', upper=NULL,  weighted=FALSE,
  X=NULL,nosym=FALSE,spobj=NULL,spdata=NULL,vario=NULL)
{



if(!is.numeric(neighb))  stop("neighb must be a numeric vector")
if(sum(neighb-floor(neighb))) stop("neighb must be a positive integer numeric vector")

estimates=best_T=best_K= NULL


if(!is.null(try(CkCorrModel(corrmodel), TRUE)))
{
    bivariate<-CheckBiv(CkCorrModel(corrmodel))
    spacetime<-CheckST(CkCorrModel(corrmodel))
    space=!spacetime&&!bivariate
}else {stop("correlation model is not valid\n")}

############## computing semivariogram  if necessary###############
if(!is.null(vario))
{
if(!inherits(vario,"GeoVariogram"))  stop("A GeoVariogram object is needed as input for vario\n")
 if( !((vario$bivariate&&bivariate)||(!is.null(vario$bint)&&spacetime)||(is.null(vario$bint)&&!bivariate)) )
       stop("The GeoVariogram object is not of the same type of the correlation model\n")

 semiv=vario
}
else stop("A GeoVariogram object is needed as input for vario\n")
######################################################################


######################################################



####################
coremax=parallel::detectCores()
if(is.na(coremax)||coremax==1) parallel=FALSE
####################


K=length(neighb); res=double(K)

P=NULL
if(spacetime){
   if(!is.numeric(maxtime))  stop("neighb must be a numeric vector")
   P=length(maxtime)
   res=double(K*P)
}
#######################################################################################
##################################### SPATIAL #########################################
#######################################################################################
if(space||bivariate) {

M=1
##################################### spatial not parallel ############################
if(!parallel)
{

progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")
pb <- progressr::progressor(along = 1:K)

for(M in 1:K) {

    pb(sprintf("k=%g", M)) 

  aa=GeoFit(data=data, coordx=coordx, coordy=coordy, coordz=coordz,coordt=coordt,coordx_dyn=coordx_dyn,copula=copula,corrmodel=corrmodel, distance=distance,
                         fixed=fixed,anisopars=anisopars,est.aniso=est.aniso, grid=grid, likelihood=likelihood, 
                         lower=lower,neighb=neighb[M],
                          maxtime=maxtime, memdist=memdist,model=model,n=n, 
                          optimizer=optimizer, 
                         radius=radius, start=start,  
                         type=type, upper=upper,  weighted=weighted,X=X,nosym=nosym,spobj=spobj,spdata=spdata)
 #res=GeoResiduals(aa) 
 #vario = GeoVariogram(data=res$data, coordx=coords,coordy=coordy, coordz=coordz,coordt=coordt,coordx_dyn=coordx_dyn,cloud=FALSE,distance=distance,
    # grid=grid,maxdist
       #maxdist=0.3) # empirical variogram 
 #semiv=GeoCovariogram(res, show.vario=TRUE, vario=vario,pch=20)



 clest=  aa$param
 estimates=rbind(estimates,unlist(clest))

 if(aa$convergence=="Successful")
 {
  ## first method slightly faster
 #cc=GeoCorrFct(semiv$centers,t=semiv$centert,corrmodel=corrmodel, model=model,distance=distance, param=c(aa$param,aa$fixed),radius=radius,n=n,covariance=TRUE,variogram=TRUE)$corr
 # res[M]=sum(cc - semiv$variograms)^2

 ## second method slighty slower but with graphics..
 cc=GeoCovariogram(fitted=aa,distance=distance,show.vario=TRUE, vario=semiv,pch=20,invisible=TRUE)
 res[M]=cc 
 }
 else { res[M]=Inf}

 }
}

##################################### end spatial not parallel ############################
##################################### spatial parallel ############################
if(parallel)
{
if(is.null(ncores)){ n.cores <- coremax - 1 }
else
{  if(!is.numeric(ncores)) stop("number of cores not valid\n")
   if(ncores>coremax||ncores<1) stop("number of cores not valid\n")
   n.cores=ncores
}
cat("Performing",K,"estimations using",n.cores,"cores...\n")

future::plan(multisession, workers = n.cores)

progressr::handlers(global = TRUE)
pb <- progressr::progressor(along = 1:K)
k=1
# Ejecutar el bucle de forma paralela usando foreach
results <- foreach(M = 1:K, .combine = rbind,
                   .options.future = list(seed = TRUE)) %dofuture% {
 
  
  # Ejecutar la función GeoFit
  aa <- GeoFit(data = data, coordx = coordx, coordy = coordy, coordz=coordz,coordt = coordt, coordx_dyn = coordx_dyn,
               copula = copula, corrmodel = corrmodel, distance = distance, fixed = fixed, anisopars = anisopars,
               est.aniso = est.aniso, grid = grid, likelihood = likelihood, lower = lower, neighb = neighb[M],
               maxtime = maxtime, memdist = memdist, model = model, n = n, optimizer = optimizer,
               radius = radius, start = start, type = type, upper = upper, weighted = weighted,
               X = X, nosym = nosym, spobj = spobj, spdata = spdata)

  # Almacenar resultados
  #clest <- aa$param
  estimates <- unlist(aa$param)

  cc <- GeoCovariogram(fitted = aa, distance = distance, show.vario = TRUE, vario = semiv, pch = 20, invisible = TRUE)
  res=ifelse(aa$convergence == "Successful",cc,Inf)
  
  # Calcular el resultado de acuerdo a la convergencia
  #if (aa$convergence == "Successful") {
#
 #   res <- cc
  #} else {
   # res <- Inf
  #}
    pb(sprintf("k=%g", M))  # Actualizar el progreso
  # Retornar los resultados como una fila
  c(estimates, res = res)
}
# Separar las columnas de resultados en 'estimates' y 'res'
estimates <- results[, -ncol(results)]  # Todas las columnas menos la última
rownames(estimates) <- NULL
res <- results[, ncol(results)]         # Última columna
future::plan(sequential)
}

}
#######################################################################################
##################################### SPATIO TEMPORAL ################################
#######################################################################################


if(spacetime){
######################################## Space time  NO parallel ############################
if(!parallel)
{

progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")
pb <- progressr::progressor(along = 1:(K*P))

  F=1
for(L in 1:P) { 
for(M in 1:K) {

  pb(sprintf("k=%g", F)) 
  aa=GeoFit(data=data, coordx=coordx, coordy=coordy, coordz=coordz,coordt=coordt, coordx_dyn=coordx_dyn,copula=copula,corrmodel=corrmodel, distance=distance,
                         fixed=fixed,anisopars=anisopars,est.aniso=est.aniso, grid=grid, likelihood=likelihood, 
                         lower=lower,neighb=neighb[M],
                          maxtime=maxtime[L], memdist=memdist,model=model,n=n, 
                          optimizer=optimizer, 
                         radius=radius, start=start,  
                         type=type, upper=upper,  weighted=weighted,X=X,nosym=nosym,spobj=spobj,spdata=spdata)
 
 estimates=rbind(estimates,unlist(aa$param))
 
 if(aa$convergence=="Successful")
 {  cc=GeoCovariogram(fitted=aa,distance=distance,show.vario=TRUE, vario=semiv,pch=20,fix.lagt=1,fix.lags=1,invisible=TRUE)
    res[F]=cc; F=F+1 
 }
 else { res[F]=Inf; F=F+1 }
}}

}
######################################## end Space time NO parallel ############################
######################################## Space time   parallel ############################
if(parallel)
{
  if(is.null(ncores)){ n.cores <- coremax - 1 }
  else
  {  if(!is.numeric(ncores)) stop("number of cores not valid\n")
    if(ncores>coremax||ncores<1) stop("number of cores not valid\n")
    n.cores=ncores
  }
  cat("Performing",K,"estimations using",n.cores,"cores...\n")
  
  future::plan(multisession, workers = n.cores)
  progressr::handlers(global = TRUE)
  progressr::handlers("txtprogressbar")
  pb <- progressr::progressor(along = 1:(K*P))
  
  iteration_count <- 1  # Cambio de F por un nombre descriptivo
  estimates <- NULL
  results <- vector("list", P)  # Lista para guardar futuros resultados

  # Crear los futuros para P paralelizados
  for (L in 1:P) {
    results[[L]] <- future::future({
      sub_results <- vector("list", K)

      for (M in 1:K) {
        progressr::with_progress({
          pb(sprintf("k=%g", iteration_count))

          # Llamada a la función GeoFit con los parámetros específicos
          aa <- GeoFit(data = data, coordx = coordx, coordy = coordy, coordz=coordz,coordt = coordt, coordx_dyn = coordx_dyn,
                       copula = copula, corrmodel = corrmodel, distance = distance, fixed = fixed, anisopars = anisopars,
                       est.aniso = est.aniso, grid = grid, likelihood = likelihood, lower = lower, neighb = neighb[M],
                       maxtime = maxtime[L], memdist = memdist, model = model, n = n, optimizer = optimizer,
                       radius = radius, start = start, type = type, upper = upper, weighted = weighted, X = X,
                       nosym = nosym, spobj = spobj, spdata = spdata)

          # Verificar si la estimación fue exitosa
          if (aa$convergence == "Successful") {
            clest <- unlist(aa$param)
            cc <- GeoCovariogram(fitted = aa, distance = distance, show.vario = TRUE, vario = semiv, pch = 20,
                                 fix.lagt = 1, fix.lags = 1, invisible = TRUE)
            result_res <- cc  # Guardar el resultado de la covariograma
          } else {
            clest <- Inf
            result_res <- Inf  # En caso de no converger
          }

          iteration_count <<- iteration_count + 1  # Incrementar el contador de iteraciones
          sub_results[[M]] <- list(estimates = clest, res = result_res)
        })
      }
      sub_results
    }, seed = TRUE)
  }

  # Recuperar y combinar resultados una vez completados
  results <- lapply(results, future::value)
  estimates <- do.call(rbind, lapply(results, function(res) do.call(rbind, lapply(res, function(x) x$estimates))))
  res <- unlist(lapply(results, function(res) unlist(lapply(res, function(x) x$res))))

  # Restablece el plan a secuencial si deseas
  future::plan(sequential)
}
######################################## end Space time   parallel ############################
}
################# end SPATIOTEMPORAL ######################################################



if(space||bivariate) {
indexmin=which.min(res)
bestK=neighb[indexmin]
}
if(spacetime)
{
  indexmin=which.min(res)
  bestKT=c(as.matrix(expand.grid(neighb,maxtime))[indexmin,])
  bestK=as.numeric(bestKT[1]);best_T=as.numeric(bestKT[2]);
}

a=list(best_neighb=bestK,best_maxtime=best_T,res=res,estimates=estimates,best_est=estimates[indexmin,])

return(a)
}