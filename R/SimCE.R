
######################## space time case: separable with Circ embeeding+ftt ##########################
CE_Space_Time <- function(coords, time.seq, param,corrmodel,distance){

  time.seq <- time.seq - time.seq[1]
  Ns <- nrow(coords); Nt <- length(time.seq)
  Nt.ext <- (2*(Nt-1)); N0 <- Ns*Ns*Nt.ext
#################################################
###### main separable models ####################
#################################################
if(corrmodel=="Matern_Matern"){
  param_s <- list(nugget = 0, sill = 1, scale = param$scale_s, smooth = param$smooth_s)
  param_t <- list(nugget = 0, sill = 1, scale = param$scale_t, smooth = param$smooth_t) 
  corrmodel="Matern"   
  }
if(corrmodel=="GenWend_GenWend"){
  param_s <- list(nugget = 0,sill = 1, scale = param$scale_s, smooth = param$smooth_s,power2=param$power2_s)
  param_t <- list(nugget = 0,sill = 1, scale = param$scale_t, smooth = param$smooth_t,power2=param$power2_t)  
  corrmodel="GenWend"                            
  }
if(corrmodel=="GenWend_Matern_GenWend_Matern"){
  param_s <- list(nugget = 0,sill = 1, scale = param$scale_s, smooth = param$smooth_s,power2=param$power2_s)
  param_t <- list(nugget = 0,sill = 1, scale = param$scale_t, smooth = param$smooth_t,power2=param$power2_t)  
  corrmodel="GenWend_Matern"                         
  }
if(corrmodel=="Kummer_Kummer"){
  param_s <- list(nugget = 0,sill = 1, scale = param$scale_s, smooth = param$smooth_s,power2=param$power2_s)
  param_t <- list(nugget = 0,sill = 1, scale = param$scale_t, smooth = param$smooth_t,power2=param$power2_t)  
  corrmodel="Kummer"                            
  }
if(corrmodel=="Kummer_Matern_Kummer_Matern"){
  param_s <- list(nugget = 0,sill = 1, scale = param$scale_s, smooth = param$smooth_s,power2=param$power2_s)
  param_t <- list(nugget = 0,sill = 1, scale = param$scale_t, smooth = param$smooth_t,power2=param$power2_t)  
  corrmodel="Kummer_Matern"                         
  }
#################################################
################################################# 
  cova.mat.s <- GeoCovmatrix(coordx = coords, corrmodel = corrmodel, param = param_s,distance=distance)$covmatrix
  cova.mat.s <- (1 - param$nugget)*cova.mat.s
  diag(cova.mat.s) <- 1
  ## cholesky decomposition
  chol.s <- t(chol(cova.mat.s))
  auxiliar.seq <- c(time.seq, rev(time.seq[-c(1,Nt)]))
  cova.mat.t <- (1 - param$nugget)*GeoCorrFct(x = auxiliar.seq, 
                                                              corrmodel = corrmodel, 
                                                              param = param_t)$corr
  cova.mat.t[auxiliar.seq == 0] <- 1
  spectrum.t <- sqrt(Re( fft(cova.mat.t) ))
  
  X <- complex(real = rep(0, Ns*Nt.ext), imaginary = rep(0, Ns*Nt.ext) ) 
  for(t in 1:Nt.ext){
    Z0 <- complex(real = rnorm(Ns), imaginary = rnorm(Ns) ) 
    X[1:Ns + Ns*(t-1)] <- (spectrum.t[t]*chol.s)%*%Z0
  }
  X <- matrix(X, ncol = Ns, byrow = T)
  X <- apply(X, 2, function(x) fft(x, T)/sqrt(Nt.ext))
  X <- sqrt(param$sill)*X[1:Nt,]
  result <-  c(t(Re(X)))
  return(result)
}

############################################################
SimCE<-function(M,N,x,y,z,corrmodel,param,mean.val,max.ext)
{
xlim=range(x);ylim=range(y)
Covariate.sim.B <- GRF.CE.Sim(nsim = 1, xrange=c(xlim[1],xlim[2]), yrange=c(ylim[1],ylim[2]), M = M, N = N, 
                              corrmodel = corrmodel, 
                              param = param,
                              mean.val = mean.val, 
                              max.ext = max.ext)
return(Covariate.sim.B)
}
##################
GRF.CE.Sim <- function(nsim,xrange, yrange, M,N, corrmodel, param,mean.val, max.ext){
  fail <- TRUE
  approx.flag <- FALSE
  pow <- seq(1, max.ext, 1)
  k0 <- 1
  ext <- 2^pow[k0]
  while(fail){
    CovDecomp <- Cov_Mat_CE(xrange, yrange, M,N, corrmodel, param, ext, approx = approx.flag) 
    nan.flag <- is.nan( sum(CovDecomp$prefix) )
    if(!nan.flag){
      grid.points <- as.matrix(CovDecomp$grid)
      GRF.sim <- Cov_Mat_CE.sim(CovDecomp, nsim = nsim)
      fail <- FALSE
    }else{
      k0 <- k0 + 1
      ext <- 2^pow[k0]
      if(k0 > max.ext){
        ext <- 4
        approx.flag <- TRUE
      }
    }
  }
  return( list(X = GRF.sim, grid.points = grid.points, k0 = k0, approx.flag = approx.flag))
}

##########################################
Cov_Mat_CE <- function(xrange, yrange, M,N, corrmodel, param, ext = 2, approx = FALSE){ #Variance = 1
  #Grid
  mygrid <- grid.prep(xrange, yrange, M = M, N = N, ext = ext)
  Rx <- mygrid$M.ext * mygrid$cell.width
  Ry <- mygrid$N.ext * mygrid$cell.height
  m.abs.diff.row1 <- abs(mygrid$mcens.ext[1] - mygrid$mcens.ext)
  m.diff.row1 <- pmin(m.abs.diff.row1, Rx - m.abs.diff.row1)
  n.abs.diff.row1 <- abs(mygrid$ncens.ext[1] - mygrid$ncens.ext)
  n.diff.row1 <- pmin(n.abs.diff.row1, Ry - n.abs.diff.row1)
  cent.ext.row1 <- expand.grid(m.diff.row1, n.diff.row1)
  D.ext.row1 <- matrix(sqrt(cent.ext.row1[, 1]^2 + cent.ext.row1[, 2]^2), mygrid$M.ext, mygrid$N.ext)
  param$mean <- 0
  SIGMA.Y.ext.row1 <- my.GeoCorrFct(x = c(D.ext.row1), corrmodel = corrmodel, param = param)
  SIGMA.Y.ext.row1 <- matrix(SIGMA.Y.ext.row1, nrow = mygrid$M.ext, ncol = mygrid$N.ext)
  ddim <- d <- dim(SIGMA.Y.ext.row1)
  dp <- prod(d)
  sdp <- sqrt(dp)
  fft.sigma <- Re(fft(SIGMA.Y.ext.row1, TRUE))
  if( approx ){
    fft.sigma[fft.sigma < 0] <- 0
    prefix <- sqrt(fft.sigma)
    print( "Negative eigenvalues detected. Truncation is done." )
  }else{
    prefix <- suppressWarnings(sqrt(fft.sigma))
  }
  result <- list(dp = dp, d = ddim, sdp = sdp, M = M, N = N, prefix = prefix, grid = expand.grid(x = mygrid$mcens, y = mygrid$ncens) )
  
  return(result)
}

##################
grid.prep <- function(xrange, yrange, M, N, ext = 2){
  cell.width <- diff(xrange)/M
  cell.height <- diff(yrange)/N
  
  mgrid <- seq(xrange[1], xrange[2], by = cell.width)
  ngrid <- seq(yrange[1], yrange[2], by = cell.height)
  mcens <- (mgrid + 0.5 * cell.width)[-(M + 1)]
  ncens <- (ngrid + 0.5 * cell.height)[-(N + 1)]
  
  if (ext <= 1) 
    mgrid.ext <- ngrid.ext <- mcens.ext <- ncens.ext <- M.ext <- N.ext <- NULL else {
      M.ext <- ext * M
      N.ext <- ext * N
      mgrid.ext <- seq(xrange[1], xrange[2] + (ext - 1) * diff(xrange), by = cell.width)
      ngrid.ext <- seq(yrange[1], yrange[2] + (ext - 1) * diff(yrange), by = cell.height)
      mcens.ext <- (mgrid.ext + 0.5 * cell.width)[-(M.ext + 1)]
      ncens.ext <- (ngrid.ext + 0.5 * cell.height)[-(N.ext + 1)]
    }
  
  return(list(M = M, N = N, mgrid = mgrid, ngrid = ngrid, mcens = mcens, ncens = ncens, 
              cell.width = cell.width, cell.height = cell.height, M.ext = M.ext, N.ext = N.ext, 
              mgrid.ext = mgrid.ext, ngrid.ext = ngrid.ext, mcens.ext = mcens.ext, ncens.ext = ncens.ext))
}
##################
my.GeoCorrFct <- function(x = x, corrmodel = corrmodel, param = param){
  if(corrmodel != "none"){
    if(param$nugget == 0){
      result <- GeoCorrFct(x = x, corrmodel = corrmodel, param = param)$corr
    }else{
      result <- GeoCorrFct(x = x, corrmodel = corrmodel, param = param)$corr
      result[x == 0] <- 1
    }
  }else{
    result <- rep(0, length(x))
    result[x == 0] <- 1
  }
  return(result)
}
##################
Cov_Mat_CE.sim <- function(ex1,nsim = 2){
  
  dp <- ex1$dp
  d <- ex1$d
  sdp <- ex1$sdp
  M <- ex1$M
  N <- ex1$N
  prefix <- ex1$prefix
  X <- matrix(0, nrow = (M*N), ncol = nsim)
  for(i in 1:nsim){
    std <- rnorm(dp)
    realz <- prefix * (fft(matrix(std, d[1], d[2]))/sdp)
    realz <- as.vector(Re(fft(realz, TRUE)/sdp)[1:M, 1:N])
    X[,i] <- realz
  }
  return(X)
}



