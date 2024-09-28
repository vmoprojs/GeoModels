#########################################################
#########################################################
############ functions for TB############################
#########################################################
#########################################################
spectral_density_1dR <- function(param, corrmodel, u_range = c(-1,1), N){

  model=CkCorrModel(corrmodel)
  av <- param$scale
  nu1v <- param$smooth
  params_other <- param$power2
  if(is.null(params_other)) params_other <- 0
  norm_u <- seq(u_range[1],u_range[2],l = N)  
  
  #result <- .C("spectral_density_1d", as.double(norm_u), 
   #                                   as.integer(N), 
    #                                  as.double(av), 
     #                                 as.double(params_other),
      #                                as.double(nu1v), 
       #                               as.integer(model), 
        #                              simu1 = as.double(numeric(N)),PACKAGE='GeoModels',DUP=TRUE,NAOK=TRUE)$simu1
   result=dotCall64::.C64("spectral_density_1d",
        SIGNATURE = c("double","integer","double", "double","double","integer","double"),
                     norm_u, N, av, params_other, nu1v, model, simu1 = dotCall64::numeric_dc(N),
        INTENT =    c("r","r","r", "r","r", "r","w"),PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$simu1

  return(result)
}

##################################################### 
######### bivarate case #############################
#####################################################
tbm2d<- function(coord, coordt, param, corrmodel,L,bivariate){

  N=1; n <- dim(coord)[1]; d <- 2
    if(corrmodel == "Matern"){
       CC=1; N=1
       a <- as.numeric(param['scale']);  nu1 <- as.numeric(param['smooth'])
       a0 = a; nu0 =nu1
       parameters <- list("CC" = CC, "a" = a,"nu1" = nu1,"nu2" = 0 )
       P <- 1;  vtype = 0}

  if (corrmodel=="Bi_matern"||corrmodel=="Bi_Matern"){
    corrmodel <- "Matern"
    a <- matrix(0,2,2)
    a[1,1] <- as.numeric(param['scale_1'])
    a[2,2] <- as.numeric(param['scale_2'])
    a[1,2] <- a[2,1] <- as.numeric(param['scale_12'])
    
    nu1 <- matrix(0,2,2)
    nu1[1,1] <- as.numeric(param['smooth_1'])
    nu1[2,2] <- as.numeric(param['smooth_2'])
    nu1[1,2] <- nu1[2,1] <- as.numeric(param['smooth_12'])
    CC <- matrix(0,2,2)
    CC[1,1] <- CC[2,2] <- 1
    CC[1,2] <- CC[2,1] <- as.numeric(param['pcol'])
    a0 <- min(a)
    nu0 <- min(nu1)
    parameters <- list("C" = CC, "a" = a,"nu1" = nu1,"nu2" =matrix(0,2,2) )
    P<-2; vtype = 0  

  }
   parametersg <- list("a" = a0,"nu1" = nu0)

      A <- matrix(0, P, P*L*N); B <- matrix(0, P, P*L*N)
      #S <- ceiling(1e7*runif(3)); set.seed(S[1])
      G <- matrix(rgamma(P*L*N, nu0, scale = 1),P*L*N,d); #set.seed(S[2])
      u <- matrix(rnorm(P*L*N*d), P*L*N, d)/sqrt(G*2)/a0/(2*pi); #set.seed(S[3])
      phi <- 2*pi*runif(P*L*N)
      sequen <- c(seq(0,n-0.5, by = ceiling(1e6/P/N)),n)

    m = c()
   for (i in 1:(length(sequen)-1)){ m1 <- sequen[i+1]-sequen[i]; m = c(m1,m)} 
simu11 = as.numeric( rep(0,N*P*sum(m)*(length(sequen)-1)))

result=dotCall64::.C64("for_c",
         SIGNATURE = c("integer","double","double", "double",
                     "double","integer","integer","integer","integer",
                     "double","double","double","double","double",
                     "integer","integer","integer",
                     "double","double","integer","integer","double",
                     "double"), 
                         d,c(a),c(nu1), c(CC), c(parameters$nu2),P,N,L,CkCorrModel(corrmodel),
                         c(u),a0,nu0,c(A),c(B),c(sequen),length(sequen),n,coord,phi,vtype,m,
                         simu1=dotCall64::numeric_dc(length(simu11)),L,
         INTENT =    c("r","r","r","r","r","r","r","r","r", "r",
                       "r","r","r","r","r","r","r","r","r", "r",
                       "r", "rw","r"),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$simu1
  simu =  matrix(result,n,P)
  return(simu)
}
#######################################################
######### univariate case #############################
#######################################################
tbm2d_uni <- function(coord, coordt, param, corrmodel, L, bivariate){
  
  N=1; n <- dim(coord)[1]; d <- 2
 ################################### 
  if(corrmodel == "Matern"){
    CC=as.double(param['sill']); N=1
    a <- as.double(param['scale']);  
    nu1 <- as.double(param['smooth'])
    other <- as.double(0)
    a0 = a; nu0 = nu1
    parameters <- list("CC" = CC, "a" = a,"nu1" = nu1, "nu2" = 0, "other" = other)
    P <- 1;  vtype = 0
  }
  
  if(corrmodel == "Kummer"){
    CC=as.double(param['sill']); N=1
    a <- as.double(param['scale']);  
    nu1 <- as.double(param['smooth']);
    other <- as.double(param['power2'])
    a0 = a; nu0 = nu1
    a0 <- a <- param$scale/sqrt(2*nu1)
    parameters <- list("CC" = CC, "a" = a,"nu1" = nu1,"nu2" = 0, "other" = other)
    P <- 1;  vtype = 0
  }
  
  if(corrmodel == "Kummer_Matern"){
    CC=as.double(param['sill']); N=1 
    a <- as.double(param['scale']);  
    nu1 <- as.double(param['smooth']);
    other <- as.double(param['power2'])
    a0 = a; nu0 = nu1
    a0 <- a <- param$scale/sqrt(2*nu1)
    corrmodel = "Kummer"
    parameters <- list("CC" = CC, "a" = a,"nu1" = nu1,"nu2" = 0, "other" = other)
    P <- 1;  vtype = 0
  }
  
  if(corrmodel == "GenWend"){
    CC <- as.double(param['sill']); N=1
    a <- as.numeric(param['scale']);  
    nu1 <- as.numeric(param['smooth'])
    other <- as.numeric(param['power2'])
    a0 = a; nu0 = nu1
    parameters <- list("CC" = CC, "a" = a,"nu1" = nu1, "nu2" = 0, "other" = other)
    P <- 1;  vtype = 0
  }
  

    if(corrmodel == "GenWend_Matern"){
    CC <- as.double(param['sill']); N=1
    a <- as.numeric(param['scale']);  
    nu1 <- as.numeric(param['smooth'])
    other <- 1/as.numeric(param['power2'])
    a0 = a; nu0 = nu1
     a0 <- a <- param$scale*(gamma(other+2*nu1+1)/gamma(other))^(-(1+2*nu1))
    corrmodel = "GenWend"
    parameters <- list("CC" = CC, "a" = a,"nu1" = nu1, "nu2" = 0, "other" = other)
    P <- 1;  vtype = 0
  }

  #parametersg <- list("a" = a0,"nu1" = nu0)11
  parametersg <- list("a" = a0, "nu1" = nu0, other = other, CC = CC)
  A <- matrix(0, P, P*L*N); B <- matrix(0, P, P*L*N)
  #S <- ceiling(1e7*runif(3)); set.seed(S[1])
  
  #Simulando frecuencias de una densidad espectral Matern
  #G <- matrix(rgamma(P*L*N, nu0, scale = 1),P*L*N,d); #set.seed(S[2])
  #u <- matrix(rnorm(P*L*N*d), P*L*N, d)/sqrt(G*2)/a0/(2*pi); #set.seed(S[3])
  #Simulando frecuencias desde g
  
  u <- frequency_sampler(P, L, N, d, parametersg, corrmodel)
  #u <- matrix(rnorm(P*L*N*d), P*L*N, d)/sqrt(G*2)/a0/(2*pi); #set.seed(S[3])
  phi <- 2*pi*runif(P*L*N)
  sequen <- c(seq(0,n-0.5, by = ceiling(1e6/P/N)),n)
  m = c()
  for (i in 1:(length(sequen)-1)){ m1 <- sequen[i+1]-sequen[i]; m = c(m1,m)}
  simu11 = as.numeric( rep(0,N*P*sum(m)*(length(sequen)-1)))
  result=dotCall64::.C64("for_c",
                         SIGNATURE = c("integer","double", "double","double",
                                       "double","integer","integer","integer","integer","double",
                                       "double","double","double","double",
                                       "integer","integer","integer",
                                       "double","double","integer","integer","double",
                                       "double","double"),
                         d_v = d, a_v = c(a), nu1_v = c(nu1), C_v = c(CC), nu2_v = c(parameters$nu2), 
                         P = P, N = N, L = L, model = CkCorrModel(corrmodel),
                         u = c(u), a0 = a0, nu0 = nu0, A = c(A), B = c(B), sequen = c(sequen),
                         largo_sequen = length(sequen), n = n, coord = coord, phi = phi, vtype = vtype,
                         m1 = m,
                         simu1=dotCall64::numeric_dc(length(simu11)), L1 = L, other,
                         INTENT = c("r", "r", "r", "r", "r","r","r","r", "r",
                                    "r", "r", "r", "r", "r", "r","r","r","r", "r",
                                    "r", "r", "rw","r","r"),
                         PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)$simu1
  simu =  matrix(sqrt(CC)*result,n,P)
  return(simu)
}


frequency_sampler <- function(P, L, N, d, parametersg, corrmodel){
  u_range <- 100*10*parametersg$CC*c(-parametersg$a, parametersg$a)
  u_range <- c( pmax(2*pi*u_range[1], -20.5), pmin(2*pi*u_range[2], 20.5)) 
  spectral.dens <- parametersg$CC*spectral_density_1dR(param = list(smooth = parametersg$nu1, 
                                                     scale = parametersg$a,
                                                     power2 = parametersg$other), 
                                        corrmodel, u_range = u_range, N = L)
  spectral.distr <- cumsum(spectral.dens)
  spectral.dens.fun <- approxfun(spectral.distr , seq(u_range[1], u_range[2], l = L) )
  u.test <- runif(P*L*N*d, min(spectral.distr), max(spectral.distr))  
  u <- spectral.dens.fun(u.test)
  return(u)
}


