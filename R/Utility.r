####################################################
### Authors:  Moreno Bevilacqua, Víctor Morales Oñate.
### Email:  moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Departamento de Estadistica
### Universidad de Valparaiso
### File name: Utility.r
### Description:
### This file contains a set of procedures
### for the set up of all the package routines.
### Last change: 28/03/2020.
####################################################

# Check if the correlation is bivariate
CheckBiv <- function(numbermodel)
{
    CheckBiv <- NULL
    if(!is.null(numbermodel)){
    if(numbermodel >= 101 & numbermodel <= 140) CheckBiv <- TRUE
    else CheckBiv <- FALSE
    }
    return(CheckBiv)
}
# Check if the correlation is spatial or spatial-temporal
CheckST <- function(numbermodel)
{
    CheckST <- NULL
     if(!is.null(numbermodel)){
    if(numbermodel >40 & numbermodel <= 100) CheckST <- TRUE
    else  CheckST <- FALSE
     }
    return(CheckST)
}


# Check the type of distances
CheckDistance<- function(distance)
{
    CheckDistance <- NULL
    CheckDistance <- switch(distance,
    eucl=0,
    Eucl=0,
    chor=1,
    Chor=1,
    geod=2,
    Geod=2)
    return(CheckDistance)
}
### Procedures are in alphabetical order.
CkCorrModel <- function(corrmodel)
  {
    CkCorrModel <- NULL
    # Correlation function are in alphabetical order
    CkCorrModel <- switch(corrmodel,
            #spatial models
                             cauchy=1,Cauchy=1,
                             Matern1=2,
                             Matern2=3,
                             exponential=4,Exponential=4,Exp=4,exp=4,Matern0=4,
                             dagum = 5, Dagum = 5,
                             GenWend_Matern=6,GenWend_matern=6,Genwend_Matern=6,genWend_Matern=6,
                             GenWend_Matern2=7,GenWend_matern2=7,Genwend_Matern2=7,genWend_Matern2=7,
                             gencauchy=8,Gencauchy=8,GenCauchy=8,genCauchy=8,
                             Shkarofski=10,shkarofski=10,
                             Wend0=11,wend0=11,
                             stable=12,Stable=12,
                             Wend1=13,wend1=13,
                             matern=14,Matern=14,
                             Wend2=15,wend2=15,
                             wave=16,Wave=16,
                             Multiquadric=17,multiquadric=17,
                             Sinpower=18,sinpower=18,
                             genWend=19,Genwend=19,GenWend=19,genwend=19,
                             smoke=20,Smoke=20,
             # spatial-temporal non-separable models
                             gneiting=42,Gneiting=42,  #ok
                             iacocesare=44,Iacocesare=44, #ok
                             porcu=46,Porcu=46,
                             stein=48,Stein=48,          #ok
                             porcu1=50,Porcu1=50,
                             gneiting_GC=52,Gneiting_GC=52, #ok
                             gneiting_GC2=54,Gneiting_GC2=54,
                             sinpower_st=56,Sinpower_st=56,    #ok
                             multiquadric_st=58,Multiquadric_st=58,   #ok
                             gneiting_mat_time=61,Gneiting_mat_time=61, #ok
                             gneiting_mat_space=62,Gneiting_mat_space=62, #ok
                             Wen0_space=63,wen0_space=63,  #ok
                             Wen0_time=64,wen0_time=64,    #ok
                             Wen1_space=65,wen1_space=65,  #ok
                             Wen1_time=66,wen1_time=66,    #ok
                             Wen2_space=67,wen2_space=67,  #ok
                             Wen2_time=68,wen2_time=68,    #ok
                             Wen_time=88,
                             Wen_space=87,
              # spatial-temporal separable models
                             Wend0_Wend0=69,wend0_wend0=69, #ok
                             Wend0_Wend1=70,wend0_wend1=70, #ok
                             Wend0_Wend2=71,wend0_wend2=71, #ok
                             Wend1_Wend0=72,wend1_wend0=72, #ok
                             Wend1_Wend1=73,wend1_wend1=73, #ok
                             Wend1_Wend2=74,wend1_wend2=74, #ok
                             Wend2_Wend0=75,wend2_wend0=75, #ok
                             Wend2_Wend1=76,wend2_wend1=76, #ok
                             Wend2_Wend2=77,wend2_wend2=77, #ok 
                             exp_cauchy=82,Exp_Cauchy=82,
                             exp_exp=84,Exp_Exp=84,
                             Matern_Matern=86, matern_matern=86,  #ok
                             #exp_gauss=86, Exp_Gauss=86,
                             #exp_cos=88,Exp_Cos=88,
                             #matern_cauchy=90,Matern_Cauchy=90,
                             #matern_exp=92, Matern_Exp=92,
                             stable_stable=94,Stable_Stable=94,
                             prove=96,
              # Bivariate models
                             Bi_wend0_sep=111,Bi_Wend0_sep=111,
                             Bi_Wend0_contr=129,Bi_wend0_contr=129,
                             Bi_wend0=112,Bi_Wend0=112,
                             
                             Bi_smoke=117,Bi_Smoke=117,
                              Bi_smoke_sep=119,Bi_Smoke_sep=119,
                              Bi_smoke_contr=121,Bi_Smoke_contr=121,
                             

                             Bi_wend1_sep=113,Bi_Wend1_sep=113,
                             Bi_Wend1_contr=131,Bi_wend1_contr=131,
                             Bi_wend1=114,Bi_Wend1=114,  
                             

                             Bi_wend2_sep=115,Bi_Wend2_sep=115,
                             Bi_wend2=116,Bi_Wend2=116,
                             Bi_Wend2_contr=120,Bi_wend2_contr=120,

                             Bi_matern_contr=118,Bi_Matern_contr=118, 
                             Bi_matern_sep=122,  Bi_Matern_sep=122, 
                             Bi_matern=128,Bi_Matern=128,

                             Bi_LMC_contr=124,
                             Bi_LMC=126,

                             Bi_GenWend_sep=130,Bi_genWend_sep=130,
                             Bi_GenWend_contr=134,Bi_genWend_contr=134,
                             Bi_GenWend=132,Bi_genWend=132,

                             Bi_Matern_Cauchy=136,  
                             Bi_GenMatern_Cauchy=137,   Bi_genMatern_Cauchy=137,
                             Bi_genCauchy=138,Bi_GenCauchy=138,
                             Bi_Stable=139,Bi_stable=139,
                             
            ########Tapers
                   ##spatial
                             Bohman=28,
                             Wendland0=30,wendland0=30,
                             Wendland1=32,wendland1=32,
                             Wendland2=34, wendland2=34,
                             unit_matrix=36,
                   ##bivariate
                             Bi_Wendland1=140,Bi_wendland1=140,
                             Bi_Wendland2=142,Bi_wendland2=142,
                             Bi_Wendland3=144,Bi_wendland3=144,
                             Bi_wendland_asy=146,Bi_wendland_asy=146,
                             unit_matrix_biv=147,
                   ##spacetime separable
                             Wendland0_Wendland0=200,
                             Wendland0_Wendland1=202,
                             Wendland0_Wendland2=204,
                             Wendland1_Wendland0=206,
                             Wendland1_Wendland1=208,
                             Wendland1_Wendland2=210,
                             Wendland2_Wendland0=212,
                             Wendland2_Wendland1=214,
                             Wendland2_Wendland2=216,

                    ##spacetime no separable
                             Wendland0_time=218,
                             Wendland1_time=220,
                             Wendland2_time=222,

                             Wendland0_space=224,
                             Wendland1_space=226,
                             Wendland2_space=228,
                             unit_matrix_st=230)
    return(CkCorrModel)
  }



# checking models valid only on the sphere
CheckSph<- function(numbermodel)
  {
    Check <- FALSE
    if(numbermodel %in% c(17,18,56,58,20,117))  Check=TRUE    
    return(Check)
  }




CkInput <- function(coordx, coordy, coordt, coordx_dyn, corrmodel, data, distance, fcall, fixed, grid,
                      likelihood, maxdist, maxtime,  model, n,  optimizer, param,
                       radius, start, taper, tapsep, type, varest, vartype, weighted,
                       X)
  {
    error <- NULL
    replicates=1
    # START Include internal functions:
    CheckParamRange <- function(param)
    { 
      #  if(!is.na(param['df'])) if(param['df'] > 1/2 || param['df'] < 0 ) return(FALSE)
        #if(!is.na(param['tail'])) if(param['tail'] >0.5) return(FALSE)
        if(!is.na(param['shape'])) if(param['shape'] <0) return(FALSE)
        if(!is.na(param['nugget'])) if(param['nugget'] < 0||param['nugget'] >= 1) return(FALSE)
        if(!is.na(param['nugget_1'])) if(param['nugget_1'] < 0) return(FALSE)
        if(!is.na(param['nugget_2'])) if(param['nugget_2'] < 0) return(FALSE)
        if(!is.na(param['power'])) if(param['power'] <=0 || param['power'] > 2) return(FALSE)
        if(!is.na(param['power_s'])) if(param['power_s'] <=0 || param['power_s'] > 2) return(FALSE)
        if(!is.na(param['power_t'])) if(param['power_t'] <=0 || param['power_t'] > 2) return(FALSE)
        if(!is.na(param['power1'])) if(param['power1'] <=0 || param['power1'] > 2) return(FALSE)
        if(!is.na(param['power2'])) if(param['power2'] <= 0) return(FALSE)
        if(!is.na(param['power2_1'])) if(param['power2_1'] <= 0) return(FALSE)
        if(!is.na(param['power2_12'])) if(param['power2_12'] <= 0) return(FALSE)
        if(!is.na(param['power2_2'])) if(param['power2_2'] <= 0) return(FALSE)
        if(!is.na(param['sep'])) if(param['sep'] < 0 || param['sep'] > 1) return(FALSE)
        if(!is.na(param['scale'])) if(param['scale'] <= 0) return(FALSE)
        if(!is.na(param['scale_s'])) if(param['scale_s'] <= 0) return(FALSE)
        if(!is.na(param['scale_t'])) if(param['scale_t'] <= 0) return(FALSE)
        if(!is.na(param['scale_1']))   if(param['scale_1'] < 0) return(FALSE)
        if(!is.na(param['scale_2']))   if(param['scale_2'] < 0) return(FALSE)
        if(!is.na(param['scale_12']))   if(param['scale_12'] <= 0) return(FALSE)
        if(!is.na(param['sill'])) if(param['sill'] <= 0) return(FALSE)
        if(!is.na(param['sill_1'])) if(param['sill_1'] <= 0) return(FALSE)
        if(!is.na(param['sill_2'])) if(param['sill_2'] <= 0) return(FALSE)
        #if(!is.na(param['smooth'])) if(param['smooth'] <= 0) return(FALSE)
        if(!is.na(param['smooth_s'])) if(param['smooth_s'] < 0) return(FALSE)
        if(!is.na(param['smooth_t'])) if(param['smooth_t'] < 0) return(FALSE)
        #if(!is.na(param['smooth_1'])) if(param['smooth_1'] < 0) return(FALSE)
        if(!is.na(param['smooth_2'])) if(param['smooth_2'] < 0) return(FALSE)
        #if(!is.na(param['smooth_12'])) if(param['smooth_12'] < 0) return(FALSE)
        if(!is.na(param['pcol'])) if(param['pcol'] > 1 || param['pcol'] < -1  ) return(FALSE)
        return(TRUE)  
    }
   
    bivariate=CheckBiv(CkCorrModel(corrmodel))
    if(!bivariate)  if(length(coordt)>0&&is.list(X)) X=X[[1]]
    if(!bivariate) {if(is.null(X))  {X=1;num_betas=1} 
                    else num_betas=ncol(X)  }
    if( bivariate) {if(is.null(X))  {X=1;num_betas=c(1,1)} 
                    else { 
                               if(is.list(X))  num_betas=c(ncol(X[[1]]),ncol(X[[2]]))
                               else  num_betas=c(ncol(X),ncol(X))
                          } }
    #bivariate<-CheckBiv(CheckCorrModel(corrmodel))
    #spacetime<-CheckST(CheckCorrModel(corrmodel))
    #if(is.null(bivariate)&&is.null(bivariate))
    # END Include internal functions
    if(fcall!="Kriging"){
    ### START Checks inserted input
    # START common check fitting and simulation
        if(missing(coordx_dyn)){
    if(missing(coordx) || !is.numeric(coordx)){
        error <- 'insert a suitable set of numeric coordinates\n'
        return(list(error=error))}

    if(!is.null(coordy) & !is.numeric(coordy)){
        error <- 'insert a suitable set of numeric coordinates\n'
        return(list(error=error))}
    }
    if(missing(corrmodel) || !is.character(corrmodel)){
        error <- 'insert the correlation model\n'
        return(list(error=error))}
    if(!is.null(grid) & !is.logical(grid)){
        error <- 'the parameter grid need to be a logic value\n'
        return(list(error=error))}
    if(!is.null(model) & !is.character(model)){
        error <- 'insert the name of the random field\n'
        return(list(error=error))}
    if(is.null(CkModel(model))){
        error <- 'the model name of the random field is not correct\n'
        return(list(error=error))}
    if(is.null(replicates) || (abs(replicates-round(replicates))>0) || replicates<1){
        error <- 'the parameter replicates need to be a positive integer\n'
        return(list(error=error))}
    if(is.null(CheckDistance(distance))){
        error <- 'the name of the distance is not correct\n'
        return(list(error=error))}
    if(is.null(CkCorrModel(corrmodel))){
        error <- 'the name of the correlation model is not correct\n'
        return(list(error=error))}
    if(radius<0){
        error <- 'the radius of the sphere must be positive\n'
        return(list(error=error))}
    if(CkModel(model)==11&&(n<1||!is.numeric(n))){
        error <- 'the parameter n for the Binomial RF is wrong \n'
        return(list(error=error))}

    # START check fitting
    if(fcall=="Fitting"){
         if(!is.null(coordx_dyn))
            {
                if(!is.list(coordx_dyn)){
                error<-"Dynamical coordinates must be a list object"
              return(list(error=error))}
              if(!is.list(data)){
                error<-"Data must be a list"
              return(list(error=error))}
          }
        else {if(missing(data) || !is.numeric(data)){
            error <- 'insert a numeric vector or matrix of data\n'
            return(list(error=error))}}

        if(!is.null(fixed) & !is.list(fixed)){
            error <- 'insert fixed values as a list of parameters\n'
            return(list(error=error))}

        if(!is.null(fixed)){
            namfixed <- names(fixed)
        if(!all(namfixed %in% c(NuisParam(model,CheckBiv(CkCorrModel(corrmodel)),num_betas),CorrelationPar(CkCorrModel(corrmodel))))){
                error <- 'some names of the fixed parameters is/are not correct\n'
                return(list(error=error))}
        if(!CheckParamRange(unlist(fixed))){
            error <- 'some fixed values are out of the range\n'
            return(list(error=error))}}

        if(!is.null(likelihood) & !is.character(likelihood)){
            error <- 'insert the type of likelihood objects\n'
            return(list(error=error))}
        for(i in 1:length(maxdist)){
        if(!is.null(maxdist[i])){
            error <- "insert a positive numeric value for the maximum spatial distance\n"
            if(!is.numeric(maxdist)) return(list(error=error))
            else if(maxdist[i]<0) return(list(error=error))}}

       #  if(!is.null(neighb)&&(!is.null(maxdist[1])||!is.null(maxtime[1])))
       #     {
       #      error <- "maxdist or neighb should be chosen\n"
       #      return(list(error=error))}
        if(!is.null(maxtime)){
            error <- "insert a positive numeric value for the maximum time interval\n"
            if(!is.numeric(maxtime)) return(list(error=error))
            else if(maxtime<0) return(list(error=error))}

        if(!is.null(optimizer) & !is.character(optimizer)){
            error <- 'insert the type of maximising algorithm\n'
            return(list(error=error))}

        if(!is.null(varest) & !is.logical(varest)){
            error <- 'the parameter std.err need to be a logical value\n'
            return(list(error=error))}

        if(type=="Tapering"){
        if(is.null(taper) || is.null(maxdist)){
          error <- 'tapering need a taper correlation model and/or a compact support\n'
          return(list(error=error))}
        if(!taper %in% c("Bohman","Wendland0","Wendland1","Wendland2",
                         "Wendland0_Wendland0","Wendland0_Wendland1","Wendland0_Wendland2",
                         "Wendland1_Wendland0","Wendland1_Wendland1","Wendland1_Wendland2",
                         "Wendland2_Wendland0","Wendland2_Wendland1","Wendland2_Wendland2",
                          "Wendland0_space","Wendland1_space","Wendland2_space",
                          "Wendland0_time","Wendland1_time","Wendland2_time",
                           "tap_st_dyn_wendland1","tap_wen_t","tap_wen_s","Bi_Wendland1",
                         "Bi_Wendland2","Bi_Wendland3","Bi_wendland_asy","Bi_Wendland_asy",
                         "unit_matrix","unit_matrix_biv","unit_matrix_st")){
         
            error <- 'insert a correct name for the taper correlation model\n'
            return(list(error=error))}}

        if(!is.null(type) & !is.character(type)){
            error <- 'insert the configuration of the likelihood objects\n'
            return(list(error=error))}
        if(is.null(CkType(type))){
            error <- 'the type name of the likelihood objects is not correct\n'
            return(list(error=error))}

        if(is.null(CkLikelihood(likelihood))){
           error <- 'the setting name of the likelihood objects is not correct\n'
           return(list(error=error))}

        if(likelihood == "Full"){
            if(!any(type == c("Restricted", "Standard", "Tapering","Tapering1","CV","Tapering2"))){
                error <- 'insert a type name of the likelihood objects compatible with the full likelihood\n'
                return(list(error=error))}}

        if(likelihood == "Marginal"){
            if(!any(type == c("Difference", "Pairwise"))){
                error <- 'insert a type name of the likelihood objects compatible with the composite-likelihood\n'
                return(list(error=error))}}

        if(varest & (likelihood == "Conditional" || likelihood == "Marginal"|| likelihood == "Marginal_2") & (!is.null(vartype) & !is.character(vartype))){
            error <- 'insert the type of estimation method for the variances\n'
            return(list(error=error))}

        if(varest & is.null(CkVarType(vartype)) & (likelihood == "Conditional" || likelihood == "Marginal"|| likelihood == "Marginal_2")){
            error <- 'the name of the estimation type for the variances is not correct\n'
            return(list(error=error))}

        if(!is.null(start)){
            if(!is.list(start)){
                error <- 'insert starting values as a list of parameters\n'
                return(list(error=error))}

       namstart <- names(start)
  if(!bivariate){
       if(num_betas>1){
       if(ncol(X)!=sum(substr(c(namstart,namfixed),1,4)=="mean"))
       {
        error <- 'number of covariates must be equal to to the regressian mean parameters\n'
                return(list(error=error)) }}}
  if(bivariate){


       if(num_betas[1]>1&&num_betas[2]>1){
  
       if(is.list(X)) 
            if(ncol(X[[1]])!=sum(substr(c(namstart,namfixed),1,6)=="mean_1")&&
          ncol(X[[2]])!=sum(substr(c(namstart,namfixed),1,6)=="mean_2"))
       { error <- 'number of covariates must be equal to to the regressian mean parameters\n'
                return(list(error=error)) }}}


                if(!all(namstart %in% c(NuisParam(model,CheckBiv(CkCorrModel(corrmodel)),num_betas), CorrelationPar(CkCorrModel(corrmodel))))){
                error <- 'some names of the starting parameters is/are not correct\n'
                return(list(error=error))}

            if(any(namstart=='mean') & (type=='Difference' || type=='Restricted')){
                error <- 'the mean parameter is not allow with the difference composite likelihood\n'
                return(list(error=error))}
    if(!CheckParamRange(unlist(start))){
                error <- 'some starting values are out of the range\n'
                return(list(error=error))}

             #   if(model %in% c("Gamma","Weibull"))
              #  {
               #     if(!is.na(as.numeric(unlist(start)['sill'])))
                #    { error <- 'sill parameter must be foxed and equal to 1\n'
                #return(list(error=error))}
                #}

            if(!is.null(fixed))
                if(any(namstart %in% namfixed)){
                    error <- 'some fixed parameter name/s is/are matching with starting parameter name/s\n'
                    return(list(error=error))}}

        if(!is.null(weighted) & !is.logical(weighted)){
            error <- 'insert if the composite likelihood need to be weighted'
            return(list(error=error))}

        # START - checks the format of the coordinates and dataset
  if(!is.null(coordx_dyn))
     {AS=1}
  else
    {
    dimdata <- dim(data) # set the data dimension
    if(is.null(coordt)) # START 1) spatial random field
      {   
        if(CheckST(CkCorrModel(corrmodel)))
          {
            error <- 'temporal coordinates are missing\n'
            return(list(error=error))}
            if(grid) # START regular grid
              {   if(CheckBiv(CkCorrModel(corrmodel)))  {
                      if(is.null(dimdata))
                      {
                          error <- c('insert an array of d x d x 2  spatial observations\n')
                          return(list(error=error))}
                  if(length(dimdata)!=3)
                      {
                          error <- c('the dimension of the data matrix is not correct\n')
                          return(list(error=error))}
                   if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2])
                      {
                          error <- c('the number of coordinates does not match with the number of spatial observations\n')
                          return(list(error=error))}}
                  else {                  
                if(is.null(dimdata))
                  {
                    error <- c('insert a matrix d x d of spatial observations\n')
                    return(list(error=error))}
                if(length(dimdata)!=2)
                  {
                    error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))}
                if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2])
                  {
                    error <- c('the number of coordinates does not match with the number of spatial observations\n')
                    return(list(error=error))
                  }}
              } # END regular grid
            else # START irregular grid
              {
                numsite <- length(data)
                if(CheckBiv(CkCorrModel(corrmodel))) numsite<- length(data)/2
                
                if(is.null(numsite))
                  {
                    error <- c('insert a vector of spatial observations\n')
                    return(list(error=error))
                  }
                if(is.null(coordy))
                  {
                    dimcoord <- dim(coordx)
                    if(is.null(dimcoord))
                      { error <- c('insert a suitable set of coordinates\n')
                        return(list(error=error))}
                    else
                      { if(dimcoord[1]!=numsite || dimcoord[2]!=2)
                          {
                            error <- c('the number of coordinates does not match with the number of spatial observations\n')
                            return(list(error=error))
                          }}}
                else
                  { if(length(coordx)!=length(coordy))
                      {error <- c('the number of the two coordinates does not match\n')
                        return(list(error=error))}
                    if(length(coordx)!=numsite)
                      {error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))}
                  }
              }
      } # END 1) spatial random field
    else # 2) case: spatial-temporal random field
      {
        if(!is.numeric(coordt))
          { error <- 'insert a numerical vector of temporal coordinates\n'
            return(list(error=error))}
        if(length(coordt)<=1)
          { error <- 'insert a numerical vector of temporal coordinates\n'
            return(list(error=error))}
    if(grid) # START regular grid
              {
                if(is.null(dimdata))
                  { error <- c('insert an array of d x d x t spatial-temporal observations\n')
                    return(list(error=error))}
                if(length(dimdata)!=3)
                  { error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))}
                if(length(coordx)!=dimdata[1] || length(coordy)!=dimdata[2])
                  { error <- c('the number of coordinates does not match with the number of spatial observations\n')
                    return(list(error=error))}
                if(dimdata[3]!=length(coordt))
                  { error <- c('the time coordinate does not match with the third dimension of the data array\n')
                    return(list(error=error))}
              } # END regular grid
            else # START irregular grid
              {
                if(is.list(coordx_dyn)) {
                 if(length(coordx_dyn)!= length(coordt))
                      { error <- c('the time coordinate does not match with spatial dynamic coordinates number\n')
                    return(list(error=error))}}
                if(is.null(dimdata))
                  { error <- c('insert a matrix of t x d spatial-temporal observations\n')
                    return(list(error=error))}
                if(length(dimdata)!=2)
                  { error <- c('the dimension of the data matrix is not correct\n')
                    return(list(error=error))}
                if(is.null(coordy))
                  {  if(dimdata[2]!=nrow(coordx) || ncol(coordx)!=2)
                      { error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))  }}
                else
                  {
                    if(length(coordx)!=length(coordy))
                      {error <- c('the number of the two spatial coordinates does not match\n')
                        return(list(error=error))}
                    if(length(coordx)!=dimdata[2])
                      { error <- c('the number of coordinates does not match with the number of spatial observations\n')
                        return(list(error=error))}
                  }
                if(dimdata[1]!=length(coordt))
                  { error <- c('the time coordinate does not match with the number of the matrix rows\n')
                    return(list(error=error))}
              } # END irregular grid
      }
      # END - check the format of the inserted coordinates and dataset
      } 
    }
    # END check fitting
    # START check simulation
    if(fcall=="Simulation"){
        if(type=="Tapering")
          {
           if(!taper %in% c("Bohman","Wendland1","Wendland2","Wendland3",
                         "Wendland1_Wendland1","Wendland1_Wendland2","Wendland1_Wendland3",
                         "Wendland2_Wendland1","Wendland2_Wendland2","Wendland2_Wendland3",
                         "Wendland3_Wendland1","Wendland3_Wendland2","Wendland3_Wendland3",
                         "tap_s_dyn_wendland1","tap_s_dyn_wendland1","tap_st_dyn_wendland1","Bi_Wendland1_sep",
                         "Bi_Wendland2_sep","Bi_Wendland3_sep","unit_matrix","unit_matrix_biv","unit_matrix_st")){
            error <- 'insert a correct name for the taper correlation model\n'
            return(list(error=error))}

   
         if(!is.null(coordx_dyn))
            {if(!is.list(coordx_dyn)){
                error<-"Dynamical coordinates mst be a list object"
              return(list(error=error))}}
          if(is.list(coordx_dyn)) {
                 if(length(coordx_dyn)!= length(coordt))
                      { error <- c('the time coordinate does not match with spatial dynamic coordinates number\n')
                    return(list(error=error))}}
            if(!is.null(maxdist)){
            error <- "insert a positive numeric value for the maximum spatial distance\n"
            if(!is.numeric(maxdist)) return(list(error=error))
            else if(maxdist<0) return(list(error=error))}
            if(!is.null(maxtime)){
            error <- "insert a positive numeric value for the maximum temporal distance\n"
            if(!is.numeric(maxtime)) return(list(error=error))
            else if(maxtime<0) return(list(error=error))}
            if(!is.null(tapsep)){
            error <- "separability parameter of spacetime taper must be between 0  and 1\n"
            if(!is.numeric(tapsep)) return(list(error=error))
            else if(tapsep<0||tapsep>1) return(list(error=error))}
          } 

        if(is.null(param) || !is.list(param)){
            error <- 'insert the parameters as a list\n'
            return(list(error=error))}
        biv<-CheckBiv(CkCorrModel(corrmodel))
       #print(length(c(unique(c(NuisParam("Gaussian",biv,num_betas),NuisParam(model,biv,num_betas))),CorrelationPar(CkCorrModel(corrmodel)))))
             if(length(param)!=length(c(unique(c(NuisParam("Gaussian",biv,num_betas),
                    NuisParam(model,biv,num_betas))),
                    CorrelationPar(CkCorrModel(corrmodel)))))
             {
            error <- "some parameters are missing or does not match with the declared model\n"
            return(list(error=error))}

        if(!all( names(param) %in% c(unique(c(NuisParam("Gaussian",biv,num_betas),NuisParam(model,biv,num_betas))),
                                      CorrelationPar(CkCorrModel(corrmodel))))){
            error <- 'some names of the parameters are not correct\n'
            return(list(error=error))}
           
        if(is.list(coordx_dyn)) {
                if(biv) coordt=c(1,2)
                 if(length(coordx_dyn)!= length(coordt)){
                      { error <- c('the number of temporal instants does not match with  the dynamic  spatial coordinates number\n')
                    return(list(error=error))}}    
        }

        if(!CheckParamRange(unlist(param))){
            error <- 'some parameters are out of the range\n'
          return(list(error=error))}
     }# END check simulation
    }
 else{   

    if(missing(coordx)) {
    error <- "spatial locations must be a matrix of dimension 2\n"
    return(list(error=error))}
    else     {
    if(is.vector(coordx)&&!length(coordx)==2)    {
           error <- "spatial locations must be a vector of dimension 2\n"
           return(list(error=error))}
    if(is.matrix(coordx)&&!ncol(coordx)==2)       {
           error <- "spatial locations must be  a matrix of dimension 2\n"
           return(list(error=error))}
    }
   if((!is.null(coordt)&&!is.numeric(coordt))  ){
           error <- "time  must be a vector\n"
           return(list(error=error))}
   if(!type %in% c("simple","ordinary","Simple","Ordinary")){
           error <-"kriging type can be  Simple or Ordinary\n"
   return(list(error=error))}
   if(missing(data) || !is.numeric(data)){
           error <- "insert a numeric vector of data\n"
           return(list(error=error))}
  if(is.list(coordx_dyn)) {
                 if(length(coordx_dyn)!= length(coordt))
                      { error <- c('the time coordinate does not match with spatial dynamic coordinates number\n')
                    return(list(error=error))}}
        }

      # END check kriging
  }

CkLikelihood <- function(likelihood)
  {
    CkLikelihood <- switch(likelihood,
                              None=0,
                              Conditional=1,
                              Full=2,
                              Marginal=3,
                              Marginal_2=1)
    return(CkLikelihood)
  }

CkModel <- function(model)
  {
    CkModel <- switch(model,
                         None=0,
                         Gaussian=1,Gauss=1,
                         Binary=2,
                         Tukeygh=9,
                         SkewGaussian=10,SkewGauss=10,
                         Binomial=11,
                         StudentT=12,
                         Wrapped=13,
                         Geom=14,Geometric=14,
                         PoisBin=15,
                         BinomialNeg=16,
                         PoisBinNeg=17,
                         SkewStudentT=18,
                         Binomial2=19,
                         SinhAsinh=20,
                         Gamma=21,
                         LogGaussian=22,LogGauss=22,
                         Gamma2=23,
                         LogLogistic=24,
                         Logistic=25,
                         Weibull=26,
                         TwoPieceStudentT=27,
                         Beta=28,
                         TwoPieceGaussian=29,TwoPieceGauss=29,
                         Poisson=30,poisson=30,
                         Binomial_TwoPieceGaussian=31,
                         Binomial_TwoPieceGauss=31,
                         BinomialNeg_TwoPieceGaussian=32,
                         BinomialNeg_TwoPieceGauss=32,
                         Kumaraswamy=33,
                         Tukeyh=34,tukeyh=34,
                         Gaussian_misp_StudentT=35,
                         Gaussian_misp_Poisson=36,
                         Gaussian_misp_SkewStudentT=37,
                         TwoPieceTukeyh=38,
                         TwoPieceBimodal=39,
                         Tukeyh2=40,tukeyh2=40,
                         Gaussian_misp_Tukeygh=41,
                         Kumaraswamy2=42,
                         PoissonZIP=43,
                         Gaussian_misp_PoissonZIP=44,
                         BinomialNegZINB=45
                         )
    return(CkModel)
  }

CkVarType <- function(type)
  {
    CkVarType <- switch(type,
                           Sampling=1,
                           SubSamp=2,
                           Theoretical=3)
    return(CkVarType)
  }
  
CkType <- function(type)
  {
    CkType <- switch(type,
                        Difference=1,
                        Pairwise=2,
                        Restricted=3,
                        Standard=4,
                        Tapering=5,
                        Tapering2=5,
                        Tapering1=6,
                        GeoWLS=7,
                        CV=8)
    return(CkType)
  }

#####  names of the correlation models ###############
CorrParam <- function(corrmodel)
{
   return(CorrelationPar(CkCorrModel(corrmodel)))
}
#####  names of the correlation models ###############
#####  names of the correlation models ###############
CorrelationPar <- function(corrmodel)
  {
    param <- NULL  
    if(is.null(corrmodel)){param <- NULL}
    else { 
    # Exponential and Gaussian and spherical and wave correlation :
     if(corrmodel %in% c(2,3,4,16)) {
      param <- c('scale')
      return(param)}
        if(corrmodel %in% c(45)) {
      param <- c('scale','pmu')
      return(param)}
        if(corrmodel %in% c(10)) {
      param <- c('scale_1','scale_2','smooth')
      return(param)}
    # Generalised Cauchy correlation model:
   if(corrmodel %in% c(8,5)) {
      param <- c('power1', 'power2','scale')
      return(param)}
    # Generalised wend correlation model abnd reparametrized version:
     if(corrmodel %in% c(19,6,7)) {
        param <- c('power2', 'scale','smooth')
        return(param)}
    # sine power on sphere 
    if(corrmodel==18){
      param <- c('power')
      return(param)}
     # multiquadric on sphere and stable
   if(corrmodel %in% c(12,17)){
      param <- c('power', 'scale')
      return(param)}     
    # Cauchy and wendx x=0,1,2 
    if(corrmodel %in% c(1,11,13,15)){
        param <- c('power2', 'scale')
        return(param)}
    # Whittle-Matern correlation model:
    if(corrmodel %in% c(14,20)){
      param <- c('scale', 'smooth')
      return(param)}
    # Gneiting or Porcu model:
    if(corrmodel %in% c(42,46,50,60,52,54)) {
      param <- c('power_s', 'power_t','scale_s','scale_t','sep')
      return(param)}
    # Iaco-Cesare model:
    if(corrmodel==44){
      param <- c('power2','power_s', 'power_t','scale_s','scale_t')
      return(param)}
    # Stein model:
    if(corrmodel==48){
      namesparam <- c('power_t','scale_s','scale_t','smooth_s')
      return(param)}

     if(corrmodel==61){
         param <- c('power_s','power2_s','scale_s','scale_t','sep','smooth_t')
         return(param)}  
    if(corrmodel==62){
         param <- c('power_t','power2_t','scale_s','scale_t','sep','smooth_s')
         return(param)}  
    # Wendland dynamic spatio temporal model:
     if(corrmodel %in% c(63,65,67)) {
         param <- c('power_t','power2_s','power2_t','scale_s','scale_t','sep')
         return(param)}  
     if(corrmodel==87){
         param <- c('power_t','power2_s','power2_t','scale_s','scale_t','sep','smooth_s')
         return(param)}   
         # Wendland dynamic spatio temporal model:
     if(corrmodel %in% c(64,66,68)) {
              param <- c('power_s','power2_s','power2_t','scale_s','scale_t','sep')
         return(param)}   
         if(corrmodel==88){
              param <- c('power_s','power2_s','power2_t','scale_s','scale_t','sep','smooth_t')
         return(param)}         
      # Wendland  separable spacetime
     if(corrmodel %in% c(69,70,71,72,73,74,75,77,96))
        {param <- c('power2_s','power2_t','scale_s','scale_t')
         return(param)}        
    
  ##########Separable spatial-temporal correlations##########################################  
    # Exponential-exponential 
    if(corrmodel==84){
      param <- c('scale_s','scale_t')
      return(param)}
    # Matern_Matern
    if(corrmodel==86){
      namesparam <- c('scale_s','scale_t','smooth_s','smooth_t')
      return(param)}  

       # sinpower_st
    if(corrmodel==56){
      namesparam <- c('scale_s','scale_t','smooth_t')
      return(param)}     
    # stable_stable  
     if(corrmodel==94||corrmodel==58){
       param <- c('power_s','power_t','scale_s','scale_t')
       return(param)}
 ###############################################################################################
    ###Multivariate models
    ## biv sep exp and sep wenhole and sep wen
  # if(corrmodel==130||corrmodel==132){
   #   param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol','scale')
   #   return(param)}
   ## biv sep matern
   if(corrmodel==122||corrmodel==119){
      param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol','scale','smooth')
      return(param)}
      ## biv sep wendland
   if(corrmodel==111||corrmodel==113|corrmodel==115){  
      param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol','power2','scale')
      return(param)}   

       if(corrmodel==130){  
      param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol','power2','scale','smooth')
      return(param)} 
       ## biv LMC parsimonious
   if(corrmodel==124){
      param <- c('a_1','a_12','a_2','nugget_1','nugget_2','scale_1','scale_2')
      return(param)}
   ## biv LMC not parsimonious
   if(corrmodel==126){
      param <- c('a_1','a_12','a_2','a_21','nugget_1','nugget_2','scale_1','scale_2')
      return(param)}
      ## biv full endland bivariate models    
   if(corrmodel==112||corrmodel==114||corrmodel==116){
     param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol',
        'power2_1', 'power2_12', 'power2_2','scale_1','scale_12','scale_2')
    return(param)}
       if(corrmodel==132){
     param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol',
        'power2_1', 'power2_12', 'power2_2','scale_1','scale_12','scale_2','smooth_1','smooth_12','smooth_2')
    return(param)}

      ## biv full matern a bivariate models
       if(corrmodel==128||corrmodel==117){
     param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol','scale_1','scale_12','scale_2',
        'smooth_1','smooth_12','smooth_2')
    return(param)}
     ## biv  matern with contrainsts
   if(corrmodel==118||corrmodel==121){
     param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol','scale_1','scale_2','smooth_1','smooth_2')
     return(param)}
      ## biv contr wend bivariate models
      if(corrmodel==129||corrmodel==131||corrmodel==120){
     param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol','power2_1','power2_2','scale_1','scale_2')
     return(param)}

        if(corrmodel==134){
     param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol','power2_1','power2_2','scale_1','scale_2','smooth_1','smooth_2')
     return(param)}
        ## biv asy
     # if(corrmodel==130){
     #param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol','scale','smooth')
     #return(param)}
    if(corrmodel==136){
       param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol','scale_1','scale_12','scale_2',
        'smooth_1','smooth_12','smooth_2','power2_2')
       return(param)
    }
     if(corrmodel==137){
       param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol','scale_1','scale_12','scale_2',
        'smooth_1','smooth_12','smooth_2','power2_1','power2_12','power2_2')
       return(param)
    }
       if(corrmodel==138){
     param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol','scale_1','scale_12','scale_2',
        'power1_1','power1_12','power1_2','power2_1','power2_12','power2_2')
    return(param)}
        if(corrmodel==139){
     param <- c('sill_1','sill_2','nugget_1','nugget_2','pcol','scale_1','scale_12','scale_2',
        'power_1','power_12','power_2')
    return(param)}
    }
    return(param)
  }
  #############################################################  
  #############################################################
NuisParam <- function(model,bivariate=FALSE,num_betas=c(1,1))
{
  param <- NULL
 if(!bivariate&&num_betas==c(1,1)) num_betas=1 
 #if(bivariate) num_betas=c(1,1)
  ############################################################# 
if(!bivariate)      {

  if(num_betas==1) {mm='mean'}
  else {mm='mean' 
        for(i in 1:(num_betas-1)) mm=c(mm,paste("mean",i,sep=""))}

  if( (model %in% c('Gaussian' ,'Gauss' ,'Binomial','Binomial2','BinomialNeg','','Poisson','Gaussian_misp_Poisson',
      'Geom','Geometric','Wrapped','PoisBin','PoisBinNeg','LogGaussian','LogGauss','Logistic')))
  {
    param <- c(mm, 'nugget', 'sill')
    return(param)}

 if( (model %in% c('PoissonZIP','Gaussian_misp_PoissonZIP','BinomialNegZINB')))
  {
    param <- c(mm, 'nugget','pmu','sill')
    return(param)
  }

      if((model %in% c("Weibull","weibull",'Gamma','gamma','LogLogistic',"Loglogistic"))){
      param <- c(mm, 'nugget', 'sill','shape')
      return(param)} 
   
      if((model %in% c('Gamma2','gamma2','Beta','Kumaraswamy','Kumaraswamy2'))) {
      #param <- c(mm, 'nugget', 'sill','shape1','shape2')
      param <- c(mm, 'nugget', 'sill','shape1','shape2','min','max')
      return(param)}     
  # Skew Gaussian univariate random field:
   if((model %in% c('SkewGaussian','SkewGauss','TwoPieceGaussian','TwoPieceGauss',
     'Binomial_TwoPieceGaussian','Binomial_TwoPieceGauss',
     'BinomialNeg_TwoPieceGaussian','BinomialNeg_TwoPieceGauss'))  ){
      param <- c(mm, 'nugget', 'sill','skew')
      return(param)}
    # T univariate ra
  # Skew T univariate random field:
  if((model %in% c('SkewStudentT',"TwoPieceStudentT","Gaussian_misp_SkewStudentT")) ){
      param <- c(mm, 'df','nugget', 'sill','skew')
      return(param)}
  if((model %in% c("TwoPieceBimodal")) ){
      param <- c(mm, 'df','nugget', 'sill','shape','skew')
      return(param)}
    # T univariate random field:
  if((model %in% c('StudentT','Gaussian_misp_StudentT')) ){
      param <- c(mm, 'df','nugget', 'sill')
      return(param)}  
     if( (model %in% c('Tukeyh','tukeyh'))){
      param <- c(mm, 'nugget', 'sill','tail')
      return(param)}
   if( (model %in% c('Tukeyh2','tukeyh2'))){
      param <- c(mm, 'nugget', 'sill','tail1','tail2')
      return(param)}
   # Tukeygh  or SinhAsinhGaussian univariate random field: 
  if( (model %in% c('Tukeygh','SinhAsinh',"TwoPieceTukeyh","Gaussian_misp_Tukeygh"))) {
      param <- c(mm, 'nugget', 'sill','skew','tail')
      return(param)}
}
  #############################################################  
  #############################################################
  ############################################################# 

    if(bivariate)     
   {   
     if(num_betas[1]==1&&num_betas[2]==1) {mm1='mean_1';mm2='mean_2'}
  else {mm1='mean_1';mm2='mean_2' 
        for(i in 1:(num_betas[1]-1)) mm1=c(mm1,paste("mean_1",i,sep=""))
        for(i in 1:(num_betas[2]-1)) mm2=c(mm2,paste("mean_2",i,sep=""))}

   mm=c(mm1,mm2)

   if( model %in% c('Gaussian' ,'Gauss' ,'Binomial','Binomial2','BinomialNeg','Geom','Geometric','Poisson',
        'Wrapped','PoisBin','PoisBinNeg','LogGaussian','LogGauss',"Logistic")){
    param <- mm
    return(param)} 
      # Skew Gaussian bivariate random field:
     if(model %in% c('SkewGaussian','SkewGauss','TwoPieceGaussian','TwoPieceGauss')){
      param <- c(mm,'skew_1','skew_2')
      return(param)}  
      }  
    if((model %in% c('Weibull','Gamma','LogGauss','LogGaussian',"LogLogistic"))){
      param <- c(mm,'shape_1','shape_2')
      return(param)}  
    if((model %in% c('StudentT','Gaussian_misp_StudentT'))) {
      param <- c(mm,'df_1','df_2')
      return(param)}  
     if((model %in% c('Tukeyh','tukeyh'))) {
      param <- c(mm,'tail_1','tail_2')
      return(param)}  
 if((model %in% c('Tukeygh','SinhAsinh',"TwoPieceTukeyh"))) {
      param <- c(mm,'skew_1','skew_2','tail_1','tail_2')
      return(param)}  
    ###################  
  return(param)
}
####################################################################################
#########################################################################################
#########################################################################################
#########################################################################################




StartParam <- function(coordx, coordy, coordt,coordx_dyn, corrmodel, data, distance, fcall, fixed, grid,
                      likelihood,  maxdist, neighb,maxtime, model, n, param, parscale,
                      paramrange, radius, start, taper, tapsep, type,
                      typereal, varest, vartype, weighted, winconst, winstp,winconst_t, winstp_t,X,memdist)
{
    ### START Includes internal functions:
    replicates=1
    # Check if the correlation is bivariate
    CheckBiv <- function(corrmodel)
    {
        CheckBiv <- NULL
        if(corrmodel >= 101 & corrmodel <= 140) CheckBiv <- TRUE
        else CheckBiv <- FALSE
        return(CheckBiv)
    }
    # Check if the correlation is spatial or spatial-temporal
    CheckST <- function(corrmodel)
    {
        CheckST <- NULL
        if(corrmodel >40 & corrmodel <= 100) CheckST <- TRUE
        else  CheckST <- FALSE
        return(CheckST)
    }
    # Check the type of distances
    CheckDistance<- function(distance)
    {
        CheckDistance <- NULL
        CheckDistance <- switch(distance,
                                eucl=0,
                                Eucl=0,
                                chor=1,
                                Chor=1,
                                geod=2,
                                Geod=2,
                                proj=3,
                                Proj=3)
        return(CheckDistance)
    }
    ### END Includes internal functions
    # Set the correlation and  if the correlation is space-time(T or F) or bivariate (T o F)  or univariate (case spacetime=F and bivariate=F)p
    corrmodel<-CkCorrModel(corrmodel)
    bivariate <- CheckBiv(corrmodel); if(bivariate) coordt=c(1,2)
    spacetime <- CheckST(corrmodel)
    isdyn=!is.null(coordx_dyn)

    if(!bivariate)
       {
        if(is.null(X))  {X=1;num_betas=1}
           else 
        {if(is.list(X))  num_betas=ncol(X[[1]])
           else  num_betas=ncol(X) }
    }
    if(bivariate){
        if(is.null(X))  {X=1;num_betas=c(1,1)}
        else
        { if(is.list(X))  num_betas=c(ncol(X[[1]]),ncol(X[[2]]))
            else  num_betas=c(ncol(X),ncol(X)) }}
    namesnuis <- NuisParam(model,bivariate,num_betas)


    ltimes=length(coordt)

    if(grid) { cc=as.matrix(expand.grid(coordx,coordy))
               coordx=cc[,1];coordy=cc[,2]; 
             }


    ### Set returning variables and initialize the model parameters:
    # Initialises the starting and fixed parameters' names
    error <- NULL
    ns<-NULL
    namesfixed <- namesstart <- namessim <- NULL
    numfixed <- numstart <- 0
    # Set the model, likelihood, correlation and the nuisance parameters:
   
    model <- CkModel(model)
    flagnuis <- NULL
    namescorr <- CorrelationPar(corrmodel)
    numparamcorr <- length(namescorr)
    paramcorr <- rep(1, numparamcorr)
    names(paramcorr) <- namescorr
    flagcorr <- NULL
    ### START settings the data structure:
    # set the coordinates sizes:



    if(is.null(coordx_dyn))
    {

      if(is.null(coordy)){coordy <- coordx[,2]
                        coordx <- coordx[,1]}

      numcoord <- numcoordx <- numcoordy <- length(coordx)
      if(bivariate && !is.null(nrow(coordx)) && !is.null(nrow(coordy))) {  #heterotopic case
        numcoordx=nrow(coordx); 
        numcoordy=nrow(coordy);
        numcoord=numcoordy+numcoordx}
      
      ns<-rep(numcoord,ltimes)
    }
    else
    {
       env <- new.env()
       coords=do.call(rbind,args=c(coordx_dyn),envir = env) 

       if(is.list(X))  X=do.call(rbind,args=c(X),envir = env)
       ns=lengths(coordx_dyn)/2 
       coordx <- coords[,1]; coordy <- coords[,2]
       numcoord <- numcoordx <- numcoordy <- length(coordx)
    }


   if((spacetime||bivariate)&&is.null(coordx_dyn)) {coordx=rep(coordx,ltimes);coordy=rep(coordy,ltimes);}
    
    NS=cumsum(ns)
    if(spacetime||bivariate)   NS=c(0,NS)[-(length(ns)+1)]


    # initialize tapering variables:
    tapering=ia=idx=ja=colidx=rowidx=integer(1)
    nozero<-NULL
    tapmodel=0
    cutoff <- FALSE
    distance<-CheckDistance(distance)

    ### END settings the data structure
    # START code for the simulation procedure
    if(fcall=="Fitting"){
        ### Parameters' settings:
        nuisance=nuisance1=nuisance2=NULL
        likelihood <- CkLikelihood(likelihood)
        vartype <- CkVarType(vartype)
        type <- CkType(type)
    
     if((!bivariate&&num_betas==1)||(bivariate&&num_betas==c(1,1)))
     {
        #if(model==1||model==10||model==18||model==9||model==20||model==12||model==13){ 
          if(model %in% c(1,10,12,18,9,20,13,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42)) 
          {
           if(!bivariate) {
                           mu <- mean(unlist(data))
                           if(any(type==c(1, 3, 7,8)))# Checks the type of likelihood
                           if(is.list(fixed)) fixed$mean <- mu# Fixs the mean
                           else fixed <- list(mean=mu)
                           nuisance <- c(mu, 0, var(c(unlist(data))))
                           if(likelihood==2 && (CkType(typereal)==5 || CkType(typereal)==7) ) tapering <- 1
                           if(model %in% c(10,29,31,32))         nuisance <- c(nuisance,0)
                           if(model %in% c(18,20,27,37,38,40,41))      nuisance <- c(0,nuisance,0)
                            if(model %in% c(39))      nuisance <- c(0,0,nuisance,0)
                           if(model %in% c(21,24,12,26,34,35))   nuisance <- c(0,nuisance)
                           #if(model %in% c(23,28,33))  nuisance <- c(0,0,0,nuisance)
                           if(model %in% c(23,28,33,42))  nuisance <- c(0,0,0,nuisance,0,0)
                       }
     if(bivariate) {
                           if(is.null(coordx_dyn)) { mu1 <- mean(data[1,]); mu2 <- mean(data[2,])}
                           else                   { mu1 <- mean(data[[1]]); mu2 <- mean(data[[2]])}
                           if(any(type==c(1, 3, 7, 8)))# Checks the type of likelihood
                           if(is.list(fixed)) {fixed$mean_1 <- mu1;fixed$mean_2<- mu2}
                           else fixed <- list(mean_1=mu1,mean_2=mu2)
                           nuisance <- c(mu1,mu2)
                           if(model %in% c(10,29,31,32))  {nuisance <- c(nuisance,0.1,0.2)}
                           if(model %in% c(26))  {nuisance <- c(nuisance,0.1,0.2)}
                           if(model %in% c(21))  {nuisance <- c(nuisance,0.1)}

                           if(likelihood==2 && (CkType(typereal)==5 || CkType(typereal)==7)) tapering <- 1
                 }
        }
        if(model %in% c(11,14,15,16,19,17,30,45)){
    
            p <- mean(unlist(data)[!is.na(unlist(data))])
            mu=0
            if(model==2||model==11) mu <- qnorm(p/n)
            if(model==14||model==16||model==19) mu <- 0
            if(model==15) mu <- -1
            if(model==17||model==30) mu <- 1
            nuisance <- c(mu, 0, 1)
            if(model==45) {nuisance <- c(mu, 0, 0,1)}
        }
        if(model %in% c(43,44)) nuisance <- c(0, 0, 0, 1)
      }
 #if(num_betas>1)
 if((!bivariate&&num_betas>1)||(bivariate&&num_betas[1]>1&&num_betas[2]>1) )
     {
        if(model %in% c(1,10,12,18,9,20,13,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42)) {
    if(!bivariate) {
         if(any(type==c(1, 3, 7,8)))# Checks the type of likelihood
            if(is.list(fixed)) {
                               mu <- mean(unlist(data));fixed$mean <- mu# Fixs the mean
                               for(i in 1:(num_betas-1)) fixed[[paste("mean",i,sep="")]]=1  # fixed$meani=1
                           }
            else  fixed <- list(mean=mu)
            for(i in 1:num_betas) nuisance=c(nuisance,1);
            nuisance=c(nuisance,0,var(c(unlist(data))))
             if(model %in% c(10,29,31,32))        nuisance=c(nuisance,1)  
             if(model %in% c(21,24,12,26,34,35))  nuisance=c(nuisance,1) 
             if(model %in% c(18,20,27,37,38,40,41))     nuisance=c(1,nuisance,1) 
             if(model %in% c(39))     nuisance=c(1,1,nuisance,1) 
          #  if(model %in% c(23,28,33))         nuisance=c(nuisance,1,1,1)  
            if(model %in% c(23,28,33,42))         nuisance=c(nuisance,1,1,1,1,1)  
             }
    if(bivariate) {
            if(any(type==c(1, 3, 7,8)))# Checks the type of likelihood
            if(is.list(fixed)) {
                if(!is.list(data))
                {
                                mu1 <- rowMeans(unlist(data))[1];fixed$mean_1 <- mu1# Fixs the mean
                                mu2 <- rowMeans(unlist(data))[2];fixed$mean_2 <- mu2# Fixs the mean
                                for(i in 1:(num_betas[1]-1)) fixed[[paste("mean_1",i,sep="")]]=1 
                                for(i in 1:(num_betas[2]-1)) fixed[[paste("mean_2",i,sep="")]]=1  # fixed$meani=1
                }
                  if(is.list(data))
                {
                                mu1 <- mean(data[[1]]);fixed$mean_1 <- mu1# Fixs the mean
                                mu2 <- mean(data[[2]]);fixed$mean_2 <- mu2# Fixs the mean
                                for(i in 1:(num_betas[1]-1)) fixed[[paste("mean_1",i,sep="")]]=1 
                                for(i in 1:(num_betas[2]-1)) fixed[[paste("mean_2",i,sep="")]]=1  # fixed$meani=1
                }
                           }
            else  fixed <- list(mean_1=mu1,mean_2=mu2)
            for(i in 1:num_betas[1]) nuisance1=c(nuisance1,1);
            for(i in 1:num_betas[2]) nuisance2=c(nuisance2,1);
            nuisance=c(nuisance1,nuisance2 )
             if(model %in% c(10,29,31,32))        nuisance=c(nuisance,1,1)  
             if(model %in% c(21,24,12,26,34,35))  nuisance=c(nuisance,1,1) 
             if(model %in% c(18,20,27,37,38,40,41))     nuisance=c(1,nuisance,1) 
              if(model %in% c(39))     nuisance=c(1,1,nuisance,1) 
            #if(model %in% c(23,28,33))         nuisance=c(nuisance,1,1,1)
            if(model %in% c(23,28,33,42))         nuisance=c(nuisance,1,1,1,1,1)  

            }
         }
     if(model %in% c(2,11,14,15,16,19,17,30)) nuisance <- c(0,rep(1,num_betas-1) ,0, 1)
     if(model %in% c(45)) nuisance <- c(0,rep(1,num_betas-1) ,0,0, 1)
     if(model %in% c(43,44)) nuisance <- c(0,rep(1,num_betas-1) ,0, 0,1)

     }
        # Update the parameter vector      
        names(nuisance) <- namesnuis
        namesparam <- sort(c(namescorr, namesnuis))
        param <- c(nuisance, paramcorr)
        param <- param[namesparam]
        numparam <- length(param)
        flag <- rep(1, numparam)
        namesflag <- namesparam
        names(flag) <- namesflag
        # Update the parameters with fixed values:
        if(!is.null(fixed)){
            fixed <- unlist(fixed)
            namesfixed <- names(fixed)
            numfixed <- length(namesfixed)
            if(numfixed==numparam){ error <- 'there are not parameters left to estimate\n';return(list(error=error))}
            flag[pmatch(namesfixed, namesflag)] <- 0
            param <- param[-pmatch(namesfixed, namesparam)]
            numparamcorr <- numparamcorr-sum(namesfixed %in% namescorr)
            namesparam <- names(param)
            numparam <- length(param)
        }
        flagcorr <- flag[namescorr]
        flagnuis <- flag[namesnuis]
        # Update the parameters with starting values:
        if(!is.null(start)){
            start <- unlist(start)
            namesstart <- names(start)
            if(any(type == c(1, 3, 7))){
                if(!bivariate) {   # univariate case
                       if(any(model==c(1,10,12,18,20,9,13,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42)))
                       if(any(namesstart == 'mean'))  start <- start[!namesstart == 'mean']
                       if(num_betas>1)
                       for(i in 1:(num_betas-1)) {  if(any(namesstart == paste("mean",i,sep="")))  {namesstart <- names(start) ; 
                       if(any(model==c(1,10,12,18,20,9,13,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42)))
                                                 start <- start[!namesstart == paste("mean",i,sep="")]}}
                }
                if(bivariate) {          
                                  if(any(namesstart == 'mean_1'))  start <- start[!namesstart == 'mean_1']        
                                  if(any(namesstart == 'mean_2'))  {namesstart <- names(start) ; start <- start[!namesstart == 'mean_2']}
                      
                       if(num_betas[1]>1)
                       for(i in 1:(num_betas[1]-1)) {  if(any(namesstart == paste("mean_1",i,sep="")))  {namesstart <- names(start) ; 
                       if(any(model==c(1,10,12,18,20,9,13,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42)))
                                                 start <- start[!namesstart == paste("mean_1",i,sep="")]}}            
                       if(num_betas[2]>1)
                       for(i in 1:(num_betas[2]-1)) {  if(any(namesstart == paste("mean_2",i,sep="")))  {namesstart <- names(start) ; 
                       if(any(model==c(1,10,12,18,20,9,13,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42)))
                                                 start <- start[!namesstart == paste("mean_2",i,sep="")]}}  
                                  }
                }
            namesstart <- names(start)
            numstart <- length(start)
            param[pmatch(namesstart,namesparam)] <- start
            }
        ### set the scale of the parameters:
        # Insert here!
        # set the range of the parameters if its the case
        paramrange=TRUE
 
        #if(paramrange) paramrange <- SetRangeParam(namesparam, numparam)
        #else 
        paramrange <- list(lower=NULL, upper=NULL)
        ### If the case set the sub-sampling parameters to the default values
        if(is.null(winconst) || !is.numeric(winconst)) winconst <- 0
        if(is.null(winstp) || !is.numeric(winstp)) winstp <- 0
        if(is.null(winconst_t) || !is.numeric(winconst_t)) winconst_t <- 0
        if(is.null(winstp_t) || !is.numeric(winstp_t)) winstp_t <- 0
        ### Set the data format:
        if(spacetime||bivariate){ # setting spam indexes
            if(spacetime) numtime <- ltimes
            if(bivariate) numtime <- 2

                if(typereal=="Tapering"||typereal=="Tapering1"||typereal=="Tapering2"){
                tapering<-1
                idx<-integer((numcoord*numtime)^2)
                ja<-integer((numcoord*numtime)^2)
                ia<-integer(numcoord*numtime+1)
                tapmodel<-CkCorrModel(taper)
                              }}
        else{              #    setting spam indexes
            numtime <- 1
            coordt <- 0
            data <- matrix(data, ncol=numcoord, nrow=replicates)
            if(typereal=="Tapering"||typereal=="Tapering1"||typereal=="Tapering2"){
                tapering<-1
                
                idx<-integer((numcoord*numtime)^2)
                ja<-integer((numcoord*numtime)^2)
                ia<-integer(numcoord*numtime+1)
                tapmodel<-CkCorrModel(taper)
                }}
    }
    # END code for the fitting procedure
    # START code for the simulation procedure
    if(fcall=="Simulation"){
       
        namesnuis <- sort(unique(c(namesnuis,NuisParam("Gaussian",bivariate,num_betas))))
        param <- unlist(param)
        numparam <- length(param)
        namesparam <- names(param)

        if(!bivariate) namessim <- c("mean","sill","nugget","scale",namescorr[!namescorr=="scale"])
        if(bivariate)  namessim <- c("mean_1","mean_2","scale",
                             namescorr[!namescorr=="scale"])  

        if(spacetime) numtime <- ltimes
        else {numtime <- 1; coordt <- 0}
        if(bivariate) numtime <- 2

        if((typereal=="Tapering"&&type=="Tapering")||(typereal=="Tapering1"&&type=="Tapering1")||(typereal=="Tapering2"&&type=="Tapering2")){
                tapering<-1
                idx<-integer((numcoord*numtime)^2)
                ja<-integer((numcoord*numtime)^2)
                ia<-integer(numcoord*numtime+1)
                tapmodel<-CkCorrModel(taper)
                }
        }
    # END code for the simulation procedure

    ### Compute the spatial and spatial-temporal distances:
    numpairs <- integer(1)
    srange <- double(1)
    trange <- double(1)
    if(is.null(maxdist)) srange<-c(srange,double(1)) else {srange<-c(srange,as.double(maxdist))}                # cutoff<-TRUE
    if(is.null(maxtime)) trange<-c(trange,double(1)) else {trange<-c(trange,as.double(maxtime))}                # cutoff<-TRUE
    isinit <- as.integer(1)
    if(is.null(tapsep))  tapsep=c(0.5,0.5)
    else  {if(length(tapsep)==1) tapsep=c(tapsep,0)}
    mem=FALSE
    if(tapering||memdist)  { mem=TRUE }   #### NB
    if(mem&&!tapering)  
      {        
                nn=numcoord*numtime
                if(spacetime&&isdyn)  nn=sum(ns)
                colidx<-rowidx<-integer(nn*(nn-1)/2)
      }
    if(bivariate) {
      if(!srange[1]&&!srange[2])  srange=c(srange,0,0)
      if(is.na(srange[3])) srange[3]=srange[2];
                  if(is.na(srange[4])) srange[4]=srange[2];}

    if(CheckSph(corrmodel))   radius=1
    aa=double(5);for(i in 1:length(tapsep)) aa[i]=tapsep[i];tapsep=aa

# gb=.C('SetGlobalVar',as.integer(bivariate),as.double(coordx),as.double(coordy),as.double(coordt),
#           as.integer(grid),ia=ia,idx=idx,
#           isinit=isinit,ja=ja,as.integer(mem),as.integer(numcoord),as.integer(numcoordx), as.integer(numcoordy),
#           numpairs=numpairs,as.double(radius),srange, as.double(tapsep), as.integer(spacetime),
#           as.integer(numtime),trange,as.integer(tapering),as.integer(tapmodel),
#           as.integer(distance),as.integer(weighted),
#           colidx=as.integer(colidx),rowidx=as.integer(rowidx),
#           as.integer(ns),as.integer(NS),as.integer(isdyn),
#           PACKAGE='GeoModels', DUP=TRUE, NAOK=TRUE)


if(is.null(neighb)){
 gb=dotCall64::.C64('SetGlobalVar',SIGNATURE = c(
         "integer","double","double","double","integer", "integer","integer",  #7
         "integer","integer","integer","integer", "integer","integer", #6
         "integer","double","double","double", "integer",  #5
         "integer","double", "integer","integer","integer","integer", #6
         "integer","integer", # 2
         "integer","integer","integer"),  # 3
     bivariate, coordx, coordy, coordt,grid,ia=ia,idx=idx,  #7
           isinit=isinit,ja=ja, mem, numcoord, numcoordx,  numcoordy, #6
           numpairs=numpairs, radius,srange,  tapsep,  spacetime, #5
            numtime,trange, tapering, tapmodel,distance, weighted, #6
           colidx= colidx,rowidx= rowidx, # 2
            ns, NS, isdyn, #3
 INTENT = c("r","r","r","r","r","w","w", #7
            "rw","w", "rw", "r", "r", "r", #6
           "rw", "r", "rw", "r", "r", #5
             "r",  "rw", "r", "r", "r", "r", #6
             "w", "w",#2
             "r", "r", "r"),
             PACKAGE='GeoModels', VERBOSE = 0, NAOK = TRUE)

  ##
rm(colidx);rm(rowidx)
if(type=="Tapering") {rm(idx);rm(ja);rm(ia)}
##
## number  of selected pairs
numpairs <- gb$numpairs
## indexes for composite 
    colidx=gb$colidx
    rowidx=gb$rowidx
    colidx <- colidx[1:numpairs]
    rowidx  <- rowidx[1:numpairs]
    idx <- gb$idx
    ja <- gb$ja
    ia <- gb$ia;
    isinit <- gb$isinit
    nozero <- numpairs/(numcoord*numtime)^2
    idx <- idx[1:numpairs]
    ja  <- ja[1:numpairs]
}
else   ######## case  with neighboord!!!
{ 
if(!spacetime&&!bivariate)   #  spatial case
   {
          ######
          indices <- function(X,Y)
          {
             res = NULL;res_d = NULL
             for(i in 2:ncol(X))
             {
                sol = cbind(X[,1],X[,i])
                res = rbind(res,sol)
                sol_d = cbind(Y[,1],Y[,i])
                res_d = rbind(res_d,sol_d)
             }
            xx=as.numeric(res[,1]); yy=as.numeric(res[,2])
            sol=xx+yy+xx*yy+xx^2+yy^2
            ids <- !duplicated(sol)
            return(list(xy = res[ids,],d = res_d[ids,][,2]))
         }
     ##########################################
         nn2Geo <- function(x, K = 2,distance)
         {
            nearest = RANN::nn2(x,k = K)
            #nearest = FNN::get.knn(k, k=K)
            #########  cases geod (2) or chordal (1) distances
            if(distance==2||distance==1){
                  agc=NULL;nnn=nrow(x)
                  for(i in 1:nnn){
                  a=fields::rdist.earth(matrix(x[i,],ncol=2), x[nearest$nn.idx[i,2:K],], miles = FALSE, R = 1)
                  agc=rbind(agc,a)
                 }
             if(distance==2)  agc=radius*agc   # geodesic
             if(distance==1)  agc=2*radius*sin(0.5*agc)   # chordal  
             nearest$nn.dists=cbind(rep(0,nnn),agc)
             }
            ########################################### 
            sol = indices(nearest$nn.idx,nearest$nn.dists)
            lags <- sol$d;rowidx <- sol$xy[,1];colidx <- sol$xy[,2]
         return(list (lags=lags, rowidx = rowidx, colidx = colidx))
         }
     ##########################################
  K=neighb
  x=cbind(coordx, coordy)
  sol = nn2Geo(x,K,distance)
  nn = length(sol$lags)
  sol$lagt=0
  gb=list(); gb$colidx=sol$colidx;
             gb$rowidx=sol$rowidx ;
             gb$numpairs=n
  ss=.C("SetGlobalVar2", as.integer(numcoord),  as.integer(numtime),  
    as.double(sol$lags),as.integer(nn),
    as.double(sol$lagt),as.integer(nn),
    as.integer(spacetime),as.integer(bivariate)) 
    ## number  of selected pairs
    numpairs <- gb$numpairs
## indexes for composite 
    colidx=gb$rowidx 
    rowidx=gb$colidx
    idx <- 0;ja <- 0;ia <- 0
    isinit <- 1
    nozero <- numpairs/(numcoord*numtime)^2
    idx <- 0;ja  <- 0
   } #### end spatial case 


}

if(is.null(coordt)) coordt=1
  
    ### Returned list of objects:
    return(list(bivariate=bivariate,coordx=coordx,coordy=coordy,coordt=coordt,corrmodel=corrmodel,
                colidx = colidx ,rowidx=rowidx,
                data=data,distance=distance,
                error=error,flagcorr=flagcorr,flagnuis=flagnuis,fixed=fixed,likelihood=likelihood,
                lower=paramrange$lower,model=model,n=n,namescorr=namescorr,namesfixed=namesfixed,
                namesnuis=namesnuis,namesparam=namesparam,namessim=namessim,namesstart=namesstart,ns=ns,NS=NS,
                num_betas=num_betas,
                numcoord=numcoord,numcoordx=numcoordx,numcoordy=numcoordy,
                numfixed=numfixed,numpairs=numpairs,numparam=numparam,numparamcorr=numparamcorr,
                numstart=numstart,numtime=numtime,param=param,setup=list(                ## setup is a list
                ia=ia,idx=idx,ja=ja,nozero=nozero,tapmodel=tapmodel,tapsep=tapsep),  radius=radius,                            ## with tapered matrix informations
                spacetime=spacetime,srange=srange,start=start,upper=paramrange$upper,type=type,
                trange=trange,vartype=vartype,weighted=weighted,winconst=winconst,winstp=winstp,
                winconst_t=winconst_t,winstp_t=winstp_t,X=X))
}

DeviceInfo <- function()
{
    .C("DeviceInfo",PACKAGE='GeoModels',DUP = TRUE, NAOK=TRUE)
}
