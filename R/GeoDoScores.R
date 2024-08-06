####################################################
### File name: GeoDoScores.r
####################################################

GeoDoScores=function(data,method="cholesky",matrix)   {


########## internal function ############################
getInv2=function(covmatrix,b){

if(!covmatrix$sparse){
               U =MatDecomp(covmatrix$covmatrix,method);Inv=0
               if(is.logical(U)){print(" Covariance matrix is not positive definite");stop()}
               vec=forwardsolve(U, b)
               Invc=forwardsolve(U, vec,transpose=T) ## R^-1 %*% c
               Inv=FastGP::rcppeigen_invert_matrix(covmatrix$covmatrix)
               diagInv=diag(Inv)

             }
if(covmatrix$sparse){
               cc=covmatrix$covmatrix
               if(spam::is.spam(cc))  U = try(spam::chol.spam(cc),silent=TRUE)
               else                    U = try(spam::chol.spam(spam::as.spam(cc)),silent=TRUE)
               if(inherits(U,"try-error")) {print(" Covariance matrix is not positive definite");stop()}#Inv=spam::chol2inv.spam(U)
             
               vec=  spam::forwardsolve(U, b)
               Invc= spam::backsolve(U, vec) ## R^-1 %*% c
               Inv= spam::solve.spam(U)
               diagInv=diag(Inv)

                #backsolve(U, forwardsolve(U, b))
             
        }
     return(list(a=Invc,b=Inv,c=diagInv))
}
################################################
if(!inherits(matrix,"GeoCovmatrix"))  stop("A GeoCovmatrix object is needed as input\n")

varcov=matrix$covmatrix
rownames(varcov)=c();colnames(varcov)=c()
if(nrow(varcov)!=length(data)) stop("The dimension of the covariance  matrix and/or the vector data are not correct  \n")
data=c(unname(data))


###########################
nsites=length(matrix$coordx)
ntime=1
if(matrix$spacetime) ntime=length(matrix$coordt)
if(matrix$bivariate) ntime=2
dime = nsites*ntime


MM=0
param=matrix$param


namesnuis=matrix$namesnuis
nuisance <- param[namesnuis]
if(length(param$mean)==1){
 sel=substr(names(nuisance),1,4)=="mean"
 mm=as.numeric(nuisance[sel])
 MM=(matrix$X)%*%mm
}
else {MM=param$mean}
########
if(matrix$model %in% c(1,35,12,34,25))  data=data-MM#Gaussian #StudentT Tukey  Logistic
if(matrix$model %in% c(10))             #SkewGauussian
                 { kk=param['skew'];
                   data=data-(MM+kk*sqrt(2/pi))}
if(matrix$model %in% c(11))     data=data-(matrix$n)*pnorm(MM) #binomial
if(matrix$model %in% c(30,36))  data=data-exp(MM) #poisson
if(matrix$model %in% c(43,44))  {p=pnorm(param['pmu']);data=data-(1-p)*exp(MM)} #poisson inlated
if(matrix$model %in% c(27)) #two piece t models
     {kk=param['skew'];dd=param['df'];ss=param['sill'];
      data=data-(MM-(2*kk*sqrt(ss*dd)*gamma((dd-1)/2))/(gamma(dd/2)*sqrt(pi)))
     }
if(matrix$model %in% c(29)) #two piece gaussian
     {kk=param['skew'];ss=param['sill'];
      data=data-(MM-(2*kk*sqrt(2*ss/pi)))
     }
if(matrix$model %in% c(38)) #two piece tukeyh
     {kk=param['skew'];ss=param['sill'];tt=param['tail'];
      data=data-(MM-(2*kk*sqrt(2*ss/pi)/(1-tt)))
     }
if(matrix$model %in% c(40))      #tukeyh2
     {ss=param['sill'];t1=param['tail1'];t2=param['tail2'];
      data=data-(MM+sqrt(ss)*(t1-t2)/(sqrt(2*pi)*(1-t1)*(1-t2)))
     }
if(matrix$model %in% c(16))     data=data-(matrix$n)*(1-pnorm(MM))/pnorm(MM) #binomialnegative
if(matrix$model %in% c(45))     data=data-(1-pnorm(param['pmu']))*(matrix$n)*(1-pnorm(MM))/pnorm(MM) #binomialnegative inflated
if(matrix$model %in% c(20))      #sas
     {ss=param['sill'];kk=param['skew'];tt=param['tail']; 
      data=data-(MM+sqrt(ss)*sinh(kk/tt)*exp(0.25)*(besselK(.25,(tt+1)/(2*tt))+besselK(.25,(1-tt)/(2*tt)))/(sqrt(8*pi)))
     }
if(matrix$model %in% c(39))      #twopiecebimodal
     {ss=param['sill'];sk=param['skew'];nu=param['df'];delta=param['shape'];alpha=2*(delta+1)/nu;nn=2^(1-alpha/2)
      data=data-(MM-sqrt(ss)*sk*2^(1/alpha+1)*gamma(nu/2+1/alpha)/(nn^(1/alpha)*gamma(nu*0.5)))
     }
#######


cc=getInv2(matrix,data)


temp=cc$a# inv%*%data
inv=cc$b   #inv
vv=cc$c    # diag inv

D=diag(1/vv,dime,dime)
DD=diag(sqrt(1/vv),dime,dime)

z=crossprod(D,temp)
zz=crossprod(DD,temp)



MAD=median(z)
RMSE=sqrt((1/dime)*crossprod(z,z))
#RMSE=sqrt((1/dime)*Rfast::Crossprod(z,z))
MAE=(1/dime)*(sum(abs(z)))


LSCORE=(1/(2*dime))*(sum(log(2*pi/vv))+sum(zz^2))
CRPS=(1/dime)*(sum((1/vv)^0.5*zz*(2*pnorm(zz)-1))+2*sum((1/vv)^0.5*pnorm(zz))+sum((1/vv)^0.5)/sqrt(pi))

###########################
scores = list(RMSE = RMSE,
               LSCORE = LSCORE,
               MAD=MAD,
               CRPS = CRPS,
               MAE=MAE)
return(scores)
}

###################################################################################################
###################################################################################################
