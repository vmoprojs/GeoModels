GeoScores=function( data_to_pred, probject=NULL,pred=NULL,mse=NULL,
   score=c("brie","crps","lscore","pit","pe"))

{

if(is.null(probject))
{if(is.null(pred)&&is.null(mse)) stop("Geokrig object or (pred and mse) are mandatory")}   


data_to_pred=c(data_to_pred)
if(inherits(probject,"GeoKrig")||inherits(probject,"GeoKrigloc")) 
{pred=c(probject$pred);mse=c(probject$mse)}

 stopifnot(is.numeric(pred))
 stopifnot(is.numeric(data_to_pred))
 if(!is.null(mse)) stopifnot(is.numeric(mse))
 N1= length(pred)
 N2= length(data_to_pred)
 if(N1!=N2) stop("length of data and predictions does not match\n")



 if(!is.null(mse)) sqrtvv=sqrt(mse)
rmse=mae=mad=lscore=crps=pit=brie=NULL
 

if(!is.null(mse)){
if(sum(grepl("pit", score))==1) 
      if(!is.null(probject)){ probject$data_to_pred=data_to_pred
                              pit = GeoPit(probject,type="Uniform")$data} 
      else                  { pit = pnorm(data_to_pred, mean = pred, sd = sqrtvv)}
if(sum(grepl("brie", score))==1) 
brie = mean(pnorm(data_to_pred, mean = pred, sd = sqrtvv)-1*I(pred<data_to_pred))
}


if(sum(grepl("pe", score))==1 )
{  err=data_to_pred-pred
   rmse=sqrt(sum(err^2)/N2)
   mae=mean(abs(err)) 
   mad=median(abs(err))
}

if(!is.null(mse)){
if(sum(grepl("lscore", score))==1||sum(grepl("crps", score))==1){
err=data_to_pred-pred ;std=err/sqrtvv
if(sum(grepl("lscore", score))==1)  lscore=0.5*sum(std^2+log(2*pi*sqrtvv))/N2 
if(sum(grepl("crps", score))==1)    crps=sum( sqrtvv*( std*(2*pnorm(std)-1 ) +2*pnorm(std)-1/sqrt(pi)))/N2
}
}

a=list(brie=brie,rmse=rmse,mae=mae,mad=mad,lscore=lscore,crps=crps,pit=pit)
a=a[!sapply(a,is.null)]
return(a)
}