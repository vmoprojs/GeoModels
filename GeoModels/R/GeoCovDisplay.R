####################################################
### File name: GeoCovDisplay.r
####################################################

GeoCovDisplay=function(covmatrix,limits=FALSE,pch=2)   {
#if(class(covmatrix)!="CovMat") stop("A CovMat object is needed as input\n")
if(!inherits(covmatrix,"GeoCovmatrix")) stop("A GeoCovmatrix object is needed as input\n")
covmat=as.matrix(covmatrix$covmatrix)
#if(!is.matrix(covmat)) stop("The function needs a covariance matrix as input\n")
opar=par(no.readonly = TRUE)

nrw=nrow(covmat)
ncl=ncol(covmat)
aa=1*(abs(covmat)>1e-150)
if(sum(aa)==nrw*ncl) return(cat("No image plot: The covariance matrix is dense")) 
image(t(apply(aa, 2, rev)),col=gray.colors(2,start =1, end = 0),pch=2,xaxt='n',yaxt='n')
if(limits){
if(covmatrix$bivariate) {abline(h=0.5);
	                     abline(v=0.5)
	                    }
if(covmatrix$spacetime) {abline(h=seq(0,1,covmatrix$numtime));
                        abline(v=seq(0,1,covmatrix$numtime))}
}
 par(opar)
}
