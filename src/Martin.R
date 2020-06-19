


###############################################################
######### Examples of  spatial bivariate RFs cauchy model ###########
###############################################################
rm(list=ls())
require(GeoModels)

# Define the spatial-coordinates of the points:
set.seed(5)
x <- runif(150, 0, 1)
y <- runif(150, 0, 1)
coords=cbind(x,y)

#########	
mean_1=0; mean_2=0
nugget_1=0; nugget_2=0
############
corrmodel="Bi_GenCauchy"

sill_1=1        ## variance 1
sill_2=2        ## variance 2

scale_1=0.01    # scale 1
scale_2=0.02    # scale 2
scale_12=0.5*(scale_1+scale_2)   # cross scale

power1_1=1.3
power1_12=1.6
power1_2=1.5

power2_1=4
power2_2=3
power2_12=3.5

pcol=0.2

# parameters
param=list(mean_1=mean_1,mean_2=mean_2,
          scale_1=scale_1,scale_2=scale_2,scale_12=scale_12,
          sill_1=sill_1,sill_2=sill_2,nugget_1=nugget_1,nugget_2=nugget_2,
          pcol=pcol,
          power1_1=power1_1,power1_12=power1_12,power1_2=power1_2,
          power2_1=power2_1,power2_12=power2_12,power2_2=power2_2)

#### covariance matrix


data <- GeoSim(coordx=coords, corrmodel=corrmodel, method="cholesky", # svd
              param=param)$data

data[1,]  # fist component
data[2,]  # second component

#CC <- GeoCovmatrix(coordx=coords, corrmodel=corrmodel, 
#              param=param)$covmatrix
#min(eigen(CC)$values)






