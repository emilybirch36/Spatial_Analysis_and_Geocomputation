library(rgdal)
library(car)
library(tmap)
library(ggplot2)
library(tidyr)
library(sp)
library(spdep)
library(knitr)
library(tmaptools)
library(sp)
library(gstat)
library(GSIF)
library(rgeos)
library(spgwr)
library(spatialreg)
tmap_mode("plot")

# not run
list.files(path = ".", pattern = NULL, all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

getwd()


fgdb <- "/Users/emilybirch/Documents/UCL/T1_Spatial _Analysis_Geocomputation/Germany_covid_shapefiles/germany/Final_Germany_files/final_file"

subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(fgdb)
print(fc_list)


# Get feature layer to process
targetfile <- readOGR(dsn = fgdb, layer = "Germany_final")

View(targetfile@data)

# Map covid deaths
tm_shape(targetfile) + tm_polygons(col = "covid_deat", palette = "YlOrRd", style = "quantile")


targetfile@data[9] <- NULL

targetfile@data[9] <- NULL


targetfile@data$SO2_mean[is.na(targetfile@data$SO2_mean)] <- mean(targetfile@data$SO2_mean)
targetfile@data$pm2.5_mean[is.na(targetfile@data$pm2.5_mean)] <- mean(targetfile@data$pm2.5_mean)
targetfile@data$alcohol_ca[is.na(targetfile@data$alcohol_ca)] <- mean(targetfile@data$alcohol_ca)
targetfile@data$dist_borde[is.na(targetfile@data$dist_borde)] <- mean(targetfile@data$dist_borde)
targetfile@data$income_25k[is.na(targetfile@data$income_25k)] <- mean(targetfile@data$income_25k)
targetfile@data$road_len[is.na(targetfile@data$road_len)] <- mean(targetfile@data$road_len)

targetfile@data$covid_deat[is.na(targetfile@data$covid_deat)] <- 0

targetfile@data$dist_borde[targetfile@data$dist_borde == 0] <- 50


targetfile@data$covid_deat <- as.numeric(targetfile@data$covid_deat)


# Run this only if you want to try removing the rows with missing data from targetfile@data
# If you change margin from 1 to 2, it removes the column instead, thus you can work only with variables without missing data
# https://gis.stackexchange.com/questions/89512/r-dealing-with-missing-data-in-spatialpolygondataframes-for-moran-test
coordinates(targetfile) <- ~x+y
# DISPLAY NA ROWS IN targetfile  
targetfile@data[!complete.cases(targetfile@data),] 

# FUNCTION TO REMOVE NA's IN sp DataFrame OBJECT
#   x           sp spatial DataFrame object
#   margin      Remove rows (1) or columns (2) 
sp.na.omit <- function(x, margin=1) {
  if (!inherits(x, "SpatialPointsDataFrame") & !inherits(x, "SpatialPolygonsDataFrame")) 
    stop("MUST BE sp SpatialPointsDataFrame OR SpatialPolygonsDataFrame CLASS OBJECT") 
  na.index <- unique(as.data.frame(which(is.na(x@data),arr.ind=TRUE))[,margin])
  if(margin == 1) {  
    cat("DELETING ROWS: ", na.index, "\n") 
    return( x[-na.index,]  ) 
  }
  if(margin == 2) {  
    cat("DELETING COLUMNS: ", na.index, "\n") 
    return( x[,-na.index]  ) 
  }
}


# DELETE NA's IN targetfile AND SHOW CHANGE IN dim
targetfile2 <- sp.na.omit(targetfile)     
dim(targetfile)
dim(targetfile2) 

# PLOT DELETED POINTS IN RED    
plot(targetfile, col="red", pch=20)
plot(targetfile2, col="black", pch=20, add=TRUE)


#Moran's I
Wl <- nb2listw(nb) # a listw object is a weights list for use in autocorrelation measures.
moran(file$norm_death, Wl, n=length(Wl$neighbours), S0=Szero(Wl)) 


 
### Ordinary regression model ###
# Extract and map residuals (not working, probably due to missing data: inc, smo susceptible)

model <- glm(covid_deat ~ log(dist_borde) + pm2.5_mean + road_len + SO2_mean + log(alcohol_ca) + X_no2mean, data = targetfile2@data, family= "poisson")



### Get VIF of model, then plot it ###
vif_values <- vif(model)
vif_values
barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue")

# extract map residuals
targetfile2$modelres <- residuals(model)
tm_shape(targetfile2)+tm_polygons("modelres", palette = "-RdBu", style = "quantile")


# Plot residuals
ggplot(data = targetfile2@data, aes(modelres)) + geom_histogram()

ggplot(data = targetfile2@data, aes(sample = modelres)) + geom_qq() + geom_qq_line()


# Spatial weight matrix 
model.W = nb2listw(poly2nb(targetfile2), zero.policy = TRUE, na.action(na.exclude))

# to determine whether a spatial regression model is needed
# first do Moran's I for residual autocorrelation, 
# then do Lagrange multiplier test for spatial lag dependence
# then do lagrange multiplier test for spatial error dependence
# Moran test, use to test for residual autocorrelation. Limited bc it doesn't reveal the type of lag/error
lm.morantest(model, model.W, zero.policy = TRUE, na.action(na.exclude))


#LM tests (lag and error)
lm.LMtests(model, model.W, test = "RLMlag", zero.policy = TRUE)
lm.LMtests(model, model.W, test = "RLMerr", zero.policy = TRUE)






### Spatial lag model ###
lagmodel <- lagsarlm(covid_deat ~ log(dist_borde) + pm2.5_mean + road_len + SO2_mean + log(alcohol_ca) + X_no2mean, data = targetfile2@data, listw = model.W, zero.policy = TRUE, na.action(na.exclude))
summary(lagmodel)

# Get residuals of lag model as a column, then plot them
targetfile2$lagmodelres <- residuals(lagmodel)
tm_shape(targetfile2)+tm_polygons("lagmodelres", palette = "-RdBu", style = "quantile")

targetfile2$lagmodelfit <- exp(fitted.values(lagmodel))

# targetfile2$err.res <- residuals(targetfile2.err)





## spatial error model 
# then fit the spatial error model to the data 
# targetfile2.err <- errorsarlm(covid_deat ~ log(dist_borde) + pm2.5_mean + road_len + SO2_mean + log(alcohol_ca) + X_no2mean, data=targetfile2@data, listw=model.W)
#summary(targetfile2.err) #

# targetfile2$err.fit <- exp(fitted.values(targetfile2.err)) #
 ## tm_shape(targetfile2)+tm_polygons("err.fit", palette = "-RdBu", style = "quantile") #



# SPATIAL ERROR MODEL

errormodel <- errorsarlm(covid_deat ~ log(dist_borde) + pm2.5_mean + road_len + SO2_mean + log(alcohol_ca) + X_no2mean, data = targetfile@data, listw = model.W, zero.policy = TRUE, na.action = na.exclude)
summary(errormodel)

# Get residuals of error model as a column, then plot them
targetfile$errormodelres <- residuals(errormodel)
# plot residuals now
tm_shape(targetfile)+tm_polygons("errormodelres", palette = "-RdBu", style = "quantile")

# Fit the spatial errormodel values for the dataset
# plot later all together 
targetfile2$errormodelfit <- exp(fitted.values(errormodel))





### Spatial Durbin model ###
durbinmodel <- lagsarlm(covid_deat ~ log(dist_borde) + pm2.5_mean + road_len + SO2_mean + log(alcohol_ca) + X_no2mean, data = targetfile2@data, listw = model.W, type = "mixed", zero.policy = TRUE, na.action = na.exclude)
summary(durbinmodel)

# Get residuals of durbin model as a column, then plot them
targetfile2$durbinmodel.res <- residuals(durbinmodel)
# plot
tm_shape(targetfile2)+tm_polygons("durbinmodel.res", palette = "-RdBu", style = "quantile")

# fit
targetfile2$durbin.fit <- exp(fitted.values(durbinmodel))
# plot later all together 
# tm_shape(targetfile2)+tm_polygons("durbinfit.fit", palette = "-RdBu", style = "quantile")







### Arrange maps of the three modelled and the observed covid death values ###
# Computationally intensive / time-consuming, can be done one-by one too

m1 <- tm_shape(targetfile2)+tm_polygons("lagmodelfit", style = "quantile")
m2 <- tm_shape(targetfile2)+tm_polygons("errormodelfit", style = "quantile")
m3 <- tm_shape(targetfile2)+tm_polygons("durbin.fit", style = "quantile")
m4 <- tm_shape(targetfile2)+tm_polygons("covid_deat", style = "quantile")
tmap_arrange(m1,m2,m3,m4)









### GEOGRAPHICALLY WEIGHTED REGRESSION ###

# You need to have the data in projected CRS to do it
# Transformation is the following: proj4string from https://epsg.io/5243

target.pr <- spTransform(targetfile2, CRS("+proj=lcc +lat_1=48.66666666666666 +lat_2=53.66666666666666 +lat_0=51 +lon_0=10.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

                        
# Get polygon centroids
target.c <- gCentroid(target.pr, byid=T)

# Get distance matrix from the polygon centroids
target.d <- spDists(target.c, longlat=F)

# Apply Gaussian kernel to distances
h <- 40000
w1 <- exp(-1/2*((target.d[1,]/h)^2))
plot(sort(w1, decreasing = T), type="l")

# Create spatial weight matrix for GWR
N <- length(target.pr)
W1 <- matrix(0, N, N)
diag(W1) <- w1
#Plot first ten rows to make sure its okay (value 1,1 has to be 0)
kable(W1[1:10,1:10], caption="First 10 rows and columns of GWR spatial weight matrix", booktabs=T)

# Create parameter for location 1 (these are not necessary for the solution just copied it there to have it from the tutorial)
# They are commented because they don't affect the program, 2:7's 7 should be changed to number of vars + 1 when you use it)
X <- as.matrix(cbind(1, model$model[,2:7]))
X <- apply(X, 2, as.numeric)
y <- as.matrix(model$model[,1], drop=F)
b1 <- solve(t(X)%*%W1%*%X)%*%t(X)%*%W1%*%y # Parameter object for location 1

# Choosing appropriate bandwidth (chooses one with lowest RMSE), twG outputs the bandwidth in metres
bwG <- gwr.sel(covid_deat ~ log(dist_borde) + pm2.5_mean + road_len + SO2_mean + log(alcohol_ca) + X_no2mean, data = target.pr, gweight = gwr.Gauss, verbose = FALSE, longlat = FALSE)
bwG

# Fitting GWR model with twG bandwidth
gwrG <- gwr(covid_deat ~ log(dist_borde) + pm2.5_mean + road_len + SO2_mean + log(alcohol_ca) + X_no2mean, data = target.pr, bandwidth = bwG, gweight = gwr.Gauss, hatmatrix = TRUE, longlat=FALSE)
gwrG

# Do moran's I test to see if there is autocorrelation remaining (probably yes)
gwr.morantest(gwrG, model.W)

# Plot a histogram of the variable estimates stored in gwrG's SDF object (not complusory, only interesting)
hist(gwrG$SDF$SO2_mean, 100)

# Create maps of coefficient values for each variable
# Mapping can be modified in order to interpret the values more precisely. For example you can go with the own breaks of variables (quantile, jenks)
coef.list <- colnames(gwrG$SDF@data[2:8]) # Extract variable names
map.list <- list()

for(i in 1:length(coef.list))
  #for(i in 1:2)# Loop through variable names and create a map for each one.
{
  if(i==1)
  {
    map <- tm_shape(gwrG$SDF)+tm_polygons(coef.list[i], style="fixed", palette="-RdBu", breaks=seq(-3,2,0.5), legend.show=TRUE)+tm_layout(title=coef.list[i])
  }
  else{map <- tm_shape(gwrG$SDF)+tm_polygons(coef.list[i], style="fixed", palette="-RdBu", breaks=seq(-3,2,0.5), legend.show=FALSE)+tm_layout(title=coef.list[i])}
  
  map.list[[coef.list[i]]] <- map
}

tmap_arrange(map.list)


### GEOGRAPHICALLY WEIGHTED REGRESSION WITH ADAPTIVE BANDWIDTH ###
# Choosing adaptive bandwidth based on the number of nearest neighbours https://crd230.github.io/gwr.html
# This accounts for spatial density of polygon centroids in contrast to the above displayed fixed bandwidth one
bwGa <- gwr.sel(covid_deat ~ log(dist_borde) + pm2.5_mean + road_len + SO2_mean + log(alcohol_ca) + X_no2mean, data = target.pr, gweight = gwr.Gauss, verbose = FALSE, longlat = FALSE, adapt = TRUE)
bwGa

# Fitting model with adaptive bandwidth
gwrGa <- gwr(covid_deat ~ log(dist_borde) + pm2.5_mean + road_len + SO2_mean + log(alcohol_ca) + X_no2mean, data = target.pr, adapt = bwGa, gweight = gwr.Gauss, hatmatrix = TRUE, longlat=FALSE)
gwrGa

#Plot gwrGa bandwidths
target.pr$adapt <- gwrGa$bandwidth
tm_shape(target.pr)+tm_polygons("adapt", palette = "Reds", style = "quantile")

# Map the coefficient values the same way as the fixed bandwidth ones
coef.list <- colnames(gwrGa$SDF@data[2:7]) # Extract variable names
map.list <- list()

for(i in 1:length(coef.list))
  #for(i in 1:2)# Loop through variable names and create a map for each one.
{
  if(i==1)
  {
    map <- tm_shape(gwrGa$SDF)+tm_polygons(coef.list[i], style="fixed", palette="-RdBu", breaks=seq(-3,2,0.5), legend.show=TRUE)+tm_layout(title=coef.list[i])
  }
  else{map <- tm_shape(gwrGa$SDF)+tm_polygons(coef.list[i], style="fixed", palette="-RdBu", breaks=seq(-3,2,0.5), legend.show=FALSE)+tm_layout(title=coef.list[i])}
  
  map.list[[coef.list[i]]] <- map
}

tmap_arrange(map.list)

# Same thing with quantile breaks
for(i in 1:length(coef.list))
  #for(i in 1:2)# Loop through variable names and create a map for each one.
{
  if(i==1)
  {
    map <- tm_shape(gwrG$SDF)+tm_polygons(coef.list[i], style="quantile", palette="-RdBu", legend.show=TRUE)+tm_layout(title=coef.list[i])
  }
  else{map <- tm_shape(gwrGa$SDF)+tm_polygons(coef.list[i], style="quantile", palette="-RdBu", legend.show=TRUE)+tm_layout(title=coef.list[i])}
  
  map.list[[coef.list[i]]] <- map
}

tmap_arrange(map.list)


