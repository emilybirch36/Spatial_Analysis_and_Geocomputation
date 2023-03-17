# Spatial_Analysis_and_Geocomputation
Spatial analysis in R to determine which factors influence Covid-19 deaths.


Group project:
My task was focusing on factors influencing Covid-19 deaths in Germany, whereas other group memebers focused on England, the Netherlands and France.


Project summary:
- Comparing 4 countries of similar governmental handling of early Covid-19 restrictions/ policy.
- Factors: age, pollution, BAME, distance to border, underlying health conditions (alcoholism, obesity), age, income, smoking, road density.
- Data for Germany sourced from: OpenStreetMap data from geofabrik, Robert Koch Institut, Copernicus Sentinel 5 and Bundesamt für Statistik.
- Exploratory data analysis (ESDA): histograms and maps of covid deaths per 100k people (normalised). Global Moran's I to test for global spatial autocorrelation (testing attribute random deistirbution by clustering). Getis and Ord's Gi* statistic to measure how high and low values are clustered across the study area, based on neighbouring values. Local Moran's I to test local spatial autocorrelation, a local indicator of spatial association (LISA). 
- Method: a linear regression model (Ordinary Least Squares) using the ratio of Covid-19 deaths and population in each spatial unit as the dependent variable. all other influencing factors are the independent variables. VIF (Variance Inflation Factor) was used to account for possinle multicolinearity between inependent variables.
- Lagrange multiplier tests were used to evaluate spatial lag and error models (if there was significant autocorrelation remaining). An adjacency-based spatial regression method (the spatial Durbin model) was also used. 
- The model fit was also evaluated by mapping model residuals to observe potential clustering. Fixed bandwidth geographically-weighted regression was used to look at spatial heterogeneity. the remaining aurocorrelation was calculated with the gwr.morantest function. Each variable's coefficient estimates were mapped to visually examine possible factors influencing covid deaths in different locations. 
- Akaike Information Criterion (AIC) was used to compare the efficiency and fit of spatial lag, error, Durbin models, and fixed and adaptive bandwidth GWRs. A lower value indicates an improvment in model fit. This determined the best model. 






