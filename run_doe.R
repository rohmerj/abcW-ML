library(clhs)
library(geodist)

rm(list=ls())

#set.seed(123456)

############################################
## FUNCTION
############################################

## scale between 0 and 1
scale <- function(x,m,M){
	y <- (x - m) / (M - m)
}

############################################
## PARAMETER of the STUDY
############################################
NIT = 25# number of iterations
prop = c(0.2,0.3,0.4,0.6,0.7,0.8,0.9)# proportion of training samples
point = c(3,1)# which point: point 3 = Point 1 in the study, point 1 = Point 2 in the study
design <- c("ORD","LHS")# which design: ORD = resampling or LHS = cLHS algorithm

############################################
## LOAD DATA
for (DESIGN in design){
for (PROP in prop){
for (POINT in point){

load(paste0("./data/Hs_extractPt",POINT,".RData"))
pt = unlist(pt)

## DOE
N <- length(Cyclone)
doe <- matrix(0,N,3)

for (ii in 1:N){

	distance.proxy <- geodist( x= Cyclone[[ii]], y = pt) / 1000
	ff <- order((distance.proxy))

	U <- max(Cyclone[[ii]][ff[1:1],"U10_max"],na.rm = T)
	R <- max(Cyclone[[ii]][ff[1:1],"Radius"],na.rm = T)

	doe[ii,] <- c(U,R,min(sqrt(distance.proxy)))
}

# scale between 0 and 1
for (ii in 1:ncol(doe)){
	doe[,ii] <- scale(doe[,ii],min(doe[,ii],na.rm=T),max(doe[,ii],na.rm=T))
}

doe <- data.frame(doe)
colnames(doe) <- c("U10","Radius","distance")

############################################
## ANALYSIS
############################################
for (IT in 1:NIT){ ## ITERATIONS

	if (DESIGN == "ORD"){
		ss <- sample(1:N,round(PROP*N),replace = FALSE)
	}
	if (DESIGN == "LHS"){
		ss <- clhs(data.frame(doe), size = round(PROP*N), use.cpp = TRUE)
	}
	X.tr <- doe[ss,]
	Y.tr <- (Hs[ss])

	X.te <- doe[-ss,]
	Y.te <- (Hs[-ss])

	save(
      		ss,Y.tr,Y.te,X.tr,X.te,doe,
      		file = paste0("./doe/DOE_",DESIGN,"Pt",POINT,"_N",PROP,"IT",IT,".RData")
	)

} ## IT

}#PROP

}#POINT

}#DESIGN


