library(ranger)
library(clhs)
library(geodist)

rm(list=ls())

set.seed(12345)

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
prop = 0.1#c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)# proportion of training samples
point = 1#c(3,1,0,-1)# which point: point 3 = Point C1 in the study, point 1 = Point C2, point 0 = Point O1, point -1 = Point O2
design = "LHS"#c("ORD","LHS")# which design: ORD = resampling or LHS = cLHS algorithm
PPP = 3 #number of predictors corresponding to the distances where the TC characteristics are extracted
############################################
## LOAD DATA
############################################
load(paste0("./data/Cyclone.RData"))
N <- length(Cyclone)

NIT = 25

for (DESIGN in design){
for (PROP in prop){
for (POINT in point){

Q2 <- Q2CV <- RMSECV <- RMSE <- MRAE <- MRAECV <- MS <- MSCV <- SS <- SSCV <- CACV <- CA <- matrix(0,PPP,NIT)
qq = 0
for (PP in PPP:PPP){
	qq = qq + 1

load(paste0("./data/Hs_extractPt",POINT,".RData"))
pt = unlist(pt)
## DOE

doe <- matrix(0,N,3*PP)

for (ii in 1:N){

	distance.proxy <- geodist( x= Cyclone[[ii]], y = pt, measure = "geodesic") / 1000
	ff <- order((distance.proxy))
	if (length(ff) >= PP){
		U <- Cyclone[[ii]][ff[1:PP],"U10_max"]
		R <- Cyclone[[ii]][ff[1:PP],"Radius"]
		D <- sqrt(distance.proxy)[ff[1:PP]]
	}else{
		U <- Cyclone[[ii]][ff[1:length(ff)],"U10_max"]
		U <- c(U,rep(U[length(ff)],PP - length(ff)))
		R <- Cyclone[[ii]][ff[1:length(ff)],"Radius"]
		R <- c(R,rep(R[length(ff)],PP - length(ff)))
		D <- sqrt(distance.proxy)[ff[1:length(ff)]]
		D <- c(D,rep(D[length(ff)],PP - length(ff)))
	}
	doe[ii,] <- c(U,R,D)
}

# scale between 0 and 1
for (ii in 1:ncol(doe)){
	doe[,ii] <- scale(doe[,ii],min(doe[,ii],na.rm=T),max(doe[,ii],na.rm=T))
}

doe <- data.frame(doe)

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

	### Uncomment for checking the validation against the true values
	'
	save(
      	ss,Y.tr,Y.te,X.tr,X.te,doe,
      	file = paste0("./doe/DOE_",DESIGN,"Pt",POINT,"_N",PROP,"IT",IT,"NumP",PP,".RData")
	)
	'

	df <-  data.frame(X.tr, Y = Y.tr)

	Y.hatB <- Q1.hatB <- Q3.hatB <- Y.teB <- list()
	ssB <- matrix(sample(1:nrow(df),nrow(df),replace = F),ncol=10)
	
	for (bb in 1:10){
		dfBB <- df[-ssB[,bb],]
		Y.teB[[bb]] <- df[ssB[,bb],"Y"]
		X.teB <- df[ssB[,bb],-ncol(df)]
		dfBB <- na.omit(dfBB )
		X.teB <- na.omit(X.teB)
		modBB<-ranger(Y ~ ., data = dfBB,keep.inbag=TRUE, quantreg = TRUE)
		test = predict(modBB,X.teB)
		Y.hatB[[bb]] <- test$predictions
		Q1.hatB[[bb]] <- unlist(predict(modBB, X.teB, type = "quantiles", quantiles = 0.25)$predictions)
		Q3.hatB[[bb]] <- unlist(predict(modBB, X.teB, type = "quantiles", quantiles = 0.75)$predictions)
	}

	(Q2CV[qq,IT] <- 1-sum((unlist(Y.teB)-unlist(Y.hatB))^2)/sum((unlist(Y.teB)-mean(unlist(Y.teB)))^2)) ## Q2
	(MRAECV[qq,IT] = (mean(abs(unlist(Y.teB)-unlist(Y.hatB)))))## MAE

	yte = unlist(Y.teB)
	hat = unlist(Y.hatB)
	shat = unlist(Q3.hatB) -  unlist(Q1.hatB)

	### Uncomment for checking the validation against the true values
	'
	plot(yte,hat,
		xlim= c(0,15),ylim= c(0,15),
		pch=15,lwd=2,
		cex=0.9,cex.main=1.5,cex.lab=1.5,xlab="Hs - True value [m]",ylab="Hs - Predicted value [m]")
	
	for (j in 1:length(yte)) lines(c(yte[j],yte[j]),c(q1[j],q3[j]),lwd=1.2)
	abline(0,1,lwd=4,lty=2,col="grey50")
	title(paste0("Point ",qq), line = 0.5,adj = 0,cex.main = 2)
	'	

	MSCV[qq,IT] = mean(shat) ## mean uncertainty widths
	SSCV[qq,IT] = sd(shat) ## sdt dev. uncertainty widths

	### PRED ###################################################
	mod <- ranger(Y~., data=df,keep.inbag=TRUE, quantreg = TRUE)
	Y.hat <- predict(mod,X.te)$predictions
	(Q2[qq,IT] <- 1-sum((unlist(Y.te)-unlist(Y.hat))^2)/sum((unlist(Y.te)-mean(unlist(Y.te)))^2)) ## Q2
	(MRAE[qq,IT] = (mean(abs(unlist(Y.te)-unlist(Y.hat)))))## MAE

	q1hat = predict(mod, X.te, type = "quantiles", quantiles = 0.25)$predictions
	q3hat = predict(mod, X.te, type = "quantiles", quantiles = 0.75)$predictions
	
	'
	plot(Y.te,Y.hat,
		xlim= c(0,15),ylim= c(0,15),
		pch=15,lwd=2,
		cex=0.9,cex.main=1.5,cex.lab=1.5,xlab="Hs - True value [m]",ylab="Hs - Predicted value [m]")
	for (j in 1:length(Y.te)) lines(c(Y.te[j],Y.te[j]),c(q1hat[j],q3hat[j]),lwd=1.2)
	abline(0,1,lwd=4,lty=2,col="grey50")
	'

	MS[qq,IT] = mean(q3hat - q1hat) ## mean uncertainty widths
	SS[qq,IT] = sd(q3hat - q1hat) ## sdt dev. uncertainty widths

} ## IT

}#PP

### SAVE
save(Q2,Q2CV,MRAE, MRAECV,MS, MSCV,SS, SSCV, 
				file = paste0("./valid/Test_",DESIGN,"Pt",POINT,"_N",PROP,".RData")
		)

}#PROP

}#POINT

}#DESIGN
