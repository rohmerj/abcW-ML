library(ranger)

rm(list = ls())

set.seed(12345)

DESIGN = "LHS" # type of design, ORD = resampling, LHS = cLHS algorithm

for (PROP in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)){ # proportion of training samples

print(PROP)

for (IT0 in 1:25){ # iteration

Q2 <- Q2CV <- RMSECV <- RMSE <- MRAE <- MRAECV <- MS <- MSCV <- SS <- SSCV <- CACV <- CA <- NULL
# XX for independent test set
# XXCV for cross validation

qq = 0
for(POINT in c(3,1)){

	qq = qq + 1
	load(paste0("./data/Hs_extractPt",POINT,".RData"))

	load(paste0("./doe/DOE_",DESIGN,"Pt",POINT,"_N",PROP,"IT",IT0,".RData"))
	X.tr <- doe[ss,]
	Y.tr <- (Hs[ss])

	X.te <- doe[-ss,]
	Y.te <- (Hs[-ss])

	df <-  data.frame(X.tr, Y = Y.tr)

	Y.hatB <- Q1.hatB <- Q3.hatB <- Y.teB <- list()
	ssB <- matrix(sample(1:nrow(df),nrow(df),replace = F),ncol=10)
	
	for (bb in 1:10){
		dfBB <- df[-ssB[,bb],]
		Y.teB[[bb]] <- df[ssB[,bb],"Y"]
		X.teB <- df[ssB[,bb],-ncol(df)]
		modBB<-ranger(Y ~ ., data = dfBB,mtry=3,keep.inbag=TRUE, quantreg = TRUE)
		test = predict(modBB,X.teB)
		Y.hatB[[bb]] <- test$predictions
		Q1.hatB[[bb]] <- unlist(predict(modBB, X.teB, type = "quantiles", quantiles = 0.25)$predictions)
		Q3.hatB[[bb]] <- unlist(predict(modBB, X.teB, type = "quantiles", quantiles = 0.75)$predictions)
	}

	(Q2CV[qq] <- 1-sum((unlist(Y.teB)-unlist(Y.hatB))^2)/sum((unlist(Y.teB)-mean(unlist(Y.teB)))^2)) ## Q2
	(MRAECV[qq] = (mean(abs(unlist(Y.teB)-unlist(Y.hatB)))))## MAE
	CA0 = NULL
	q1 <- unlist(Q1.hatB)
	q3 <- unlist(Q3.hatB)
	te <- unlist(Y.teB)
	for (kk in 1:length(te)){
		CA0[kk] = ifelse(te[kk] >= q1[kk] & te[kk] <= q3[kk], 1, 0)
	}
	CACV[qq] = mean(CA0)

	yte = unlist(Y.teB)
	hat = unlist(Y.hatB)
	shat = unlist(Q3.hatB) -  unlist(Q1.hatB)

	MSCV[qq] = mean(shat) ## mean uncertainty widths
	SSCV[qq] = sd(shat) ## sdt dev. uncertainty widths

	### PRED ###################################################
	mod <- ranger(Y~., data=df, mtry = 3,keep.inbag=TRUE, quantreg = TRUE)
	Y.hat <- predict(mod,X.te)$predictions
	(Q2[qq] <- 1-sum((unlist(Y.te)-unlist(Y.hat))^2)/sum((unlist(Y.te)-mean(unlist(Y.te)))^2)) ## Q2
	(MRAE[qq] = (mean(abs(unlist(Y.te)-unlist(Y.hat)))))## MAE

	q1hat = predict(mod, X.te, type = "quantiles", quantiles = 0.25)$predictions
	q3hat = predict(mod, X.te, type = "quantiles", quantiles = 0.75)$predictions

	CA0 = NULL
	for (kk in 1:length(Y.te)){
		CA0[kk] = ifelse(Y.te[kk] >= q1hat[kk] & Y.te[kk] <= q3hat[kk], 1, 0)
	}
	CA[qq] = mean(CA0)

	MS[qq] = mean(q3hat - q1hat) ## mean uncertainty widths
	SS[qq] = sd(q3hat - q1hat) ## sdt dev. uncertainty widths

}

write.table(data.frame(Q2,Q2CV,MRAE, MRAECV,CA, CACV, MS, MSCV,SS, SSCV), 
				file=paste0("./valid/","Test_N",PROP,"_IT",IT0,".txt")

}

}

