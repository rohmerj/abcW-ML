library(POT)
library(ranger)
library(transport)

rm(list= ls())

############################################
## FUNCTION
############################################
retlev <- function(us, scale,shape,tau,npy,prob){
	# us: POT threshold
	# tau: occurrence probability of excesses
	# shape,scale: GPD parameters
	# npy: number of excesses  per year
	# probability level
	us + scale * ((tau*npy/prob)^shape - 1)/ shape
}

############################################
## PARAMETER of the STUDY
############################################
prop = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)# proportion of training samples

n = 685 # number of cyclones
npy = n/1000 # rate of TC occurrence

T = c(100,500) # targeted return period
prob = 1/T # corresponding probability level

point = c(3,1)# which point: point 3 = Point 1 in the study, point 1 = Point 2 in the study
DESIGN <- "LHS" # NAME of the DESIGN

lim = c(0.05,0.075,0.125,0.15) # ABC threshold
N.sim = 2500 #  number of ABC random samples

US = c(5.04,NA,6.16,NA,6.12) # POT threshold

NIT = 25 # number of iterations

############################################
## COMPUTE
############################################
for (LIM in lim){

for (PROB in prob){

for (PROP in prop){

for (POINT in point){

for (IT0 in 1:NIT){

############################################
## DATA
############################################
us = US[POINT]
load(paste0("./data/Hs_extractPt",POINT,".RData"))
pt = unlist(pt)

load(paste0("./doe/DOE_",DESIGN,"Pt",POINT,"_N",PROP,"IT",IT0,".RData"))

Hs1 = Hs[ss]
data.s <- Hs1[which(Hs1 > us)]
data0 <- Hs[which(Hs > us)]

############################################
## REF - MLE
############################################
mle0 <- fitgpd(data0,us)
tau = length(data0)/685 
RL00 = unlist(retlev(us,coef(mle0)["scale"],coef(mle0)["shape"],tau,npy,PROB))

############################################
## SIM - wasserstein - NO error - REF
############################################
DD = scale.sim = shape.sim = NULL
for (i in 1:N.sim){
	#print(i)
	scale.sim[i] = rnorm(1,coef(mle0)["scale"],coef(mle0)["scale"]*0.5)
	if (scale.sim[i] < 0) scale.sim[i] = 0.01   
	shape.sim[i] = rnorm(1,coef(mle0)["shape"],abs(coef(mle0)["shape"])*0.25)
	data.sim = rgpd(n, us, scale.sim[i], shape.sim[i])
	DD[i] = wasserstein1d(data0,data.sim)
}
ff = which(DD < quantile(DD,LIM))
fit = cbind(scale.sim,shape.sim)

SC0 = quantile(fit[ff,1],c(0.05,0.5,0.95))
SH0 = quantile(fit[ff,2],c(0.05,0.5,0.95))
RL0 = quantile(unlist(retlev(us,fit[ff,1],fit[ff,2],tau,npy,PROB)),c(0.05,0.5,0.95))

############################################
## SIM - wasserstein - NO error - SUBSET
############################################
mle.s <- fitgpd(data.s,us)
n.s = length(data.s)
tau = length(data.s)/685
RL0.s = unlist(retlev(us,coef(mle.s)["scale"],coef(mle.s)["shape"],tau,npy,PROB))

DD = scale.sim = shape.sim = NULL
for (i in 1:N.sim){
	#print(i)
	scale.sim[i] = rnorm(1,coef(mle.s)["scale"],coef(mle.s)["scale"]*0.5)
	if (scale.sim[i] < 0) scale.sim[i] = 0.01   
	shape.sim[i] = rnorm(1,coef(mle.s)["shape"],abs(coef(mle.s)["shape"])*0.25)
	data.sim = rgpd(n.s, us, scale.sim[i], shape.sim[i])
	DD[i] = wasserstein1d(data.s,data.sim)
}
ff = which(DD < quantile(DD,LIM))
fit = cbind(scale.sim,shape.sim)

SC.s = quantile(fit[ff,1],c(0.05,0.5,0.95))
SH.s = quantile(fit[ff,2],c(0.05,0.5,0.95)) 

RL.s = quantile(unlist(retlev(us,fit[ff,1],fit[ff,2],tau,npy,PROB)),c(0.05,0.5,0.95))

############################################
## ML fitting
############################################
X.tr <- doe[ss,]
Y.tr <- (Hs[ss])

X.te <- doe[-ss,]
Y.te <- (Hs[-ss])

mod <- ranger(Y~., data= data.frame(X.tr, Y = Y.tr),mtry = 3,keep.inbag=TRUE, quantreg = TRUE)
y.hat <- unlist(predict(mod, X.te)$predictions)

############################################
## SIM - wasserstein - NO error - ML
############################################
Hs2 = c(Hs[ss],y.hat)
data.m <- Hs2[which(Hs2 > us)]

mle.m <- fitgpd(data.m,us)
tau = length(data.m)/685 
RL0.m = unlist(retlev(us,coef(mle.m)["scale"],coef(mle.m)["shape"],tau,npy,PROB))

DD = scale.sim = shape.sim = NULL
for (i in 1:N.sim){
	#print(i)
	scale.sim[i] = rnorm(1,coef(mle.s)["scale"],coef(mle.s)["scale"]*0.5)
	if (scale.sim[i] < 0) scale.sim[i] = 0.01   
	shape.sim[i] = rnorm(1,coef(mle.s)["shape"],abs(coef(mle.s)["shape"])*0.25)
	data.sim = rgpd(n, us, scale.sim[i], shape.sim[i])
	DD[i] = wasserstein1d(data.m,data.sim)
}
ff = which(DD < quantile(DD,LIM))
fit = cbind(scale.sim,shape.sim)

SC.m = quantile(fit[ff,1],c(0.05,0.5,0.95))
SH.m = quantile(fit[ff,2],c(0.05,0.5,0.95)) 
RL.m = quantile(unlist(retlev(us,fit[ff,1],fit[ff,2],tau,npy,PROB)),c(0.05,0.5,0.95))

############################################
## SIM - wasserstein - ML + error
############################################
n.te = nrow(X.te)
qq = seq(0,1,by=0.005)
q.hatr <- unlist(predict(mod, X.te, type = "quantiles", quantiles = qq)$predictions)
qapprox = function(xout,qy=q.hatr,x=qq){
	qa = NULL
	for (i in 1:nrow(qy)){
		qa[i] = approx(qq,qy[i,],xout[i])$y
	}
	return(qa)
}

DD = scale.sim = shape.sim = NULL
for (i in 1:N.sim){
	scale.sim[i] = rnorm(1,coef(mle.m)["scale"],coef(mle.m)["scale"]*0.5)
	if (scale.sim[i] < 0) scale.sim[i] = 0.01   
	shape.sim[i] = rnorm(1,coef(mle.m)["shape"],abs(coef(mle.m)["shape"])*0.25)
	data.sim = rgpd(n, us, scale.sim[i], shape.sim[i])
	qr = runif(n.te)
	#q.hatr <- unlist(predict(mod, X.te, type = "quantiles", quantiles = qr)$predictions)
	q.hat = qapprox(qr)

	Hs2r = c(Hs[ss],q.hat)
	data.mr <- Hs2r[which(Hs2r > us)]
	DD[i] = wasserstein1d(data.mr,data.sim)
}

ff = which(DD < quantile(DD,LIM))
fit = cbind(scale.sim,shape.sim)

SC.mr = quantile(fit[ff,1],c(0.05,0.5,0.95))
SH.mr = quantile(fit[ff,2],c(0.05,0.5,0.95)) 
RL.mr = quantile(unlist(retlev(us,fit[ff,1],fit[ff,2],tau,npy,PROB)),c(0.05,0.5,0.95))

### store results for RL
df.rl <- data.frame(
	 med = c(RL.mr[2],RL.m[2],RL.s[2],RL0[2]),
	 q1 = c(RL.mr[1],RL.m[1],RL.s[1],RL0[1]),
	 q3 = c(RL.mr[3],RL.m[3],RL.s[3],RL0[3]),
	 unc = c(RL.mr[3] - RL.mr[1],RL.m[3] - RL.m[1],RL.s[3] - RL.s[1],RL0[3] - RL0[1]),
	 bias = c(RL.mr[2] - RL0[2],RL.m[2] - RL0[2],RL.s[2] - RL0[2],0),
	 uncr = c(RL.mr[3] - RL.mr[1],RL.m[3] - RL.m[1],RL.s[3] - RL.s[1],RL0[3] - RL0[1]) / (RL0[3] - RL0[1])
	)

### store results for SCALE
df.sc <- data.frame(
	 med = c(SC.mr[2],SC.m[2],SC.s[2],SC0[2]),
	 q1 = c(SC.mr[1],SC.m[1],SC.s[1],SC0[1]),
	 q3 = c(SC.mr[3],SC.m[3],SC.s[3],SC0[3]),
	 unc = c(SC.mr[3] - SC.mr[1],SC.m[3] - SC.m[1],SC.s[3] - SC.s[1], SC0[3] - SC0[1]),
	 bias = c(SC.mr[2] - SC0[2],SC.m[2] - SC0[2],SC.s[2] - SC0[2],0),
	 uncr = c(SC.mr[3] - SC.mr[1],SC.m[3] - SC.m[1],SC.s[3] - SC.s[1],SC0[3] - SC0[1]) / (SC0[3] - SC0[1])
	)

### store results for SHAPE
df.sh <- data.frame(
	 med = c(SH.mr[2],SH.m[2],SH.s[2],SH0[2]),
	 q1 = c(SH.mr[1],SH.m[1],SH.s[1],SH0[1]),
	 q3 = c(SH.mr[3],SH.m[3],SH.s[3],SH0[3]),
	 unc = c(SH.mr[3] - SH.mr[1],SH.m[3] - SH.m[1],SH.s[3] - SH.s[1], SH0[3] - SH0[1]),
	 bias = c(SH.mr[2] - SH0[2],SH.m[2] - SH0[2],SH.s[2] - SH0[2],0),
	 uncr = c(SH.mr[3] - SH.mr[1],SH.m[3] - SH.m[1],SH.s[3] - SH.s[1],SH0[3] - SH0[1]) / (SH0[3] - SH0[1])
	)

### save
save(	df.rl,
	df.sc,
	df.sh,
	mle0, mle.m, mle.s,
	RL00, RL0.m, RL0.s,
	file=paste0("./results/Test_Pt",POINT,"_N", PROP, "Design", DESIGN, "T", 1/PROB, "IT", IT0,"abc",LIM, ".RData")
	)

}# IT0

}# POINT

}# PROP

}# RP

}# ABC lim

