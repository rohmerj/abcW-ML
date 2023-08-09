library(POT)
library(ranger)
library(transport)

rm(list= ls())

set.seed(12345)

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
prop = c(0.1,0.5)# proportion of training samples

n = 685 # number of cyclones
npy = n/1000 # rate of TC occurrence

T = c(100,500) # targeted return period
prob = 1/T # corresponding probability level

point = 1#c(3,1,0,-1)# which point: point 3 = Point C1 in the study, point 1 = Point C2, point 0 = Point O1, point -1 = Point O2
DESIGN <- "LHS" # NAME of the DESIGN

lim = c(0.10) # ABC threshold
N.sim = 2500 #  number of ABC random samples

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
# POT threshold
if (POINT == 3) us = 6.16
if (POINT == 1) us = 5.04
if (POINT == 0) us = 9.20
if (POINT <0) us = 8.85

load(paste0("./data/Hs_extractPt",POINT,".RData"))
pt = unlist(pt)

load(
	file = paste0("./doe/DOE_",DESIGN,"Pt",POINT,"_N",PROP,"IT",IT0,"NumP",3,".RData")
)

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

SC0 = fit[ff,1]
SH0 = fit[ff,2]
RL0 = unlist(retlev(us,fit[ff,1],fit[ff,2],tau,npy,PROB))

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

SC.s = fit[ff,1]
SH.s = fit[ff,2]

RL.s = unlist(retlev(us,fit[ff,1],fit[ff,2],tau,npy,PROB))

############################################
## ML fitting
############################################
X.tr <- doe[ss,]
Y.tr <- (Hs[ss])

X.te <- doe[-ss,]
Y.te <- (Hs[-ss])

mod <- ranger(Y~., data= data.frame(X.tr, Y = Y.tr),keep.inbag=TRUE, quantreg = TRUE)
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

SC.m = fit[ff,1]
SH.m = fit[ff,2]
RL.m = unlist(retlev(us,fit[ff,1],fit[ff,2],tau,npy,PROB))

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

SC.mr = fit[ff,1]
SH.mr = fit[ff,2]
RL.mr = unlist(retlev(us,fit[ff,1],fit[ff,2],tau,npy,PROB))


### save

save(	SC0,SH0,RL0,
	SC.m,SH.m,RL.m,
	SC.s,SH.s,RL.s,
	SC.mr,SH.mr,RL.mr,
	mle0, mle.m, mle.s,
	RL00, RL0.m, RL0.s,
	file=paste0("./results/seuilABC/Test_Pt",POINT,"_N", PROP, "Design", DESIGN, "T", 1/PROB, "IT", IT0,"abc",LIM, ".RData")
	)

}# IT0

}# POINT

}# PROP

}# RP

}# ABC lim

