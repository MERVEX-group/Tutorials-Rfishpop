library(Rfishpop)

############################## EX 1 ##################################

ctrPop<-list(years=seq(1980,2020,by=1),niter=1,N0=100000,ages=0:20,
             minFage=4,maxFage=10,tc=0.5,seed=NULL)
number_ages<-length(ctrPop$ages);number_years<-length(ctrPop$years)
M<-matrix(rep(0.2,number_ages*number_years),ncol = number_years)
colnames(M)<-ctrPop$years
rownames(M)<-ctrPop$ages
ctrBio<-list(M=M,CV_M=0, L_inf=100, t0=0, k=0.15, CV_L=0, CV_LC=0, 
             a=10^(-5), b=3, a50_Mat=3, ad_Mat=-0.5,CV_Mat=0)
ctrSEL<-list(type="Logistic", par=list(a50_Sel=3, ad_Sel=-1),CV_SEL=0)
f=matrix(rep(0.5,number_years),ncol=number_years,nrow=1,byrow=TRUE)
ctrFish<-list(f=f,ctrSEL=ctrSEL)
a_BH=100000; b_BH=10000; CV_REC_BH=0
SR<-list(type="BH",par=c(a_BH,b_BH,CV_REC_BH))

Pop.Mod<-Population.Modeling(ctrPop=ctrPop,ctrBio=ctrBio,ctrFish=ctrFish,SR=SR)


############################## EX 2 ##################################

# population size for each age and year
N<-Pop.Mod$Matrices$N; round(N[,seq(1,41,by=10),1],2)

# instantaneous fishing mortality for each age and year
F<-Pop.Mod$Matrices$F; F[,seq(1,41,by=10),1]

# proportion of mature at each age and year
Mat<-Pop.Mod$Matrices$Mat; Mat[,seq(1,41,by=10),1]

# number of catches for each age and year
C_N<-Pop.Mod$Matrices$C_N; C_N[,seq(1,41,by=10),1]

# weight of catches for each age and year
C_W<-Pop.Mod$Matrices$C_W; C_W[,seq(1,41,by=10),1]

# weight of the catches for each year 
C<-Sum.Pop.Mod(Pop.Mod,c("C")); C$C[,,1]

# total biomass for each year
BIO<-Sum.Pop.Mod(Pop.Mod,c("BIO")); BIO$BIO[,,1]

# maturity biomass for each year (spawning stock biomass)
SSB<-Sum.Pop.Mod(Pop.Mod,c("SSB")); SSB$SSB[,,1]


############################## EX 3 ##################################

a_BH=100000; b_BH=10000
R<-RBH(SSB$SSB[,,1],a_BH,b_BH)
plot(SSB$SSB[,,1],R,type="b", pch=19, col="red", xlab="spawning biomass", ylab="recruitment")
steepness_value<-steepness(Pop.Mod,Fish.years=3,Bio.years=3,type="steepness",
                           Method="mean",par=NULL)
steepness_value[,,1]


############################## EX 4 ##################################

f.grid<-seq(0.00,2,by=0.01)
bpr<-BPR(Pop.Mod,f.grid,Bio.years=3,Fish.years=3,plot=c(TRUE,1),Method="mean",
         par=NULL)
head(bpr[,,1])

ypr<-YPR(Pop.Mod,f.grid,3,3,plot=c(TRUE,1), Method="mean",par=NULL)
head(ypr[,,1])


############################## EX 5 ##################################

RE0<-BYR.eq(Pop.Mod,0,3,3,FALSE,Method="mean",par=NULL)
N_eq0<-RE0$N; N_eq0[,,1]
DPM0<-RE0$DPM; DPM0[,,1]

RE0.5<-BYR.eq(Pop.Mod,0.5,3,3,FALSE,Method="mean",par=NULL)
N_eq0.5<-RE0.5$N; round(N_eq0.5[,,1],2)
DPM0.5<-RE0.5$DPM; round(DPM0.5[,,1],5)


############################## EX 6 ##################################

ma<-cbind(N[,41,],N_eq0.5)
colnames(ma)<-c("itself","effort"); round(ma,5)
round(ma[,1]-ma[,2],5)

age<-seq(0,20,by=1)
plot(age,ma[,1],type="b",pch=19,col="red",ylim=c(min(ma),max(ma)),
     xlab="age",ylab="N")
lines(age,ma[,2],pch=19,col="blue",type="b",lty=2)
legend(13, max(ma), legend=c("by itself", "fishing effort = 0.5"),
       col=c("red", "blue"), lty=1:2, cex=0.8)


############################## EX 7 ##################################

rf<-RF(Pop.Mod, 3,3,Method="mean",par=NULL,FM_type="F_msy",iters=1,plot=TRUE); rf$F_msy[,,1]
Fmsy<-rf$F_msy[2]; Fmsy


############################## EX 8 ##################################

ctrPop2<-list(years=seq(1980,2020,by=1),niter=1,N0=N_eq0,ages=0:20,
             minFage=4,maxFage=10,tc=0.5,seed=NULL)

fmsy<-rf$F_msy[1]
f=matrix(c(rep(fmsy*2,25),rep(fmsy,number_years-25)),ncol=number_years,nrow=1,byrow=TRUE)
ctrFish2<-list(f=f,ctrSEL=ctrSEL)

Pop.Mod2<-Population.Modeling(ctrPop=ctrPop2,ctrBio=ctrBio,ctrFish=ctrFish2,SR=SR)


############################## EX 9 ##################################

SEL<-Sum.Pop.Mod(Pop.Mod2,c("SEL"))
q_A<-SEL[["SEL"]][,,1]
gamma<-1; CV_A=0.1

BIO2<-Sum.Pop.Mod(Pop.Mod2,c("BIO"))
BIO2<-Sum.Pop.Mod(Pop.Mod2,c("BIO"))$BIO[,,1]; BIO2


IA0.1<-Sampling_Survey(Pop.Mod=Pop.Mod2,type="biomass",q_A=q_A,gamma=gamma,CV_A=CV_A,tsampling=0)
IA0.1<-IA0.1$biomass[,,1]; IA0.1

mat<-cbind(IA0.1,BIO2); mat

years<-seq(1980,2020,by=1)
plot(years,IA0.1,type="b",pch=19,col="red",ylim=c(min(mat),max(mat)),
                xlab="years",ylab="biomass")
lines(years,BIO2,pch=19,col="blue",type="b")
legend(1998, max(mat), legend=c("sample (CV = 0.1)", "population"),
       col=c("red", "blue"), lty=1, cex=0.8)


CV_A=0.2
IA0.2<-Sampling_Survey(Pop.Mod=Pop.Mod2,type="biomass",q_A=q_A,gamma=gamma,CV_A=CV_A,tsampling=0)
IA0.2<-IA0.2$biomass[,,1]; IA0.2

mat<-cbind(IA0.2,BIO2); mat

plot(years,IA0.2,type="b",pch=19,col="red",ylim=c(min(mat),max(mat)),
     xlab="years",ylab="biomass")
lines(years,BIO2,pch=19,col="blue",type="b")
legend(1998, max(mat), legend=c("sample (CV = 0.2)", "population"),
       col=c("red", "blue"), lty=1, cex=0.8)


C2<-Sum.Pop.Mod(Pop.Mod2,c("C"))$C[,,1]; C2

IC0.1<-Sampling_Catch(Pop.Mod=Pop.Mod2,type="catch weight",CV_CN=0.1)
IC0.1<-IC0.1$weight[,,1]; IC0.1

mat<-cbind(IC0.1,C2); mat

plot(years,IC0.1,type="b",pch=19,col="red",ylim=c(min(mat),max(mat)),
     xlab="years",ylab="weight")
lines(years,C2,pch=19,col="blue",type="b")
legend(1998, max(mat), legend=c("sample (CV = 0.1)", "population"),
       col=c("red", "blue"), lty=1, cex=0.8)


IC0.2<-Sampling_Catch(Pop.Mod=Pop.Mod2,type="catch weight",CV_CN=0.2)
IC0.2<-IC0.2$weight[,,1]; IC0.2

mat<-cbind(IC0.2,C2); mat

plot(years,IC0.2,type="b",pch=19,col="red",ylim=c(min(mat),max(mat)),
     xlab="years",ylab="weight")
lines(years,C2,pch=19,col="blue",type="b")
legend(1998, max(mat), legend=c("sample (CV = 0.1)", "population"),
       col=c("red", "blue"), lty=1, cex=0.8)
