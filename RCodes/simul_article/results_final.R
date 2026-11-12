K<-3;H<-4
setwd("~/Documents/GitHub/BayesianAnalysisMultiblockData/")
library(dplyr)
library(ggplot2)
library(tidyr)
true_value<-c(c(.7, .5, 0,0), c(.5, .7, 0,0), c(.3, 0, .5,0))
ESTIM<-NULL
sigma<-c(0.02,0.04,0.06,0.08,0.12,0.16,0.18)

#### CHOICE BETWEEN Lambda and G #####

B<-250

PIP<-NULL
for(i in 1:length(sigma)){
  load(paste("RCodes/simul_article/res",i,".l_vs_g.RData",sep=""))
  for(b in 1:B){
    PIP<-rbind(PIP,cbind(res$gsomcoa1[[b]]$summary[,c("block","dim","pip")],choice=1,sigma=sigma[i]))
    PIP<-rbind(PIP,cbind(res$gsomcoa2[[b]]$summary[,c("block","dim","pip")],choice=2,sigma=sigma[i]))
    PIP<-rbind(PIP,cbind(res$gsomcoa3[[b]]$summary[,c("block","dim","pip")],choice=3,sigma=sigma[i]))
  }
}

PIP$statut<-rep(factor(true_value>0,labels=c("inactive","active")),dim(PIP)[1]/length(true_value))
PIP$statut_id<-factor(PIP$pip>=0.5,labels=c("inactive","active"))

PIP$choice<-factor(PIP$choice)
PIP$block<-factor(PIP$block,labels=paste("block",1:K))
PIP$dim<-factor(PIP$dim,labels=paste("dim",1:H))

res1 <- PIP |>
  group_by(sigma, block, dim, choice) |>
  summarise(n = sum(statut_id == "active"), .groups = "drop")

print(res1, n=32)


p<-ggplot(res1,aes(x=sigma,y=proportion,group=choice,color=choice))+
  geom_point()+
  geom_line()+
  facet_grid(block~dim)+
  theme_bw()+
  theme(legend.position = "bottom")+
  scale_colour_manual(
    values = c("1" = "red", "2" = "lightgreen","3"="lightblue"),
    name = "choice",
    # Intégration d'expressions mathématiques dans les étiquettes de la légende
    labels = c(
      expression(paste("norm ",lambda,"=1",sep="")),
      expression(paste("noise ",g[k],"=",(N-1)*p[k],sep="")),
      expression(paste("noise ",g[k],"=10",(N-1)*p[k],sep=""))
    ))
p

ggsave(filename = "article/figs/choice_l_vs_g.eps", plot = p, width = 12, height = 6, units = "in")


#### CHOICE of LAMBDA #####
PIP<-NULL

for(k in 1:length(sigma)){
  load(paste("RCodes/simul_article/res",k,".RData",sep=""))
  B<-length(res$mcoa)
  for(b in 1:B){
    PIP<-rbind(PIP,cbind(res$gsomcoa1[[b]]$summary[,c("block","dim","pip")],lambda=.1,sigma=sigma[k]))
    PIP<-rbind(PIP,cbind(res$gsomcoa2[[b]]$summary[,c("block","dim","pip")],lambda=1,sigma=sigma[k]))
    PIP<-rbind(PIP,cbind(res$gsomcoa3[[b]]$summary[,c("block","dim","pip")],lambda=10,sigma=sigma[k]))
  }
}

PIP$statut<-rep(factor(true_value>0,labels=c("inactive","active")),dim(PIP)[1]/length(true_value))
PIP$statut_id<-factor(PIP$pip>=0.5,labels=c("inactive","active"))

PIP$lambda<-factor(PIP$lambda)
PIP$sigma<-factor(PIP$sigma)
PIP$block<-factor(PIP$block,labels=paste("block",1:K))
PIP$dim<-factor(PIP$dim,labels=paste("dim",1:H))

res1<-PIP |>
  group_by(sigma,block,dim,lambda) |>
  summarise(n = sum(statut_id == "active"), .groups = "drop")
res1$proba<-res1$n/B



p<-ggplot(res1,aes(x=sigma,y=proba,group=lambda,color=lambda))+
  geom_point()+
  geom_line()+
  facet_grid(block~dim)+
  theme_bw()+
  theme(legend.position = "bottom")
p
ggsave(filename = "article/figs/choice_lambda.eps", plot = p, width = 12, height = 6, units = "in")

#### CHOICE OF SPIKE #####

B<-250

PIP<-NULL
for(i in 1:length(sigma)){
  load(paste("RCodes/simul_article/res",i,"AB.RData",sep=""))
  for(b in 1:B){
    PIP<-rbind(PIP,cbind(res$gsomcoa1[[b]]$summary[,c("block","dim","pip")],a=1,b=1,sigma=sigma[i]))
    PIP<-rbind(PIP,cbind(res$gsomcoa2[[b]]$summary[,c("block","dim","pip")],a=1,b=3,sigma=sigma[i]))
    PIP<-rbind(PIP,cbind(res$gsomcoa3[[b]]$summary[,c("block","dim","pip")],a=.5,b=.5,sigma=sigma[i]))
  }
}

PIP$statut<-rep(factor(true_value>0,labels=c("inactive","active")),dim(PIP)[1]/length(true_value))
PIP$statut_id<-factor(PIP$pip>=0.5,labels=c("inactive","active"))

PIP$spike<-factor(paste("(",PIP$a,",",PIP$b,")",sep=""))
PIP$block<-factor(PIP$block,labels=paste("block",1:K))
PIP$dim<-factor(PIP$dim,labels=paste("dim",1:H))

res1 <- PIP |>
  group_by(sigma, block, dim, spike) |>
  summarise(n = sum(statut_id == "active"), .groups = "drop")
res1$proportion<-res1$n/B




p<-ggplot(res1,aes(x=sigma,y=proportion,group=spike,color=spike))+
  geom_point()+
  geom_line()+
  facet_grid(block~dim)+
  theme_bw()+
  theme(legend.position = "bottom")

p
ggsave(filename = "article/figs/choice_spike.eps", plot = p, width = 12, height = 6, units = "in")


#### RESULTS of ESTIMATION of SIGMA 

ESTIM<-NULL
resSIG<-NULL
for(k in 1:length(sigma)){
load(paste("RCodes/simul_article/res",k,".RData",sep=""))
mcoa<-NULL
omcoa<-NULL
gomcoa<-NULL
gsomcoa<-NULL
#gsomcoa2<-NULL
#gsomcoa3<-NULL
#simul<-NULL

B<-length(res$mcoa)
for(b in 1:B){
mcoa<-rbind(mcoa,res$mcoa[[b]])
omcoa<-rbind(omcoa,res$omcoa[[b]])
gomcoa<-rbind(gomcoa,res$gomcoa[[b]])
gsomcoa<-rbind(gsomcoa,res$gsomcoa2[[b]]$summary)
}
mcoa$model<-"mcoa"
mcoa$sigma<-sigma[k]
omcoa$model<-"omcoa"
omcoa$sigma<-sigma[k]
gomcoa$model<-"Gibbs_omcoa"
gomcoa$sigma<-sigma[k]
gsomcoa$model<-"Gibbs_sparse_omcoa"
gsomcoa$sigma<-sigma[k]

colnames(gomcoa)[3]<-"estimate"
colnames(gsomcoa)[3]<-"estimate"

selectCol<-c("block","dim","estimate","lower","upper","model","sigma")
resSIG<-rbind(resSIG,gomcoa[,selectCol],gsomcoa[,selectCol])

ESTIM=rbind(ESTIM,mcoa |>
  group_by(sigma,model,block,dim)|>
  summarise(estim=mean(estimate),sd=sd(estimate))|>
  as_tibble())

ESTIM=rbind(ESTIM,
  omcoa |>
  group_by(sigma,model,block,dim)|>
  summarise(estim=mean(estimate),sd=sd(estimate))|>
  as_tibble())

ESTIM=rbind(ESTIM,
  gomcoa |>
  group_by(sigma,model,block,dim)|>
  summarise(estim=mean(estimate),sd=sd(estimate))|>
    as_tibble())

ESTIM=rbind(ESTIM,
  gsomcoa |>
  group_by(sigma,model,block,dim)|>
  summarise(estim=mean(estimate),sd=sd(estimate))|>
    as_tibble())

}
resSIG$true_value<-rep(true_value,dim(resSIG)[1]/length(true_value))

resSIG$width<-resSIG$upper-resSIG$lower


width<-resSIG |> 
  group_by(sigma,dim,block,model) |>
  summarise(M=mean(width))

width|> filter(sigma>=0.12,model=="Gibbs_sparse_omcoa") |> group_by(sigma) |> summarise(min(M))

isInCI<-function(x){
  x<-as.vector(x)
return(as.numeric(abs(x$true_value-(x$lower+x$upper)/2)<=(x$upper-x$lower)/2))
}
I<-NULL
for(b in 1:dim(resSIG)[1]){
x<-resSIG[b,]
I<-c(I,isInCI(x))
}
resSIG$isInCI<-I

cov_percent<-resSIG |>
  #filter(dim %in% 1:3) |>
  group_by(block,dim,model,sigma) |>
  summarise(P=sum(isInCI)/B)

cov_percent$model<-factor(cov_percent$model)
cov_percent$dim<-factor(cov_percent$dim,labels=paste("dim",1:H))
cov_percent$block<-factor(cov_percent$block,labels=paste("block",1:K))
p<-ggplot(cov_percent|>filter(dim!="dim 4"),aes(x=sigma,y=P,color=model,group=model))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = .95)+
  facet_grid(block~dim)+
  theme_bw()+
  theme(legend.position = "bottom")
p
ggsave(filename = "article/figs/cov_percent.eps", plot = p, width = 10, height = 5, units = "in")

ESTIM$true_value<-rep(true_value,dim(ESTIM)[1]/length(true_value))
ESTIM$sigma<-factor(ESTIM$sigma)
ESTIM$model<-factor(ESTIM$model)
ESTIM$block=factor(ESTIM$block,labels=paste("block",1:K))
ESTIM$dim=factor(ESTIM$dim,labels=paste("dim",1:H))
pd <- position_dodge(width = 0.4)
ESTIM2=ESTIM|>filter(dim!="dim 4")
TRUE_DF <- distinct(ESTIM2, block, dim, true_value)

p<-ggplot(ESTIM2,aes(x = sigma, y = estim, color = model)) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = estim - 2*sd, ymax = estim + 2*sd),
                width = 0.05, position = pd) +
  geom_hline(data = TRUE_DF, aes(yintercept = true_value), linetype = "dashed") +
  theme_bw() +ylim(c(-0.5,2.0))+
  facet_grid(block ~ dim)+
  theme(legend.position = "bottom")
p
ggsave(filename = "article/figs/bias_estimation.eps", plot = p, width = 10, height = 5, units = "in")


#### RANK STUDY ####

RANK<-NULL
for(k in 1:length(sigma)){
  load(paste("RCodes/simul_article/res",k,".RData",sep=""))
  for(b in 1:B){
    RANK<-rbind(RANK,data.frame(effective_rank=names(which.max(res$gsomcoa2[[b]]$rank)),sigma=sigma[k],model="Gibbs_sparse_omcoa"))
  }
}
RANK$effective_rank<-as.numeric(RANK$effective_rank)

rank <- RANK |>
  count(sigma,effective_rank) |>
  complete(sigma, effective_rank = 0:H, fill = list(n = 0))

rank$effective_rank<-factor(rank$effective_rank,labels=paste("Reff=",0:H,sep=""))

rank <- rank |>
  group_by(sigma) |>
  mutate(proportion = n / sum(n)) |>
  ungroup()

rank$model="Gibbs_sparse_omcoa"

load("RCodes/simul_article/resGD.RData")
resGD<-as.data.frame(resGD)
resGD$effective_rank<-as.numeric(resGD$effective_rank)
resGD$sigma<-as.numeric(resGD$sigma)

rank2 <- resGD |>
  count(sigma,effective_rank) |>
  complete(sigma, effective_rank = 0:H, fill = list(n = 0))

rank2$effective_rank<-factor(rank$effective_rank,labels=paste("Reff=",0:H,sep=""))

rank2 <- rank2 |>
  group_by(sigma) |>
  mutate(proportion = n / sum(n)) |>
  ungroup()

rank2$model<-"Gavish–Donoho"

load("RCodes/simul_article/resHoff.RData")
resH<-data.frame(matrix(NA,nrow=dim(resHoff)[1],ncol=2))
for(b in 1:dim(resHoff)[1]) resH[b,]<-c(names(which.max(resHoff[b,1:H])),as.numeric(resHoff[b,6]))
colnames(resH)<-c("effective_rank","sigma")

resH$effective_rank<-as.numeric(resH$effective_rank)
resH$sigma<-as.numeric(resH$sigma)

rank3 <- resH |>
  count(sigma,effective_rank) |>
  complete(sigma, effective_rank = 0:H, fill = list(n = 0))

rank3$effective_rank<-factor(rank$effective_rank,labels=paste("Reff=",0:H,sep=""))

rank3 <- rank3 |>
  group_by(sigma) |>
  mutate(proportion = n / sum(n)) |>
  ungroup()

rank3$model<-"Hoff"




ColSelec<-c("sigma","effective_rank","proportion","model")
rank4<-rbind(rank[,ColSelec],rank2[,ColSelec],rank3[,ColSelec])
p<-ggplot(rank4,aes(x=sigma,y=proportion,color=model,group=model))+
  geom_point() +
  geom_line()+
  #geom_errorbar(aes(ymin = q1, ymax = q2),
  #              width = 0.05, position = pd) +
  theme_bw() +
  ylim(c(0,1))+
  facet_grid(.~effective_rank)+
  theme(legend.position = "bottom")

p

ggsave(filename = "article/figs/Reff_comparison.eps", plot = p, width = 10, height = 5, units = "in")



##############################


RANK<-NULL

for(i in 1:length(sigma)){
  load(paste("RCodes/simul_article/res",i,"AB.RData",sep=""))
  for(b in 1:B){
    RANK<-rbind(RANK,data.frame(effective_rank=names(which.max(res$gsomcoa1[[b]]$rank)),a=1,b=1,sigma=sigma[i]))
    RANK<-rbind(RANK,data.frame(effective_rank=names(which.max(res$gsomcoa2[[b]]$rank)),a=1,b=3,sigma=sigma[i]))
    RANK<-rbind(RANK,data.frame(effective_rank=names(which.max(res$gsomcoa3[[b]]$rank)),a=0.5,b=0.5,sigma=sigma[i]))
  }
}
RANK$spike<-factor(paste("(",RANK$a,",",RANK$b,")",sep=""))

rank<-RANK |>
  group_by(spike,sigma,effective_rank) |>
  count()
rank$prop<-rank$n/B
rank$effective_rank<-factor(rank$effective_rank,levels=0:H,labels=paste("Reff=",0:H,sep=""))
p<-ggplot(rank,aes(x=sigma,y=prop,group=spike,color=spike))+
  geom_point()+
  geom_line()+
  facet_grid(~effective_rank)+
  theme_bw()+
  theme(legend.position = "bottom")

p

ggsave(filename = "article/figs/figure4a.eps", plot = p, width = 12, height = 6, units = "in")

#### SENSITIVITY to lambda against noise



