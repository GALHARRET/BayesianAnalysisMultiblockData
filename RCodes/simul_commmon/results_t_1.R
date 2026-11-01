K<-3;H<-5

library(dplyr)
true_value<-c(c(.7, .5, 0,0,0), c(.5, .7, 0,0,0), c(.3, 0, .5,0,0))
ESTIM<-NULL
for(k in 1:7){
load(paste("RCodes/simul_commmon/res",k,".RData",sep=""))
mcoa<-NULL
#mcoa2<-NULL
omcoa<-NULL
gomcoa<-NULL
gsomcoa1<-NULL
gsomcoa2<-NULL
gsomcoa3<-NULL
simul<-NULL

B<-length(res$mcoa) 
for(b in 1:B){
simul<-c(simul,rep(b,K*H))
mcoa<-rbind(mcoa,res$mcoa[[b]])
omcoa<-rbind(omcoa,res$omcoa[[b]])
gomcoa<-rbind(gomcoa,res$gomcoa[[b]])
gsomcoa1<-rbind(gsomcoa1,res$gsomcoa1[[b]]$summary)
gsomcoa2<-rbind(gsomcoa2,res$gsomcoa2[[b]]$summary)
gsomcoa3<-rbind(gsomcoa3,res$gsomcoa3[[b]]$summary)  
}
mcoa$simul<-simul
mcoa$method<-"mcoa"
mcoa$sigma<-k/10
omcoa$simul<-simul
omcoa$method<-"omcoa"
omcoa$sigma<-k/10
gomcoa$simul<-simul
gomcoa$method<-"gomcoa"
gomcoa$sigma<-k/10
gsomcoa1$simul<-simul
gsomcoa1$method<-"gsomcoa1"
gsomcoa1$sigma<-k/10
gsomcoa2$simul<-simul
gsomcoa2$method<-"gsomcoa2"
gsomcoa2$sigma<-k/10
gsomcoa3$simul<-simul
gsomcoa3$method<-"gsomcoa3"
gsomcoa3$sigma<-k/10

ESTIM=rbind(ESTIM,mcoa |>
  group_by(sigma,method,block,dim)|>
  summarise(estim=mean(estimate),sd=sd(estimate))|>
  as_tibble())

ESTIM=rbind(ESTIM,
  omcoa |>
  group_by(sigma,method,block,dim)|>
  summarise(estim=mean(estimate),sd=sd(estimate))|>
  as_tibble())

ESTIM=rbind(ESTIM,
  gomcoa |>
  group_by(sigma,method,block,dim)|>
  summarise(estim=mean(mean),sd=sd(mean))|>
    as_tibble())

ESTIM=rbind(ESTIM,
  gsomcoa1 |>
  group_by(sigma,method,block,dim)|>
  summarise(estim=mean(mean),sd=sd(mean))|>
    as_tibble())

ESTIM=rbind(ESTIM,
  gsomcoa2 |>
  group_by(sigma,method,block,dim)|>
  summarise(estim=mean(mean),sd=sd(mean))|>
    as_tibble())

ESTIM=rbind(ESTIM,
  gsomcoa3 |>
  group_by(sigma,method,block,dim)|>
  summarise(estim=mean(mean),sd=sd(mean))|>
    as_tibble())
}





dim(ESTIM)
ESTIM$true_value<-rep(true_value,630/length(true_value))

library(ggplot2)

ESTIM2 =ESTIM|>filter(method %in% c("mcoa","omcoa","gomcoa","gsomcoa1"))
ESTIM2$method<-factor(ESTIM2$method,levels=c("mcoa","omcoa","gomcoa","gsomcoa1"),
labels=c("mcoa","omcoa","gibbs_omcoa","gibbs_sparse_omcoa"))
ESTIM2$block<-factor(ESTIM2$block,labels = paste("block",1:K))
ESTIM2$dim<-factor(ESTIM2$dim,labels = paste("dim",1:H))
ESTIM2$sigma<-factor(ESTIM2$sigma)
pd <- position_dodge(width = 0.4)
TRUE_DF <- distinct(ESTIM2, block, dim, true_value)

ESTIM3<-ESTIM2 |>
  filter(dim %in% c("dim 1","dim 2","dim 3"))|>
  mutate(dim=droplevels(dim))
TRUE_DF2<-TRUE_DF|>
  filter(dim %in% c("dim 1","dim 2","dim 3"))|>
  mutate(dim=droplevels(dim))

ggplot(ESTIM3,aes(x = sigma, y = estim, color = method)) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = estim - 2*sd, ymax = estim + 2*sd),
                width = 0.05, position = pd) +
  geom_hline(data = TRUE_DF2, aes(yintercept = true_value), linetype = "dashed") +
  theme_bw() +
  facet_grid(block ~ dim)+ 
  theme(legend.position = "bottom")



PIP<-NULL
for(k in 1:7){
load(paste("RCodes/simul_commmon/res",k,".RData",sep=""))
for(b in 1:B){
PIP<-rbind(PIP,cbind(res$gsomcoa1[[b]]$summary[,c("block","dim","pip")],lambda=1,sigma=k/10))
PIP<-rbind(PIP,cbind(res$gsomcoa2[[b]]$summary[,c("block","dim","pip")],lambda=0.1,sigma=k/10))
PIP<-rbind(PIP,cbind(res$gsomcoa3[[b]]$summary[,c("block","dim","pip")],lambda=10,sigma=k/10))
}  
}
PIP$lambda<-factor(PIP$lambda)
PIP$sigma<-factor(PIP$sigma)
PIP$block<-factor(PIP$block,labels=paste("block",1:K))
PIP$dim<-factor(PIP$dim,labels=paste("dim",1:H))
ggplot(PIP,aes(x=sigma,y=pip,color=lambda))+
geom_boxplot()+
  facet_grid(block~dim)+
  theme_bw()

PIP |> group_by(block,dim,lambda) |> summarise(mean(pip),sd(pip))

RANK<-NULL
for(k in 1:7){
load(paste("RCodes/simul_commmon/res",k,".RData",sep=""))
names(res$gsomcoa1[[1]]$rank)
for(b in 1:B){
RANK<-rbind(RANK,data.frame(res$gsomcoa1[[b]]$rank,lambda=1,sigma=k/10))
RANK<-rbind(RANK,data.frame(res$gsomcoa2[[b]]$rank,lambda=0.1,sigma=k/10))
RANK<-rbind(RANK,data.frame(res$gsomcoa3[[b]]$rank,lambda=10,sigma=k/10))
}  
}
RANK$lambda<-factor(RANK$lambda)
RANK$sigma<-factor(RANK$sigma)

colnames(RANK)[1:2]<-c("rank","proba")

ggplot(RANK,aes(x=sigma,y=proba,color=rank))+
geom_boxplot()+
  theme_bw()+
  facet_grid(~lambda)

RANK |> group_by(lambda,sigma,rank) |> summarise(mean(proba),sd(proba))


