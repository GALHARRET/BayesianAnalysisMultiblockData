K<-3;H<-5
setwd("~/Documents/GitHub/BayesianAnalysisMultiblockData/")
library(dplyr)
library(ggplot2)
true_value<-c(c(.7, .5, 0,0,0), c(.5, .7, 0,0,0), c(.3, 0, .5,0,0))
ESTIM<-NULL
sigma<-c(0.02,0.04,0.06,0.08,0.12,0.16,0.18)  
resSIG<-NULL

for(k in 1:length(sigma)){
load(paste("RCodes/simul_article/res",k,".RData",sep=""))
mcoa<-NULL
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
mcoa$model<-"mcoa"
mcoa$sigma<-sigma[k]
omcoa$simul<-simul
omcoa$model<-"omcoa"
omcoa$sigma<-sigma[k]
gomcoa$simul<-simul
gomcoa$model<-"gomcoa"
gomcoa$sigma<-sigma[k]
gsomcoa1$simul<-simul
gsomcoa1$model<-"gsomcoa1"
gsomcoa1$sigma<-sigma[k]
gsomcoa2$simul<-simul
gsomcoa2$model<-"gsomcoa2"
gsomcoa2$sigma<-sigma[k]
gsomcoa3$simul<-simul
gsomcoa3$model<-"gsomcoa3"
gsomcoa3$sigma<-sigma[k]

selectCol<-c("block","dim","mean","lower","upper","model","sigma")
resSIG<-rbind(resSIG,gomcoa[,selectCol],gsomcoa1[,selectCol],gsomcoa2[,selectCol],gsomcoa3[,selectCol])

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
  summarise(estim=mean(mean),sd=sd(mean))|>
    as_tibble())

ESTIM=rbind(ESTIM,
  gsomcoa1 |>
  group_by(sigma,model,block,dim)|>
  summarise(estim=mean(mean),sd=sd(mean))|>
    as_tibble())

ESTIM=rbind(ESTIM,
  gsomcoa2 |>
  group_by(sigma,model,block,dim)|>
  summarise(estim=mean(mean),sd=sd(mean))|>
    as_tibble())

ESTIM=rbind(ESTIM,
  gsomcoa3 |>
  group_by(sigma,model,block,dim)|>
  summarise(estim=mean(mean),sd=sd(mean))|>
    as_tibble())
}
head(resSIG)
resSIG$true_value<-rep(true_value,B*length(sigma)*4)

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
  filter(dim %in% 1:3) |>
  group_by(block,dim,model,sigma) |>
  summarise(P=sum(isInCI)/500)
cov_percent$model<-factor(cov_percent$model,levels=c("gomcoa","gsomcoa1","gsomcoa2","gsomcoa3"),labels=c("gibbs_omcoa","gibbs_sparse_omcoa(l=1)","gibbs_sparse_omcoa(l=.1)","gibbs_sparse_omcoa(l=10)"))
cov_percent$dim<-factor(cov_percent$dim,labels=paste("dim",1:3))
cov_percent$block<-factor(cov_percent$block,labels=paste("block",1:3))
p<-ggplot(cov_percent,aes(x=sigma,y=P,color=model,group=model))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = .95)+
  facet_grid(block~dim)+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave(filename = "article/figs/figure1.eps", plot = p, width = 10, height = 5, units = "in")

ESTIM$true_value<-rep(true_value,630/length(true_value))


ESTIM2 =ESTIM|>filter(model %in% c("mcoa","omcoa","gomcoa","gsomcoa1"))
ESTIM2$model<-factor(ESTIM2$model,levels=c("mcoa","omcoa","gomcoa","gsomcoa1"),
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

p<-ggplot(ESTIM3,aes(x = sigma, y = estim, color = model)) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = estim - 2*sd, ymax = estim + 2*sd),
                width = 0.05, position = pd) +
  geom_hline(data = TRUE_DF2, aes(yintercept = true_value), linetype = "dashed") +
  theme_bw() +ylim(c(-0.5,2.0))+
  facet_grid(block ~ dim)+ 
  theme(legend.position = "bottom")

ggsave(filename = "article/figs/figure2.eps", plot = p, width = 10, height = 5, units = "in")

ESTIM2 =ESTIM|>filter(model %in% c("gsomcoa1","gsomcoa2","gsomcoa3"))
ESTIM2$model<-factor(ESTIM2$model,labels=c("l=1","l=0.1","l=10"))
ESTIM2$block<-factor(ESTIM2$block,labels = paste("block",1:K))
ESTIM2$dim<-factor(ESTIM2$dim,labels = paste("dim",1:H))
ESTIM2$sigma<-factor(ESTIM2$sigma)
pd <- position_dodge(width = 0.4)
TRUE_DF <- distinct(ESTIM2, block, dim, true_value)

ESTIM3<-ESTIM2 |>
  filter(dim %in% c("dim 4","dim 5"))|>
  mutate(dim=droplevels(dim))
TRUE_DF2<-TRUE_DF|>
  filter(dim %in% c("dim 4","dim 5"))|>
  mutate(dim=droplevels(dim))

p<-ggplot(ESTIM3,aes(x = sigma, y = estim, color = model)) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = estim - 2*sd, ymax = estim + 2*sd),
                width = 0.05, position = pd) +
  geom_hline(data = TRUE_DF2, aes(yintercept = true_value), linetype = "dashed") +
  theme_bw() +
  facet_grid(block ~ dim)+ 
  theme(legend.position = "bottom")

ggsave(filename = "article/figs/figure3.eps", plot = p, width = 10, height = 5, units = "in")



PIP<-NULL
B<-500
for(k in 1:length(sigma)){
load(paste("RCodes/simul_article/res",k,".RData",sep=""))
for(b in 1:B){
PIP<-rbind(PIP,cbind(res$gsomcoa1[[b]]$summary[,c("block","dim","pip")],lambda=1,sigma=sigma[k]))
PIP<-rbind(PIP,cbind(res$gsomcoa2[[b]]$summary[,c("block","dim","pip")],lambda=0.1,sigma=sigma[k]))
PIP<-rbind(PIP,cbind(res$gsomcoa3[[b]]$summary[,c("block","dim","pip")],lambda=10,sigma=sigma[k]))
}
}

PIP$statut<-rep(factor(true_value>0,labels=c("inactive","active")),3*500*7)
PIP$statut_id<-factor(PIP$pip>=0.5,labels=c("inactive","active"))

PIP$lambda<-factor(PIP$lambda)
PIP$sigma<-factor(PIP$sigma)
PIP$block<-factor(PIP$block,labels=paste("block",1:K))
PIP$dim<-factor(PIP$dim,labels=paste("dim",1:H))

ggplot(PIP,aes(x=sigma,y=statut_id,col=lambda,group=lambda))+
  geom_point()+
  geom_line()+
  theme_bw()+
  facet_grid(block~dim)

res1<-PIP |> 
  group_by(sigma,block,dim,lambda) |>
  count(statut_id)
res1$proba<-res1$n/500


p<-ggplot(res1|>filter(statut_id=="active"),aes(x=sigma,y=proba,group=lambda,color=lambda))+
  geom_point()+
  geom_line()+
  facet_grid(block~dim)+
  theme_bw()+ 
  theme(legend.position = "bottom")

ggsave(filename = "article/figs/figure4.eps", plot = p, width = 12, height = 6, units = "in")


#### RANK STUDY ####

RANK<-NULL
for(k in 1:7){
load(paste("RCodes/simul_article/res",k,".RData",sep=""))
names(res$gsomcoa1[[1]]$rank)
for(b in 1:B){
RANK<-rbind(RANK,data.frame(res$gsomcoa1[[b]]$rank,lambda=1,sigma=sigma[k]))
RANK<-rbind(RANK,data.frame(res$gsomcoa2[[b]]$rank,lambda=0.1,sigma=sigma[k]))
RANK<-rbind(RANK,data.frame(res$gsomcoa3[[b]]$rank,lambda=10,sigma=sigma[k]))
}  
}
RANK$lambda<-factor(RANK$lambda)
RANK$sigma<-factor(RANK$sigma)

colnames(RANK)[1:2]<-c("rank","proba")


RK2<-RANK |> group_by(lambda,sigma,rank) |> summarise(m=mean(proba),
q1=quantile(proba,probs=0.025),q2=quantile(proba,probs=0.975))

pd <- position_dodge(width = 0.4)
RK2$lambda<-factor(RK2$lambda,labels=c("lambda=0.1","lambda=1","lambda=10"))
p<-ggplot(RK2,aes(x=sigma,y=m,color=rank,group=rank))+
  geom_point() +
  geom_line()+
  theme_bw() +
  facet_grid( ~ lambda)+ 
  theme(legend.position = "bottom")

ggsave(filename = "figure5.eps", plot = p, width = 10, height = 5, units = "in")

RANK3<-NULL
B<-500
for(k in 1:length(sigma)){
load(paste("RCodes/simul_article/res",k,".RData",sep=""))
for(b in 1:B){
RANK3<-rbind(RANK3,data.frame(choice=names(which.max(res$gsomcoa1[[b]]$rank)),lambda=1,sigma=sigma[k]))
RANK3<-rbind(RANK3,data.frame(choice=names(which.max(res$gsomcoa2[[b]]$rank)),lambda=0.1,sigma=sigma[k]))
RANK3<-rbind(RANK3,data.frame(choice=names(which.max(res$gsomcoa3[[b]]$rank)),lambda=10,sigma=sigma[k]))
}  
}
rank4<-RANK3 |>
  filter(lambda==1)|>
  group_by(lambda,sigma,choice) |>
  count()
rank4$prop<-rank4$n/B
rank4$model<-"gibbs_sparse_omcoa(l=1)"
rank4$effective_rank<-factor(rank4$choice,labels=c(paste("Reff=",0:5,sep="")))
#rank4$sigma<-factor(rank4$sigma)

load("RCodes/simul_article/resGD.RData")
rank5<-resGD |>
  as.data.frame()|>
  group_by(sigma,effective_rank)|>
  count()
B<-1000
rank5$prop<-rank5$n/B
#rank5$sigma<-factor(rank5$sigma)
rank5$effective_rank<-factor(rank5$effective_rank,levels=0:5,labels=paste("Reff=",0:5,sep=""))
rank5$model<-"Gavish–Donoho"

load("RCodes/simul_article/resHoff.RData")

resH<-data.frame(matrix(NA,nrow=dim(resHoff)[1],ncol=2))
for(b in 1:dim(resHoff)[1]) resH[b,]<-c(names(which.max(resHoff[b,1:6])),as.numeric(resHoff[b,7]))
colnames(resH)<-c("effective_rank","sigma")

rank6<-resH |>
  as.data.frame() |>
  group_by(sigma,effective_rank)|>
  count()
B<-300
rank6$prop<-rank6$n/B
#rank6$sigma<-factor(rank6$sigma)
rank6$effective_rank<-factor(rank6$effective_rank,levels=0:5,labels=paste("Reff=",0:5,sep=""))
rank6$model<-"Hoff"
rank6$sigma<-as.numeric(rank6$sigma)
print(rank6)

rank<-rbind(rank4[,c("sigma","effective_rank","prop","model")],rank5[,c("sigma","effective_rank","prop","model")],rank6[,c("sigma","effective_rank","prop","model")])
table(rank$model)
p<-ggplot(rank|>filter(model %in% c("gibbs_sparse_omcoa(l=1)","Gavish–Donoho","Hoff")),aes(x=sigma,y=prop,color=model,group=model))+
  geom_point() +
  geom_line()+
  #geom_errorbar(aes(ymin = q1, ymax = q2),
  #              width = 0.05, position = pd) +
  theme_bw() +
  ylim(c(0,1))+
  facet_grid(.~effective_rank)+ 
  theme(legend.position = "bottom")

p

ggsave(filename = "article/figs/figure6.eps", plot = p, width = 10, height = 5, units = "in")

x<-RANK3 |>
  group_by(lambda,sigma,choice) |>
  count()

x$prop<-x$n/sum(x$n[1:3])

library(xtable)
xtable(x[,-3])
