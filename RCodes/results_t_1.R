K<-3;H<-5
B<-1000
library(dplyr)
true_value<-c(c(.7, .5, .3,0,0), c(.5, .3, .7,0,0), c(.3, .7, .5,0,0))
ESTIM<-NULL
for(k in c(1,4)){
load(paste("RCodes/res",k,".RData",sep=""))
mcoa<-NULL
#mcoa2<-NULL
omcoa<-NULL
gomcoa<-NULL
gsomcoa1<-NULL
gsomcoa2<-NULL
gsomcoa3<-NULL
simul<-NULL

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
ESTIM$true_value<-rep(true_value,12)

library(ggplot2)

ESTIM2 =ESTIM|>filter(method %in% c("mcoa","omcoa","gomcoa","gsomcoa1"))
ESTIM2$method<-factor(ESTIM2$method,levels=c("mcoa","omcoa","gomcoa","gsomcoa1"),
labels=c("mcoa","omcoa","gibbs_omcoa","gibbs_sparse_omcoa"))
ESTIM2$block<-factor(ESTIM2$block,labels = paste("block",1:K))
ESTIM2$dim<-factor(ESTIM2$dim,labels = paste("dim",1:H))
pd <- position_dodge(width = 0.1)
TRUE_DF <- distinct(ESTIM2, block, dim, true_value)
ggplot(ESTIM2, aes(x = sigma, y = estim, color = method)) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = estim - 2*sd, ymax = estim + 2*sd),
                width = 0.1, position = pd) +
  geom_hline(data = TRUE_DF, aes(yintercept = true_value), linetype = "dashed") +
  theme_bw() +
  facet_grid(block ~ dim)
