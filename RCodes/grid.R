K<-3;R<-3
N<-100
dim<-rep(1:R,K)
block<-c(rep(1,R),rep(2,R),rep(3,R))
w <- c(c(.7, .5,0), c(.5, .7,0), c(.3,0,.5))
P <- c(10, 5, 15)
pk<-c(rep(P[1],K),rep(P[2],K),rep(P[3],K))
tab<-data.frame(dim=dim,block=block,weight=w,threshold=w/(N*pk)^.25)

tab1<-tab[sort(tab$threshold,decreasing = F,index=T)$ix,]

library(xtable)
objet_xtable<-xtable(tab1[tab1$threshold>0,],digits = c(0,0,0,2,2))
print(objet_xtable, include.rownames = FALSE)
