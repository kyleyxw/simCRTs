################################################################################
# example function to simulate data for model-based type I error rate
# beta: log hazard ratio under the null scenario i.e. log(1)
# clstr_info: STRIDE study cluster size file. You can create your own cluster size file. 
# clstr_size: a k*1 vector with the size of each cluster
# clstr: number of clusters k (k=86 for SRIDE study)
# lambda: annual event rate
# q: follow-up in years
# cmptrsk: a vector of competing risk rate
# dropout: annual dropout rate
# tau_f0: a vector of frailty tau
# tau_c0: a vector of copula tau
# m: number of simulations
################################################################################
############################ Example Parameter Settings ########################
beta<-log(1)
clstr_info <- read.delim('STRIDE_cluster_sizes.txt', header = TRUE, sep = "\t", dec = ".")
clstr_size <- clstr_info$clustersize
clstr<-100
lambda<-0.08
q<-40/12
cmptrsk<-c(0.02, 0.04, 0.08, 0.12)
dropout<-0.03
tau_f0<-c(0.001)
tau_c0<-c(0.001)
m = 2
################################################################################
############################# Package Settings #################################
#Check the existence of required packages
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

#Load the required packages
usePackage("survival")
usePackage("cmprsk")
usePackage("crrSC")
usePackage("dplyr")
#################################################################################

modelbased_typeIerror <- function(beta=log(1), clstr_size, clstr,lambda, q, cmptrsk, dropout, tau_f, tau_c, m){
  
  dir.create("./modelbased_typeIerror")
  theta<-1/(1-exp(q*log(1-dropout)))
  
  p3<-length(cmptrsk)
  p2<-length(tau_f0)
  p1<-length(tau_c0)
  
  smry<-NULL
  for (o3 in 1:p3){ 
    for (o2 in 1:p2){
      for (o1 in 1:p1){
        
        pb1<-txtProgressBar(min=0,max=m,style=2)
        result <- NULL
        gamma<-cmptrsk[o3]
        
        for (j in 1:m) {
          new_clstr_size <- sample(clstr_size,clstr,replace=T)
          n <- sum(new_clstr_size)
          clusters <- rep(1:clstr,new_clstr_size)
          tau_c<-tau_c0[o1]
          delta<-1/(1-tau_c)
          temp<-matrix(,n,4)
          
          for (i in 1:n) {
            u<-runif(1,min=0,max=1)
            w<-runif(1,min=0,max=1)
            f<-function(z) z+(delta-1)*log(z)-(-log(u)+(delta-1)*log(-log(u))-log(w))
            z<-uniroot(f,lower=-log(u),upper=1000000)$root
            v<-exp(-(z^delta-(-log(u))^delta)^(1/delta))
            temp[i,]<-c(u,w,z,v)
          }
          
          data<-data.frame(temp)
          colnames(data)<-c("u","w","z","v")
          data$id <- seq(1:n)
          data$clusters <- clusters
          
          #Kendall's tau, which is used to measure the dependence within cluster. 
          tau_f<-tau_f0[o2]
          
          #If a is the shape parameter of gamma distribution, tau=1/(2*a+1)
          a<-0.5*(1/tau_f-1)
          
          #Let rate parameter=shape parameter in order to make mean of frailty
          b<-a 
          
          #Generate the frailty for each cluster
          alpha<-rgamma(clstr,a,b) 
        
          #Re-assign treatment indicator to make it a cluster randomized design
          data$alpha <- alpha[data$clusters] 
          data$x <- data$clusters%%2
          
          #Simulate event time, competing risk time and censoring time
          data$time_event<--log(data$u)/(lambda*exp(beta*data$x)*data$alpha)
          data$time_cmpt<--log(data$v)/(gamma*data$alpha)
          data$time_event[data$time_event==0] <- 0.01
          data$time_cmpt[data$time_cmpt==0] <- 0.01
          data$time_c<-sapply(runif(n,0,theta),function(x){min(x,q)})
          data$time_t<-apply(data[,c("time_event","time_cmpt","time_c")],1,min,na.rm = T)
          data$status<- as.numeric(data$time_t==data$time_event)
          data$status[data$time_t==data$time_cmpt] <- 2
          
          #Simulate multi-state data
          data$status0<-0
          data$start1<-0
          data$end1<-data$time_c
          data$status1<-0
          data$start2<-NA
          data$end2<-NA
          data$status2<-NA
          
          data$status0[data$status==1]<-1
          data$start1[data$status==1]<-0
          data$end1[data$status==1]<-data$time_event[data$status==1]
          data$status1[data$status==1]<-1
          data$start2[data$status==1]<-data$time_event[data$status==1]
          data$end2[data$status==1 & data$time_cmpt<=data$time_c]<-data$time_cmpt[data$status==1 & data$time_cmpt<=data$time_c]
          data$end2[data$status==1 & data$time_cmpt>data$time_c]<-data$time_c[data$status==1 & data$time_cmpt>data$time_c]
          data$status2[data$status==1 & data$time_cmpt<=data$time_c]<-2
          data$status2[data$status==1 & data$time_cmpt>data$time_c]<-0
          
          data$status0[data$status==2]<-0
          data$start1[data$status==2]<-NA
          data$end1[data$status==2]<-NA
          data$status1[data$status==2]<-NA
          data$start2[data$status==2]<-0
          data$end2[data$status==2]<-data$time_cmpt[data$status==2]
          data$status2[data$status==2]<-2
          
          #Create column labels and transfrom wide format to long format
          data_w<-data[,c('id','clusters','x','time_t','status','status0')]
          data_l1<-data[,c('id','clusters','x','start1','end1','status1')]
          data_l2<-data[,c('id','clusters','x','start2','end2','status2')]
          colnames(data_l1)<-c('id','clusters','x','start','end','status')
          colnames(data_l2)<-c('id','clusters','x','start','end','status')
          data_l1$seq<-1
          data_l2$seq<-2
          data_l0<-rbind(data_l1,data_l2)
          data_l<-data_l0[complete.cases(data_l0),]
          data_l<-data_l[order(data_l$cluster,data_l$id,data_l$seq),]
          data_l<-data_l[!(data_l$end-data_l$start<=0.001),]
          data_l$ind <- as.numeric(data_l$end<=data_l$start)
          
          
          #Fit Different models
          test.ms<-coxph(formula = Surv(start,end,as.factor(status))~x,id=id, data = data_l)
          test.ms_cluster<-coxph(formula = Surv(start,end,as.factor(status))~x,id=id,cluster = clusters, data = data_l)
          test.cox<-coxph(Surv(data_w$time_t,data_w$status0)~data_w$x)
          test.cox_cluster<-coxph(Surv(data_w$time_t,data_w$status0)~data_w$x+cluster(data_w$clusters))
          test.crr<-crr(data_w$time_t,data_w$status,data_w$x)
          test.zhou<-crrc(ftime=data_w$time_t,fstatus=data_w$status,cov1=data_w$x,cluster=data_w$clusters,maxiter = 100)
          
          
          #Summarize model
          q.95 <- -qnorm((1-0.95)/2)  
          ms.coef<- data.frame(Coefficient = test.ms$coefficients[1],
                               SE = sqrt(test.ms$var[1,1]), 
                               coverage = test.ms$coefficients[1]-q.95*sqrt(test.ms$var[1,1])<beta & test.ms$coefficients[1]+q.95*sqrt(test.ms$var[1,1])>beta, 
                               modelName = "ms")
          ms_cluster.coef <- data.frame(Coefficient = test.ms_cluster$coefficients[1],
                                        SE = sqrt(test.ms_cluster$var[1,1]), 
                                        coverage = test.ms_cluster$coefficients[1]-q.95*sqrt(test.ms_cluster$var[1,1])<beta & test.ms_cluster$coefficients[1]+q.95*sqrt(test.ms_cluster$var[1,1])>beta,
                                        modelName = "ms_cluster")
          cox.coef <- data.frame(Coefficient = summary(test.cox)$coefficients[1],
                                 SE = summary(test.cox)$coefficients[3], 
                                 coverage = confint(test.cox)[1]<beta & confint(test.cox)[2]>beta, 
                                 modelName = "Cox")
          cox_cluster.coef <- data.frame(Coefficient = summary(test.cox_cluster)$coefficients[1],
                                         SE = summary(test.cox_cluster)$coefficients[4], 
                                         coverage = confint(test.cox_cluster)[1]<beta & confint(test.cox_cluster)[2] >beta, 
                                         modelName = "Marginal Cox")
          crr.coef <- data.frame(Coefficient = test.crr$coef,
                                 SE = sqrt(test.crr$var),
                                 coverage = test.crr$coef -q.95*sqrt(test.crr$var)<beta & test.crr$coef +q.95*sqrt(test.crr$var) >beta,
                                 modelName = "Fine and Gray")
          zhou.coef <- data.frame(Coefficient = test.zhou$coef,
                                  SE = sqrt(test.zhou$var),
                                  coverage = test.zhou$coef -q.95*sqrt(test.zhou$var)<beta & test.zhou$coef +q.95*sqrt(test.zhou$var) >beta,
                                  modelName = "Zhou")
          
          full.coef <- data.frame(rbind(ms.coef, ms_cluster.coef, cox.coef, cox_cluster.coef,crr.coef, zhou.coef),row.names = NULL)
          full.coef$tau1 = tau_c
          full.coef$tau2 = tau_f
          full.coef$lambda1 = lambda
          full.coef$lambda2 = gamma
          full.coef$dropout = dropout
          full.coef$HR = exp(beta)
          full.coef$mean_cluster_size = mean(new_clstr_size)
          full.coef$CV_cluster_size = sqrt(var(new_clstr_size))/mean(new_clstr_size)
          full.coef$coe.true = beta
          
          result <- data.frame(rbind(result,full.coef))
          setTxtProgressBar(pb1,j)
          
        }
        colnames(result)<-c("coefficient","SE","coverage","modelName","tau1","tau2","lambda1","lambda2","dropout","HR","mean.cluster.size","CV.cluster.size","TrueBeta")
        write.csv(result,file = paste0('./modelbased_typeIerror/tau_f_',tau_f,"_tau_c_",tau_c,"_cmptrsk_",gamma,".csv"))
        close(pb1)
      }
    }  
  }
  warnings()
}

#################################################################################
#Call the function
modelbased_typeIerror(beta=beta, clstr_size=clstr_size, clstr=clstr,lambda =lambda, 
                 q = q, cmptrsk=cmptrsk, dropout=dropout, tau_f=tau_f, tau_c=tau_c, m=m)
#Combine the results into a single file
csv.list <- list.files(path="./modelbased_typeIerror/", pattern=".csv$", full.names=TRUE)
full<-NULL
for (i in 1:length(csv.list)) {
  update<-read.csv(csv.list[i])
  full<-rbind(full,update)
}
system("rm ./modelbased_typeIerror/*csv")
write.csv(full,file = paste0('./modelbased_typeIerror/[m]full_cluster',clstr,"_hr",exp(beta),".csv"))

