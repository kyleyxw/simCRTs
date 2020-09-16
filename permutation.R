################################################################################
# example function for permutation tests with small number of clusters e.g.30
# beta: log hazard ratio, under the null scenario -> type I error
#                         under the alternative -> power
# clstr_info: STRIDE study cluster size file. You can create your own cluster size file. 
# clstr_size: a k*1 vector with the size of each cluster
# clstr: small number of clusters k e.g.k=30
# lambda: annual event rate
# q: follow-up in years
# cmptrsk: a number of competing risk rate
# dropout: annual dropout rate
# tau_f0: a vector of frailty tau
# tau_c0: a vector of copula tau
# m: number of simulations
################################################################################
############################ Example Parameter Settings ########################
beta<-log(1)
clstr_info <- read.delim('STRIDE_cluster_sizes.txt', header = TRUE, sep = "\t", dec = ".")
clstr_size <- clstr_info$clustersize
clstr<-30
lambda<-0.08
q<-40/12
cmptrsk<-0.08
dropout<-0.03
tau_f0<-c(0.001)
tau_c0<-c(0.001,0.003)
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
usePackage("parallel")
usePackage("numbers")
usePackage("mvtnorm")
usePackage("survival")
usePackage("cmprsk")
usePackage("crrSC")
usePackage("dplyr")
usePackage("matrixcalc")
usePackage("prodlim")
#################################################################################
############################# Permutation Settings ##############################

# Simulation of Permutations
S=500   

pmt <- matrix(0, S, clstr) 

for (s in 1:S){
  trt <- sample(1:clstr, clstr / 2)
  pmt[s, trt] <- 1
}
pmt <- unique(pmt) # indicator matrix
R <- dim(pmt)[1]

#################################################################################

permutation <- function(beta, clstr_size, clstr=200,lambda, q, cmptrsk, dropout, tau_f, tau_c, m){
  dir.create("./permutations")
  theta<-1/(1-exp(q*log(1-dropout)))
  p2<-length(tau_f0)
  p1<-length(tau_c0)
  
  #Initialize the result
  smry<-NULL 
  smry2<-NULL
  fulldata<-NULL
  gamma<-cmptrsk
  
  #Start the simulation
  for (o2 in 1:p2){
    for (o1 in 1:p1){
      
      pb1<-txtProgressBar(min=0,max=m,style=2)
      for (j in 1:m) {
        result <- NULL
        
        #Setting all the parameters
        new_clstr_size <- sample(clstr_size,clstr,replace=T)
        n <- sum(new_clstr_size)
        
        #Adding the randomization for the treament, 
        #each cluster will randomly be given a cluster_id from 1 to clstr, 
        #and odd number will be given treatment
        ran_clstr <-sample(1:clstr)
        ran_x <- pmt[1,]
        (obs_trt <- row.match(ran_x, pmt))
        
        clusters <- rep(ran_clstr,new_clstr_size)
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
        
        #Start to construct the dataframe
        data<-data.frame(temp)
        colnames(data)<-c("u","w","z","v")
        data$id <- seq(1:n)
        data$clusters <- clusters
        
        #Kendall's tau, which is used to measure the dependence within cluster
        tau_f<-tau_f0[o2] 
        #If a is the shape parameter of gamma distribution, tau=1/(2*a+1)
        a<-0.5*(1/tau_f-1)
        #Let rate parameter=shape parameter in order to make mean of frailty 
        b<-a 
        #Generate the frailty for each cluster
        alpha<-rgamma(clstr,a,b) 
        
        #Re-assign treatment indicator to make it a cluster randomized design
        data$alpha <- alpha[data$clusters] 
        #Re-=assign the treatment group
        data$x <- ran_x[data$clusters] 
        
        #Simulate event time, competing risk time and censoring time
        data$time_event<--log(data$u)/(lambda*exp(beta*data$x)*data$alpha)
        data$time_cmpt<--log(data$v)/(gamma*data$alpha)
        data$time_event[data$time_event==0] <- 0.01
        data$time_cmpt[data$time_cmpt==0] <- 0.01
        data$time_c<-sapply(runif(n,0,theta),function(x){min(x,q)})
        data$time_t<-apply(data[,c("time_event","time_cmpt","time_c")],1,min, na.rm = T)
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
        
        
        #Permutate 500 times 
        for(p in 1:S){
          #Extract one scenario of treatment from the permutation 
          trt <- pmt[p,]
          #Re-assign the treatment group
          data$x <- trt[data$clusters] 
          
          #Change the wide format data to long format for the multistate models
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
          test.ms_cluster<-coxph(formula = Surv(start,end,as.factor(status))~x,id=id,cluster = clusters, data = data_l)
          test.cox_cluster<-coxph(Surv(data_w$time_t,data_w$status0)~data_w$x+cluster(data_w$clusters))
          test.zhou<-crrc(ftime=data_w$time_t,fstatus=data_w$status,cov1=data_w$x,cluster=data_w$clusters,maxiter = 100)
          
          
          #Summarize model
          q.95 <- -qnorm((1-0.95)/2)  # 95% multiplier
          ms_cluster.coef <- data.frame(Coefficient = test.ms_cluster$coefficients[1],
                                        SE = sqrt(test.ms_cluster$var[1,1]), 
                                        z_s = test.ms_cluster$coefficients[1]/sqrt(test.ms_cluster$var[1,1]),
                                        modelName = "ms_cluster")
          cox_cluster.coef <- data.frame(Coefficient = summary(test.cox_cluster)$coefficients[1],
                                         SE = summary(test.cox_cluster)$coefficients[4], 
                                         z_s = summary(test.cox_cluster)$coefficients[1]/summary(test.cox_cluster)$coefficients[4],
                                         modelName = "Marginal Cox")
          zhou.coef <- data.frame(Coefficient = test.zhou$coef,
                                  SE = sqrt(test.zhou$var),
                                  z_s = test.zhou$coef/sqrt(test.zhou$var),
                                  modelName = "Zhou")
          
          full.coef <- data.frame(rbind(ms_cluster.coef,cox_cluster.coef,zhou.coef),row.names = NULL)
          full.coef$tau_c = tau_c
          full.coef$tau_f = tau_f
          full.coef$lambda2 = gamma
          full.coef$HR = exp(beta)
          
          result <- data.frame(rbind(result,full.coef))
          
          setTxtProgressBar(pb1,j)
          
        }
        
        colnames(result)<-c("coefficient","SE","z_s","ModelName","tau_c","tau_f","lambda2","HR")
        close(pb1)
        
        #Combine every permutation from each simulation
        fulldata <- rbind(fulldata, result)
        
        #Summarize the combined full data result
        sum_info <- result %>%group_by(ModelName,tau_c,tau_f,lambda2,HR)%>%
          summarize(P_beta = mean(abs(coefficient)>=abs(coefficient[obs_trt])),
                    P_zscore = mean(abs(z_s)>=abs(z_s[obs_trt])))
        
        #Combine the smry
        smry<-rbind(smry,sum_info)
      }
    }
  }
  #Further summarize the result for each scenario
  sum_info2 <- smry %>% group_by(ModelName,tau_c,tau_f,lambda2,HR)%>%
    summarize(typeIerror_beta = mean(P_beta<0.05),
              typeIerror_z = mean(P_zscore<0.05))
  
  smry2 <- rbind(smry2,sum_info2)
  
  #output the data
  write.csv(smry2,file = paste0('./permutations/[p]full_cluster',clstr,"_hr",exp(beta),".csv"))
}


#################################################################################
#Call the function
permutation(beta=beta, clstr_size=clstr_size, clstr=clstr,lambda =lambda, 
          q = q, cmptrsk=cmptrsk, dropout=dropout, tau_f=tau_f, tau_c=tau_c, m=m)

