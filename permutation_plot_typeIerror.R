# Packages Setup 
packages <- c('dplyr','extrafont','ggplot2','readr')
for(p in packages){
  if(!require(p,character.only = T)){
    install.packages(p)
  }
}


# I/O Setup
args <- commandArgs(TRUE)
full_Permutation <- args[1]
full_Modelbased <- args[2]
output_path <- args[3]

full_P <- read.csv(full_Permutation)
full_M <- read.csv(full_Modelbased)

if(any(full_P$HR)!=1){stop("This Script is suitable for situation where Hazard Ratio equals to 1")}

# Label functions
label_tau1 <- function(x){
  paste("tau[w]==",x)
}
label_tau2 <- function(x){
  paste("tau[b]==",x)
}

label_parseall <- function(variable, value) {
  plyr::llply(value, function(x) parse(text = paste(x, sep = "==")))
}

# Type I Error-Permutaion: Beta vs. Z
sum_info <- full_P%>% group_by(ModelName,tau_c,tau_f,lambda2,HR)%>%
  summarize(typeIerror_beta = mean(typeIerror_beta),
            typeIerror_z = mean(typeIerror_z))
sum_info$tau2label <- sapply(sum_info$tau_f,label_tau2)
sum_info$tau1label <- sapply(sum_info$tau_c,label_tau1)
sum_info1 <- sum_info[,-7]
sum_info2 <- sum_info[,-6]
colnames(sum_info1)[6]<- "typeIerror"
sum_info1$label<- "Beta"
colnames(sum_info2)[6]<- "typeIerror"
sum_info2$label<- "Z"
sum_info0<-rbind(sum_info1,sum_info2)

ggplot(sum_info0)+
  geom_point(aes(colour = factor(ModelName),x = lambda2,y=typeIerror,fill = factor(ModelName),pch=factor(label)),size = 4,position ="jitter")+
  geom_hline(yintercept =c(0.036,0.064), colour = "black", lty = 2,lwd=0.5)+
  geom_vline(xintercept = c(0.02,0.04,0.08,0.12), colour = "black", lty = 2,lwd=0.5)+
  theme_test() +
  scale_y_continuous(breaks = c(0.036,0.05,0.064))+
  scale_x_continuous(breaks = c(0.02,0.04,0.08,0.12))+
  labs(x="Competing Risk Rate",y="Type I Error")+
  scale_colour_manual(name = "Method",
                      labels = c("Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                      values = c("red", "blue", "black")) +  
  scale_fill_manual(name = "Method",
                    labels =c( "Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                    values = c("red", "blue", "black")) +
  scale_shape_manual(name = "Type I Error",values = c("Beta"=19,"Z"=17))+
  theme(plot.title = element_text(hjust=0.5,color="Black",size=20,face="bold",family="Calibri"),
        axis.title = element_text(color="Black",size=20,face="bold",family="Calibri"),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(color="Black",size=17,face="bold",family="Calibri"),
        legend.text = element_text(color="Black",size=17,family="Calibri"),
        legend.position="top",legend.box = "vertical",
        strip.text = element_text(color="Black",size=14,family="Calibri"))+
  guides()+
  facet_grid(tau1label~tau2label, scale="fixed",labeller = label_parseall)

ggsave(paste(output_path,"type_I_error-permutation.png"),width = 10,height = 10)

# Type I Error: Comparison with Model-based Test
sum_info1$label<- "Permutation"
sum_info1<-sum_info1[,-5]

result<-full_M[,-c(1,2)]
colnames(result)[9]<-"dropout"
result$coverage1[result$coverage=="TRUE"]<-1
result$coverage1[result$coverage=="FALSE"]<-0

result<-result[,-3]
colnames(result)[12]<-"coverage"
result$coverage<-as.numeric(result$coverage)
mBase <- result %>%group_by(modelName,tau1,tau2,lambda1,lambda2,dropout,HR)%>%
  summarize(avg_bias = mean(coefficient), 
            avg_SE =mean(SE),coverage = mean(coverage), 
            sample_SE = sqrt(var(coefficient)/length(coefficient)))


mBase0 <- mBase[mBase$modelName%in%c("ms_cluster","Marginal Cox","Zhou"),]
mBase1 <- mBase0[,c(1,2,3,5,10)]
mBase1$tau2label <- sapply(mBase1$tau1,label_tau2)
mBase1$tau1label <- sapply(mBase1$tau2,label_tau1)
mBase1$label <- "Model-Based"
mBase1$coverage <- 1-mBase1$coverage
colnames(mBase1) <- colnames(sum_info1)
mBase1$ModelName <- rep(c("Marginal Cox", "ms_cluster", "Zhou"),nrow(mBase1)/3)
mBase1$ModelName <- factor(mBase1$ModelName)
mBase1 <- as.data.frame(mBase1)
sum_info1 <- as.data.frame(sum_info1)
sum_info0 <- rbind(sum_info1,mBase1)

ggplot(sum_info0)+
  geom_point(aes(colour = factor(ModelName),x = lambda2,y=typeIerror,fill = factor(ModelName),pch=factor(label)),size = 4,position ="jitter")+
  #geom_hline(yintercept =c(0.036,0.064), colour = "black", lty = 2,lwd=0.5)+
  geom_vline(xintercept = c(0.02,0.04,0.08,0.12), colour = "black", lty = 2,lwd=0.5)+
  theme_test() +
  scale_y_continuous(breaks = c(0.036,0.05,0.064,0.1,0.15))+
  scale_x_continuous(breaks = c(0.02,0.04,0.08,0.12))+
  labs(x="Competing Risk Rate",y="Type I Error")+
  scale_colour_manual(name = "Method",
                      labels = c("Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                      values = c("red", "blue", "black")) +  
  scale_fill_manual(name = "Method",
                    labels =c( "Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                    values = c("red", "blue", "black")) +
  scale_shape_manual(name = "Type I Error",values = c("Permutation"=19,"Model-Based"=17))+
  theme(plot.title = element_text(hjust=0.5,color="Black",size=20,face="bold",family="Calibri"),
        axis.title = element_text(color="Black",size=20,face="bold",family="Calibri"),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(color="Black",size=17,face="bold",family="Calibri"),
        legend.text = element_text(color="Black",size=17,family="Calibri"),
        legend.position="top",legend.box = "vertical",
        strip.text = element_text(color="Black",size=14,family="Calibri"))+
  guides()+
  facet_grid(tau1label~tau2label, scale="fixed",labeller = label_parseall)

ggsave(paste(output_path,"type_I_error-permutationVSmodelbased.png"),width = 10,height = 10)