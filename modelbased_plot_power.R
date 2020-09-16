# Packages Setup 
packages <- c('dplyr','extrafont','ggplot2','readr')
for(p in packages){
  if(!require(p,character.only = T)){
    install.packages(p)
  }
}



# I/O Setup
args <- commandArgs(TRUE)
full_data <- args[1]
true_data <- args[2]
output_path <- args[3]

full0 <- read.csv(full_data)
true0 <- read.csv(true_data)

if(any(full0$HR)==1){stop("This Script is suitable for situation where Hazard Ratio is not equal to 1")}


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


# Data Manipulation
q.95 <- -qnorm((1-0.95)/2)  # 95% multiplier

# True Data
truebeta<-true0[,-c(1,5,7,9)]
truebeta1<-truebeta
truebeta1$modelName1[truebeta$modelName=="Cox"]<-"Marginal Cox"
truebeta1$modelName1[truebeta$modelName=="Fine and Gray"]<-"Zhou"
truebeta1$modelName1[truebeta$modelName=="ms"]<-"ms_cluster"
truebeta1$modelName<-truebeta1$modelName1
truebeta1<-truebeta1[,-6]
truebeta0<-rbind(truebeta,truebeta1)
full1 <- left_join(full0,truebeta0)
full1$coverage <- full1$coefficient -q.95*full1$SE<full1$beta & full1$coefficient +q.95*full1$SE>full1$beta

result<-full0[,3:15]
result$reject<-as.numeric(result$pvalue<=0.05)
sum_info <- result %>%group_by(modelName,tau1,tau2,lambda2)%>%summarise(mean_beta = mean(coefficient), 
                                                                             avg_SE =mean(SE),coverage = mean(coverage), 
                                                                             sample_SE = sd(coefficient),
                                                                             power = mean(reject))


sub <- left_join(sum_info,truebeta0)
sub$rela_bias<-(sub$mean_beta-sub$beta)/sub$beta
sub$SE_ratio<-sub$sample_SE/sub$avg_SE
sum_info<-sub
sub$modelName <- as.factor(sub$modelName)

levels(sub$modelName) <- c("Cox Model","Fine and Gray", "Marginal Cox Model","Multi-State Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray")

sub$cluster <- sub$modelName%in%c("Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray")
sub$cluster <- as.integer(sub$cluster)
sub$cluster <- factor(sub$cluster)
levels(sub$cluster) <- c("Without Clustering","With Clustering")
sub$labels <- rep(c(1,3,4,2,5,6),each=nrow(sub)/6)
sub <- sub[order(sub$labels),]
sub$tau2label <- sapply(sub$tau2,label_tau2)
sub$tau1label <- sapply(sub$tau1,label_tau1)

# Relative Bias
ggplot(sub,aes(colour = factor(labels),x = lambda2,y=rela_bias,fill = factor(labels),pch=factor(labels)))+
  geom_point(size = 4,position = "jitter")+
  geom_hline(yintercept =0, colour = "black", lty = 2,lwd=0.5)+
  geom_vline(xintercept = c(0.02,0.04,0.08,0.12), colour = "black", lty = 2,lwd=0.5)+ 
  theme_test() +
  scale_x_continuous(breaks = c(0.02,0.04,0.08,0.12))+
  #theme_gray(base_size = 20)+
  labs(x="Competing Risk Rate",y="Relative Bias")+
  scale_colour_manual(name = "Method",
                      labels = c("Cox Model","Multi-State Cox Model","Fine and Gray", "Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                      values = c("#FF0000", "Black", "#0000FF","#FF0000", "Black", "#0000FF")) +  
  scale_fill_manual(name = "Method",
                    labels =c("Cox Model","Multi-State Cox Model","Fine and Gray", "Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                    values = c("#FF0000", "Black", "#0000FF","#FF0000", "Black", "#0000FF")) +
  scale_shape_manual(name = "Method",
                     labels = c("Cox Model","Multi-State Cox Model","Fine and Gray", "Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                     values = c(19,19,19,17,17,17))+
  theme(plot.title = element_text(hjust=0.5,color="Black",size=20,face="bold",family="Calibri"),
        axis.title = element_text(color="Black",size=20,face="bold",family="Calibri"),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(color="Black",size=17,face="bold",family="Calibri"),
        legend.text = element_text(color="Black",size=17,family="Calibri"),
        legend.position="top",legend.box = "horizontal",
        strip.text = element_text(color="Black",size=14,family="Calibri"))+
  guides(colour=guide_legend(nrow=2,byrow=TRUE))+
  facet_grid(tau1label~tau2label, scale="fixed",labeller = label_parseall)

ggsave(paste(output_path,"relative_bias.png"),width = 10,height = 10)


# Standard Error Ratio
ggplot(sub,aes(colour = factor(labels),x = lambda2,y=SE_ratio,fill = factor(labels),pch=factor(labels)))+
  geom_point(size = 4,position = "jitter")+
  #geom_hline(yintercept =0.95, colour = "black", lty = 2,lwd=0.5)+
  geom_vline(xintercept = c(0.02,0.04,0.08,0.12), colour = "black", lty = 2,lwd=0.5)+
  theme_test() +
  #scale_y_continuous(breaks = c(0,0.50,0.60,0.70,0.80,0.90,0.95))+
  scale_x_continuous(breaks = c(0.02,0.04,0.08,0.12))+
  labs(x="Competing Risk Rate",y="Standard Error Ratio")+
  scale_colour_manual(name = "Method",
                      labels = c("Cox Model","Multi-State Cox Model","Fine and Gray", "Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                      values = c("#FF0000", "Black", "#0000FF","#FF0000", "Black", "#0000FF")) +  
  scale_fill_manual(name = "Method",
                    labels =c("Cox Model","Multi-State Cox Model","Fine and Gray", "Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                    values = c("#FF0000", "Black", "#0000FF","#FF0000", "Black", "#0000FF")) +
  scale_shape_manual(name = "Method",
                     labels = c("Cox Model","Multi-State Cox Model","Fine and Gray", "Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                     values = c(19,19,19,17,17,17))+
  theme(plot.title = element_text(hjust=0.5,color="Black",size=20,face="bold",family="Calibri"),
        axis.title = element_text(color="Black",size=20,face="bold",family="Calibri"),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(color="Black",size=17,face="bold",family="Calibri"),
        legend.text = element_text(color="Black",size=17,family="Calibri"),
        legend.position="top",legend.box = "horizontal",
        strip.text = element_text(color="Black",size=14,family="Calibri"))+
  guides(colour=guide_legend(nrow=2,byrow=TRUE))+
  facet_grid(tau1label~tau2label, scale="fixed",labeller = label_parseall)

ggsave(paste(output_path,"standard_error_ratio.png"),width = 10,height = 10)

# Power
ggplot(sub,aes(colour = factor(labels),x = lambda2,y=power,fill = factor(labels),pch=factor(labels)))+
  geom_point(size = 4,position = "jitter")+
  #geom_hline(yintercept =0.95, colour = "black", lty = 2,lwd=0.5)+
  geom_vline(xintercept = c(0.02,0.04,0.08,0.12), colour = "black", lty = 2,lwd=0.5)+
  theme_test() +scale_y_continuous(breaks = c(0.90,0.95,1))+
  scale_x_continuous(breaks = c(0.02,0.04,0.08,0.12))+
  labs(x="Competing Risk Rate",y="Power")+
  scale_colour_manual(name = "Method",
                      labels = c("Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                      values = c("#FF0000", "Black", "#0000FF","#FF0000", "Black", "#0000FF")) +  
  scale_fill_manual(name = "Method",
                    labels =c( "Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                    values = c("#FF0000", "Black", "#0000FF","#FF0000", "Black", "#0000FF")) +
  scale_shape_manual(name = "Method",
                     labels = c( "Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                     values = c(17,17,17))+
  theme(plot.title = element_text(hjust=0.5,color="Black",size=20,face="bold",family="Calibri"),
        axis.title = element_text(color="Black",size=20,face="bold",family="Calibri"),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(color="Black",size=17,face="bold",family="Calibri"),
        legend.text = element_text(color="Black",size=17,family="Calibri"),
        legend.position="top",legend.box = "horizontal",
        strip.text = element_text(color="Black",size=14,family="Calibri"))+
  guides(colour=guide_legend(nrow=1,byrow=TRUE))+
  facet_grid(tau1label~tau2label, scale="fixed",labeller = label_parseall)

ggsave(paste(output_path,"power.png"),width = 10,height = 10)

# Coverage
ggplot(sub,aes(colour = factor(labels),x = lambda2,y=coverage,fill = factor(labels),pch=factor(labels)))+
  geom_point(size = 4,position = "jitter")+
  geom_hline(yintercept =0.95, colour = "black", lty = 2,lwd=0.5)+
  geom_vline(xintercept = c(0.02,0.04,0.08,0.12), colour = "black", lty = 2,lwd=0.5)+
  theme_test() +scale_y_continuous(breaks = c(0,0.50,0.60,0.70,0.80,0.90,0.95))+
  scale_x_continuous(breaks = c(0.02,0.04,0.08,0.12))+
  labs(x="Competing Risk Rate",y="Coverage")+
  scale_colour_manual(name = "Method",
                      labels = c("Cox Model","Multi-State Cox Model","Fine and Gray", "Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                      values = c("#FF0000", "Black", "#0000FF","#FF0000", "Black", "#0000FF")) +  
  scale_fill_manual(name = "Method",
                    labels =c("Cox Model","Multi-State Cox Model","Fine and Gray", "Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                    values = c("#FF0000", "Black", "#0000FF","#FF0000", "Black", "#0000FF")) +
  scale_shape_manual(name = "Method",
                     labels = c("Cox Model","Multi-State Cox Model","Fine and Gray", "Marginal Cox Model","Marginal Multi-State Cox Model","Marginal Fine and Gray"),
                     values = c(19,19,19,17,17,17))+
  theme(plot.title = element_text(hjust=0.5,color="Black",size=20,face="bold",family="Calibri"),
        axis.title = element_text(color="Black",size=20,face="bold",family="Calibri"),
        axis.text.x = element_text(angle = 90),
        legend.title = element_text(color="Black",size=17,face="bold",family="Calibri"),
        legend.text = element_text(color="Black",size=17,family="Calibri"),
        legend.position="top",legend.box = "horizontal",
        strip.text = element_text(color="Black",size=14,family="Calibri"))+
  guides(colour=guide_legend(nrow=2,byrow=TRUE))+
  facet_grid(tau1label~tau2label, scale="free_y",labeller = label_parseall)

ggsave(paste(output_path,"coverage.png"),width = 10,height = 10)
