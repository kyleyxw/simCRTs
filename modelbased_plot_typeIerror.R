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
output_path <- args[2]

full <- read.csv(full_data)

if(any(full0$HR)!=1){stop("This Script is suitable for situation where Hazard Ratio equals to 1")}


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
result<-full[,-c(1,2)]
colnames(result)[9]<-"dropout"
result$coverage1[result$coverage=="TRUE"]<-1
result$coverage1[result$coverage=="FALSE"]<-0
result$coverage1[result$coverage=="FALSE"]<-0
result<-result[,-3]
colnames(result)[12]<-"coverage"
result$coverage<-as.numeric(result$coverage)
sum_info <- result %>%group_by(modelName,tau1,tau2,lambda1,lambda2,dropout,HR)%>%
  summarize(avg_bias = mean(coefficient), 
            avg_SE =mean(SE),coverage = mean(coverage), 
            sample_SE = sqrt(var(coefficient)/length(coefficient)))

sub <- sum_info
sub$rela_bias<-sub$mean_beta
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

# Type I Error
ggplot(sub,aes(colour = factor(labels),x = lambda2,y=1-coverage,fill = factor(labels),pch=factor(labels)))+
  geom_point(size = 4,position = "jitter")+
  geom_hline(yintercept =0.05, colour = "black", lty = 2,lwd=0.5)+
  geom_vline(xintercept = c(0.02,0.04,0.08,0.12), colour = "black", lty = 2,lwd=0.5)+
  theme_test() +scale_y_continuous(breaks = c(0,0.05,0.10,0.20,0.30,0.40,0.50))+
  scale_x_continuous(breaks = c(0.02,0.04,0.08,0.12))+
  labs(x="Competing Risk Rate",y="Type I Error")+
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

ggsave(paste(output_path,"type_I_error.png"),width = 10,height = 10)