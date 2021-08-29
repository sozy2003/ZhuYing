# Ying Zhu, so_zy2003@126.com
# 2021-07


# clean environment variables
rm(list=ls()) 

# load related packages
if(!require(ggplot2))install.packages("ggplot2")
library("ggplot2") 

# Set ggplot2 ploting parameter
main_theme <-  theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=7),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=7),
                   text=element_text(family="sans", size=7))


#  Design of experiment
design <-  read.table("metadata.txt", header=T, row.names= 1, sep="\t") 

#bray_curtis distance
betadiv <-  read.table("bray_curtis.txt", header=T, row.names= 1, sep="\t")

#reorde the row 
betadiv = betadiv[rownames(design), rownames(design)]

betadis <- as.matrix(betadiv)

betadis[upper.tri(betadis,diag = T)] <- NA


# subset for each time point

day <- unique(sort(design$Day))

p=list() #save the dataframe in the variable p
for (i in day){
  p[[i]] <- betadis[design$Day==i, design$Day==i]
  nrow <- dim(p[[i]])[1]
  ncol <- dim(p[[i]])[2]
  row <- rep(rownames(p[[i]]), ncol)
  col <- rep(colnames(p[[i]]),each=nrow)
  p[[i]] <- data.frame(row, col, value=as.numeric(p[[i]]))
  p[[i]]$day <- rep(i,nrow(p[[i]]))
  
}
p

# combine the data frame from each time point.

frame <- rbind(p[[1]],p[[3]],p[[6]],p[[9]],p[[12]],p[[15]],
               p[[18]],p[[21]], p[[24]], p[[27]],p[[30]],p[[37]],
               p[[44]])

# delete missing data
frame <- na.omit(frame)

#boxplot to check pattern
p=ggplot(frame,aes(x=as.factor(day),y=value))+
  geom_boxplot()+
  labs(x="Day", y="bray_curtis") + theme_classic()
p=p+geom_jitter(position=position_jitter(0.17), size=0.001, alpha=0.3,color='blue')+
  geom_smooth(method = "lm", formula = y ~ poly(x,3))
p=p+main_theme
p.boxplot=p
p.boxplot
# ggsave(paste("beta_box_bray_curtis",".pdf", sep=""), p.boxplot, width = 4, height =2.5 )



#scatter point
p=ggplot(frame,aes(x=day,y=value))+
  labs(x="Day", y="Bray-Curtis Distance") + 
  geom_point(size=1,alpha=1)+main_theme
p=p+geom_smooth(method = "lm", formula = y ~ poly(x,3))

p.point=p
p.point
# ggsave(paste("beta_scatter_bray_curtis",".pdf", sep=""), p.point, width = 4, height =2.5 )



#GAM
if(!require("mgcv"))install.packages("mgcv")

#使用以下模型进行拟合
#安装R包
# package <- c("segmented", "splines", "Hmisc", "rms", "mgcv", "caret")

# if(!require("rms"))install.packages("rms")
# if(!require("caret"))install.packages("caret")
# library(segmented)
# library(splines)
# library(Hmisc)
# library(rms)
library(mgcv)
# library(caret)

#参考
#https://mp.weixin.qq.com/s?__biz=MzI2OTQyMzc5MA==&mid=2247497698&idx=2&sn=8245fed6dc19f8653cc115420b36bdf6&chksm=eae23223dd95bb35f72774371d1f6cf66da3fc5f9fd5acf9e858f19fbdfaad270e772696b05b&mpshare=1&scene=1&srcid=1209XoHstr8nzT5tIqis8VF1&sharer_sharetime=1607513547446&sharer_shareid=9361991bf2e30b686dcc15ce09f07200&key=2f88c2a11d638eeafdd76a9dd2d50c3822d093e373d2a6efef8bf5aa265836fc4344a0d255738c76158124dd2c1ac6aa2d667dcfa241cf3d4b73f8bf229248a5c33edf55b069733dee8a68ef91debacadf3aa74aa0307a388b46e6b49a4d2d82f305127fd4b8e75b8d5229a25c8d96aa0e2016bfdeef9a6179bd95f4e008880c&ascene=1&uin=NTMwNzU3NTAw&devicetype=Windows+10+x64&version=6300002f&lang=zh_CN&exportkey=A5RvloGdI%2FJwpR21OCce5So%3D&pass_ticket=eZl8prOEAdh9VZz0W3qBFHhPF4fEnLnJpfErZ5SmVWzglgalirHTwMScQm3%2FTCCj&wx_header=0
#线性拟合
model.lm <- lm(value~day, data=frame)
summary(model.lm)
#0.1959, -0.0003201 

#曲线方程
model.log <- lm(value~log(day), data=frame)
summary(model.log)
#0.1953, 0.006377

#建立分段回归
model.segmented <- segmented(model.lm)
summary(model.segmented)
#0.1843,0.1151

#样条回归
model.spline <- lm(value~rcs(day,c(10,20,30)),data=frame)
summary(model.spline)
#0.1898, 0.06094


#GAM
model.gam <- gam(value~s(day),data=frame)
summary(model.gam) 

#没有给出残差标准误，用预测值和实际值进行比较
pr.gam <- predict(model.gam,frame)#生成预测值

#计算RMSE和R方
data.frame(RMSE=RMSE(pr.gam,
                     frame$value),
           R2=R2(pr.gam, frame$value))
#0.1813905, 0.1415742
#尝试散点拟合
ggplot(frame, aes(day, log(value)))+
  geom_point()+
  stat_smooth(method=gam, formula = y~s(x))


# fit with mean value 
library(dplyr)
temp <- frame%>%
  group_by(day)%>%
  summarize(mean(value))

temp <- as.data.frame(temp)
colnames(temp) <- c("day", "mean")

#scatter point of mean value
p=ggplot(temp,aes(x=day,y=mean))+
  labs(x="Day", y="mean values of Bray_curtis distance") + theme_classic()+
  geom_point(size=1,alpha=1)
p1=p+geom_smooth(method = "lm", formula = y ~ x)
p2=p+geom_smooth(method = "lm", formula = y ~ poly(x,3))
p1
p2
p.mean.point=p # good 

# ggsave(paste("beta_mean_scatter_bray_curtis_linear",".pdf", sep=""), 
#        p1, width = 5, height =4 )
# 
#  ggsave(paste("beta_mean_scatter_bray_curtis",".pdf", sep=""), 
#         p2, width = 5, height =4 )
