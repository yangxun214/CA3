#1.载入R包
library(survival)
library(plyr)
#2.清理工作环境
rm(list = ls()) 
#3.读入数据


setwd("E:/lab/CA3")
aa<- read.csv('clinical_meta-374.csv')

#4.查看数据数据性质
str(aa)
#5.查看结局，0=复发，1未复发
aa$event<-factor(aa$event)
summary(aa$event)


#=========单因素=========================
#1.构建模型的y，只需修改这一行代码和第3-(2)步代码
y<- Surv(time=aa$time,event=aa$event==1)#0为复发

#2.批量单因素回归模型建立：Uni_cox_model
Uni_cox_model<-function(x){
  FML<-as.formula(paste0("y~",x))
  Cox<-coxph(FML,data = aa) 
  Sum<-summary(Cox)
  CI<-paste0(round(Sum$conf.int[,3:4],3),collapse = "-") 
  Pvalue<-round(Sum$coefficients[,5],3)
  HR<-round(Sum$coefficients[,2],3)
  Unicox<-data.frame("Characteristics"=x,
                     "Hazard Ratio"=HR,
                     "CI95"=CI,
                     "P value"=Pvalue)
  return(Unicox)
}

#3.将想要进行的单因素回归变量输入模型
#3-(1)查看变量的名字和序号
names(aa)
#3-(2)输入变量序号
variable.names<- colnames(aa)[c(5:15)] #例：这里选择了3-10号变量

#4.输出结果
Uni_cox <- lapply(variable.names, Uni_cox_model)
Uni_cox<- ldply(Uni_cox,data.frame)


#5.优化表格，这里举例HR+95% CI+P 风格
Uni_cox$CI<-paste(Uni_cox$CI5,'-',Uni_cox$CI95)
Uni_cox<-Uni_cox[,-5]

#查看单因素cox表格
View(Uni_cox)


#=========多因素=========================

#1.提取单因素p<0.05变量
Uni_cox$Characteristics[Uni_cox$P.value<0.05]

#2.多因素模型建立
mul_cox_model<- as.formula(paste0 ("y~",
                                   paste0(Uni_cox$Characteristics[Uni_cox$P.value<0.05],
                                          collapse = "+")))
mul_cox<-coxph(mul_cox_model,data=aa)
cox4<-summary(mul_cox) 

#3.提取多因素回归的信息
mul_HR<- round(cox4$coefficients[,2],3) 
mul_PValue<- round(cox4$coefficients[,5],3) 
mul_CI1<-round(cox4$conf.int[,3],3)
mul_CI2<-round(cox4$conf.int[,4],3)

#4.多因素结果优化并成表：mul_cox1
mul_CI<-paste(mul_CI1,'-',mul_CI2)
mul_cox1<- data.frame("HR"=mul_HR,"CI"=mul_CI, "P"=mul_PValue)

write.csv(Uni_cox,'CA3-cli-Uni_cox.csv') 
write.csv(mul_cox1,'CA3-cli-mul_cox1.csv') 

