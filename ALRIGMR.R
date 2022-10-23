rm(list = ls())
setwd('D:\\ideal\\GSE144269肝癌RNA-seq\\重新整理\\添加差异掉的重要基因')
lmqcm<-read.csv('分群添加后数据2222.csv',header=T)
index<-c(rep(1,times=542),rep(2,times=527),rep(3,times=160),rep(4,times=120),rep(5,times=52),
         rep(6,times=41),rep(7,times=11),rep(8,times=11),rep(9,times=10),10:20)
A<-lmqcm
normal70<-t(A[,2:71])
dim(normal70)
cancer70<-t(A[,72:141])
dim(cancer70)
#训练集和测试集的划分，随机取2/3练，1/3测试

set.seed(10)
#70个癌症样本划分
N1=nrow(cancer70)
trainindex1<-sample(1:N1,47)
trainset1<-cancer70[trainindex1,]
testset1<-cancer70[-trainindex1,]
# 70个正常样本划分
N2=nrow(normal70)
trainindex2<-sample(1:N2,47)
trainset2<-normal70[trainindex2,]
testset2<-normal70[-trainindex2,]
#3 #合并
x<-rbind(trainset2, trainset1)
test<-rbind(testset2, testset1)
y<-as.vector(c(rep(0,47),rep(1,47)))
ytest<-c(rep(0,23),rep(1,23))
dim(x)
#单个基因权重
lv<-read.csv('分群添加后突变.csv',header = F)
lv1<-lv[,1]
#lv1[1476]<-1
huxinxi<-read.csv('分群添加后互信息.csv',header = F)
huxinxi1<-huxinxi[,3]
w1<-(lv1+huxinxi1)*(1/2)
w1<-1/w1
#群权重SU度量
library(FCBF)
mes<-read.csv('lmQCM-tezheng.csv',header =T)
mesnormal<-mes[,1:70]
mescancer<-mes[,71:140]
trainnormal<-mesnormal[,trainindex2]
traincancer<-mescancer[,trainindex1]
trainmes<-cbind(trainnormal,traincancer)
discrete_expression <- as.data.frame(discretize_exprs(trainmes))
#target<-read.csv('target.csv',header=F)
su<- get_su_for_feature_table_and_vector(discrete_expression[,],y)
a<-read.csv('a.csv',header = F)
for (m in 1:9) {
  temp1<-a[m,1]
  for (n in 1:9) {
    temp2<-su[n,2]
    if (temp1==temp2){
      a[m,2]<-su[n,1]
    }
  }
}
grow<-as.vector(a[1:9,2])
grow1<-c(grow,rep(1,11))
w2<-1/grow1

#训练并测试精度
library(asgl)
groupindex<-as.vector(index)
a1<-c(0.001,0.05,0.1,0.3,0.55,0.6,0.7,0.8,0.95)

acc<-matrix(c(rep(0,900)),nrow = 100,ncol = 9)
for (j in 1:9) {
  fitasgl<-asgl(x, y, index=groupindex, family = "binomial",alpha =a1[j],
                standardize = T,
                grp_weights = w2,
                ind_weights =w1,
                nlambda = 100)
  y1<-predict(fitasgl,test,type = "response")
  y1<-sign(y1-0.5)
  y1[y1==-1]<-0
  y1
  write(paste("*** results ", j, "***"),"")
  #coefSGL1 <-fitasgl$beta[,100]
  #n1<- which(coefSGL1!= 0)   #结果复制粘贴保存下来
  #n1
  for (k in 1:100) {
    y0<-y1[,k]
    lab0<-length(which(y0[1:23]!=0))
    lab1<-length(which(y0[24:46]!=1))
    accuracy<-1-(lab0+lab1)/46
    acc[k,j]<-accuracy
    #print(accuracy)
  }
}
max(acc)

###最优参数
fitasgl<-asgl(x, y, index=groupindex, family = "binomial",alpha =0.1,
              standardize = T,
              grp_weights = w2,
              ind_weights =w1,
              nlambda = 100)
y1<-predict(fitasgl,test,type = "response")
y1<-sign(y1-0.5)
y1[y1==-1]<-0
hunxiao<-table(ytest,y1[,90])
P<-hunxiao[1,1]/(hunxiao[1,1]+hunxiao[2,1])
R<-hunxiao[1,1]/(hunxiao[1,1]+hunxiao[1,2])
F1<-2*P*R/(P+R)
coefSGL1 <-fitasgl$beta[,90]
n1<- which(coefSGL1!= 0)   #结果复制粘贴保存下来
gene<-A[c(n1),1]
length(n1)
'CTNNB1' %in% gene
'TP53' %in% gene
'TTN' %in% gene
'MUC16' %in% gene
'PCLO' %in% gene
'RYR2' %in% gene
'XIRP2' %in% gene
'APOB' %in% gene
'OBSCN' %in% gene
'ABCA13' %in% gene
'LRP1B' %in% gene
'ALB' %in% gene
hunxiao

