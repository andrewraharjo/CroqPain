file<-read.table("fin.data.txt",header=TRUE)

#See the correlation
cor(file)

#print variance accounted for
#comp1 for RHP, comp2 for NCR comp3 for NG comp 4 for ETW comp5 for RAFE
fit=princomp(file,cor=TRUE)
summary(fit)
loadings(fit) #principal component loadings

#---------------------------------------------------------PCA-------------#
my.classes = read.table("RHPvsRAFE.txt",header=TRUE)
plot(my.classes,cex=0.9,col="blue",main="RAFE vs RHP")
options(digits=3)
par(mfrow=c(1,1))

# Scale the data

standardize <- function(x) {(x - mean(x))}
my.scaled.classes = apply(my.classes,2,function(x) (x-mean(x)))
plot(my.scaled.classes,cex=0.9,col="blue",main="Plot of RHP vs. RAFE",sub="Mean Scaled",xlim=c(-200,500))

# Find Eigen values of covariance matrix

my.cov = cov(my.scaled.classes)
my.eigen = eigen(my.cov)
rownames(my.eigen$vectors)=c("RHP","RAFE")
colnames(my.eigen$vectors)=c("PC1","PC")
# Note that the sum of the eigen values equals the total variance of the data

sum(my.eigen$values)
var(my.scaled.classes[,1]) + var(my.scaled.classes[,2])

# The Eigen vectors are the principal components. We see to what extent each variable contributes

loadings = my.eigen$vectors

# Let's plot them 

pc1.slope = my.eigen$vectors[1,1]/my.eigen$vectors[2,1]
pc2.slope = my.eigen$vectors[1,2]/my.eigen$vectors[2,2]

abline(0,pc1.slope,col="red")
abline(0,pc2.slope,col="green")

textxy(12,10,"(-0.710,-0.695)")
textxy(-12,10,"(0.695,-0.719)")

# See how much variation each eigenvector accounts for

pc1.var = 100*round(my.eigen$values[1]/sum(my.eigen$values),digits=2)
pc2.var = 100*round(my.eigen$values[2]/sum(my.eigen$values),digits=2)
xlab=paste("PC1 - ",pc1.var," % of variation",sep="")
ylab=paste("PC2 - ",pc2.var," % of variation",sep="")

# Multiply the scaled data by the eigen vectors (principal components)

scores = my.scaled.classes %*% loadings
sd = sqrt(my.eigen$values)
rownames(loadings) = colnames(my.classes)

plot(scores,ylim=c(-10,10),main="Data in terms of EigenVectors / PCs",xlab=xlab,ylab=ylab)
abline(0,0,col="red")
abline(0,90,col="green")

# Correlation BiPlot

scores.min = min(scores[,1:2])
scores.max = max(scores[,1:2])

plot(scores[,1]/sd[1],scores[,2]/sd[2],main="BiPlot",xlab=xlab,ylab=ylab,type="n")
rownames(scores)=seq(1:nrow(scores))
abline(0,0,col="red")
abline(0,90,col="green")

# This is to make the size of the lines more apparent
factor = 5

# First plot the variables as vectors
arrows(0,0,loadings[,1]*sd[1]/factor,loadings[,2]*sd[2]/factor,length=0.1, lwd=2,angle=20, col="red")
text(loadings[,1]*sd[1]/factor*1.2,loadings[,2]*sd[2]/factor*1.2,rownames(loadings), col="red", cex=1.2)

# Second plot the scores as points
text(scores[,1]/sd[1],scores[,2]/sd[2], rownames(scores),col="blue", cex=0.7)

somelabs = paste(round(my.classes[,1],digits=1),round(my.classes[,2],digits=1),sep=" , ")
identify(scores[,1]/sd[1],scores[,2]/sd[2],labels=somelabs,cex=0.8)
#----------------------------------------PCAdone----------------------------------------


#------------------------------------------------Varcomp

library(lattice)
my.file <- read.table("fin.data.txt", header=TRUE)
library(gclus)
my.abs     <- abs(cor(my.file))
my.colors  <- dmat.color(my.abs)
my.ordered <- order.single(cor(my.file))
cpairs(my.file, my.ordered, panel.colors=my.colors, gap=0.5)

# Do the PCA 

my.prc <- prcomp(my.file, center=TRUE, scale=TRUE)
screeplot(my.prc, main="Scree Plot", xlab="Components")
screeplot(my.prc, main="Scree Plot", type="line" )

# DotPlot PC1

load    <- my.prc$rotation
sorted.loadings <- load[order(load[, 1]), 1]
myTitle <- "Loadings Plot for PC1" 
myXlab  <- "Variable Loadings"
dotplot(sorted.loadings, main=myTitle, xlab=myXlab, cex=1.5, col="red")

# DotPlot PC2

sorted.loadings <- load[order(load[, 2]), 2]
myTitle <- "Loadings Plot for PC2"
myXlab  <- "Variable Loadings"
dotplot(sorted.loadings, main=myTitle, xlab=myXlab, cex=1.5, col="red")

# Now draw the BiPlot
biplot(my.prc, cex=c(1, 0.7))

# Apply the Varimax Rotation
my.var <- varimax(my.prc$rotation)


#--------------------------------------------------

#compare 2 variables
x<-as.vector(as.matrix(file$RHP))
x1<-as.vector(as.matrix(file$NCR))
x2<-as.vector(as.matrix(file$NG))
x3<-as.vector(as.matrix(file$ETW))
y<-as.vector(as.matrix(file$RAFE))

prediction<-lm(y~x , data=file)
best_fitline=data.frame(file , fitted.value=fitted(prediction), residual=resid(prediction))
# Obtaining the confidence bands:
confidence_band=predict(prediction, interval="confidence")
# Obtaining the prediction bands:
predict_band=predict(prediction, interval="prediction")
# Create a new data frame containing the values of X at which we want the predictions to be made
pred.frame <- data.frame(x= seq(6, 16,by=1))
# Confidence bands
pc <- predict(predictio, int="c", newdata=pred.frame)
pp <- predict(prediction, int="p", newdata=pred.frame)
require(graphics)
plot(x , y , ylim=range(x , pp, na.rm=T))
pred.Size <- pred.frame$x
matlines(pred.Size , pc, lty=c(1,2,2), lwd=1.5, col=2)
matlines(pred.Size , pp, lty=c(1,3,3), lwd=1.5, col=1)

#See Confidence Intervals
#y=data$RAFE;
#x=data$RHP;
lm1 <- lm(y~x)
p_conf2 <- predict(lm1,interval="confidence")
plot(file$RAFE~file$ETW,RAFElim=c(0,60),RHPlim=c(0,650));## data

#Compare between RHP vs RAFE
final=function(x,y){
  n=length(x)
  yhat<-predict(lm(y~x))
  result<-lm(y~x)
  x1=summary.lm(result)
  mean_y=(sum(y)/n)
  
  #qtable
  t2=qt(.975,n-1)
  t1=qt(.975,n-2)
  
  #confidence interval
  confin<-confint(result,level=0.95)
  
  #Find beta0cap beta1cap
  beta0cap <- x1$coefficients[1, "Estimate"]
  beta1cap <- x1$coefficients[2, "Estimate"]
  
  #Find Sbeta0 and Sbeta1
  Sbeta0 <- x1$coefficients[1, "Std. Error"]
  Sbeta1 <- x1$coefficients[2, "Std. Error"]
  
  #Test Statistic
  #Beta1(0)
  T_1= (beta1cap-0)/Sbeta1
  
  #Beta0 (0)
  T_0=(beta0cap-0)/Sbeta0
  
  #Beta0= Ybar
  T_mean=(beta0cap-mean_y)/Sbeta0
  
  if(abs(T_1)>= t2){
    print("we reject the null hypotesis for beta1")
  } else {
    print("we accept null hypotesis for beta1")
  }

  if(abs(T_0)>=t2){    
    print("we reject the null hypotesis for beta0")
  } else {
    print("we accept null hypotesis for beta0")
  }
  
  if(abs(T_mean)>=t2){
    print("we reject the null hypotesis for beta0 = Ybar")
  } else {
    print("we accept null hypotesis for beta0 = Ybar")
  }
  
  final_result=list(beta1cap,beta0cap,confin,mean_y)
  return(final_result)
}
plot(x,y)
final(x,y)

#ANOVA
model1<-lm(y~x+x1+x2+x3)
model2<-lm(y~x+x1+x2)
model3<-lm(y~x+x1+x3)
model4<-lm(y~x3+x1+x)
model5<-lm(y~x1+x3+x)

summary(model1)
anova(model1)

summary(model2)
anova(model2)

summary(model3)
anova(model3)

summary(model4) 
anova(model4)

summary(model5)#best model so far
anova(model5)

confidence=function (x,y){
  result=lm(y~x)
  ci_interval=confint(result,level=0.95)
  ci_interval
}
