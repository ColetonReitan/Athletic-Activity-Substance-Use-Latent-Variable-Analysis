library(descr)
library(Hmisc)
library(gmodels) # crosstabs
library(psych) #EFA, etc
library(lavaan) #SEM
library(stringr) # diagramming
library(DiagrammeR) # diagramming
library(dplyr) # data management
library(semPlot) # plots and diagrams
library(plotrix)
library(GPArotation)

# Import a csv file
mydataset <- read.csv("C:/Users/colet/Documents/Wesleyan QAC/QAC 313 Latent Variable Analysis/Athletic Involvement Data.csv") 


# select a subset of variables to exclude observations with missing data
colnames(mydataset) <- tolower(colnames(mydataset))
names(mydataset)

var.keep <- c("age", "excur", "pacur", "spcur", "cesd10", "cesd09", "cesd05", "dai01", "sai06", "dai04")
data2 <- mydataset[ , var.keep]

data2[data2==-9] <- NA
data2[data2==-8] <- NA
data2[data2==-7] <- NA
data2$spcur[data2$spcur == 1] <- NA
data3 <- na.omit(data2)
data3$sai06[data3$sai06 == 8] <- 0
data3$dai01[data3$dai01 == 8] <- 0
data3$dai02[data3$dai02 == 8] <- 0
data3$dai04[data3$dai04 == 8] <- 0

names(data3)

# Random 50/50 split data set into development and validation sets
set.seed(12345)
samp = sort(sample(nrow(data3), nrow(data3)*.5))
data10<-data3[samp,]
#data10$cesd05r<-4-data10$cesd05
data20<-data3[-samp,]

# Frequency tables
#freq(data10$mpsq06, plot=F)
#freq(data10$mpsq11, plot=F)
freq(data10$pshro, plot=F)
freq(data10$psyears, plot=F)
freq(data10$spcur, plot=F)
freq(data10$cesd10, plot=F)
freq(data10$cesd09, plot=F)
freq(data10$cesd05, plot=F)
freq(data10$dai01, plot=F)
freq(data10$sai06, plot=F)
freq(data10$dai04, plot=F)
freq(data10$dai02, plot=F)

# Means for all variables
colMeans(data10,na.rm=TRUE)


# correlation heatmap
cor.plot(data10,numbers=TRUE)

fa4<-fa(data10, nfactors=4, rotate="promax")
fa5<-fa(data10, nfactors=5, rotate="promax")
fa6<-fa(data10, nfactors=6, rotate="promax")
fa4
fa5
fa6
fa3<-fa(data10, nfactors=3, rotate="promax")
fa3
fa2<-fa(data10, nfactors=2, rotate="promax")
fa2



# eigenvalues & scree plot
fa3$e.values
plot(fa3$e.values)
plot(fa3$e.values, pch=16, col = "blue", xlim=c(1,10),
     xlab="Factor", ylab="Eigenvalue", main="Athletic Involvement EFA Scree Plot")
lines(fa3$e.values, col = "blue")
abline(h=1, col="red")

# view sorted loadings for 5 factor solution
print.psych(fa4, cut=.3, digits=2, sort=TRUE)
# view sorted loadings for 3 factor solution
print.psych(fa3, cut=.3, digits=2, sort=TRUE)
# confidence intervals for factor loadings and correlations
fa3c<-fa(data10, nfactors=3, rotate="promax", n.iter=100)
fa3c

# coefficient alpha
sub1a <- subset(data10, select = excur:spcur)
sub1b <- subset(data10, select = cesd10:cesd05)
sub1c <- subset(data10, select = dai01:dai04)

alpha(sub1a, check.keys=TRUE)
alpha(sub1b, check.keys=TRUE)
alpha(sub1c, check.keys=TRUE)




# define factor analysis model (enclose model definition in 
# single quotes)
model <- '
# latent variable definitions
 
substance_use =~ dai01 + dai04 + sai06
depressed =~ cesd10 + cesd09 + cesd05
athletic_activity =~  excur + pacur + spcur

# covariance between factors
substance_use ~~ athletic_activity
depressed ~~ athletic_activity
depressed ~~ substance_use
'
# fit the model (default estimator = ML (maximum likelihood); 
# MLR=robust maximum likelihood)
# MLR is best when observed variables are nonnormal quantitative 
# variables) 
fit <- cfa(model, data=data10, meanstructure=T, estimator="MLR")
summary(fit, fit.measures=T, standardized=T)

# extract factor scores 
fscores<-predict(fit)
head(fscores, n=20)

# draw CFA diagram
cfaplot<-semPaths(fit, intercept = FALSE, whatLabel = "Std.all", 
                  residuals = FALSE)



# reverse score negative loading variables so all loading are 
# positive (reflect higher score on construct)
data10$cesd05r<-4-data10$cesd05

# create composite variables
data10$depression<-rowMeans(subset(data10, select = c(cesd10,cesd09,cesd05r)), 
                            na.rm = FALSE)
data10$physical_activ<-rowMeans(subset(data10, select = c(excur,pacur,spcur)), 
                                na.rm = FALSE)
data10$substance_use<-rowMeans(subset(data10, select = c(dai01,sai06,dai04)), 
                               na.rm = FALSE)

describe(data10$depression)
describe(data10$physical_activ)
describe(data10$substance_use)
#describe(data10$age)

#data1$female[data1$gender==2]<-1
#data1$female[data1$gender==1]<-0
#freq(data1$female)

# standardize quantitative variables
(data10$depression <- scale(data10$depression, center = TRUE, 
                            scale = TRUE))
(data10$physical_activ <- scale(data10$physical_activ, center = TRUE, 
                                scale = TRUE))
(data10$substance_use <- scale(data10$substance_use, center = TRUE, 
                               scale = TRUE))
#(data1$age<- scale(data1$age, center = TRUE, 
#        scale = TRUE))


# define model (enclose model definition in single quotes) . 
model1 <- ' 

# structural regression model
# direct effects to outcome
             substance_use ~ s1*physical_activ + s2*depression +a2*age
            
           # direct effects to intervening/confounding variable
             depression ~ t1*physical_activ +a3*age

           # indirect effects (a*b)
           # indirect effect of physical activity and age on substance use        
             t1s2 := t1*s2
             a3s2 := a3*s2
             

           # total effects
             tot_t1s2 := s1+(t1*s2)
             tot_a3s2 := a2+(a3*s2)
             
'
fit <- cfa(model1, data=data10, meanstructure=T, estimator="MLR")
summary(fit, fit.measures=T, standardized=T)



path <- sem(model1, data = data10, estimator="MLR")
summary(path, fit.measures=T, standardized=T)

# draw CFA diagram
cfaplot<-semPaths(fit, intercept = FALSE, whatLabel = "Std.all", 
                  residuals = FALSE)

# get r-square values for endogeneous (dependent) variables
inspect(path, 'r2')


sem1 <- sem(model1, data = data10, estimator="MLR" )

summary(sem1, fit.measures=T, standardized=T)

# get r-square values for endogeneous (dependent) variables
inspect(sem1, 'r2')

# assessing local fit 

# print residual correlation matrix
resid(sem1)$cov

# print modification indices (decrease in chi-square from adding/deleting 
# a parameter)
mi1<-modindices(sem1)
mi1
# print only modification indices > 10.83 (change is significant at p<.001)
modsub<-subset(mi1, mi > 10.83)
# sort by size of modification indices (descending)
modsort <- modsub[order(-modsub$mi),] 
modsort