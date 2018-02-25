setwd("~/Stats")
file1="Data_example.txt"
data <-read.delim(file1, dec=",") 

names(data)
head(data)
summary(data)

data$Nom_de_colonne

rm(list=ls())

data2=t(scale(t(data[,-c(1,2)]),center=TRUE,scale=TRUE)) 

#moyenne et erreur standard
aggregate(data[,-c(1:2)],by=list(data$Genotype,data$Traitement),mean,na.rm=T)

SE<-function(data){
  sterr<- sqrt(var(data,na.rm=TRUE)/ length(data[is.na(data)!=TRUE]))
  return(sterr)
}
aggregate(data[,-c(1)],by=list(data$Traitement),SE) 

MEAN<-aggregate(data[,-c(1)],by=list(data$Traitement),mean,na.rm=T)
MEAN<-as.data.frame(MEAN)
write.table(MEAN,"MEAN_data.txt",dec=".",quote=FALSE,row.names=FALSE)

# graphiques
hist(data$Diam_moyen,c="grey",main="Histogramme Var1",ylab="Frequence",xlab="valeurs Var1")
boxplot(data) 
boxplot(data$Diam_moyen)
plot(data$Diam_moyen,data$Hauteur_moyenne) 
pairs(data) 

options("contrasts")
options(contrasts=c("contr.sum","contr.sum"))

R1 <- aov(Saccharose~Traitement, data= data)
res<-residuals(R1)
shapiro.test(res)

bartlett.test(data,residus) #si un seul facteur

library(car)
leveneTest(res~data$Genotype*data$Traitement) #si plusieurs facteurs


t.test(data$Diam_moyen,y=c(1,2)) 

data2<-data.frame(data2, data$Genotype, data$Traitement)

R1<-aov(Saccharose~Traitement, data= data) 
R2<-aov(Saccharose ~Genotype + Traitement+ Traitement%in%Genotype, data= data) 
summary(R2)

data$ID<-interaction(data$Genotype, data$Traitement)
kruskal.test(data$Diam_moyen~data$ID) 


R1 <- aov(Saccharose~ID, data= data)
require(laercio) # package spécifique
LDuncan(R1) 

TukeyHSD(R1,ordered=TRUE) # test de Tukey (toutes les comparaisons 2à 2)
require(graphics)
plot(TukeyHSD(R1)) 

library(multcomp) #pour des comparaisons spécifiques
ht <- glht(R1, linfct = mcp ("ID" = c("Geno1.Stress-Geno1.Control=0", "Geno2.Stress-Geno2.Control=0")))
confint(ht, level = 0.95)
summary(ht)

pairwise.wilcox.test(data$Diam_moyen, data$ID, p.adjust.method="bonf")

library(pgirmess)
kruskalmc(data$Diam_moyen~data$ID, data=data, probs=0.05)


cor(data[,-c(1,2)], method = "spearman") # matrice de corrélations
cor.test(data$Diam_moyen,data$Glucose,
         alternative = c("two.sided", "less", "greater"),
         method = c("pearson", "kendall", "spearman"), # test à choisir
         exact = NULL, conf.level = 0.95, continuity = FALSE) # corrélation deux à deux


rcorr(as.matrix(data2), type = "pearson")  
library(corrplot)
M<-cor(data2)
corrplot(M, type="upper", order="hclust", tl.col="black", tl.srt=45) 

library(ppcor)
pcor(data2, method = c("kendall")) 


library(nlme)
library(lme4)
library(MASS)
X<-lm(X.MS ~ Glucose + Fructose + Saccharose + Ac.citrique + Ac.malique + Ac.quinique + luteine + lycopene + beta.carotene + phytoene + VitCtot + VitCred,data=data)
Y<-stepAIC(X,direction=c("both"))
summary(Y)


data$ID<-interaction(data$Genotype,data$Traitement)
i=1
for (i in c(3:22))  # mettre le nombre de colonne du fichier
{ x<-data[,i]
R1<-aov(x~data$ID,data=data)
aov1<-anova(R1)
print(aov1)}

R1<-aov(Diam_moyen~data$ID,data=data)
res<-residuals(R1)
shapiro.test(res)
p.value<-shapiro.test(res)$p.value
Leven<-leveneTest(res~data$Genotype*data$Traitement)
p.value2<-Leven$`Pr(>F)`[1]
if (p.value>0.05 & p.value2>0.05) { LDuncan(R1)
} else { Kruskal.test(data$Diam_moyen~data$ID,data=data) }


data <-read.delim(file1, dec=".")
library(laercio)
library(pgirmess)
options(contrasts=c("contr.sum","contr.sum"))
sink(file=("Resultats_analyse.txt", append=F)
     names(data)
     i=1
     for (i in c(3:10))
     {
       x<-data[ ,i]
       print(i)
       R1 <- aov(x ~ Genotype+Traitement+Traitement%in%Genotype, data= data)
       aov1<-anova(R1)
       print(aov1)
       res<-residuals(R1)
       P<-shapiro.test(res)$p.value
       Leven<-leveneTest(res~data$Genotype*data$Traitement)
       p.value2<-Leven$`Pr(>F)`[1]
       if (p.value2>0.05 & P>0.05) { print(LDuncan(R1))
       }  else { data$ID<-interaction(data$Genotype,data$Traitement)
       KW<-kruskal.test(x~data$ID)
       p.value3<-KW$p.value
       if (p.value3<0.05) { print(kruskalmc(x~data$ID, data=NULL, probs=0.05))
       } else { print("Aucune difference stat !!!")}
       }}
     sink()
     
