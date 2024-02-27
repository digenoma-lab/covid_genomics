# Checkear los resultados de Barbara
# 15/08/23
# Pastor
# https://github.com/digenoma-lab/covid_genomics

library("ggpubr")
library('ggplot2')
library("dplyr")

mydata <- read.table('../Barbara/heterogeneity_groups.csv', h=TRUE,sep=",", row.names = 1)
summary(as.factor(mydata$label))
mydata$X <- as.factor(mydata$X)
mydata$variable <- as.factor(mydata$variable)
mydata$variant <- as.factor(mydata$variant)
mydata$location <- as.factor(mydata$location)
mydata$label <- as.factor(mydata$label)
mydata$window <- as.factor(1:30)

summary(mydata)

# Data visualization (Distribucion normal)
ggdensity(mydata$value, main = "Density plot",
          xlab = "Value")

# QQ plot (En los extremos no se ve bien)
ggqqplot(mydata$value)

# Normality test
# La hipotesis nula es que "la distribucion de la muestra es normal".
# Si el test sale significativo, significa que la distribucion NO es normal.

shapiro.test(mydata$value)
#W = 0.98803, p-value = 9.104e-15
# No tenemos distribucion normal

# Normalization
# 1. Log

mydata$logscale <- log(mydata$value)
shapiro.test(mydata$logscale) # nop

# 2. Mean = 0  

mydata$scaled <- (mydata$value - mean(mydata$value)) / sd(mydata$value)
shapiro.test(mydata$scaled) # nop
summary(mydata)

# Modelo lineal 

# Valores sin ajustar
fm4 <- 
  mydata %>% group_by(window) %>% 
  do(lm = summary(lm(value ~ label,data = .)))

fm4$lm
# label significativos 16, 20, 21, 22 and 23

# Usando datos ajustados
fm5 <- 
  mydata %>% group_by(X) %>% 
  do(lm = summary(lm(value ~ label,data = .)))

fm5$lm
# label significativos 16, 20,21,22 and 23
# Las mismas ventanas son significativas

# Adjustando p-value mediante FDR

#                      ***Atencion***
# * Aqui no estoy seguro cual deberia ser el valor de n                          ***
# * n es el numero de comparaciones (length(p)) del vector de p-values           ***   
# * Aqui supongo que el numero de comparaciones es 3, label_1, label_2 y label_3 ***
# **********************************************************************************

# venatana 16
fm4$lm[[16]] # 0.0327
p.adjust(p= 0.0327 , "fdr") #0.0981 -> no significativo

# venatana 20
fm4$lm[[20]] # 0.0253
p.adjust(p= 0.0253 , "fdr") # 0.0759 -> no significativo

# venatana 21
fm4$lm[[21]] # 0.0177
p.adjust(p= 0.0177 , "fdr",n=  3) # 0.0531 -> no significativo

# venatana 22
fm4$lm[[22]] # 0.0114
p.adjust(0.0114 , "fdr",n=  3) # 0.0342 -> significativo

# venatana 23
fm4$lm[[23]] # 0.0119, 0.0349 
p.adjust(p= 0.0119 , "fdr",n=  3) # 0.0357 -> significativo
p.adjust(p= 0.0349 , "fdr",n=  3) # 0.1047 -> no significativo

# Significativos son 22 y 23

# Plot 
ggplot(mydata, aes(x = window, y = value, col = label))+
  geom_boxplot(width=0.7, position=position_dodge(width = 0.7))

# Agregando 'location' como una covarainte
# solo ventaba 23 es significativa
# sin embargo no significa que el modelo sea mejor

fm6 <- 
  mydata %>% group_by(X) %>% 
  do(lm = summary(lm(value ~ label + location,data = .)))

fm6$lm

# 16
fm6$lm[[16]] # 0.03793, 0.02407
p.adjust(p= 0.03793 , "fdr",n= 4) # 0.15172 -> no significativo
p.adjust(p= 0.02407 , "fdr",n= 4) # 0.09628 -> no significativo

# 21
fm6$lm[[21]] # 0.018859, 0.012779 
p.adjust(p= 0.018859 , "fdr",n= 4) # 0.075436 -> no significativo
p.adjust(p= 0.012779 , "fdr",n= 4) # 0.051116 -> no significativo

# 23
fm6$lm[[23]] # 0.00770, 0.00595 
p.adjust(p= 0.00770 , "fdr",n= 4) # 0.0308 -> significativo
p.adjust(p= 0.00595 , "fdr",n= 4) # 0.0238 -> significativo

# Agregando 'variant' como una covarainte
# solo 21 es significativo, pero NO despues de calcular p-value adjusted.
# El modelo no es mejor agregando las variantes
# El Standard error aumenta y el R^2 adjusted no mejora

fm7 <- 
  mydata %>% group_by(X) %>% 
  do(lm = summary(lm(value ~ label + variant,data = .)))

fm7$lm

fm7$lm[[23]] # 0.0475
p.adjust(p= 0.0475 , "fdr",n= 4) # 0.19



