#############################################################################################################################
# Code for the manuscript:
#   "Environmental and socio-economic determinants of the occurrence of spacio-temporal clusters of malaria in Colombia"
#    
#   Gutiérrez and Tapias 
#
#   April 2022
#   
#
#   
#
#############################################################################################################################

library(RCurl)
library(dplyr)
library(rpql)
library(mvtnorm)
library(lme4)
library(glmmTMB)
library(rpql)
library(parallel)
library(scales)
library(ggeffects)

#read dataset
url_path = "https://raw.githubusercontent.com/juandavidgutier/malaria-spatio-temporal/master/data_final.csv"
malaria <- read.csv(url_path)

malaria_na <- na.omit(malaria)
malaria_first <- malaria_na[,-c(1:9)]
malaria_last <- malaria_first[,-c(116:125)]


#exclude variables with correlation > 0.8, it prevents rank errors in glmer
malaria_cor <- cor(malaria_last)
malaria_cor[upper.tri(malaria_cor)] <- 0
diag(malaria_cor) <- 0
subdata_malaria <- malaria_last[,!apply(malaria_cor,2,function(x) any(x > 0.80))]


response <- malaria_na$response
cod_munici <- malaria_na$Code_DANE
total_population <- malaria_na$population


all_data <- cbind(response, cod_munici, total_population, subdata_malaria)

str(all_data)


sub_socio_economic <- all_data[,4:length(all_data)]


## Continous responses are scaled so as to avoid large variances and easier intepretation
n_cols <- ncol(sub_socio_economic)
for (i in 1:n_cols){
  #sub_socio_economic[,i] = sub_socio_economic[,i] <- scale(sub_socio_economic[,i])
  sub_socio_economic[,i] = sub_socio_economic[,i] <- rescale(sub_socio_economic[,i], to = c(0, 1), from = range(sub_socio_economic[,i], na.rm = TRUE, finite = TRUE))
  
}


sub_socio_economic <- as.matrix(sub_socio_economic)

ini <- all_data[1:3]
socio_economic_scaled <- cbind(ini, sub_socio_economic, malaria_first[,c(116:125)])

socio_economic_scaled$Rurality <-as.factor(socio_economic_scaled$Rurality)
socio_economic_scaled$Subregion <-as.factor(socio_economic_scaled$Subregion)
socio_economic_scaled$Tipology.DNP <-as.factor(socio_economic_scaled$Tipology.DNP)
socio_economic_scaled$cod_munici <-as.factor(socio_economic_scaled$cod_munici)



########univariated glm
n_cols <- ncol(socio_economic_scaled)
for (i in 4:88){ 
  variable_name <- names(socio_economic_scaled[i])
  print(variable_name)
  
  glmmTMB_m = glmmTMB(response ~ socio_economic_scaled[,i] + (1 | cod_munici) + offset(log(total_population)), 
                      family = binomial(link = "logit"), data = socio_economic_scaled)
  print(summary(glmmTMB_m))
  
}

selected <- c("response", "cod_munici", "total_population", "forest", "grasses", "min.rainfall",
               "max.rainfall", "Population.density", "Percentage.of.ethnic.population", 
               "Percentage.of.afrocolombian.population", "Coverage.of.aqueduct", 
               "Sewer.coverage", 
               "Quantitative.housing.deficit", 
               "Average.of.prenatal.controls",  
               "Mortality.rate..x.per.1.000.inhabitants", 
               "Infant.mortality.rate.in.children.under.1.year.of.x.per.1.000.live.births.", 
               "Component.of.results.Education",  "Health.Results.Component",
               "Accumulated.annual.precipitation", "MPI..Rural", 
               "Population.in.condition.of.misery.in.the.rural.area")

socio_economic_scaled_selected <- select(socio_economic_scaled,  selected)

####
#LASSO
malaria01 <- socio_economic_scaled_selected

#set1
str(malaria01)


#fixed effects including intercept
#intercept
b <- rep.int(1, length(malaria01$response))

#fixed effects without intecept
X <- model.matrix(~ ., data = malaria01[,-(1:3)]) 
colnames(X) <- c("Intercept", paste0("x", 1:(ncol(X)-1)))

#unpenmodel
unpenfit <- glmer(response ~ X - 1 + (1|cod_munici), data = malaria01, 
                  family = "binomial")



summary(unpenfit)

isSingular(unpenfit, tol = 1e-4) #FALSE

#random effects
Z <- list(cod_munici = X[,1,drop=FALSE]) 

#id
id <- list(cod_munici = as.integer(factor(malaria01$cod_munici)))

#rpql
lambda_seq <- lseq(1e-10,1,length=1000) %>%
  cbind(., 0) ## So that random effects are not penalized
fit_sat <- build.start.fit(unpenfit, gamma = 1)


fitseq <- rpqlseq(y = malaria01$response, X = X, Z = Z, id = id, family = nb2(), 
                  lambda = lambda_seq, start = fit_sat,
                  pen.type = "lasso", offset = log(malaria01$total_population))


summary(fitseq$best.fit[[1]])


cols01 <- c(2, 4, 5,  8, 9, 18, 19)#vars of X selected by rpql

variables_inlasso <- cbind(b, socio_economic_scaled_selected[4:21])

selected_by_lasso <- dplyr::select(malaria_na, Code_DANE, population, response, forest, min.rainfall, max.rainfall, Percentage.of.afrocolombian.population,
                                   Coverage.of.aqueduct, MPI..Rural, Population.in.condition.of.misery.in.the.rural.area)
selected_by_lasso$Code_DANE <- as.factor(selected_by_lasso$Code_DANE)


#############################################################
#MULTIVARIATED MODEL WITH SELECTED VARAIBALES
#best model
fit_glmmTMB8 <- glmmTMB(response ~ forest + min.rainfall + 
                          Coverage.of.aqueduct + 
                          offset(log(population)) + 
                          (1|Code_DANE), family = binomial(link = "logit"),
                          data=selected_by_lasso)
summary(fit_glmmTMB8)

#aOR's
exp(0.041214)#aOR bosques
exp(0.001123)#aOR min.precipitacion
exp(-0.043726)#aOR Cobertura.de.acueducto..Censo
#exp(0.0720745)#aOR IPM...Rural

#ci
confint(fit_glmmTMB8, method="Wald")


#OR's of variables in best model = fit_glmmTMB8
fit_glmmTMB_for <- glmmTMB(response ~ forest +
                             offset(log(population)) + 
                             (1|Code_DANE), family = binomial(link = "logit"),
                             data=selected_by_lasso)
summary(fit_glmmTMB_for)

#OR
exp(0.08617)#univariate OR forest

#ci
confint(fit_glmmTMB_bos, method="Wald")


fit_glmmTMB_min <- glmmTMB(response ~ min.rainfall + 
                             offset(log(population)) + 
                             (1|Code_DANE), family = binomial(link = "logit"),
                             data=selected_by_lasso)
summary(fit_glmmTMB_min)

#OR
exp(0.0021370)#univariate OR min.rainfall

#ci
confint(fit_glmmTMB_min, method="Wald")


fit_glmmTMB_aque <- glmmTMB(response ~ Coverage.of.aqueduct +  
                             offset(log(population)) + 
                             (1|Code_DANE), family = binomial(link = "logit"),
                             data=selected_by_lasso)
summary(fit_glmmTMB_aque)

#OR
exp(-0.07170)#univariate OR Coverage.of.aqueduct

#ci
confint(fit_glmmTMB_aque, method="Wald")





#Figure 2
library(ggplot2)
mydf <- ggpredict(fit_glmmTMB8, terms = "forest [all]")
p <- ggplot(mydf, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
p <- p + ggtitle("A")
p <- p + labs(x = expression(paste("Forest coverage (%)")))
p <- p + labs(y = expression(paste("Probability to belongs to a cluster")))
p <- p +   theme (plot.title = element_text(size=rel(2), 
                                            vjust=2, 
                                            hjust = 0.5,
                                            face="bold"))
p <- p + theme(axis.title.y = element_text(size = rel(1.2), angle = 90))
p <- p + theme(axis.title.x = element_text(size = rel(1.2), angle = 00))
p




mydf <- ggpredict(fit_glmmTMB8, terms = "min.rainfall [all]")
p <- ggplot(mydf, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
p <- p + ggtitle("B") 
p <- p + labs(x = expression(paste("Min. multi-year average rainfall (mm)")))
p <- p + labs(y = expression(paste("Probability to belongs to a cluster")))
p <- p +   theme (plot.title = element_text(size=rel(2), 
                                            vjust=2, 
                                            hjust = 0.5,
                                            face="bold"))
p <- p + theme(axis.title.y = element_text(size = rel(1.2), angle = 90))
p <- p + theme(axis.title.x = element_text(size = rel(1.2), angle = 00))
p


mydf <- ggpredict(fit_glmmTMB8, terms = "Coverage.of.aqueduct [all]")
p <- ggplot(mydf, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)
p <- p + ggtitle("C")
p <- p + labs(x = expression(paste("Aqueduct coverage (%)")))
p <- p + labs(y = expression(paste("Probability to belongs to a cluster")))
p <- p +   theme (plot.title = element_text(size=rel(2), 
                                            vjust=2, 
                                            hjust = 0.5,
                                            face="bold"))
p <- p + theme(axis.title.y = element_text(size = rel(1.2), angle = 90))
p <- p + theme(axis.title.x = element_text(size = rel(1.2), angle = 00))
p

