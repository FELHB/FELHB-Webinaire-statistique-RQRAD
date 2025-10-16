########### Code pour le webinaire du RQRAD sur les statistiques ##############
################### Par Félix L'Heureux Bilodeau ##############################
library(MASS)
library(ggplot2) #faire des graphique
library(Matrix) # créer une matrice de covariace positive
library(gridExtra) # grouper des graphiques
library(readxl) # importer les données
library(nlme) # modèle mixte
library(dplyr) # sélection des colonnes
library(lme4) # modèle mixte
library(AICcmodavg) # Valeur prédite du modèle
library(multcomp) #comparaison multiple
library(emmeans) # comparaison multiple
library(MuMIn) #IC GLMM
library(hett) # modèle robuste
library(glmmTMB) #modèle mixte généralisé

############## # Diapos 3 - refgression linéaire simple ###################

# Définir la taille de l'échantillon et le R^2 souhaité
n <- 50
r <- 0.9

# Créer la matrice de covariance avec la corrélation entre les deux variables
sigma <- matrix(c(20, r * sqrt(20 * 5),
                  r * sqrt(20 * 5), 5),
                nrow = 2)
# Générer des données suivant une distribution normale bivariée
data <- mvrnorm(n, mu = c(301, 5), Sigma = sigma)
#Transformer les données en dataframe et ajouter des noms de colonnes
df <- as.data.frame(data)
names(df) <- c("GES", "Temperature")

ggplot(df, aes(x = Temperature, y = GES)) +
  geom_point(color = "#0b3a57", size = 2) +  # Nuage de points
  geom_smooth(method = "lm", color = "black", fill = "gray", alpha = 0.3) +  # Droite de régression avec intervalle de confiance
  labs(x = "Température", y = "GES") +
  theme_minimal()


############## # Diapos 4 - regression linéaire multiple ###################

# Définir les corrélations
n <- 100  # Nombre d'observations
cor_GES_Temperature <- 0.6   # Corrélation positive avec GES
cor_GES_pH <- 0.3            # Corrélation positive avec GES
cor_GES_TEE <- -0.2          # Corrélation négative avec GES
cor_GES_N <- 0               # Aucune corrélation avec GES
# Définir la matrice de covariance en fonction des corrélations et variances souhaitées
cov_matrix <- matrix(c(
  1, cor_GES_Temperature, cor_GES_TEE, cor_GES_pH, cor_GES_N,
  cor_GES_Temperature, 1, 0.18, 0.25, 0,
  cor_GES_TEE, 0.18, 1, 0, 0,
  cor_GES_pH, 0.25, 0, 1, 0,
  cor_GES_N, 0, 0, 0, 1
), nrow = 5)

# Utiliser la méthode nearPD pour obtenir la matrice définie positive la plus proche
cov_matrix_pd <- as.matrix(nearPD(cov_matrix)$mat)

#Générer les données bivariées

data <- mvrnorm(n, mu = rep(0, 5), Sigma = cov_matrix_pd)
df <- as.data.frame(data)
names(df) <- c("GES", "Temperature", "TEE", "pH", "N")

# Redimensionner les variables pour correspondre aux plages de valeurs souhaitées
df$GES <- scales::rescale(df$GES, to = c(270, 350))
df$Temperature <- scales::rescale(df$Temperature, to = c(14, 27))
df$TEE <- scales::rescale(df$TEE, to = c(0.2, 0.5))
df$pH <- scales::rescale(df$pH, to = c(4.5, 6.5))
df$N <- scales::rescale(df$N, to = c(1.2, 2))

mod <- lm(GES ~ Temperature + TEE + pH + N, data = df)
summary(mod)

plot(df$N,df$GES)
plot(df$TEE,df$GES)

# Tracer les graphiques
g1 <-ggplot(df, aes(x = Temperature, y = GES)) +
  geom_point(color = "#0b3a57", size = 2) +  # Nuage de points
  geom_smooth(method = "lm", color = "black", fill = "gray", alpha = 0.3) +  # Droite de régression avec intervalle de confiance
  labs(x = "Température", y = "GES") +
  theme_minimal()

g2 <-ggplot(df, aes(x = TEE, y = GES)) +
  geom_point(color = "#0b3a57", size = 2) +  # Nuage de points
  geom_smooth(method = "lm", color = "black", fill = "gray", alpha = 0.3) +  # Droite de régression avec intervalle de confiance
  labs(x = "TEE", y = "GES") +
  theme_minimal()
g3 <-ggplot(df, aes(x = pH, y = GES)) +
  geom_point(color = "#0b3a57", size = 2) +  # Nuage de points
  geom_smooth(method = "lm", color = "black", fill = "gray", alpha = 0.3) +  # Droite de régression avec intervalle de confiance
  labs(x = "pH", y = "GES") +
  theme_minimal()
g4 <-ggplot(df, aes(x = N, y = GES)) +
  geom_point(color = "#0b3a57", size = 2) +  # Nuage de points
  geom_smooth(method = "lm", color = "black", fill = "gray", alpha = 0.3) +  # Droite de régression avec intervalle de confiance
  labs(x = "N", y = "GES") +
  theme_minimal()
grid.arrange(g1, g2, g3, g4, ncol = 2)


######################## # Exemple en Bloc -diapo 10 # ##########################
# Importer le jeux de donnée
# Données non disponible

Sophie$Bloc <- as.factor(Sophie$Bloc)
Sophie$Trt <- as.factor(Sophie$Trt)
Sophie$Saule <- as.factor(Sophie$Saule)
# Variables indicatrices
Sophie$Trt <- relevel(Sophie$Trt, ref = "Témoin")
Sophie$Saule <- relevel(Sophie$Saule, ref = "Sm")
Sophie <- Sophie[-c(35,36, 37), ]
Sophie <- Sophie[,-c(1,3,5,6,7,8,10,11,13,14,15,16,17) ]

mod1 <- lm(Rdt ~ Trt + Saule + Bloc, data = Sophie)
summary(mod1)

plot(mod1)

modM <- lme(Rdt ~ Trt + Saule,~1|Bloc, data = Sophie)
plot(modM)

residus<-residuals(modM,type = 'pearson')
plot(qqnorm(residus), main = 'Vérification de la normalité', ylab = 'Sample quantile', xlab = 'Theoretical Quantiles')
qqline(residus)

summary(modM)

modM2 <- lme(Rdt ~ Trt + Saule,~1|Bloc,weights=varConstPower(),
             data = Sophie)
plot(modM2)
summary(modM2)

anova(modM2)
# Obtenir les moyennes marginales estimées pour chaque niveau de Trt
emmeans_modM2 <- emmeans(modM2, ~ Trt)

# Comparaisons multiples avec ajustement de p-valeurs
comparisons <- contrast(emmeans_modM2, method = "pairwise", adjust = "tukey")

# Afficher les résultats
summary(comparisons)

######## Illustration des résultats
nv <- expand.grid(Bloc = factor(seq(from = 1, to = 4, by = 1)),
                  Trt = factor(c('Témoin', 'AB', 'AB-BRF','AB-MRF1-BRF'
                                 ,'AB-MRF1-CI','AB-MRF2-BRF')), 
                  Saule = factor(c('Pr','Sm')))

pred1 <- predictSE.lme(modM2, newdata = nv, se.fit = TRUE)

nv$Pv <- pred1$fit
nv$se <- pred1$se.fit


nv$min <- (nv$Pv - 1.96*nv$se)
nv$max <- (nv$Pv + 1.96*nv$se)

ggplot(nv, aes(x = Trt, y = Pv, fill = Saule)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = min, ymax = max), position = position_dodge(width = 0.7), width = 0.2) +
  labs(x = "Traitements",
       y = "Rendement moyen (T/ha)") +
  theme_minimal()+
  theme(
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5))+   
  scale_fill_manual(values = c("#92D050","#3C642E"))

#################### # Mesures répétées - diapo 24 # ###############################
#Données non disponible

# mettre les niveaux du traitement en facteur
donnees$Traitement<-as.factor(donnees$Traitement)
donnees$BRF<-as.factor(donnees$BRF)
donnees$Sol<-as.factor(donnees$Sol)
donnees$Plante<-as.factor(donnees$Plante)
donnees$Colonne<-as.factor(donnees$Colonne)
donnees$Semaine<-as.double(donnees$Semaine)

# suppression des colonnes abandonnées
Condition <- donnees$Colonne %in% c(1,23,33,71,72) & donnees$An ==1
donnees = subset(donnees, !Condition)

Condition2 <- donnees$Colonne %in% c(3,8,9,14,31,38,42,51,65) & donnees$An ==2
donnees = subset(donnees, !Condition2)

# suppression des colonne inutile
Flux = donnees %>%
  select(Semaine,Colonne,Sol,BRF,Plante,Tair,Tsol,TEE,Rt)

#Flux$Rtt <- sqrt(Flux$Rt)

#Modèle de base mixte
modR <- lme(Rt~Semaine+Sol+Plante+BRF+Tair+TEE,
          random =~1|Colonne,
          data=Flux,na.action=na.exclude)

plot(modR)
residusR<-residuals(modR,type = 'pearson')
plot(qqnorm(residusR), main = 'Vérification de la normalité', ylab = 'Sample quantile', xlab = 'Theoretical Quantiles')
qqline(residusR)

coef(modR, level = 1)

## grahique de l'effet de la température sur Rt
plot(x = 0, y = 0, xlim = c(24, 32),
     ylim = c(0, 400), type = "n", xlab = "Température",
     ylab = "Rt",
     main = "Effet de la température sur Rt")
abline(a= fixef(modR)["(Intercept)"],b= coef(modR, level = 1)["2",7],
       col = "black")

### Graphique de l'effet de la température sur Rt avec les pentes de chaque colonne
par(mar = c(5, 4, 4, 8))
plot(x = 0, y = 0, xlim = c(24, 32),
     ylim = c(0, 400), type = "n", xlab = "Température",
     ylab = "Rt",
     main = "Effet de la température sur Rt")
abline(a= coef(modR, level = 1)["1",1],b= coef(modR, level = 1)["2",7],
       col = "blue")
abline(a= coef(modR, level = 1)["2",1],b= coef(modR, level = 1)["2",7],
       col = "red")
abline(a= coef(modR, level = 1)["31",1],b= coef(modR, level = 1)["2",7],
       col = "green")
abline(a= coef(modR, level = 1)["72",1],b= coef(modR, level = 1)["2",7],
       col = "magenta")
abline(a= fixef(modR)["(Intercept)"],b= coef(modR, level = 1)["2",7],
       col = "black",lwd = 3)
legend(x = "topright", inset = c(-0.5, 0),
       col = c("blue", "red", "green", "magenta", "black"),
       lty = 1,
       legend = c("Colonne 1", "Colonne 2", "Colonne 3", "Colonne 4", "Moyenne"),
       xpd = TRUE)

### Modèle avec terme pour la varaince
modR <- lme(Rtt~Semaine+Sol+Plante+BRF+Tair+TEE,
            random =~1|Colonne,weights=varConstPower(),
            data=Flux,na.action=na.exclude)
plot(modR)
residusR<-residuals(modR,type = 'pearson')
plot(qqnorm(residusR), main = 'Vérification de la normalité', ylab = 'Sample quantile', xlab = 'Theoretical Quantiles')
qqline(residusR)

##### modèle avec autocorrélation ##
plot(ACF(modR), alpha=0.01, main = 'Terme d\'autocorrélation des données')

modR <- lme(Rtt~Semaine+Sol+Plante+BRF+Tair+TEE,
            random =~1|Colonne,weights=varConstPower(),
            correlation=corAR1(form = ~ Semaine|Colonne),
            data=Flux,na.action=na.exclude)

summary(modR)

########################## # Dispositif en tirroir - diapo 35 # ############################

modT <- lme(y ~ Semance + Fertilisant,
            random =~1|Bloc/Fertilisant,
            data=Ferti,na.action=na.exclude)


################# # Carré Latin - diapo36 # ###############

# Nombre de traitements et dimensions du carré latin
traitements <- c("A", "B", "C", "D", "E")
n <- length(traitements)

# Générer le carré latin de base (une matrice n x n)
carre_latin <- matrix(nrow = n, ncol = n)

# Remplissage du carré latin en décalant les traitements à chaque ligne
for (i in 1:n) {
  carre_latin[i, ] <- traitements[((1:n) + (i - 1) - 1) %% n + 1]
}

# Conversion en data frame avec noms des colonnes en français
carre_latin_df <- data.frame(
  Ligne = rep(1:n, each = n),
  Colonne = rep(1:n, times = n),
  Traitement = as.vector(carre_latin)
)

# Génération des valeurs de y_i en fonction des traitements
set.seed(42)  # Pour reproductibilité
carre_latin_df$y_i <- with(carre_latin_df, 
                           ifelse(Traitement == "A", rnorm(n, mean = 10, sd = 1),
                                  ifelse(Traitement == "B", rnorm(n, mean = 15, sd = 1),
                                         ifelse(Traitement == "C", rnorm(n, mean = 20, sd = 1),
                                                ifelse(Traitement == "D", rnorm(n, mean = 25, sd = 1),
                                                       rnorm(n, mean = 30, sd = 1))))))



modele_mixe <- lmer(y_i ~ Traitement + (1 | Ligne) + (1 | Colonne),
                    data = carre_latin_df)


################################################################################
############## Distribution généralisé #########################

######### #  Poisson  # ##########

lambdas <- c(3, 10, 15, 30)  # Par exemple, des valeurs variées pour lambda
x_vals <- 0:40  # Plage de valeurs pour x (nombre d'événements)

# Initialiser le graphique
plot(x_vals, dpois(x_vals, lambdas[1]), type = "p", pch = 16, col = "blue",
     xlab = "z", ylab = "Probabilité de z",
     main = "Poisson",
     ylim = c(0, 0.25))


# Ajouter les courbes pour chaque lambda
cols <- c("blue", "red", "green", "purple")  # Couleurs pour chaque lambda
for (i in 1:length(lambdas)) {
  points(x_vals, dpois(x_vals, lambdas[i]), pch = 16, col = cols[i])  # Utilisation de points
}

# Ajouter la légende
legend("topright", legend = paste("lambda =", lambdas),
       col = cols, lty = 1, lwd = 2, bty = "n")

################ # Binomial # ##############
n <- 20  # Nombre d'essais (constant pour tous les exemples)
probabilites <- c(0.15, 0.5, 0.7, 0.85)  # Différentes probabilités de succès

x_vals <- 0:30  # Plage de valeurs pour x (nombre de succès)

# Initialiser le graphique
plot(x_vals, dbinom(x_vals, size = n, prob = probabilites[1]), type = "p", pch = 16, col = "blue",
     xlab = "z", ylab = "Probabilitéde z",
     main = "Binomiale",
     ylim = c(0, 0.25))

# Ajouter les points pour chaque probabilité de succès
cols <- c("blue", "red", "green", "purple")  # Couleurs pour chaque probabilité
for (i in 1:length(probabilites)) {
  points(x_vals, dbinom(x_vals, size = n, prob = probabilites[i]), pch = 16, col = cols[i])
}

# Ajouter la légende
legend("topright", legend = paste("p =", probabilites),
       col = cols, pch = 16, bty = "n")


#################### #  Binomiale négative # ##########################

ks <- c(5, 15)  # Valeurs de k (nombre de succès requis)
lambdas <- c(0.3, 0.7)  # Probabilités de succès

x_vals <- 0:40  # Plage de valeurs pour x (nombre d'échecs avant d'atteindre les succès requis)

# Initialiser le graphique avec les paramètres de la première distribution
plot(x_vals, dnbinom(x_vals, size = ks[1], prob = lambdas[1]), type = "p", pch = 16, col = "blue",
     xlab = "z", ylab = "Probabilité",
     main = "Binomiale négative",
     ylim = c(0, 0.2))

# Définir les couleurs pour chaque courbe et tracer les points pour chaque combinaison de k et lambda
cols <- c("blue", "red", "green", "purple")  # Couleurs pour chaque combinaison

index <- 1  # Index pour suivre les couleurs

for (k in ks) {
  for (lambda in lambdas) {
    points(x_vals, dnbinom(x_vals, size = k, prob = lambda), pch = 16, col = cols[index])
    index <- index + 1
  }
}

# Ajouter la légende avec les combinaisons choisies de k et p
legend("topright", legend = paste("k =", rep(ks, each = length(lambdas)), ", p =", lambdas),
       col = cols, pch = 16, bty = "n")

################# # Normale # ###############################
moyennes <- c(2, 15)  # Moyennes des distributions normales
ecarts_types <- c(2, 5)  # Ecart-types des distributions normales

x_vals <- seq(-10, 40, length.out = 100)  # Plage de valeurs pour x

# Initialiser le graphique avec la première distribution normale
plot(x_vals, dnorm(x_vals, mean = moyennes[1], sd = ecarts_types[1]), type = "l", col = "blue",
     xlab = "z", ylab = "Densité",
     main = "Normale",
     ylim = c(0, 0.25))

# Définir les couleurs pour chaque courbe et tracer les lignes pour chaque combinaison de moyenne et écart-type
cols <- c("blue", "red", "green", "purple")  # Couleurs pour chaque combinaison

index <- 1  # Index pour suivre les couleurs

for (mu in moyennes) {
  for (sigma in ecarts_types) {
    lines(x_vals, dnorm(x_vals, mean = mu, sd = sigma), col = cols[index], lwd = 2)  # Tracer la courbe avec une ligne
    index <- index + 1
  }
}

# Ajouter la légende avec les combinaisons choisies de moyenne et écart-type
legend("topright", legend = paste("μ =", rep(moyennes, each = length(ecarts_types)), ", σ =", rep(ecarts_types, times = length(moyennes))),
       col = cols, lwd = 2, bty = "n")


################# # Gamma # ###################################
shapes <- c(2, 5)  # Valeurs pour shape (forme de la distribution)
rates <- c(0.5, 1.5)  # Valeurs pour rate (taux, inverse de l'échelle)

x_vals <- seq(0, 20, length.out = 100)  # Plage de valeurs pour x

# Initialiser le graphique avec la première distribution gamma
plot(x_vals, dgamma(x_vals, shape = shapes[1], rate = rates[1]), type = "l", col = "blue",
     xlab = "z", ylab = "Densité",
     main = "Gamma",
     ylim = c(0, 0.6))

# Définir les couleurs pour chaque courbe et tracer les lignes pour chaque combinaison de shape et rate
cols <- c("blue", "red", "green", "purple")  # Couleurs pour chaque combinaison

index <- 1  # Index pour suivre les couleurs

for (shape in shapes) {
  for (rate in rates) {
    lines(x_vals, dgamma(x_vals, shape = shape, rate = rate), col = cols[index], lwd = 2)  # Tracer la courbe avec une ligne
    index <- index + 1
  }
}

# Ajouter la légende avec les combinaisons choisies de shape et rate
legend("topright", legend = paste("alpha =", rep(shapes, each = length(rates)), ", beta =", rep(rates, times = length(shapes))),
       col = cols, lwd = 2, bty = "n")

######################## # exponentiel # ##############################
lambda_values <- c(0.1, 0.5, 1, 2)

# Créer un vecteur x pour l'axe des abscisses
x <- seq(0, 10, length.out = 100)

# Créer le graphique avec une mise en place initiale
plot(x, dexp(x, rate = lambda_values[1]), type = "l", col = "blue", lwd = 2, 
     xlim = c(0, 6), ylim = c(0, 2), 
     xlab = "z", ylab = "Densité",
     main = "Exponentielle")

# Ajouter les courbes pour les autres valeurs de lambda
lines(x, dexp(x, rate = lambda_values[2]), col = "red", lwd = 2)
lines(x, dexp(x, rate = lambda_values[3]), col = "green", lwd = 2)
lines(x, dexp(x, rate = lambda_values[4]), col = "purple", lwd = 2)

# Ajouter une légende
legend("topright", legend = c("lambda = 0.1", "lambda = 0.5", "lambda = 1", "lambda = 2"),
       col = c("blue", "red", "green", "purple"), lwd = 2, bty = "n")

###################### # T de Student # #############################
# Créer une séquence de valeurs pour l'axe des x
x <- seq(-5, 5, length = 100)

# Valeurs des degrés de liberté pour différentes courbes
df_values <- c(1, 3, 10, 30)
colors <- c("blue", "red", "green", "purple")
# Créer une fenêtre graphique vide
plot(x, dt(x, df = 1), type = "n", ylim = c(0, 0.4), xlab = "z", ylab = "Densité", main = "Student")

# Tracer les courbes de densité pour chaque valeur de df
for (i in seq_along(df_values)) {
  lines(x, dt(x, df = df_values[i]), col = colors[i], lty = 1)
}

# Ajouter une légende pour les différents degrés de liberté
legend("topright", legend = paste("df =", df_values), col = colors, lty = 1)


##################################################################################
################### # Modèle Gamma - diapo 42 # ######################
# Importer le jeux de donnée
# Données non disonible

Sophie$Bloc <- as.factor(Sophie$Bloc)
Sophie$Trt <- as.factor(Sophie$Trt)
Sophie$Saule <- as.factor(Sophie$Saule)
# Variables indicatrices
Sophie$Trt <- relevel(Sophie$Trt, ref = "Témoin")
Sophie$Saule <- relevel(Sophie$Saule, ref = "Sm")
Sophie <- Sophie[-c(35,36, 37), ]
Sophie <- Sophie[,-c(1,3,5,6,7,8,10,11,13,14,15,16,17) ]

modG <- glmer(Rdt ~ Trt + Saule + (1 | Bloc), 
                      data = Sophie, 
                      family = Gamma(link = "log"))

summary(modG)

###### # Test d'ajustement # ####################

# Création de la matrice de modèle (design matrix) pour les effets fixes
X.cond <- model.matrix(~ Trt + Saule, data = Sophie)

# Prédiction des effets fixes conditionnels
beta.cond <- fixef(modG)
pred.cond <- X.cond %*% beta.cond

# Prédictions sur l'échelle de la réponse
pred.count <- exp(pred.cond)  # Transformation log -> scale de la réponse
pred.count <- data.frame(pred.count)

# Simulation des paramètres des effets fixes conditionnels
set.seed(101)
pred.condpar.psim <- mvrnorm(10000, mu = beta.cond, Sigma = vcov(modG))
pred.cond.psim <- X.cond %*% t(pred.condpar.psim)
# Calcul des prédictions sur l'échelle de la réponse avec simulations postérieures
pred.ucount.psim <- exp(pred.cond.psim)

# Calcul des intervalles de confiance à 95%
ci <- t(apply(pred.ucount.psim, 1, quantile, c(0.025, 0.975)))
ci <- data.frame(ci)
colnames(ci) <- c("CI_2.5", "CI_97.5")

# Ajout des intervalles de confiance aux prédictions
pred.ucount <- data.frame(pred.count, ci)

# Simulation des Chi-carrés pour chaque observation
ChiNew <- data.frame(matrix(NA, nrow = nrow(Sophie), ncol = 10000))
for (i in 1:nrow(Sophie)) {
  for (j in 1:10000) {
    if (pred.ucount.psim[i, j] != 0) {
      ChiNew[i, j] <- ((Sophie$Rdt[i]  - pred.ucount.psim[i, j])^2) / (pred.ucount.psim[i, j])
    }
  }
}

# Calcul des Chi-carrés simulés
ChiNewSum <- colSums(ChiNew, na.rm = TRUE)

# Calcul du Chi-carré observé
Sophie$nrow <- c(1:nrow(Sophie))
pred.count$nrow <- c(1:nrow(Sophie))
CHIsq <- merge(Sophie, pred.count, by = "nrow")
CHIsq <- CHIsq %>%
  mutate(chi = ((Rdt - pred.count)^2) / pred.count)

chisq_observed <- sum(CHIsq$chi)

# Affichage de l'histogramme des Chi-carrés simulés
par(mfrow = c(1, 1))
hist(ChiNewSum, main = "Distribution du Chi-carré simulé", xlab = "Chi-carré simulé")
abline(v = chisq_observed, col = "red", lwd = 2)

#################### # Poisson - diapo 49 # ##########################################################

# Nombre d'unités expérimentales et de traitements
n_trt <- 3  # Nombre de traitements
n_ue <- 3   # Nombre d'UE par traitement
n_temps <- 4  # Nombre de temps

# Créer les facteurs pour le traitement, l'UE et le temps
Trt <- rep(paste("Trt", 1:n_trt), each = n_ue * n_temps)
UE <- rep(rep(paste("UE", 1:n_ue), each = n_temps), times = n_trt)
Temps <- rep(paste("Temps", 1:n_temps), times = n_trt * n_ue)

# Créer un facteur pour les combinaisons
data <- data.frame(Trt = Trt, UE = UE, Temps = Temps)

# Paramètres d'effet
# 1. Deux premiers traitements ont des effets différents sur le nombre d'insectes
# 2. Effet croissant du temps
set.seed(124)

# Définir une fonction pour le lambda de Poisson en fonction du traitement, du temps et de l'UE
data$lambda <- 10  # Valeur de base de lambda

# Appliquer un effet sur le lambda en fonction du traitement
data$lambda[data$Trt == "Trt1"] <- data$lambda[data$Trt == "Trt1"] + 6  # Traitement 1 a un lambda plus élevé
data$lambda[data$Trt == "Trt2"] <- data$lambda[data$Trt == "Trt2"] +2  # Traitement 2 a un lambda plus élevé


# Appliquer un effet croissant du temps (lambda augmente avec le temps)
data$lambda <- data$lambda + as.numeric(factor(data$Temps))  # Temps croissant ajoute un effet à lambda

data$UE <- rep(1:9, each = 4)
# Ajouter un effet aléatoire pour chaque UE
effet_UE <- rnorm(9, mean = 0, sd = 2)  # Effet aléatoire par UE
names(effet_UE) <- paste("UE", 1:9)  # Nommer les effets des UE
data$lambda <- data$lambda + 3*effet_UE[as.numeric(factor(data$UE))]  # Ajouter l'effet aléatoire par UE

# Générer les données de nombre d'insectes suivant une distribution de Poisson
data$y <- rpois(nrow(data), lambda = data$lambda)

data <- data %>%
  select(-lambda)
data <- data %>%
  rename(Insecte = y)

data$Trt <- as.factor(data$Trt)
data$UE <- as.factor(data$UE)
data$Temps <- as.numeric(factor(data$Temps))

modP <- glmer(Insecte ~ Trt + Temps + (1|UE), 
                        family = poisson(), 
                        data = data)

summary(modP)

#Calcul de la statistique chi-carré (somme des résidus de Pearson au carré)
chi2 <- sum(residuals(modP, type = "pearson")^2)

# Calcul de c_hat (surdispersion) en divisant par les degrés de liberté résiduels
ddl_residuels <- df.residual(modP)
c_hat <- chi2 / ddl_residuels

# Affichage de c_hat
c_hat

plot(modP)

# Comparaison des niveaux de Trt dans un modèle glm
glht_res <- glht(modP, linfct = mcp(Trt = "Tukey"))

# Résultats des comparaisons multiples
summary(glht_res)

# Calcul des moyennes marginales estimées pour Trt
emmeans_res <- emmeans(modP, ~ Trt)

# Comparaisons multiples avec ajustement de p-value
contrast(emmeans_res, method = "pairwise", adjust = "Scheffe")

####################################################

nv <- expand.grid(Temps = seq(from = 1, to = 4, by = 0.5),
                  Trt = factor(c('Trt 1', 'Trt 2', 'Trt 3')),
                  UE = factor(seq(from = 1, to = 9, by = 1)))

# Utiliser predict() pour obtenir les prédictions
pred1 <- predict(modP, newdata = nv, level = 0)  # level = 0 pour prédictions fixes

# Utiliser 'predictSE' du package 'MuMIn' pour obtenir les erreurs standard des prédictions

pred_with_se <- predictSE(modP, newdata = nv, se.fit = TRUE)

# Extraire les prédictions et leurs erreurs standard
predictions <- pred_with_se$fit
se <- pred_with_se$se.fit

nv$pred <- predictions
nv$se <- se

nv$min <- (nv$pred - 1.96*nv$se)
nv$max <- (nv$pred + 1.96*nv$se)


graph1 <- ggplot(data = nv)+
  geom_line(aes(x = Temps, y = pred, color = Trt, group = Trt),lwd=1)+
  geom_line(aes(x = Temps, y = min, color = Trt, group = Trt),lty = 6,lwd=0.65)+
  geom_line(aes(x = Temps, y = max, color = Trt, group = Trt),lty = 6,lwd=0.65)+
  labs(y= ('insecte'))+
  theme_minimal()+labs(caption = "Barre d'erreur : IC 95%")

graph1

##################### # Robuste (student) - diapo 45 # ##############################################
#https://pmarchand1.github.io/ECL8202/notes_cours/04-Regression_robuste.html#r%C3%A9gression_(t)
# Définir les paramètres
set.seed(123)
n <- 150  # Nombre d'observations
df <- 1.3   # Degrés de liberté pour la distribution de Student

# Générer la variable explicative continue TEE
TEE <- runif(n, min = 10, max = 50)  # TEE varie de 0 à 10

# Générer la variable de traitement avec deux niveaux
Traitement <- rep(c("A", "B"), each = n / 2)

# Générer la variable réponse y avec une distribution de Student
# On introduit un effet de TEE et de Traitement dans y
effet_TEE <- 0.2 * TEE
effet_Traitement <- ifelse(Traitement == "A", 1, -1)
y <- effet_TEE + effet_Traitement + rt(n, df = df)  # Ajout d'une erreur de Student

# Créer le data frame final
data <- data.frame(y = y, TEE = TEE, Traitement = Traitement)
data$Traitement <- as.factor(data$Traitement)

# modèle linéaire
mod <- lm(y ~ TEE + Traitement, data =data)

summary(mod)


plot(mod)

plot(data$TEE,data$y)
boxplot(data$Traitement,data$y, xlab = "Traitement", 
        ylab = "Valeur de y")

## modèle robuste avec distribution de student
modS <- tlm(y ~ TEE + Traitement, data = data,
            estDof = TRUE)
summary(modS)


plot(modS[["loc.fit"]][["fitted.values"]], modS[["loc.fit"]][["residuals"]],
     xlab = "Valeurs ajustées", ylab = "Résidus",
     main = "Graphique des résidus en fonction des valeurs ajustées")


plot(data$TEE, modS[["loc.fit"]][["residuals"]], 
     xlab = "TEE", ylab = "Résidus", 
     main = "Graphique des résidus en fonction de TEE")

hist(modS[["loc.fit"]][["residuals"]], breaks = 20, main = "Histogramme des résidus",
    xlab = "Résidus", col = "lightblue")

############################# # Binomiale - diapo 52 # ############################################

# Importer le jeux de donnée
Sophie<- read_excel("C:/Users/felhb/OneDrive - Université Laval/Doctorat/1. Ste-Sophie/Resultats_Saule/Saules_2023.xlsx")

Sophie$Trt <- as.factor(Sophie$Trt)
Sophie$Saule <- as.factor(Sophie$Saule)
Sophie$Bloc <- as.factor(Sophie$Bloc)
Sophie$Trt <- relevel(Sophie$Trt, ref = "Témoin")
Sophie$Saule <- relevel(Sophie$Saule, ref = "Sm")

Sophie <- Sophie[,-c(1,3,5,6,7,8,10,11,12,13,14,15) ]

#modèle
modB <- glm(Survie ~ Trt+Bloc+Saule,
            family = binomial(link = logit),
            weights = Nb, data = Sophie)

summary(modB)

plot(rstudent(modB) ~ fitted(modB), ylab = "Résidus de Student",
     xlab = "Valeurs prédites", main = "Résidus vs valeurs prédites")

plot(cooks.distance(modB)) #Valeurs extrêmes

#surdispersion
c_hat(modB) #library(AICcmodavg)

## binomiale négative https://pmarchand1.github.io/ECL8202/notes_cours/06-Modeles_generalises_mixtes2.html

#réarangement des données
Sophie$Vivant <- abs(Sophie$Survie*Sophie$Nb)
Sophie <- Sophie[-2,]

#modèle 
library(MASS)
modBN <- glm.nb(Vivant ~ Trt+Bloc+Saule, data = Sophie,
                control = glm.control(maxit = 5000))

summary(modBN, dispersion = 1.29)
cor(Sophie$Vivant,fitted(modBN))^2
plot(Sophie$Vivant,fitted(modBN))

# Vérification de la surdispersion
chi2 <- sum(residuals(modBN, type = "pearson")^2)
chi2 / df.residual(modBN)

#valeurs prédite du modèle
sim_nb <- simulate(modBN, nsim = 1000, re.form = NULL, newdata = Sophie)
sim_pred <- mutate(Sophie, pred = predict(modBN, type = "response"),
                    q025 = apply(sim_nb, 1, quantile, probs = 0.025),
                    q975 = apply(sim_nb, 1, quantile, probs = 0.975)) %>%
  arrange(pred)

#graphique
ggplot(sim_pred, aes(x = 1:nrow(sim_pred), y = pred, ymin = q025, ymax = q975)) +
  geom_ribbon(alpha = 0.3) +
  geom_line() +
  geom_point(aes(y = Vivant))+
  xlab(" # observation") + ylab("nb de survie")


########################### #Beta - diapo 46# #############################
set.seed(123)

# Paramètres
n_fermes <- 5    # Nombre de fermes
n_obs_per_ferme <- 30 # Nombre d'observations par ferme
n_total <- n_fermes * n_obs_per_ferme

# Facteurs
Ferme <- factor(rep(1:n_fermes, each = n_obs_per_ferme))
Travail <- factor(rep(c("Labour", "Non-Labour"), times = n_total / 2))

# Effet aléatoire des fermes
effet_ferme <- rnorm(n_fermes, mean = 0, sd = 0.2)
effet_ferme_rep <- rep(effet_ferme, each = n_obs_per_ferme)

# Effets fixes
effet_travail <- c(Labour = -0.5, `Non-Labour` = 0.5)
effet_fixe_travail <- effet_travail[Travail]

# Génération des valeurs moyennes (mu) dans l'intervalle (0, 1)
mu <- plogis(0.5 + effet_fixe_travail + effet_ferme_rep) # Lien logit

# Génération des données bêta
shape1 <- mu * 10
shape2 <- (1 - mu) * 10
beta_values <- rbeta(n_total, shape1 = shape1, shape2 = shape2)

# Transformation pour l'intervalle [0.01, 0.12]
Response <- 0.01 + (0.12 - 0.01) * beta_values

# Construction du jeu de données
data <- data.frame(
  Ferme = Ferme,
  Travail = Travail,
  MO = Response
)


# Visualisation
boxplot(MO ~ Travail, data = data, col = c("lightblue", "lightgreen"),
        main = "Pourcentage de MO par type de travail",
        ylab = "% Matière Organique", xlab = "Type de Travail")

# Modèle linéaire
modNorm <- lme(MO ~ Travail, ~1|Ferme, data = data)
summary(modNorm)

plot(modNorm)
residusN<-residuals(modNorm,type = 'pearson')
plot(qqnorm(residusN), main = 'Vérification de la normalité', ylab = 'Sample quantile', xlab = 'Theoretical Quantiles')
qqline(residusN)

# Modèle mixte bêta
modB <- glmmTMB(
  MO ~ Travail + (1 | Ferme), 
  data = data,
  family = beta_family(link = "logit")
)

# Résidus
residuals <- residuals(modB, type = "pearson")
plot(residuals ~ fitted(modB), main = "Résidus vs Valeurs ajustées")

summary(modB)

# Récupérer tous les objets dans l'environnement global
objects <- ls()

# Conserver seulement 'data' et 'modB'
objects_to_remove <- objects[!(objects %in% c("data", "modB"))]

# Supprimer tous les autres objets
rm(list = objects_to_remove)


############# Test d'ajustement #################
# Création de la matrice de modèle (design matrix) pour les effets fixes
X.cond <- model.matrix(~ Travail, data = data)

# Prédiction des effets fixes conditionnels
beta.c <- fixef(modB)
beta.cond <- beta.c[["cond"]]
pred.cond <- X.cond %*% beta.cond

# Prédictions sur l'échelle de la réponse
pred.mu <- plogis(pred.cond)  # Transformation logit -> échelle de la réponse

# Simulation des paramètres des effets fixes conditionnels
set.seed(101)
pred.condpar.psim <- mvrnorm(10000, mu = beta.cond, Sigma = Sigma <- vcov(modB)$cond)
pred.cond.psim <- X.cond %*% t(pred.condpar.psim)

# Calcul des prédictions sur l'échelle de la réponse avec simulations postérieures
pred.ucount.psim <- plogis(pred.cond.psim)

# Calcul des variances spécifiques à la distribution bêta
phi <- 310  # Dispersion (phi) du modèle bêta
variances <- pred.mu * (1 - pred.mu) / (1 + phi)

# Simulation des Chi-carrés pour chaque observation
ChiNew <- data.frame(matrix(NA, nrow = nrow(data), ncol = 100))
for (i in 1:nrow(data)) {
  for (j in 1:100) {
    if (!is.na(pred.ucount.psim[i, j])) {
      ChiNew[i, j] <- ((data$MO[i] - pred.ucount.psim[i, j])^2) / variances[i]
    }
  }
}

# Calcul des Chi-carrés simulés
ChiNewSum <- colSums(ChiNew, na.rm = TRUE)

# Calcul du Chi-carré observé
data$nrow <- c(1:nrow(data))
pred.mu <- data.frame(pred.mu, nrow = 1:nrow(data))
CHIsq <- merge(data, pred.mu, by = "nrow")
CHIsq <- CHIsq %>%
  mutate(chi = ((MO - pred.mu)^2) / variances)

chisq_observed <- sum(CHIsq$chi)

# Affichage de l'histogramme des Chi-carrés simulés
par(mfrow = c(1, 1))
hist(ChiNewSum, main = "Distribution du Chi-carré simulé", xlab = "Chi-carré simulé")
abline(v = chisq_observed, col = "red", lwd = 2)

# Résultat du test
p_value <- mean(ChiNewSum >= chisq_observed)
cat("Chi-carré observé :", chisq_observed, "\n")
cat("p-value :", p_value, "\n")

