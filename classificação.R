# Load required libraries
library(cluster)
library(fclust)
library(factoextra)
library(tidyverse)
library(vegan)
library(entropy)
library(poorman)

dir <- "/home/mcure/Documents/Dexter/anomalia"

# Carrega os dados
completoRN1 <- read_rds("/home/mcure/Documents/Dexter/anomalia/completo.rds")

# Carrega o ciclo anual do EVI2 e adiciona uma coluna com a latitude
# Filtra os pontos ao norte da latitude 20º S
mon.df <- read_csv(paste0(dir, "/media_mensal_evi2.csv")) %>% 
  add_column(lat = completoRN1$lat) %>% 
  add_column(plot_id = completoRN1$plot_id)%>% 
  filter(lat >= -20) %>%
  dplyr::filter(plot_id != 346) %>% 
  dplyr::filter( plot_id != 3761) %>% 
  dplyr::select(-plot_id)


# Carrega o ciclo anual da anomalia do EVI2 e adiciona uma coluna com a latitude
# Filtra os pontos ao norte da latitude 20º S
mon_anom <- read_csv(paste0(dir, "/mon_anom.csv")) %>% 
  add_column(lat = completoRN1$lat)%>% 
  filter(lat >= -20) %>%
  dplyr::filter(plot_id != 346) %>% 
  dplyr::filter( plot_id != 3761) # remove quem tem EVI2 médio abaixo de zero.

# Filtra os pontos ao norte da latitude 20º S
completoRN1 <- completoRN1 %>% 
  filter(lat >= -20) %>%
  dplyr::filter(mean_ndvi > 0)

# Read data frame
#df <- read_csv(paste0(dir,"/media_mensal_evi2.csv"))

dados_marcio <- readRDS(paste0(dir,"/dados_mcure_capII.rds"))  %>% 
  filter(lat >= -20) %>%
  dplyr::filter(mean_ndvi > 0)

dados_kmeans <- data.frame(plots = dados_marcio$plot_id,
                           mean_evi = NA,
                           max_evi = NA,
                           range_evi = NA,
                           min_evi = NA,
                           sd_evi = NA,
                           entropy_evi = NA,
                           missing = NA,
                           negative = NA,
                           iqrr = NA,
                           iqrr_anom = NA)

for(i in 1:length(dados_marcio$dados1)){
  dadosR <- dados_marcio$dados1[[i]]
  dadosR$year <- sapply(as.character(dadosR[,2]), function(x) strsplit(x, "-")[[1]][1])
  dadosR <- dadosR[dadosR$year %in% 2019:2021,]
  dados_kmeans[i,2] <- mean(as.numeric(as.vector(dadosR[,4])$evi2))
  dados_kmeans[i,3] <- mean(aggregate(as.numeric(as.vector(dadosR[,4])$evi2), list(dadosR[,5]), max)$x)
  dados_kmeans[i,4] <- mean(aggregate(as.numeric(as.vector(dadosR[,4])$evi2), list(dadosR[,5]), function(x) abs(diff(range(x))))$x)
  dados_kmeans[i,5] <- mean(aggregate(as.numeric(as.vector(dadosR[,4])$evi2), list(dadosR[,5]), min)$x)
  dados_kmeans[i,6] <- mean(aggregate(as.numeric(as.vector(dadosR[,4])$evi2), list(dadosR[,5]), sd)$x)
  
  psd <- spec.pgram(as.numeric(as.vector(dadosR[,4])$evi2), plot = F)$spec
  normalized_psd <- psd / sum(psd)
  dados_kmeans[i,7] <- entropy(normalized_psd)
  
  dados_kmeans[i,8] <- sum(is.na(dadosR[,3]))/nrow(dadosR)
  dados_kmeans[i,9] <- sum(dadosR[,4]<0)/nrow(dadosR)
  
  dados_kmeans[i,10] <- mean(aggregate(as.numeric(as.vector(dadosR[,4])$evi2), list(dadosR[,5]), IQR)$x)
  
  dados_kmeans[i,11] <- IQR(mon_anom[i,-1] %>% as.numeric)
}

dados_kmeans$sd_evi <- log(dados_kmeans$sd_evi)

dados_kmeans$Leaf_flush <- completoRN1$Leaf_flush
dados_kmeans$coupling <- completoRN1$coupling

dados_kmeans1 <- dados_kmeans %>% add_column(mon.df)

dados_kmeans_scaled <- dados_kmeans %>%
  mutate(across(c(mean_evi, max_evi, range_evi, min_evi, sd_evi, entropy_evi, iqrr, iqrr_anom), scale))

# dados_std <- dados_kmeans %>%
# #  select(-1) %>% 
#   decostand(method = "standardize")



# Para determinar o número ótimo de clusters
set.seed(1234)
cluster <- fviz_nbclust(dados_kmeans_scaled[,-c(1,8,9,12,13)], FUNcluster = kmeans, method = "wss")+
  geom_vline(xintercept = 2, "dashed", color="darkgrey")

#cluster

# Realizando a clusterização por k-means

set.seed(1234)
km.res <- kmeans(dados_kmeans_scaled[,-c(1,8,9,12,13)], 2, nstart = 25)

#print(km.res) %>% 
#  broom::tidy()

dd <- cbind(dados_kmeans_scaled[,-c(1,8,9,12,13)], cluster = km.res$cluster)
#head(dd)

clusclus <- fviz_cluster(km.res,
                         dados_kmeans_scaled[,-c(1,8,9,12,13)],
                         choose.vars = c("mean_evi", "sd_evi", "max_evi", "entropy_evi", "min_evi", "range_evi", "iqrr", "iqrr_anom"),
                         ellipse.type = "norm",
                         geom = "point",
                         outlier.color = "black")+
  scale_color_manual(values=c("orange2","blue")) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = 12*PC1, yend = 12*PC2, alpha = 0.8), arrow = arrow(length = unit(0.7, "cm")), show.legend = F, color = "red") +
  geom_text(data = loadings, aes(x = 10*PC1, y = 10*PC2, label = label), color = "black", size = 5, nudge_x = 0.03, nudge_y = 0.04)+
  theme(text = element_text(size = 12))


#clusclus
# PCA

PCA <- dados_kmeans_scaled[,-c(1,8,9,12,13)] %>%
  prcomp 

#PCA %>% summary

# Eigenvalues

PC1 <- PCA %>% 
  broom::tidy() %>%
  dplyr::filter(PC==1) %>%
  dplyr::select("value")

PC2 <- PCA %>% 
  broom::tidy() %>%
  dplyr::filter(PC==2) %>%
  dplyr::select("value")

dados_pca <- data.frame(PC1 = PC1$value, PC2 = PC2$value)

# Eigenvetors
eixo1 <- PCA$rotation[,1] %>%
  broom::tidy() %>% 
  rename(variables = names, PC1 = x) %>% 
  knitr::kable()

eixo2 <- PCA$rotation[,2] %>%
  broom::tidy() %>% 
  rename(variables = names, PC2 = x) %>% 
  knitr::kable() 

# Loadings
#PCA %>% plot

# Biplot

dados_pca<- dados_pca %>% add_column(leaf = dados_marcio$Leaf_flush)
loadings <- as.data.frame(PCA$rotation)
scores <- as.data.frame(PCA$x)
label <- rownames(loadings)

# Faz um biplot
pca_plot <- ggplot()+
  geom_jitter(data = dados_pca,
              aes(x = PC1, y = PC2, color = factor(leaf), shape = factor(dd$cluster)),
              show.legend = T) +
  scale_color_manual(values = c("purple2", "green3", "orange2")) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = 10*PC1, yend = 10*PC2, alpha = 0.8),
               arrow = arrow(length = unit(0.7, "cm")), show.legend = F, color = "red") +
  geom_text(data = loadings, aes(x = 10*PC1, y = 10*PC2, label = label), color = "black", size = 4, nudge_x = 0.03, nudge_y = 0.04) +
  xlab("PC1 (39.89 %)") +
  ylab("PC2 (27.4 %)") +
  ggtitle("PCA of vegetation index metrics")+
  labs(color = "Cluster",
       shape = "Leaf flush",
       alpha = NULL)# +
#  theme(legend.position = c(0.92, 0.23))


# png("/home/mcure/Documents/Dexter/pca_metricas.png", res = 300, width = 2500, height = 1990)
#pca_plot
# dev.off()


# Random Forest

library(randomForest)
library(caret)

# Dividir dados em treino e teste (70% treino, 30% teste)

set.seed(1234)
dados_kmeans1 <- dados_kmeans %>% add_column(mon.df)
indices <- createDataPartition(as.factor(dados_kmeans1$Leaf_flush), p = 0.7, list = FALSE)
treino <- dados_kmeans1[indices, c(2:7, 10:12, 14:25)]
teste <- dados_kmeans1[-indices, c(2:7, 10:12, 14:25)]

predictors <- colnames(treino[,-9])
target <- "Leaf_flush"

# Treinar o modelo
set.seed(1234)
#?randomForest
# Instanciar modelo de Random Forest
modelo_rf <- randomForest(as.factor(Leaf_flush) ~ ., data = treino, ntree = 5000, mtry = 3)

# Mostra erro
plot(modelo_rf)
# Fazer predições com o modelo treinado
predicoes <- predict(modelo_rf, newdata = teste)

# Avalia o modelo: taxa de acerto do modelo para cada classe
cm <- confusionMatrix(predicoes, teste$Leaf_flush %>% as.factor)

# Avaliar a precisão do modelo
table(predicoes, teste$Leaf_flush)

# Recall: Revocação (Recall) ou Sensibilidade (Sensitivity): É a proporção de observações corretamente classificadas como positivas (TP) em relação ao total de observações reais que são positivas (TP + FN). Mede a habilidade do modelo em detectar corretamente observações positivas.

recall <- cm$byClass[,"Recall"]

# Precisão: representa a proporção de previsões positivas corretas em relação ao total de previsões positivas. É uma métrica útil quando o objetivo é minimizar os falsos positivos.

precision <- cm$byClass[,"Precision"]

# F1-score: é uma medida de precisão e recall combinadas, que leva em conta tanto os falsos positivos quanto os falsos negativos.

f_measure <- cm$byClass[,"F1"]

# library(pROC)
# 
# roc_curve <- multiclass.roc(teste$Leaf_flush %>% as.factor, as.numeric(predicoes), plot = TRUE, print.auc=TRUE, legacy.axes=TRUE)

