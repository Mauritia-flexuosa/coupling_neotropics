library(tidyverse)

### Calcula anomalia
# Calcula o ciclo anual da anomalia no EVI2 =====

# Carrega os dados
dados <- read_rds("./dados1_annual_cycle.rds")%>%
  filter(prop_na <= .35) %>%
  filter(lat >= -20)
### Ciclo anual do EVI2
# ANNUAL CYCLE STORED FOR MEANS
#mon <- matrix(NA, nrow = dim(exp.leaf)[1], ncol = 13)
mon.df <- dados[,39:50]

# Cria uma lista vazia
lista <- vector("list", length = dim(dados)[1])

# para todos os pontos...
for(i in 1:dim(dados)[1]){

  # Seleciona cada um dos data.frames na coluna dados1
  df <- dados$dados1 %>%
    pluck(i)

  # Separa a coluna com a data em mês e ano
  df$Data <- as.Date(df$date, format = "%Y-%m-%d")
  df$Ano <- format(df$date, "%Y")
  df$mês <- format(df$date, "%m")

  # adiciona uma coluna com o plot_id, agora chamado de Ponto
  df$Ponto <- dados$plot_id[i]

  # Calcula o EVI2 mensal e mantém a coluna Ano
  df_mensal <- df %>%
    group_by(mês, Ano) %>%
    summarize(EVI2_mensal_por_ano = mean(without_missing$evi2))

  # Calcula a média anual de EVI2
  df_media_anual <- df %>%
    group_by(Ano) %>%
    summarize(EVI2_media_anual = mean(without_missing$evi2))

  # Para todos os valores no data.frame...
  for (j in 1:dim(df_mensal)[1]) {

    # Calcula a anomalia no EVI2 para cada mês
    anomalia_aux <- df_media_anual[which(df_media_anual[,1] %>%
                                           unlist %>%
                                           as.numeric == df_mensal[j,2] %>%
                                           as.numeric),2] - df_mensal[j,3]

    # Salva os valores em um vetor
    if(j == 1){
      anomalia <- anomalia_aux %>% as.numeric
    }else{
      anomalia <- c(anomalia, anomalia_aux %>% as.numeric)
    }

    # Transforma em data.frame
    anomalia %>%
      as.data.frame
  }

  # Junta os resultados no data.frame que continha os dados usados no cálculo
  df$anomalia <- anomalia

  # Calcula a anomalia media mensal e salva um  data.frame para cada ponto em uma lista
  lista[[i]] <- df %>%
    add_column(plot_id = rep(dados$plot_id[i], length(anomalia))) %>%
    group_by(mês,plot_id) %>%
    summarise(anomalia_mensal = mean(anomalia))

}

# Junta os data.frames da lista em apenas um data.frame
df_combinado <- do.call(rbind, lista) %>% group_by(mês)
df_combinado %>% head
# Insere acoluna Leaf_flush nos dados
df_anom <- left_join(df_combinado, dados[,c(1,16)], by = "plot_id")
df_anom %>% head
df_anom1 <- df_anom[,-4]
mon_anom <- spread(df_anom1, key = "mês", value = "anomalia_mensal")



####
library(entropy)

dados_kmeans <- data.frame(plots = dados$plot_id,
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

for(i in 1:length(dados$dados1)){
  dadosR <- dados$dados1[[i]]
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

dados_kmeans$Leaf_flush <- dados$Leaf_flush
dados_kmeans$coupling <- dados$coupling
dados_kmeans$prop_na <- dados$prop_na

# dados_kmeans <- dados_kmeans[-1389,]
#
# mon.df <- mon.df[-1389,]
#
# dados <- dados[-1389,]

dados_kmeans1 <- dados_kmeans %>%
  add_column(mon.df) #%>%
#  add_column(lat = dados$lat)

dados_kmeans <- dados_kmeans

dados_kmeans_scaled <- dados_kmeans %>%
  mutate(across(c(mean_evi, max_evi, range_evi, min_evi, sd_evi, entropy_evi, iqrr, iqrr_anom), scale))


### Classificação
for (i in 1:dim(dados_kmeans_scaled)[1]) {
  if(sum(!sapply(dados_kmeans_scaled[i,-c(1,8,9,12,13,14,15)], is.finite))>0){
    print(i)
  }
}


# PCA
library(vegan)
PCA <- dados_kmeans_scaled[,-c(1,8,9,12,13,14)] %>%
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

dados_pca<- dados_pca %>% add_column(leaf = dados$Leaf_flush)
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
       alpha = NULL)

# K-means
library(factoextra)
# Para determinar o número ótimo de clusters
set.seed(1234)
cluster <- fviz_nbclust(dados_kmeans_scaled[,-c(1,8,9,12,13,14)], FUNcluster = kmeans, method = "wss")+
  geom_vline(xintercept = 2, "dashed", color="darkgrey")

#cluster

# Realizando a clusterização por k-means

set.seed(1234)
km.res <- kmeans(dados_kmeans_scaled[,-c(1,8,9,12,13,14)], 2, nstart = 25)

#print(km.res) %>%
#  broom::tidy()

dd <- cbind(dados_kmeans_scaled[,-c(1,8,9,12,13,14)], cluster = km.res$cluster)
#head(dd)

clusclus <- fviz_cluster(km.res,
                         dados_kmeans_scaled[,-c(1,8,9,12,13,14)],
                         choose.vars = c("mean_evi", "sd_evi", "max_evi", "entropy_evi", "min_evi", "range_evi", "iqrr", "iqrr_anom"),
                         ellipse.type = "norm",
                         geom = "point",
                         outlier.color = "black")+
  scale_color_manual(values=c("orange2","blue")) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = 12*PC1, yend = 12*PC2, alpha = 0.8), arrow = arrow(length = unit(0.7, "cm")), show.legend = F, color = "red") +
  geom_text(data = loadings, aes(x = 10*PC1, y = 10*PC2, label = label), color = "black", size = 5, nudge_x = 0.03, nudge_y = 0.04)+
  theme(text = element_text(size = 12))


#clusclus

# Random Forest

library(randomForest)
library(caret)


# Dividir dados em treino e teste (70% treino, 30% teste)

set.seed(1234)
dados_kmeans1 <- dados_kmeans %>% add_column(mon.df)

# Realizar a subamostragem das classes majoritárias
# dados_balanceados <- downSample(x = dados_kmeans1[, c(2:7, 10:12, 15:26)], y = as.factor(dados_kmeans1$Leaf_flush),yname = "Leaf_flush", list = FALSE)
#
# # Verificar o balanceamento das classes
# treino <- dados_balanceados

# ev_80 <- dados_kmeans1 %>%
#   filter(Leaf_flush == "Evergreen") %>%
#   sample_n(103)
#
# sd_70 <- dados_kmeans1 %>%
#   filter(Leaf_flush == "Semideciduous") %>%
#   sample_n(439)
#
# d_70 <- dados_kmeans1 %>%
#   filter(Leaf_flush == "Deciduous") %>%
#   sample_n(222)

# dados_kmeans1 %>%
#   filter(Leaf_flush == "Deciduous") %>%
#   nrow * 0.7



indices <- createDataPartition(as.factor(dados_kmeans1$Leaf_flush), p = 0.7, list = FALSE)
treino <- dados_kmeans1[indices, c(2:7, 10:12, 15:26)]
teste <- dados_kmeans1[-indices, c(2:7, 10:12, 15:26)]


predictors <- colnames(treino[,-9])
target <- "Leaf_flush"

# Treinar o modelo
set.seed(1234)
#?randomForest
# Instanciar modelo de Random Forest
# modelo_rf <- randomForest(as.factor(Leaf_flush) ~ ., data = treino, ntree = 5000, mtry = 3, importance = T, proximity = T)
#
# set.seed(1234)
# modelo_rf1 <- randomForest(as.factor(Leaf_flush) ~ ., data = treino, ntree = 5000, mtry = 4, importance = T, proximity = T)


modelo_rf <- randomForest(as.factor(Leaf_flush) ~ ., data = treino, ntree = 5000, mtry = 3, importance = T, proximity = T)

# 0.65 <- modelo_rf
# 0.7 <- modelo_rf2


# Mostra erro

#plot(modelo_rf, ylim = c(0,0.9))


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

# data_fuz %>%
#   add_column(lat = dados$lat) %>%
#   add_column(lon = dados$lon) %>%
#   filter(Leaf_flush == "Deciduous") %>%
#    filter(group == 1)

# 2336
#
# 2341
# 2342
# 3254
# 2366
# 2367

# data_fuz %>%
#   add_column(lat = dados$lat) %>%
#   add_column(lon = dados$lon) %>%
#   filter(Leaf_flush == "Evergreen") %>%
# filter(plots == 3764)
#
# tres254 <- dados %>%
#   filter(plot_id == 3254) %>%
#   dplyr::select(dados1) %>%
#   unnest(cols = c(dados1))
#
#
# treze23 <- dados %>%
#   filter(plot_id == 1323) %>%
#   dplyr::select(dados1) %>%
#   unnest(cols = c(dados1))
#
# quinze41 <- dados %>%
#   filter(plot_id == 1541) %>%
#   dplyr::select(dados1) %>%
#   unnest(cols = c(dados1))
#
#
# test <- dados$dados1 %>% pluck(62) %>%
#   as.data.frame()
#
# test$mês <- format(test$date, "%m")
#
# test %>% group_by(mês) %>%
#   summarise(mean_evi = mean(evi2, na.rm = T)) %>%
# ggplot(aes(x = as.factor(mês), y = mean_evi))+
#   geom_line()+
#   geom_point(aes(y = mean_evi))+
#   ylab("EVI2")+
#   xlab("Date")+
#   ggtitle("Evergreen from Amazon", subtitle = "Classified as group 2")
