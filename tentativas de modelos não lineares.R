# Carrega o pacote minpack.lm
library(lme4)
library(lmerTest)


ggplot(data_fuz, aes(sample = abs(coupling)))+
  stat_qq()+
  stat_qq_line()+
  facet_wrap(~group)

# Define o modelo


m1 <- lm(abs(coupling) ~ poly(log(map),1), data = completoRN1)
plot(m1, col = factor(completoRN1$Leaf_flush))
summary(m1)
plot(completoRN1$map, predict(m1))
points(y = abs(completoRN1$coupling), x = completoRN1$map)

m2 <- lm(abs(coupling) ~ poly(log(map),2), data = completoRN1)
plot(m2)
summary(m2)

plot(y = abs(completoRN1$coupling), x = log(completoRN1$map))
points(log(completoRN1$map), predict(m2), col = factor(completoRN1$Leaf_flush))

m3 <- lm(abs(coupling) ~ poly(log(map),3), data = completoRN1)
plot(m3)
summary(m3)

plot(y = abs(completoRN1$coupling), x = log(completoRN1$map))
points(log(completoRN1$map), predict(m3), col = factor(completoRN1$Leaf_flush))

anova(m1, m2, m3)
AIC(m1, m2, m3)
summary(m2)
# Melhor é o modelos 2!

# Usa os parâmetros do m2 nesta função
função <- function(map, a, b, c){
  return(a*I(map**2) + b*map + c)
}


modelo <- nls(abs(coupling) ~ função(log(map), a, b, c), data = completoRN1, start = list(b = -2.64, a = -0.8074, c = 0.3936))
AIC(m1, m2, m3, modelo)

# Avalia o ajuste
summary(modelo)
plot(modelo, col = factor(completoRN1$Leaf_flush))


    fit <- nls(abs(coupling) ~ função(log(map), b, a, c), start = list(b = -2.64, a = -0.8074, c = 0.3936), data = completoRN1)

    ggplot(aes(I(map), abs(coupling)), data = broom::augment(fit)) +
      geom_line(aes(y = .fitted, color = "red", size = 3))+
      geom_point(aes(x = map, y = abs(coupling)))+
      guides(size = "none",
             colour = "none")
    



    
(I(x**2)/I(a**2)) - (I(y**2)/I(b**2))=1

plot(completoRN1$map, abs(completoRN1$coupling))
lines(completoRN1$map, predict(m2), col = "red")

lines(completoRN1$map, predict(modelo), col = "blue")
text(x = 1900, y = 0.7, "azul: nls(abs(coupling) ~ a*(map**2) + b*map + c, b = -2.64, a = -0.8074,
     c = 0.3936)", family="sans")
text(x = 1900, y = 0.76, "vermelho: lm(formula = abs(coupling) ~ poly(factor(map), 1))")



library(lme4)

# Definição da função inversa da Equação de Michaelis-Menten
inverse_michaelis_menten <- function(map, Vmax, Km) {
  return((map * Km) / (Vmax - map))
}
michaelis_menten <- function(map, Vmax, Km) {
  return(Vmax * map / (Km + map))
}
dados <- completoRN1 %>% dplyr::select(coupling | map) %>% scale


dados <- dados[,1:2] %>% data.frame%>% add_column(Leaf_flush = factor(completoRN1$Leaf_flush))

my_model <- function(map, a, b, c){
  return(a*I(map**2) + b*map + c)
}
#b = -2.64, a = -0.8074, c = 0.3936

modelo1 <- lmer(abs(coupling) ~ inverse_michaelis_menten(map, -2.64, 0.8)+(1|Leaf_flush),
               data = dados)

modelo2 <- lmer(abs(coupling) ~ michaelis_menten(map, -2.64, 0.8)+(1|Leaf_flush),
               data = dados,
               REML = F)


# Ajuste do modelo linear misto invertido
modelo3 <- lmer(abs(coupling) ~ my_model(map, -2.64, 0.8, 0.39)+(1|Leaf_flush),
               data = dados,
               REML = F)

AIC(modelo1, modelo2, modelo3)
# Sumário do modelo
summary(modelo1)
plot(modelo)
?lme4::nlmer

plot(completoRN1$map, abs(completoRN1$coupling))
lines(completoRN1$map, predict(modelo) %>% as.numeric, col = "green")
plot(completoRN1$map, predict(modelo), col = "blue", type = "l")

plot(y = inverse_michaelis_menten(completoRN1$map, -2.64, 0.8), x = completoRN1$map)

plot(y = inverse_michaelis_menten(completoRN1$map, -2.64, 0.8), x = completoRN1$coupling %>% scale %>% abs)
SSmicmen(dados$map, -2.64, 0.8)


plot(log(SSmicmen(dados$map, -2.64, 0.8)), log(abs(dados$coupling)))


plot(completoRN1$map, abs(completoRN1$coupling))

plot(inverse_michaelis_menten(dados$map, -.64, -4), col="blue", lwd = 3, add = T)
curve(SSmicmen(dados$map, -2.64, 0.8), col="blue", lwd = 3, add = T)


gmm <- nls(formula = coupling ~ -2*p*sqrt(map),
           data = completoRN1,
           start=c(p = -0.18),
           trace = T, 
           control = nls.control(maxiter = 3500, minFactor = 4.0e-5))



#nls.control(printEval = T)
plot(gmm)
summary (gmm)
p <- coef(gmm)

plot (coupling~map, xlab = "MAP", ylab = "Coupling", data=completoRN1)
#points(x, predict(gmm), col = "red")

curve(expr = expression(dados$map * p["a"])/(p["b"] - dados$map), col="blue", lwd = 3, add = T, data = dados)

plot(coupling~map, xlab = "MAP", ylab = "Coupling", data=completoRN1)
lines((-2*p["p"]*sqrt(completoRN1$map)), col="blue", lwd = 3, add = T)

plot(density(abs(dados$coupling), na.rm = T))

plot(abs(dados$coupling), predict(gmm), col="blue", lwd = 3, add = T)


plot(density(inverse_michaelis_menten(dados$map, -2.64,0.8)))





parabola <- function(map, p){
  return(sqrt(-2*p*map))
}

plot(parabola(completoRN1$map, p = -0.00003))

modelo4 <- nls(coupling ~ parabola(map,  p),
    trace = T, 
      start = list(p = -0.00003),
    data = completoRN1,
    control = nls.control(maxiter = 2000, minFactor = 4.0e-5))


AIC(m1, m2, m3,modelo, modelo1, modelo2, modelo3, modelo4,modelo5)



modelo5 <- nls(coupling ~ SSmicmen(log(map), Vm, K),
    trace = T, 
    start = list(Vm = .64, K= -0.8),
    data = completoRN1,
    control = nls.control(maxiter = 2000, minFactor = 4.0e-5))
plot(y = completoRN1$coupling, predict(modelo5))
plot(SSmicmen(log(completoRN1$map), 2.64, -.8))
summary(modelo5)



 # dispersão

x <- completoRN1$map %>% as.numeric

plot(completoRN1$coupling ~ completoRN1$map, ylim = c(-1,1));abline(h = 0, col = "red", lty = 2)

curve((0.21+x)/(31*x)-.5, add=TRUE, lwd = 2, col = "blue2")
curve(1-(0.21+x)/(21*x)-.5, add=TRUE, lwd = 2, col = "blue2")
