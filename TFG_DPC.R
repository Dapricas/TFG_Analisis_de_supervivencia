# Carga de paquetes y librerías
#------------------------------
#
# Se recomienda sólo ejecutar aquellos paquetes que no se tengan instalados
#en el equipo, motivo por el cual se encuentran comentados
#
# install.packages(c("survival", "survminer"))
# 
# install.packages("dplyr")
# 
# install.packages("writexl")
# 
# install.packages("volcano3D")
#
# install.packages("pec")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("qvalue")

require(survival)
require(xtable)
library(pec)
library(lava)
library(qvalue)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(plotly)
library(gridExtra)
library(dplyr)
library(volcano3D)
library(writexl)
library(readxl)

# Extracción de los datos
#------------------------

# Dado que los datos vienen de manera "traspuesta" a la forma en que el 
#software los requiere, además de que el tipo de algunos campos es 
#modificado en el proceso de carga, se deben realizar una extracción y 
#transformación previa

Datos_limpios <- read.csv("C:/Users/dapri/OneDrive/Escritorio/Mates_uni/Cuarto/TFG/Datos_limpios_bien.csv", header=FALSE, sep=";", na.strings=" N/A")
matriz_t = t(as.matrix(Datos_limpios))
data_full <- as.data.frame(matriz_t[-1,])
colnames(data_full) = matriz_t[1,]

write.csv(data_full,"C:/Users/dapri/OneDrive/Escritorio/Mates_uni/Cuarto/TFG/data_full.csv")
data_full <- read.csv("C:/Users/dapri/OneDrive/Escritorio/Mates_uni/Cuarto/TFG/data_full.csv", header = TRUE)

data = data.frame(data_full)

# Debemos asegurarnos de que el dataframe tiene la estructura deseada

data$X <- NULL

str(data)

# Observamos que los campos de los factores genéticos son de tipo factor, de modo que tendremos que transformarlos a num

len = length(data)
for (index in c(33:len)) {
  data[,index] = as.numeric(gsub(',','.',data[,index]))
}

# Se actualizan los valores "MX" de la variable tnm.m por NA's

data$tnm.m = as.factor(ifelse(data$tnm.m == ' MX', NA, data$tnm.m)-1)

data$tnm.stage = as.factor(data$tnm.stage)

# data_save = data.frame(data)

# Preprocesamiento de los datos
#------------------------------

# Visualización previa

str(data)

# Seleccionamos los campos a estudiar

data = data[,c(1,2,6,10:13,33:54707)]

# Resumen de los campos especiales

summary(data[c(1,2,3)])
summary(data[c(11,13)])
summary(as.factor(data$rfs.event))
summary(as.factor(data$os.event))

# Eliminamos los datos con rfs.delay igual a 0

data = filter(data, ifelse(is.na(rfs.delay),1,rfs.delay) != 0)

# Generaión de datasets
#----------------------

# data_lde
#---------

# Creación del dataset

data_lde = data.frame(data)

# Comprobamos el número de NA's de tiempos de vida hasta la recaída, eliminamos los datos que no tengan ese dato por ser pocos

sum(is.na(data_lde$rfs.delay), na.rm = TRUE)

if (sum(is.na(data$rfs.delay), na.rm = TRUE) > 0){
  data_lde = filter(data_lde, is.na(rfs.delay) == FALSE)
}

# Comprobamos el número de NA's de recaídas conocidas, si hay NA's los remplazamos por 0 porque los consideramos datos censurados

sum(is.na(data_lde$rfs.event), na.rm = TRUE)

if (sum(is.na(data_lde$rfs.event), na.rm = TRUE) > 0){
  data_lde$rfs.event = ifelse(is.na(data_lde$rfs.event),0,data_lde$rfs.event)
}

# data_tot
#---------

data_tot = data.frame(data)

# Comprobamos el número de NA's de tiempos de vida hasta la muerte, eliminamos los datos que no tengan ese dato por ser pocos

sum(is.na(data_tot$os.delay), na.rm = TRUE)

if (sum(is.na(data_tot$os.delay), na.rm = TRUE) > 0){
  data_tot = filter(data_tot, is.na(os.delay) == FALSE)
}

# Comprobamos el número de NA's de tiempos de vida hasta la recaída, eliminamos los datos que no tengan ese dato por ser pocos

sum(is.na(data_tot$rfs.delay), na.rm = TRUE)

if (sum(is.na(data_tot$rfs.delay), na.rm = TRUE) > 0){
  data_tot = filter(data_tot, is.na(rfs.delay) == FALSE)
}

# Comprobamos el número de NA's de muertes, si hay NA's los remplazamos por 0 porque los consideramos datos censurados

sum(is.na(data_tot$os.event), na.rm = TRUE)

if (sum(is.na(data_tot$os.event), na.rm = TRUE) > 0){
  data_tot$os.event = ifelse(is.na(data_tot$os.event),0,data_tot$os.event)
}


########################################################################################################################
#                                                                                                                      #
#                                             SUPERVIVENCIA HASTA RECIDIVA                                             #
#                                                                                                                      #
########################################################################################################################

# Regresiones de Cox independientes
#----------------------------------

# Creamos el elemento Surv

S_lde = Surv(data_lde$rfs.delay, data_lde$rfs.event)

# Generamos las regresiones de Cox independientes para cada variable y reunimos la información en el dataset hr_lde

nombres = names(data_lde)[-c(4:8)]
len = length(nombres)
hr = numeric(len)
ic_min = numeric(len)
ic_max = numeric(len)
pvalor = numeric(len)
tam = numeric(len)
resid = numeric(len)
i = 1

for (name in nombres){
  form = as.formula(paste("S_lde","~",name))
  fit_cox = coxph(form, data_lde, ties = 'efron')
  p_schoenfeld = cox.zph(fit_cox)
  coefs = summary(fit_cox)$coefficients
  ic = summary(fit_cox)$conf.int
  hr[i] = coefs[2]
  ic_min[i] = ic[3]
  ic_max[i] = ic[4]
  pvalor[i] = coefs[5]
  resid[i] = p_schoenfeld$table[1,3]
  tam[i] = summary(fit_cox)$n
  i = i + 1
}

hr_lde = data.frame(nombres,hr,ic_min,ic_max,pvalor,resid,tam)

hr_lde = hr_lde[order(hr_lde$pvalor),]

# write_xlsx(hr_lde,'C:/Users/dapri/OneDrive/Escritorio/Mates_uni/Cuarto/TFG/hr_lde.xlsx')
# 
# hr_lde = read_excel('C:/Users/dapri/OneDrive/Escritorio/Mates_uni/Cuarto/TFG/hr_lde.xlsx')

# Obtención del q-valor
#----------------------

qvalor = qvalue(hr_lde$pvalor, fdr.level = 0.05)

hr_lde$qvalor <- qvalor$qvalues

hr_lde$fdr <- qvalor$significant

# Representación de las distribuciones del p-valor y el q-valor

p1 <- ggplot(data=data.frame(hr_lde$pvalor), mapping=aes(hr_lde$pvalor)) +
  geom_histogram(fill = "cyan3") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_bw()

p2 <- ggplot(data=data.frame(hr_lde$qvalor), mapping=aes(hr_lde$qvalor)) +
  geom_histogram(fill = "cyan3") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_bw()

grid.arrange(p1, p2)

# Volcano plot
#-------------

hr_lde$logFC = log2(hr_lde$hr)
hr_lde$expression = ifelse(hr_lde$fdr == TRUE & abs(hr_lde$logFC) >= 0.5, 
                           ifelse(hr_lde$logFC> 0.5 ,'Up','Down'),
                           'Stable')

volcano_lde <- ggplot(data = hr_lde, 
                      aes(x = logFC, 
                          y = -log10(qvalor), 
                          colour=expression)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-4.5, 4.5)) +
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  labs(x="log2(hazard ratio)",
       y="-log10 (q-valor)",
       title="Volcano plot")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
volcano_lde

# Primer filtro de seleccion de variables
#----------------------------------------

hr_lde_fdr = filter(hr_lde, fdr == TRUE)
dim(hr_lde_fdr)
hr_lde_new = filter(hr_lde_fdr, resid < 0.05)
dim(hr_lde_new)

pos = match(hr_lde_new$nombres,hr_lde_fdr$nombres)
resultados_lde = hr_lde_fdr[pos,1:8]
resultados_lde
xtable(resultados_lde, digits=c(0,0,5,3,3,7,7,0,7))

# Aplicación de un método de bakward
#-----------------------------------

data_lde_bw = data_lde %>% select(hr_lde_new$nombres)
data_lde_bw$rfs.delay = data_lde$rfs.delay
data_lde_bw$rfs.event = data_lde$rfs.event

form = as.formula(paste("Surv(rfs.delay,rfs.event)~",paste(names(data_lde_bw)[1:(length(names(data_lde_bw))-2)], collapse= "+")))
bw_lde = selectCox(form, data_lde_bw)

bw_lde$In

# Representacion por curvas K-M
#------------------------------

data_lde_km = data_lde %>% select(c("rfs.delay","rfs.event"))
m.X231606_at = ifelse(data_lde$X231606_at < median(data_lde$X231606_at),1,0)
m.X1558369_at = ifelse(data_lde$X1558369_at < median(data_lde$X1558369_at),1,0)
m.X206411_s_at = ifelse(data_lde$X206411_s_at < median(data_lde$X206411_s_at),0,1)
m.X205327_s_at = ifelse(data_lde$X205327_s_at < median(data_lde$X205327_s_at),1,0)
suma_lde = m.X231606_at + m.X1558369_at + m.X206411_s_at + m.X205327_s_at
data_lde_km$conjunto = floor(suma_lde/2)

curva.km = survfit(S_lde ~ conjunto, data_lde_km, type = "kaplan-meier")

plot(curva.km, xmax=100, ylim = c(0,1), mark.time=FALSE, col=c(3,'darkgoldenrod1',2), lwd=3, bty='n', las=1,
     xlab='Tiempo de estudio (meses)', ylab='Probabilidad (Metástasis)')

survdiff(S_lde ~ conjunto, data_lde_km)





########################################################################################################################
#                                                                                                                      #
#                                                  SUPERVIVENCIA TOTAL                                                 #
#                                                                                                                      #
########################################################################################################################

# Regresiones de Cox independientes
#----------------------------------

# Creamos el elemento Surv

S_tot = Surv(data_tot$os.delay, data_tot$os.event)

# Generamos las regresiones de Cox independientes para cada variable y reunimos la información en el dataset hr_tot

nombres = names(data_tot)[-c(4:8)]
len = length(nombres)
hr = numeric(len)
ic_min = numeric(len)
ic_max = numeric(len)
pvalor = numeric(len)
tam = numeric(len)
resid = numeric(len)
i = 1

for (name in nombres){
  form = as.formula(paste("S_tot","~",name))
  fit_cox = coxph(form, data_tot, ties = 'efron')
  p_schoenfeld = cox.zph(fit_cox)
  coefs = summary(fit_cox)$coefficients
  ic = summary(fit_cox)$conf.int
  hr[i] = coefs[2]
  ic_min[i] = ic[3]
  ic_max[i] = ic[4]
  pvalor[i] = coefs[5]
  resid[i] = p_schoenfeld$table[1,3]
  tam[i] = summary(fit_cox)$n
  i = i + 1
}

hr_tot = data.frame(nombres,hr,ic_min,ic_max,pvalor,resid,tam)

hr_tot = hr_tot[order(hr_tot$pvalor),]

# write_xlsx(hr_tot,'C:/Users/dapri/OneDrive/Escritorio/Mates_uni/Cuarto/TFG/hr_tot.xlsx')
# 
# hr_tot = read_excel('C:/Users/dapri/OneDrive/Escritorio/Mates_uni/Cuarto/TFG/hr_tot.xlsx')

# Obtención del q-valor
#----------------------

qvalor = qvalue(hr_tot$pvalor, fdr.level = 0.05)

hr_tot$qvalor <- qvalor$qvalues

hr_tot$fdr <- qvalor$significant

# Representación de las distribuciones del p-valor y el q-valor

p1 <- ggplot(data=data.frame(hr_tot$pvalor), mapping=aes(hr_tot$pvalor)) +
  geom_histogram(fill = "cyan3") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_bw()

p2 <- ggplot(data=data.frame(hr_tot$qvalor), mapping=aes(hr_tot$qvalor)) +
  geom_histogram(fill = "cyan3") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_bw()

grid.arrange(p1, p2)

# Volcano plot
#-------------

hr_tot$logFC = log2(hr_tot$hr)
hr_tot$expression = ifelse(hr_tot$fdr == TRUE & abs(hr_tot$logFC) >= 0.5, 
                           ifelse(hr_tot$logFC > 0.5 ,'Up','Down'),
                           'Stable')

volcano_tot <- ggplot(data = hr_tot, 
                      aes(x = logFC, 
                          y = -log10(qvalor), 
                          colour=expression)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("blue","grey","red"))+
  xlim(c(-4.5, 4.5)) +
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  labs(x="log2 (hazard ratio)",
       y="-log10 (q-valor)",
       title="Volcano plot")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
volcano_tot

# Primer filtro de seleccion de variables
#----------------------------------------

hr_tot_fdr = filter(hr_tot, fdr == TRUE)
dim(hr_tot_fdr)
hr_tot_new = filter(hr_tot_fdr, resid < 0.05)
dim(hr_tot_new)

pos = match(hr_tot_new$nombres,hr_tot_fdr$nombres)
resultados_tot = hr_tot_fdr[pos,1:8]
resultados_tot
xtable(resultados_tot, digits=c(0,0,5,3,3,7,7,0,7))

# Aplicación de un método de bakward
#-----------------------------------

data_tot_bw = data_tot %>% select(hr_tot_new$nombres)
data_tot_bw$os.delay = data_tot$os.delay
data_tot_bw$os.event = data_tot$os.event

form = as.formula(paste("Surv(os.delay,os.event)~",paste(names(data_tot_bw)[1:(length(names(data_tot_bw))-2)], collapse= "+")))
bw_tot = selectCox(form, data_tot_bw)

bw_tot$In

# Representacion por curvas K-M
#------------------------------

data_tot_km = data_tot %>% select(c("os.delay","os.event"))
m.X204798_at = ifelse(data_tot$X204798_at < median(data_tot$X204798_at),1,0)
m.X1566163_at = ifelse(data_tot$X1566163_at < median(data_tot$X1566163_at),0,1)
m.X226652_at = ifelse(data_tot$X226652_at < median(data_tot$X226652_at),1,0)
m.X218350_s_at = ifelse(data_tot$X218350_s_at < median(data_tot$X218350_s_at),1,0)
suma_tot = m.X204798_at + m.X1566163_at + m.X226652_at + m.X218350_s_at
data_tot_km$conjunto = floor(suma_tot/2)

curva.km = survfit(S_tot ~ conjunto, data_tot_km, type = "kaplan-meier")

plot(curva.km, xmax=100, ylim = c(0,1), mark.time=FALSE, col=c(3,'darkgoldenrod1',2), lwd=3, bty='n', las=1,
     xlab='Tiempo de estudio (meses)', ylab='Probabilidad (Metástasis)')

survdiff(S_tot ~ conjunto, data_tot_km)

