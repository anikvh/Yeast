#INSTALAR Y LLAMAR LAS SIGUIENTES LIBRERIAS
library(VIM)
library(missForest)
library(arules)
library(hydroGOF)
library(class)

#FUNCION PARA LEER LOS DATOS
Readfile<-function(datos){
  path<-"C:/Users/Valeria/Documents/Escuela/2018-2/MIND(M)/PROYECTO FINAL"
  setwd(path)
  data<-read.table(file="YeastOriginal.data",sep="")
  return(data)
}

#FUNCION PARA INTRODUCIR MV
DataMV<-function(datos,porc){
  #Adicionar el porcentaje de MV de manera aleatoria
  return(prodNA(datos,noNA=porc))
}

#FUNCION PARA ANALIZAR MV
InfoMV<-function(datos){
  #Encontrar las columnas con MV
  cmv=which(colSums(is.na(datos))!=0)
  cat("\n\nCant de columnas con MV:",cmv)
  #Detectar las filas que tienen MV
  rmv=which(rowSums(is.na(datos))!=0,arr.ind=T)
  cat("\n\nCant de filas con MV:",rmv)
  #Para hallar el porcentaje de filas con MV
  pfmv=length(rmv)*100/dim(datos)[1]
  cat("\n\nPorcentaje de filas con Mv:",pfmv)
}

#FUNCION PARA GRAFICAR EL ATASET CON MV
PlotValues<-function(datos, ini_var, fin_var){
  mice_plot<-aggr(datos[ini_var:fin_var], col=c("gray77","lightskyblue2"),numbers=TRUE, sortVars=TRUE, labels=names(datos[ini_var:fin_var]), gap=3, ylab=c("Missing Values","Pattern"))
}

#FUNCION PARA IMPUTAR POR LA MODA
Mode<-function(datos){
  ux<-unique(datos)
  ux[which.max(tabulate(match(datos,ux)))]
}

###################################################################

data<-Readfile("YeastOriginal.data") #Traer los datos
InfoMV(data) #Verificar que no hay mv 
summary(data)

names(data)#Ver los nombres de las variables
boxplot(data$erl)
boxplot(data$pox)
yeast_sin<-subset(data, select=-c(Sequence_Name,erl,pox))#Quitar las variables que no cambian V6 y V7

write.table(yeast_sin, file="C:/Users/Valeria/Documents/Escuela/2018-2/MIND(M)/PROYECTO FINAL/Yeast_sin.data")#Exportar datos para comparar en knime

#iMPUTACION DE MV

data.MV=DataMV(yeast_sin,0.1) #Introducir 10% de mv
InfoMV(data.MV) #Verificar que ahora hay mv

names(yeast_sin)#Mostrar los nombres de las variables
data.MV<-subset(data.MV, select=-c(Localization_site))#Quitar la variable de la clase

PlotValues(data.MV,1,6)#Graficar MV
summary(data.MV)#summary antes de imputar

data.ImptMode<-data.MV#Crear una copia del dataset con Mv para la moda
data.ImptMode[is.na(data.ImptMode)]<-Mode(data.ImptMode[!is.na(data.ImptMode)])#Imputar la moda
summary(data.ImptMode)#Summary depues de imputar la moda

data.ImptMean<-data.MV#Crear una copia del dataset con Mv para la media
data.ImptMean[is.na(data.ImptMean)]<-mean(data.ImptMean[!is.na(data.ImptMean)])#imputar media
summary(data.ImptMean)#Summary después de imputar la media

data.ImptMedian<-data.MV#Crear una copia del dataset con Mv para la mediana
data.ImptMedian[is.na(data.ImptMedian)]<-median(data.ImptMedian[!is.na(data.ImptMedian)])#Imputar media
summary(data.ImptMedian)#summary despues de imputar la mediana

#Imputar con varios valores de k para ver cual tiene el menor error
data.ImptKnn7<-kNN(data.MV,k=7)#Imputar por knn con k=7
summary(data.ImptKnn7)#Summary para verificar la imputacion con k=7
data.ImptKnn7<-subset(data.ImptKnn7, select=-c(mcg_imp:nuc_imp))#Datos con knn sin las columnas de True/ False para poder calcular el MSE

data.ImptKnn8<-kNN(data.MV,k=8)#Imputar por knn con k=8
summary(data.ImptKnn8)#Summary para verificar la imputacion con k=8
data.ImptKnn8<-subset(data.ImptKnn8, select=-c(mcg_imp:nuc_imp))#Datos con knn sin las columnas de True/ False para poder calcular el MSE

data.ImptKnn9<-kNN(data.MV,k=9)#Imputar por knn con k=9
summary(data.ImptKnn9)#Summary para verificar la imputacion con knn
data.ImptKnn9<-subset(data.ImptKnn9, select=-c(mcg_imp:nuc_imp))#Los datos deben tener la misma dimension

datamse<-data#Agregar los datos originales a una nueva variable para calcular el mse
datamse<-subset(datamse,select=-c(Sequence_Name,erl,pox,Localization_site))#Quitar variables categoricas para que todos los datos tengan la misma dimension

Errork7<-mse(data.ImptKnn7,datamse)#MSE k=7
Errork8<-mse(data.ImptKnn8,datamse)#MSE k=8
Errork9<-mse(data.ImptKnn9,datamse)#MSE k=9

table(Errork7<Errork8)#Verificar si el errror de k=7 es menor a k=8
table(Errork7<Errork9)#Verificar si el errror de k=7 es menor a k=9
table(Errork8<Errork7)#Verificar si el errror de k=8 es menor a k=7
table(Errork8<Errork9)#Verificar si el errror de k=8 es menor a k=9
table(Errork9<Errork7)#Verificar si el errror de k=9 es menor a k=7
table(Errork9<Errork8)#Verificar si el errror de k=9 es menor a k=8

ErrorMean<-mse(data.ImptMean,datamse)
ErrorMedian<-mse(data.ImptMedian,datamse) 
ErrorMode<-mse(data.ImptMode,datamse)

table(Errork8<ErrorMean)#Verificar si el errror de k=8 es menor a mean
table(Errork8<ErrorMedian)#Verificar si el errror de k=8 es menor a median
table(Errork8<ErrorMode)#Verificar si el errror de k=8 es menor a mode
table(ErrorMean<Errork8)#Verificar si el errror de mean es menor a k=8
table(ErrorMean<ErrorMedian)#Verificar si el errror de mean es menor a median
table(ErrorMean<ErrorMode)#Verificar si el errror de mean es menor a mode
table(ErrorMedian<Errork8)#Verificar si el errror de median es menor a k=8
table(ErrorMedian<ErrorMean)#Verificar si el errror de median es menor a mean
table(ErrorMedian<ErrorMode)#Verificar si el errror de median menor a mode
table(ErrorMode<Errork8)#Verificar si el errror de mode es menor a k=8
table(ErrorMode<ErrorMean)#Verificar si el errror de mode es menor a mean
table(ErrorMode<ErrorMedian)#Verificar si el errror de mode es menor a median

BestImput<-data.frame(data.ImptKnn8,yeast_sin$Localization_site)
write.table(BestImput, file="C:/Users/Valeria/Documents/Escuela/2018-2/MIND(M)/PROYECTO FINAL/best_imput.data")#Exportar datos para hacer comparación en knime

#DECISION TREE

yeastDisc<- discretizeDF(BestImput)#Discretizar el dataset

write.table(yeastDisc, file="C:/Users/Valeria/Documents/Escuela/2018-2/MIND(M)/PROYECTO FINAL/YeastDiscretizado.data")#Exportar datos para el arbol de decision en knime

#KNN CLASSIFICATION

DataNorm<-read.csv(file="Yeast_Normalizado_kNN.csv",head=TRUE,sep=",",dec = ".")#Importar datos normalizados desde knime
DataKnn<-subset(DataNorm,select=-c(yeast_sin.Localization_site))

set.seed(101)#semilla
test<-1:446
train.yeast<-DataKnn[-test,]
test.yeast<-DataKnn[test,]

dim(train.yeast)#verificar la dimension del training set
dim(test.yeast)#verificar la dimension del test set

train.site<-DataNorm$yeast_sin.Localization_site[-test]
test.site<-DataNorm$yeast_sin.Localization_site[test]

knn.7<-knn(train.yeast,test.yeast,train.site, k=7)
knn.21<-knn(train.yeast,test.yeast,train.site, k=21)

sum(test.site==knn.7)
100*sum(test.site==knn.7)/446 #Porcentaje de clasificacion correcta k=7
table(knn.7,test.site)#matriz de confusión k=7

sum(test.site==knn.21)
100*sum(test.site==knn.21)/446 #Porcentaje de clasificacion correcta k=7
table(knn.21,test.site)#matriz de confusión k=21

#CLUSTERING
yeast_clustering<-subset(yeast_sin, select=-c(Localization_site))
write.table(yeast_clustering, file="C:/Users/Valeria/Documents/Escuela/2018-2/MIND(M)/PROYECTO FINAL/YeastClustering.data")#Exportar datos para el arbol de decision en knime





