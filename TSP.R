#install.packages("igraph")
library(igraph)

#Creamos 10 valores enteros
coordX=sample(1:100,10);coordX
coordY=coordX

#Los vertices de nuestro grafo
matrizCoord=as.matrix(cbind(coordX,coordY));matrizCoord

#Calculamos la matriz de distancias de los valores
matrizDist=as.matrix(dist(coordX)); matrizDist

#Mostramos el grafo que hemos creado
conjunto=c(integer(0),integer(0),integer(0))
for(i in 1:10){
  for(j in 1:10){
    if(i!=j) conjunto=cbind(conjunto,c(i,j,matrizDist[i,j]))
  }
}
grafo=as.data.frame(t(conjunto))

g=graph.data.frame(grafo,directed=FALSE)
lay=layout_in_circle(g, order = V(g))
tkplot(g,edge.label=paste(E(g)$V3,sep=""),vertex.color="yellow",layout=lay)

#Resolvemos el problema TSP exacto con el algoritmo de Held-Karp
minDist=matrix(nrow=bitwShiftL(1,11),ncol=10,rep(1000,bitwShiftL(1,11)*10))
path=matrix(nrow=bitwShiftL(1,11),ncol=10,rep(0,bitwShiftL(1,11)*10))

#Conjunto de nodos posibles que tenemos
S=c(2,3,4,5,6,7,8,9,10)

#Funcion para pasar los elementos de un conjunto a un valor suma de potencias de 2
#Esto se utiliza para que cada conjunto tenga un indice unico
#No se tiene en cuenta el orden de sus elementos, solo los que pertenecen al conjunto
indConj=function(conj){
  sum=0
  for(i in conj){
    sum=sum+2^i
  }
  return(sum)
}

#Algoritmo de Held-Karp TSP
n=10

#Arista minima del primer nodo al resto
for(k in 2:n){
    indI=2^k
    minDist[indI,k]=matrizDist[1,k]
    path[indI,k]=1
}

#Recorrido minimo de todos los nodos menos el primero
len=n-1
for(s in 2:len){
  auxSet=combn(c(S),s) #S in {2,...,n} con |S|=s

for(i in 1:dim(auxSet)[2]){
    auxConj=auxSet[,i] #conjunto de dimension s
    sumSet=indConj(auxConj) #El indice de ese conjunto
    for(k in auxConj){ #Para cada elemento del conjunto
      minimo=1000
      indMin=0
      conjSink=setdiff(auxConj,k) #Conjunto sin ese elemento
      sumSink=indConj(conjSink) #Indice del conjunto sin ese elemento
      for(m in conjSink){ #Para cada elemento del conjunto sin k
        camino=minDist[sumSink,m]+matrizDist[k,m] #Coste de ir de ese conjunto a k
        if(minimo>camino){
          minimo=camino
          indMin=m
        }
      }
      minDist[sumSet,k]=minimo #Para cada conjunto, su minimo
      path[sumSet,k]=indMin
    }
  }
}

#Calculo de la arista que vuelve al nodo 1 tras el recorrido
#Da el ciclo hamiltoniano minimo
conj=combn(c(S),n-1)
sumTodos=indConj(conj)
minimo=1000
indMin=0
for(k in 2:n){
  camino=minDist[sumTodos,k]+matrizDist[k,1]
  if(minimo>camino){
    minimo=camino
    indMin=k
  }
}

sumTodos=sumTodos+2^1
path[sumTodos,1]=indMin

#Valor optimo
opt=minimo

#Camino encontrado
camino=rep(0,n+1)
sumFaltan=sumTodos
nodo=1
camino[n+1]=nodo
ind=n
while(sumFaltan>0){
  resta=2^nodo
  nodo=path[sumFaltan,nodo]
  sumFaltan=sumFaltan-resta
  camino[ind]=nodo
  ind=ind-1
}

#Resultados
print(paste("El minimo es ", opt))
paste(c("El ciclo hamiltoniano es", camino),collapse=" ")

#Comprobacion de que es el camino correcto
suma=0
for(i in 2:length(camino)){
  suma=suma+matrizDist[camino[i-1],camino[i]]
}
suma==minimo

#Comprobacion con un solver de TSP (heuristica)
#install.packages("TSP")
#install.packages("tspmeta")
library("TSP")
library("tspmeta")

#Construimos la instancia TSP
instTSP=tsp_instance(matrizCoord,matrizDist)

#Resolvemos
tour=run_solver(instTSP,method="2-opt")

#Resultados
print(paste("El minimo con una heuristica es ",tour_length(tour)))
i=which(labels(tour)==1)
paste(c("El camino con una heuristica es", labels(tour)[i:10],labels(tour)[1:i]),collapse=" ")
