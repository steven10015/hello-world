#Se cargan las librerias descargadas
library(ape)
library(phangorn)
library(phytools)
#Secuencia, alineamiento, matriz, arbol
#Menor número de pasos, menos cambios, mejor arbol

fraxatin <- read.phyDat(file = "fraxatine_aligned.fasta", format = "FASTA", type = "AA")
fraxatin
#para poder crear árboles de distancia o parsimonia a través de ella. En este caso la clase debe ser AAbin (Amino Acid Sequences), por eso transformamos el objeto fraxatin en este tipo de clase.
matrizdist <- as.AAbin(fraxatin)
#La función dist.aa (método o función de objeto o clase) calcula una matriz de distancias por pares de las secuencias de aminoácidos
matrizdist <- dist.aa(matrizdist)
matrizdist
#Ahora creamos un árbol con el método de grupo de pares no ponderados con media aritmética (UPGMA) usando la matriz de distancia que acabamos de calcular.
arbolUPGMA <- upgma(matrizdist)
plot(arbolUPGMA)
#Si la longitud de dos ramas es indéntica, significa que las secuancias también son indénticas y en la matriz de distancia la diferencia es de 0.
arbolNJ <- nj(matrizdist)
plot(arbolNJ)
#Personalizar arbol
#Type: es para la forma del arbol, cuadrado, en forma c, etc
#Cex: tamaño
#Font tipo de letra
plot(arbolUPGMA, type= "p", cex=0.8, edge.width=2, edge.color="red", font=3)
plot(arbolUPGMA, type= "c", cex=0.8, edge.width=2, edge.color="blue", font=3)
plot(arbolUPGMA, type= "p", label.offset=0.0005, edge.lty=1, node.pos=2, cex=0.8, edge.width=2, edge.color="black", font=3)
#Formato phytools
plotTree(arbolNJ)
plotTree(arbolNJ, ftype="b", fsize=0.8, offset=1, color="red", lwd=2)
plotTree(ladderize(arbolNJ))
#Guardar arbol
write.tree(arbolNJ, file = "file_name.nex")
#Leerlo
read.tree(file = "file_name.nex")
#Enraizar
arbolNJraiz <-root(arbolNJ, outgroup = "Ornitorrinco", r = TRUE)
plot(arbolNJraiz)
arbolNJraiz <-root(arbolNJ, outgroup = "Ornitorrinco", r = TRUE)
plot(arbolNJraiz)
#Enraizado con UPGMA
arbolUPGMAraiz <-root(arbolUPGMA, outgroup = "Ornitorrinco", r=TRUE)
plot(arbolUPGMAraiz)
#Visualizar dos arboles a la vez
layout(matrix(c(1,2)), height=c(10,10))
par(mar=c(1,1,1,1))
plot(arbolUPGMAraiz, label.offset=0.0005, main="ARBOL UPGMA", cex=0.4)
plot(arbolNJraiz, label.offset=0.0005, main="ARBOL NJ", cex=0.4)

#Parsimonia permite comparar dos arboles al mismo tiempo
#Se estima el número de pasos de un arbol
parsimony(arbolUPGMAraiz, fraxatin)
#Ahora sin raiz
parsimony(arbolUPGMA, fraxatin)
#Con el método optim.parsimony se obtiene el árbol con mejor parsimonia. 
#Este método permite encontrar árboles bajo máxima parsimonia usando árboles de distancia de inicio.
mejorUPGMA <- optim.parsimony(arbolUPGMAraiz, fraxatin)
#Lo mismo con el arbol NJ
mejorNJ <- optim.parsimony(arbolNJraiz, fraxatin)
#Intercambio de ramas
fraxatin_parsimonia <- pratchet(fraxatin, all = TRUE)
#Vemos la cantdad de arboles con ese número de pasos
fraxatin_parsimonia
#Se los enraiza para compararlos
fraxatin_parsimoniaR <- root(phy = fraxatin_parsimonia, outgroup = "Ornitorrinco")
plot(fraxatin_parsimoniaR, cex = 0.6)
#Para hacer un árbol de consenso estricto 
#podemos usar el método ape con parámetro p de 1, que corresponde a un 100% de consenso entre ramas.
estrictode100 <- consensus(fraxatin_parsimoniaR, p = 1)
plot(estrictode100, cex = .6)
#Arbol menos estricto
estrictode30 <- consensus(fraxatin_parsimoniaR, p = 0.3)
plot(estrictode30, cex = .6)
#BookstrappingEl soporte de los árboles se puede evaluar mediante bootstrapping
#Consiste en crear varias seudoréplicas de una matriz con remplazamiento. 
#Los cambios en las réplicas se basan en el uso diferente de los caracteres
#generando árboles en los que algunos caracteres se repiten o no se usan.
arbolesbootstrap <- bootstrap.phyDat(fraxatin, FUN = pratchet, bs = 10)
#La rutina anterior genera entonces 10 árboles pseudoréplicas.
plot(arbolesbootstrap, cex = .6)
#Se genera conseno del 60%
estricto60 <- consensus(arbolesbootstrap, p = 0.6)
plot(estricto60, cex = .6)

#Determinar verosimilitud
arbolazar <- rtree(n = 11, tip.label = names(fraxatin))
plot(arbolazar, cex = .5)
#Se enraiza
arbolazarR <- root(phy = arbolazar, outgroup = "Ornitorrinco")
plot(ladderize(arbolazarR), cex = .5); add.scale.bar()
#Cañcula verosimilitud
ajustado <- pml(arbolazarR, fraxatin)
ajustado

ajustadoconDay <- optim.pml(object = ajustado, model = "Dayhoff", rearrangement = "ratchet")

ajustadoconDay$tree

ajustadoconDayraíz <- root(ajustadoconDay$tree, outgroup = "Ornitorrinco")
plot(ladderize(ajustadoconDayraíz), cex = .5); add.scale.bar()
#Ahora con un modelo diferente a Dayhoff
ajustadoconBlo <- optim.pml(object = ajustado, model = "Blosum62", rearrangement = "ratchet")

ajustadoconJTT <- optim.pml(object = ajustado, model = "JTT", rearrangement = "ratchet")

#Podemos comparar los modelos calculando el Criterio de información de Akaike AIC:
AIC(ajustadoconDay, ajustadoconBlo, ajustadoconJTT)
# Se usa el modelo Jones-Taylor-Thornton para evaluar la distancia entre secuencias de proteínas y optimiza la verosimilitud.
mejorarbol <- optim.pml(
  object = ajustadoconDay, 
  model = "JTT", 
  rearrangement = "ratchet")

mejorarbol
mejorarbolR <- root(mejorarbol$tree, outgroup = "Ornitorrinco")
plot(ladderize(mejorarbolR), cex = 0.5); add.scale.bar()

#Este es el arbol final con mayor verosimilitud
