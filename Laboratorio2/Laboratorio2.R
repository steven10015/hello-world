library(tidyverse)
data(starwars)
#Abres starwars
starwars
#Seleccionas las columnas excepto "-" el nombre
starwars  %>% select(-name)
#Seleccionas las columnas cuyos nombres contengan un "_"
starwars %>% select(contains("_"))
#Seleccionas columnas que empiecen con "S"
starwars %>% select(starts_with("s"))
#Creas un date frame con nombres y planetas de origen (homeworld)
#Recuerda que con <- creas un dataframe
homeworld <- starwars %>% select(name, homeworld)
#En estas filtras la columna especies, seleccionas solo los que sean humanos
human <- starwars %>% filter(species == "Human")
#Filtras humanos de Tatooine
starwars %>% filter(species == "Human", homeworld == "Tatooine")
#Crear un nuevo datframe con todas las especies menos los Droides
starwars_nodroids <- starwars %>% filter(species != "Droid")

#¿Cuántos registros cumplen las condiciones finales?
#Teniendo en cuenta solo los dataframes, en el caso de humanos solo 35 objetos cumplen con el requerimentos
#En nodroids cumplen las condiciones 77 objetos
#Si contaramos humans y tatooine hubo 8 objetos que cumplieron

#El comando group_by sirve para agrupar
#El comando tally cuenta los objetos de la variable que pongamos
starwars %>% group_by(species) %>% tally()
starwars %>% group_by(species, gender) %>% tally()
#Creamos un dataframe
table_gender <- starwars %>% group_by(species, gender) %>% tally()
#na.rm=T quiere decir que elima los NA (valores No Asignados o sin datos)
starwars %>% group_by(species) %>% summarise(mean_height = mean(height, na.rm = T),mean_mass = mean(mass,na.rm = T))
#¿Cómo calcularías la desviación estándar (sd) de esos parámetros?
starwars %>% group_by(species) %>% summarise(sd_height = sd(height, na.rm = T),sd_mass = sd(mass, na.rm = T) )

#Grafica de altura vs masa
#geom_point sirve para insertar los puntos en la grafica
ggplot(starwars, aes(height, mass)) + geom_point()

#Puedes modificar el color 
ggplot(starwars, aes(height, mass)) + geom_point(colour = "red")

#Modificando el color y el punto
ggplot(starwars, aes(height, mass)) + geom_point(colour = "purple", pch = 3)

#Modificando el color y el fondo 
ggplot(starwars, aes(height, mass)) + geom_point(colour = "red") + theme_light()

#Ejercicio, nuevo gráfico sin el personaje super pesado
#Primero cree un df con la variable peso y nombre para identificarlo
pesos <- starwars %>% select(name, mass)
#Luego ya se crea un dataset sin este personaje
no_heavy <- starwars %>% filter(name != "Jabba Desilijic Tiure")
#Se crea el gráfico
ggplot(no_heavy, aes(height, mass)) + geom_point(colour = "red") + theme_light()


tabla_toy <- toy %>% group_by(Sex) %>% summarise(mean_height = mean(Height_cm, na.rm = T),mean_weight = mean(Weight_Kg,na.rm = T), mean_imc = mean(IMC,na.rm = T), mean_ias = mean(IAS,na.rm = T), mean_ccintura = mean(Ccintura,na.rm = T))
tabla_femeina <- toy %>% filter(Sex == "Women") 
women_overweight <- tabla_femeina %>% filter(IMC_clas == "Overweight")

ggplot(toy, aes(IMC, Weight_Kg)) + geom_point(colour = "blue")+ theme_light()
obesity_overweight <- toy %>% filter(IMC_clas == "Overweight" | IMC_clas == "Obesity") 
ggplot(obesity_overweight, aes(IMC, Weight_Kg)) + geom_point(colour = "orange")+ theme_light()

install.packages("ape")
install.packages("phangorn")
install.packages("phytools")