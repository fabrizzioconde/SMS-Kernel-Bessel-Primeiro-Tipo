#Codigo SMS

# Instalando e carregando pacotes
# Nome dos pacotes
packages <- c("scrime","e1071","kernlab","randomForest", 
              "doParallel", "GA", "ggplot2")

# Instalando os pacotes que ainda nao foram instalados
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Carregando os pacotes
invisible(lapply(packages, library, character.only = TRUE))


# Definindo a trilha de dados
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Parâmetros para simulação
dados<-list();                       #Cria lista de dados.
genotipo<-list();                    #Cria lista de genótipos.
num_individuos<-250;                 #Número de indivíduos na amostra.
num_snp<-1000;                       #Total de marcadores simulados.
list.snp<-list(1,10,20,30,40,50,60)  #Indica quais são os marcadores causais.
list.ia<-list(2,1,3,3,3,2,2)         #Somente efeitos aditivos.


#Executa a simulação do genótipo e do fenótipo simultaneamente
simulacao<-simulateSNPglm(n.obs=num_individuos,
                          n.snp=num_snp,
                          list.ia=list.ia,
                          list.snp=list.snp,
                          beta0=0,
                          beta=c(200,200,200,900,200,200,200),
                          maf=c(0.1,0.4),
                          err.fun=rnorm,
                          rand=123)

#Definindo o genótipo
genotipo[[1]]<-as.data.frame(simulacao$x)

#Definindo o fenotipo
fenotipo<-as.data.frame(simulacao$y)
names(fenotipo)<-"fenotipo"
pdf(file="Histograma_fenotipo_6.pdf",height=5,width=9)
hist(fenotipo$fenotipo,xlab="Fenotipo simulado",ylab="Numero de touros",col="gray",main="")
dev.off()
pdf(file="Boxplot_fenotipo_6.pdf",height=5,width=9)
boxplot(fenotipo$fenotipo, ylab="Fenotipo simulado")
dev.off()

#Testando a normalidade do Fenotipo e imprimindo o resultado
(teste_normal<-shapiro.test(fenotipo$fenotipo))
ifelse(teste_normal$p.value >= 0.05, 
       "O fenótipo segue distribuição normal", 
       "O fenótipo NÂO segue distribuição normal")

#Definindo o dataframe genotipo com fenotipo
dados[[1]]<-as.data.frame(cbind(genotipo[[1]],fenotipo))
colnames(dados[[1]])[ncol(dados[[1]])]<-"fenotipo"

save(dados,file="Simulacao6.RData")

# Implementando a função validacao_cruzada para o Kernel Radial e Linear.

validacao_cruzada<-function(data,folds,gamma,cost,epsilon,kernel) 
{
  datalength = nrow(data)
  index <- 1:datalength
  size = trunc(datalength/folds)
  set.seed(123)
  geral = matrix(sample(index), ncol=folds, byrow=TRUE)
  
  mse = vector('double', folds)
  mape = vector('double', folds)
  corr = vector('double', folds)
  r2 = vector('double', folds)
  r2_adj = vector('double', folds)
  
  #  svm.pred<-foreach(i = 1:folds, 
  #                    .combine=rbind,
  #                    .inorder=TRUE,
  #                    .packages=c("e1071","rpart")) %dopar% { 
  
  for(i in 1:folds)
  {
    testv  = geral[,i];
    trainv = index[-testv]
    testset<-data.frame()
    trainset<-data.frame()
    
    testset  = as.data.frame(na.omit(data[testv,]))
    trainset = as.data.frame(na.omit(data[trainv,]))
    
    svmR.model <- svm(fenotipo~., data = trainset, kernel=kernel, gamma=gamma, cost=cost,epsilon=epsilon)
    #svmR.model <- svm(fenotipo~., data = trainset, kernel="linear", gamma=gamma, cost=cost,epsilon=epsilon)
    seltestset = as.data.frame(testset[,-ncol(testset)]);
    names(seltestset) = names(testset)[-ncol(testset)];
    svmR.pred <- predict(svmR.model, seltestset) 
    
    mse[i] = sum((svmR.pred-testset$fenotipo)^2)/dim(testset)[1]
    mape[i] = (sum( abs( (testset$fenotipo - svmR.pred) / testset$fenotipo) ) / length(testset$fenotipo)) * 100;
    corr[i]<-cor(svmR.pred,testset$fenotipo,method="pearson")
    r2[i] <- cor(svmR.pred,testset$fenotipo,method="pearson")^2
    r2_adj[i] <- 1-(1-r2[i])*((nrow(testset)-1)/(nrow(testset)-ncol(testset)-1))
  }
  
  
  
  resultado<-c(mean(mse),sd(mse),mean(mape),sd(mape),mean(corr),sd(corr),mean(r2),sd(r2),mean(r2_adj),sd(r2_adj))
  
  return(resultado)
}


########


# Implementando a função validacao_cruzada para o Kernel Bessel.

validacao_cruzada_bessel<-function(data,folds,sigma,order,degree,cost,epsilon) 
{
  datalength = nrow(data)
  index <- 1:datalength
  size = trunc(datalength/folds)
  set.seed(123)
  geral = matrix(sample(index), ncol=folds, byrow=TRUE)
  
  mse = vector('double', folds)
  mape = vector('double', folds)
  corr = vector('double', folds)
  r2 = vector('double', folds)
  r2_adj = vector('double', folds)
  
  bessel <- besseldot(sigma = sigma, 
                      order = order, 
                      degree = degree)
  
  for(i in 1:folds)
  {
    testv  = geral[,i];
    trainv = index[-testv]
    testset<-data.frame()
    trainset<-data.frame()
    
    testset  = as.data.frame(na.omit(data[testv,]))
    trainset = as.data.frame(na.omit(data[trainv,]))
    
    svmR.model <- ksvm(fenotipo~., data = trainset, kernel=bessel, gamma=gamma, cost=cost,epsilon=epsilon)
    #svmR.model <- svm(fenotipo~., data = trainset, kernel="linear", gamma=gamma, cost=cost,epsilon=epsilon)
    seltestset = as.data.frame(testset[,-ncol(testset)]);
    names(seltestset) = names(testset)[-ncol(testset)];
    svmR.pred <- predict(svmR.model, seltestset) 
    
    mse[i] = sum((svmR.pred-testset$fenotipo)^2)/dim(testset)[1]
    mape[i] = (sum( abs( (testset$fenotipo - svmR.pred) / testset$fenotipo) ) / length(testset$fenotipo)) * 100;
    corr[i]<-cor(svmR.pred,testset$fenotipo,method="pearson")
    r2[i] <- cor(svmR.pred,testset$fenotipo,method="pearson")^2
    r2_adj[i] <- 1-(1-r2[i])*((nrow(testset)-1)/(nrow(testset)-ncol(testset)-1))
  }
  
  
  
  resultado<-c(mean(mse),sd(mse),mean(mape),sd(mape),mean(corr),sd(corr),mean(r2),sd(r2),mean(r2_adj),sd(r2_adj))
  
  return(resultado)
}


#Funcao para o calculo do valor p bruto e ajustado pela correcao de Bonferroni

valor.p<-function(genotipo_fenotipo)
{
  valor_p_bruto<-vector();
  valor_p_ajustado<-vector();
  y<-genotipo_fenotipo[,ncol(genotipo_fenotipo)];
  m<-ncol(genotipo_fenotipo)-1;
  model_regression<-list();
  saida<-data.frame();
  
  for (i in 1:(ncol(genotipo_fenotipo)-1))
  {
    x<-genotipo_fenotipo[,i];
    model_regression[[i]]<-lm(y~x,data=genotipo_fenotipo);
    valor_p_bruto[i]<-ifelse(all(x==3)|all(x==2)|all(x==1),1,summary(model_regression[[i]])[[4]][2,4]);
    
  }
  valor_p_ajustado<- m*valor_p_bruto;  
  valor_p<-data.frame();
  valor_p<-cbind(valor_p_bruto,valor_p_ajustado);
  colnames(valor_p)<-c("Valor p bruto","Valor p ajustado")
  rownames(valor_p)<-names(genotipo_fenotipo)[1:(ncol(genotipo_fenotipo)-1)];
  return(valor_p)
}

################ Inicio do SMS##################


#Definicao de listas para cada kernel do SVR
#svr_test<-list()
mean_svr_RF_list<-list() #Cria lista para a media do SVR sobre o rank da RF
GA<-list()
minimo<-list()
corte<-list()
snps_selec_corte<-list()
snps_selec_ref<-list()
percentual_snps<-0.95;

i=1 #Contador do kernel do SVR utilizado

#Random Forest
ntree<-4000           #Numero de arvores dentro da floresta aleatoria.

data_temp<-as.data.frame(dados[[1]]) #Transforma a base de dados em dataframe.

set.seed(1) #Semente aleatoria para floresta aleatoria.

RF <- randomForest(fenotipo~., data=data_temp,
                   ntree=ntree,
                   mtry=ncol(data_temp)-1,
                   importance=TRUE)    

rank_RF<-sort(importance(RF)[,1],decreasing=TRUE)
#View(rank_RF)


#Erro MSE do SVR em relacao ao rank da RF

#Inicializacao dos parâmetros do SVR
gamma = 0.01
cost = 1.0
epsilon = 0.1
folds = 10
kernel = "linear"

#Etapa de Corte com o SVR sobre o Rank da RF
mean_svr_RF_list[[i]]<-vector()
passo<-10;

limite<-floor((length(names(dados[[1]][,-dim(dados[[1]])[2]]))/passo)*percentual_snps);


#Constroi a sequencia crescente de SNPs mais importantes a partir do rank da RF
for (cont in 1:(limite+1))
{
  
  print(cont)
  j<-cont*passo
  #Ajuste quando j==1, pois o R demonstra um erro quando j==1.
  if (j==1) {
    
    var_sel<-names(rank_RF)[1]
    
  }else{
    
    var_sel<-names(rank_RF)[1:j]
    
  }
  
  #Avalia os 10-fold para o SVR na sequencia crescente de SNPs a partir do rank da RF
  svr_test<- validacao_cruzada(data = dados[[1]][c(var_sel,"fenotipo")],
                               folds = folds, 
                               gamma = gamma,
                               cost = cost,
                               epsilon = epsilon, 
                               kernel = kernel)
  mean_svr_RF_list[[i]][cont]<-svr_test[1]
} 

#Constroi o Grafico do MSE do SVR sobre a RF
pdf(file="Grafico_MSE_SVR_linear.pdf",height=5,width=9)
plot(seq(passo,limite*(passo)+10,by=passo),
     mean_svr_RF_list[[i]],
     type="o",
     lwd=2,
     xlab="Grupo de Marcadores",
     ylab="MSE do SVR")

#ESCOLHE O PRIMEIRO PONTO DE CORTE
minimo[[i]]<-which.min(mean_svr_RF_list[[i]])
corte[[i]]<-(minimo[[i]]+1)*passo #Adotou-se o segundo menor MSE do SVR, pois o menor foi o referente ao primeiro marcador SNP2
abline(v=corte[[i]],col="black",lty=2)
dev.off()


snps_selec_corte[[i]]<-names(rank_RF[1:corte[[i]]])
snps_selec_corte[[i]]


#Construindo o dataframe com rank_RF juntamente com o p_valor dos SNPs selecionados na etapa de corte
rank_global_selec<-list()
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]][c(snps_selec_corte[[i]],'fenotipo')])
rank_global_selec[[i]]<-cbind(rank_RF[1:corte[[i]]],valor_p)
#View(rank_global_selec[[i]])

#Construindo o dataframe com rank_RF juntamente com o p_valor de todos os SNPs da base de dados inicial.
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]])
rank_global<-cbind(importance(RF)[,1],valor_p)
#View(rank_global)

#Selecionando somente o genotipo apos o corte
genotipo[[2]]<-dados[[1]][,names(rank_RF[1:corte[[i]]])]

#Selecionando o genotipo e o fenotipo apos o corte
dados[[2]]<-cbind(genotipo[[2]],dados[[1]]$fenotipo)
names(dados[[2]])[ncol(dados[[2]])]<-"fenotipo"


###AG para a SEGUNDA selecao dos marcadores na base de dados inicial (REFINAMENTO)

#Declara o objeto GA como uma lista

f<-function(x)
{
  inc<-which(x==1)
  
  cat("inc = ",inc, "\n")
  
  dados_validacao<- dados[[2]][,c(inc,ncol(dados[[2]]))]
  model<-validacao_cruzada(dados_validacao,folds,gamma,cost,epsilon,kernel);
  media<-model[5];
  return(media);       
}
fitness<-function(x) {f(x)}


#Parametros do GA
run = 30            #Numero maximo de geracoes do GA sem melhoria.
maxiter = 10        #Numero maximo de geracoes do GA.
pcross = 0.8        #Probabilidade de crossing over
pmut = 0.1          #Probabilidade de mutacao.
elitism = 5         #Numero de melhores indiviuos do GA que irao para proxima geracao sem alteracao alguma.
popSize=100         #Tamanho da populacaocada geracao do GA

GA[[i]]<-ga(type="binary",
            fitness=fitness, 
            nBits=ncol(genotipo[[2]]),
            popSize=popSize,
            names=colnames(genotipo[[2]]),
            maxiter=maxiter,
            seed=i,
            parallel=TRUE,
            run=run, 
            suggestions=matrix(rep(1,ncol(genotipo[[2]])),
                               ncol=ncol(genotipo[[2]])))

snps_selec_ref[[i]]<-names(GA[[i]]@solution[,which(GA[[i]]@solution==1)])

snps_selec_SVR_linear_GA_corr<-snps_selec_ref[[i]]
snps_selec_SVR_linear_GA_corr

#Relatorio do AG
summary(GA[[i]])

#Grafico do AG
plot(GA[[i]])

#Salvando o grafico do AG com ggplot2
pdf(file="Grafico_GA_SVR_linear.pdf",height=5,width=9)
geracao<-seq(1,GA[[1]]@iter,by=1)
mean_fitness   <- GA[[i]]@summary[,2]
median_fitness <- GA[[i]]@summary[,4]
best_fitness   <- GA[[i]]@summary[,1]
Estatisticas <- c(rep("Mediana", length(geracao)), rep("Media", length(geracao)),rep("Melhor", length(geracao)) )
data_grafico_1 <- data.frame(
  Geracao = rep(geracao, 3),
  Aptidao = c(mean_fitness, median_fitness, best_fitness),
  Estatisticas = Estatisticas
)

ggplot(data_grafico_1, aes(x = Geracao, y = Aptidao, group = Estatisticas)) +
  geom_line(aes(colour = Estatisticas, linetype = Estatisticas),size=2) + 
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Geracao),
                                  max(data_grafico_1$Geracao), by = 1))
dev.off()


############ SMS com SVR com kernel bessel 1 #############

i=2 #Contador do kernel "bessel_1" do SVR utilizado



#Erro MSE do SVR em relacao ao rank da RF

#Inicializacao dos parametros do SVR com kernel Bessel 1
folds = 10
sigma = 0.10
order = 1.00
degree = 1
cost = 1.0
epsilon = 0.1

#Etapa de Corte com o SVR sobre o Rank da RF
mean_svr_RF_list[[i]]<-vector()
passo<-10;

limite<-floor((length(names(dados[[1]][,-dim(dados[[1]])[2]]))/passo)*percentual_snps);


#Constroi a sequencia crescente de SNPs mais importantes a partir do rank da RF
for (cont in 1:(limite+1))
{
  
  print(cont)
  j<-cont*passo
  #Ajuste quando j==1, pois o R demonstra um erro quando j==1.
  if (j==1) {
    
    var_sel<-names(rank_RF)[1]
    
  }else{
    
    var_sel<-names(rank_RF)[1:j]
    
  }
  
  #Avalia os 10-fold para o SVR na sequencia crescente de SNPs a partir do rank da RF
  svr_test<- validacao_cruzada_bessel(data = dados[[1]][c(var_sel,"fenotipo")],
                                      folds = folds, 
                                      sigma = sigma,
                                      order = order,
                                      degree = degree,
                                      cost = cost,
                                      epsilon = epsilon)
  mean_svr_RF_list[[i]][cont]<-svr_test[1]
} 

#Constroi o Grafico do MSE do SVR sobre a RF
pdf(file="Grafico_MSE_SVR_bessel_1.pdf",height=5,width=9)
plot(seq(passo,limite*(passo)+10,by=passo),
     mean_svr_RF_list[[i]],
     type="o",
     lwd=2,
     xlab="Grupo de Marcadores",
     ylab="MSE do SVR")

#ESCOLHE O PRIMEIRO PONTO DE CORTE
minimo[[i]]<-which.min(mean_svr_RF_list[[i]])
corte[[i]]<-(minimo[[i]]+1)*passo #Adotou-se o segundo menor MSE do SVR, pois o menor foi o referente ao primeiro marcador SNP2
abline(v=corte[[i]],col="black",lty=2)
dev.off()


snps_selec_corte[[i]]<-names(rank_RF[1:corte[[i]]])
snps_selec_corte[[i]]

#Construindo o dataframe com rank_RF juntamente com o p_valor dos SNPs selecionados na etapa de corte
rank_global_selec<-list()
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]][c(snps_selec_corte[[i]],'fenotipo')])
rank_global_selec[[i]]<-cbind(rank_RF[1:corte[[i]]],valor_p)
#View(rank_global_selec[[i]])

#Construindo o dataframe com rank_RF juntamente com o p_valor de todos os SNPs da base de dados inicial.
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]])
rank_global<-cbind(importance(RF)[,1],valor_p)
#View(rank_global)

#Selecionando somente o genotipo apos o corte
genotipo[[2]]<-dados[[1]][,names(rank_RF[1:corte[[i]]])]

#Selecionando o genotipo e o fenotipo apos o corte
dados[[2]]<-cbind(genotipo[[2]],dados[[1]]$fenotipo)
names(dados[[2]])[ncol(dados[[2]])]<-"fenotipo"


###AG para a SEGUNDA selecao dos marcadores na base de dados inicial (REFINAMENTO)

#Declara o objeto GA como uma lista

f<-function(x)
{
  inc<-which(x==1)
  
  cat("inc = ",inc, "\n")
  
  dados_validacao<- dados[[2]][,c(inc,ncol(dados[[2]]))]
  model<-validacao_cruzada_bessel(dados_validacao,
                                  folds,
                                  sigma,
                                  order,
                                  degree,
                                  cost,
                                  epsilon);
  media<-model[5];
  return(media);       
}
fitness<-function(x) {f(x)}


#Parametros do GA
run = 30            #Numero maximo de geracoes do GA sem melhoria.
maxiter = 10        #Numero maximo de geracoeso GA.
pcross = 0.8        #Probabilidade de crossing over
pmut = 0.1          #Probabilidade de mutacao.
elitism = 5         #Numero de melhores individuos do GA que irao para proxima geracaos em alteracao alguma.
popSize=100         #Tamanho da populacao em cada geracao do GA.

GA[[i]]<-ga(type="binary",
            fitness=fitness, 
            nBits=ncol(genotipo[[2]]),
            popSize=popSize,
            names=colnames(genotipo[[2]]),
            maxiter=maxiter,
            seed=i,
            parallel=TRUE,
            run=run, 
            suggestions=matrix(rep(1,ncol(genotipo[[2]])),
                               ncol=ncol(genotipo[[2]])))

snps_selec_ref[[i]]<-names(GA[[i]]@solution[,which(GA[[i]]@solution==1)])

snps_selec_SVR_bessel_1_GA_corr<-snps_selec_ref[[i]]
snps_selec_SVR_bessel_1_GA_corr

#Relatorio do AG
summary(GA[[i]])

#Grafico do AG
plot(GA[[i]])

#Salvando o grafico do AG com ggplot2
pdf(file="Grafico_GA_SVR_bessel_1.pdf",height=5,width=9)
geracao<-seq(1,GA[[1]]@iter,by=1)
mean_fitness   <-  GA[[i]]@summary[,2]
median_fitness <- GA[[i]]@summary[,4]
best_fitness   <- GA[[i]]@summary[,1]
Estatisticas <- c(rep("Mediana", length(geracao)), rep("Media", length(geracao)),rep("Melhor", length(geracao)) )
data_grafico_1 <- data.frame(
  Geracao = rep(geracao, 3),
  Aptidao = c(mean_fitness, median_fitness, best_fitness),
  Estatisticas = Estatisticas
)

ggplot(data_grafico_1, aes(x = Geracao, y = Aptidao, group = Estatisticas)) +
  geom_line(aes(colour = Estatisticas, linetype = Estatisticas),size=2) + 
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Geracao), max(data_grafico_1$Geracao), by = 1))
dev.off()



############ SMS com SVR com kernel bessel 2 #############

i=3 #Contador do kernel "bessel_2" do SVR utilizado



#Erro MSE do SVR em relacao ao rank da RF

#Inicializacao dos parametros do SVR com kernel Bessel 2
folds = 10
sigma = 0.10
order = 1.50
degree = 1
cost = 1.0
epsilon = 0.1

#Etapa de Corte com o SVR sobre o Rank da RF
mean_svr_RF_list[[i]]<-vector()
passo<-10;

limite<-floor((length(names(dados[[1]][,-dim(dados[[1]])[2]]))/passo)*percentual_snps);


#Constroi a sequencia crescente de SNPs mais importantes a partir do rank da RF
for (cont in 1:(limite+1))
{
  
  print(cont)
  j<-cont*passo
  #Ajuste quando j==1, pois o R demonstra um erro quando j==1.
  if (j==1) {
    
    var_sel<-names(rank_RF)[1]
    
  }else{
    
    var_sel<-names(rank_RF)[1:j]
    
  }
  
  #Avalia os 10-fold para o SVR na sequencia crescente de SNPs a partir do rank da RF
  svr_test<- validacao_cruzada_bessel(data = dados[[1]][c(var_sel,"fenotipo")],
                                      folds = folds, 
                                      sigma = sigma,
                                      order = order,
                                      degree = degree,
                                      cost = cost,
                                      epsilon = epsilon)
  mean_svr_RF_list[[i]][cont]<-svr_test[1]
} 

#Constroi o Grafico do MSE do SVR sobre a RF
pdf(file="Grafico_MSE_SVR_bessel_2.pdf",height=5,width=9)
plot(seq(passo,limite*(passo)+10,by=passo),
     mean_svr_RF_list[[i]],
     type="o",
     lwd=2,
     xlab="Grupo de Marcadores",
     ylab="MSE do SVR")

#ESCOLHE O PRIMEIRO PONTO DE CORTE
minimo[[i]]<-which.min(mean_svr_RF_list[[i]])
corte[[i]]<-(minimo[[i]]+1)*passo #Adotou-se o segundo menor MSE do SVR, pois o menor foi o referente ao primeiro marcador SNP2
abline(v=corte[[i]],col="black",lty=2)
dev.off()


snps_selec_corte[[i]]<-names(rank_RF[1:corte[[i]]])
snps_selec_corte[[i]]

#Construindo o dataframe com rank_RF juntamente com o p_valor dos SNPs selecionados na etapa de corte
rank_global_selec<-list()
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]][c(snps_selec_corte[[i]],'fenotipo')])
rank_global_selec[[i]]<-cbind(rank_RF[1:corte[[i]]],valor_p)
#View(rank_global_selec[[i]])

#Construindo o dataframe com rank_RF juntamente com o p_valor de todos os SNPs da base de dados inicial.
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]])
rank_global<-cbind(importance(RF)[,1],valor_p)
#View(rank_global)

#Selecionando somente o genotipo apos o corte
genotipo[[2]]<-dados[[1]][,names(rank_RF[1:corte[[i]]])]

#Selecionando o genotipo e o fenotipo apos o corte
dados[[2]]<-cbind(genotipo[[2]],dados[[1]]$fenotipo)
names(dados[[2]])[ncol(dados[[2]])]<-"fenotipo"


###AG para a SEGUNDA selecao dos marcadores na base de dados inicial (REFINAMENTO)

#Declara o objeto GA como uma lista

f<-function(x)
{
  inc<-which(x==1)
  
  cat("inc = ",inc, "\n")
  
  dados_validacao<- dados[[2]][,c(inc,ncol(dados[[2]]))]
  model<-validacao_cruzada_bessel(dados_validacao,
                                  folds,
                                  sigma,
                                  order,
                                  degree,
                                  cost,
                                  epsilon);
  media<-model[5];
  return(media);       
}
fitness<-function(x) {f(x)}


#Parametros do GA
run = 30            #Numero maximo de geracoes do GA sem melhoria.
maxiter = 10        #Numero maximo de geracoeso GA.
pcross = 0.8        #Probabilidade de crossing over
pmut = 0.1          #Probabilidade de mutacao.
elitism = 5         #Numero de melhores individuos do GA que irao para proxima geracaos em alteracao alguma.
popSize=100         #Tamanho da populacao em cada geracao do GA.

GA[[i]]<-ga(type="binary",
            fitness=fitness, 
            nBits=ncol(genotipo[[2]]),
            popSize=popSize,
            names=colnames(genotipo[[2]]),
            maxiter=maxiter,
            seed=i,
            parallel=TRUE,
            run=run, 
            suggestions=matrix(rep(1,ncol(genotipo[[2]])),
                               ncol=ncol(genotipo[[2]])))

snps_selec_ref[[i]]<-names(GA[[i]]@solution[,which(GA[[i]]@solution==1)])

snps_selec_SVR_bessel_1_GA_corr<-snps_selec_ref[[i]]
snps_selec_SVR_bessel_1_GA_corr

#Relatorio do AG
summary(GA[[i]])

#Grafico do AG
plot(GA[[i]])

#Salvando o grafico do AG com ggplot2
pdf(file="Grafico_GA_SVR_bessel_2.pdf",height=5,width=9)
geracao<-seq(1,GA[[1]]@iter,by=1)
mean_fitness   <-  GA[[i]]@summary[,2]
median_fitness <- GA[[i]]@summary[,4]
best_fitness   <- GA[[i]]@summary[,1]
Estatisticas <- c(rep("Mediana", length(geracao)), rep("Media", length(geracao)),rep("Melhor", length(geracao)) )
data_grafico_1 <- data.frame(
  Geracao = rep(geracao, 3),
  Aptidao = c(mean_fitness, median_fitness, best_fitness),
  Estatisticas = Estatisticas
)

ggplot(data_grafico_1, aes(x = Geracao, y = Aptidao, group = Estatisticas)) +
  geom_line(aes(colour = Estatisticas, linetype = Estatisticas),size=2) + 
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Geracao), max(data_grafico_1$Geracao), by = 1))
dev.off()

############ SMS com SVR com kernel bessel 3 #############

i=4 #Contador do kernel "bessel_3" do SVR utilizado



#Erro MSE do SVR em relacao ao rank da RF

#Inicializacao dos parametros do SVR com kernel Bessel 3
folds = 10
sigma = 0.30
order = 1.00
degree = 1
cost = 1.0
epsilon = 0.1

#Etapa de Corte com o SVR sobre o Rank da RF
mean_svr_RF_list[[i]]<-vector()
passo<-10;

limite<-floor((length(names(dados[[1]][,-dim(dados[[1]])[2]]))/passo)*percentual_snps);


#Constroi a sequencia crescente de SNPs mais importantes a partir do rank da RF
for (cont in 1:(limite+1))
{
  
  print(cont)
  j<-cont*passo
  #Ajuste quando j==1, pois o R demonstra um erro quando j==1.
  if (j==1) {
    
    var_sel<-names(rank_RF)[1]
    
  }else{
    
    var_sel<-names(rank_RF)[1:j]
    
  }
  
  #Avalia os 10-fold para o SVR na sequencia crescente de SNPs a partir do rank da RF
  svr_test<- validacao_cruzada_bessel(data = dados[[1]][c(var_sel,"fenotipo")],
                                      folds = folds, 
                                      sigma = sigma,
                                      order = order,
                                      degree = degree,
                                      cost = cost,
                                      epsilon = epsilon)
  mean_svr_RF_list[[i]][cont]<-svr_test[1]
} 

#Constroi o Grafico do MSE do SVR sobre a RF
pdf(file="Grafico_MSE_SVR_bessel_3.pdf",height=5,width=9)
plot(seq(passo,limite*(passo)+10,by=passo),
     mean_svr_RF_list[[i]],
     type="o",
     lwd=2,
     xlab="Grupo de Marcadores",
     ylab="MSE do SVR")

#ESCOLHE O PRIMEIRO PONTO DE CORTE
minimo[[i]]<-minimo[[i]]<-which.min(mean_svr_RF_list[[i]][1]) ######## Tive que escolher manualmente o segundo menor MSE, pois o menor valor foi a base de dados com os 100 SNPs.
corte[[i]]<-(minimo[[i]]+1)*passo #Adotou-se o segundo menor MSE do SVR, pois o menor foi o referente ao primeiro marcador SNP2
abline(v=corte[[i]],col="black",lty=2)
dev.off()


snps_selec_corte[[i]]<-names(rank_RF[1:corte[[i]]])
snps_selec_corte[[i]]

#Construindo o dataframe com rank_RF juntamente com o p_valor dos SNPs selecionados na etapa de corte
rank_global_selec<-list()
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]][c(snps_selec_corte[[i]],'fenotipo')])
rank_global_selec[[i]]<-cbind(rank_RF[1:corte[[i]]],valor_p)
#View(rank_global_selec[[i]])

#Construindo o dataframe com rank_RF juntamente com o p_valor de todos os SNPs da base de dados inicial.
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]])
rank_global<-cbind(importance(RF)[,1],valor_p)
#View(rank_global)

#Selecionando somente o genotipo apos o corte
genotipo[[2]]<-dados[[1]][,names(rank_RF[1:corte[[i]]])]

#Selecionando o genotipo e o fenotipo apos o corte
dados[[2]]<-cbind(genotipo[[2]],dados[[1]]$fenotipo)
names(dados[[2]])[ncol(dados[[2]])]<-"fenotipo"


###AG para a SEGUNDA selecao dos marcadores na base de dados inicial (REFINAMENTO)

#Declara o objeto GA como uma lista

f<-function(x)
{
  inc<-which(x==1)
  
  cat("inc = ",inc, "\n")
  
  dados_validacao<- dados[[2]][,c(inc,ncol(dados[[2]]))]
  model<-validacao_cruzada_bessel(dados_validacao,
                                  folds,
                                  sigma,
                                  order,
                                  degree,
                                  cost,
                                  epsilon);
  media<-model[5];
  return(media);       
}
fitness<-function(x) {f(x)}


#Parametros do GA
run = 30            #Numero maximo de geracoes do GA sem melhoria.
maxiter = 10        #Numero maximo de geracoeso GA.
pcross = 0.8        #Probabilidade de crossing over
pmut = 0.1          #Probabilidade de mutacao.
elitism = 5         #Numero de melhores individuos do GA que irao para proxima geracaos em alteracao alguma.
popSize = 100         #Tamanho da populacao em cada geracao do GA.

GA[[i]]<-ga(type="binary",
            fitness=fitness, 
            nBits=ncol(genotipo[[2]]),
            popSize=popSize,
            names=colnames(genotipo[[2]]),
            maxiter=maxiter,
            seed=i,
            parallel=TRUE,
            run=run, 
            suggestions=matrix(rep(1,ncol(genotipo[[2]])),
                               ncol=ncol(genotipo[[2]])))

snps_selec_ref[[i]]<-names(GA[[i]]@solution[,which(GA[[i]]@solution==1)])

snps_selec_SVR_bessel_3_GA_corr<-snps_selec_ref[[i]]
snps_selec_SVR_bessel_3_GA_corr

#Relatorio do AG
summary(GA[[i]])

#Grafico do AG
plot(GA[[i]])

#Salvando o grafico do AG com ggplot2
pdf(file="Grafico_GA_SVR_bessel_3.pdf",height=5,width=9)
geracao<-seq(1,GA[[1]]@iter,by=1)
mean_fitness   <-  GA[[i]]@summary[,2]
median_fitness <- GA[[i]]@summary[,4]
best_fitness   <- GA[[i]]@summary[,1]
Estatisticas <- c(rep("Mediana", length(geracao)), rep("Media", length(geracao)),rep("Melhor", length(geracao)) )
data_grafico_1 <- data.frame(
  Geracao = rep(geracao, 3),
  Aptidao = c(mean_fitness, median_fitness, best_fitness),
  Estatisticas = Estatisticas
)

ggplot(data_grafico_1, aes(x = Geracao, y = Aptidao, group = Estatisticas)) +
  geom_line(aes(colour = Estatisticas, linetype = Estatisticas),size=2) + 
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Geracao), max(data_grafico_1$Geracao), by = 1))
dev.off()


############ SMS com SVR com kernel bessel 4 #############

i=5 #Contador do kernel "bessel_4" do SVR utilizado



#Erro MSE do SVR em relacao ao rank da RF

#Inicializacao dos parametros do SVR com kernel Bessel 4
folds = 10
sigma = 0.30
order = 1.50
degree = 1
cost = 1.0
epsilon = 0.1

#Etapa de Corte com o SVR sobre o Rank da RF
mean_svr_RF_list[[i]]<-vector()
passo<-10;

limite<-floor((length(names(dados[[1]][,-dim(dados[[1]])[2]]))/passo)*percentual_snps);


#Constroi a sequencia crescente de SNPs mais importantes a partir do rank da RF
for (cont in 1:(limite+1))
{
  
  print(cont)
  j<-cont*passo
  #Ajuste quando j==1, pois o R demonstra um erro quando j==1.
  if (j==1) {
    
    var_sel<-names(rank_RF)[1]
    
  }else{
    
    var_sel<-names(rank_RF)[1:j]
    
  }
  
  #Avalia os 10-fold para o SVR na sequencia crescente de SNPs a partir do rank da RF
  svr_test<- validacao_cruzada_bessel(data = dados[[1]][c(var_sel,"fenotipo")],
                                      folds = folds, 
                                      sigma = sigma,
                                      order = order,
                                      degree = degree,
                                      cost = cost,
                                      epsilon = epsilon)
  mean_svr_RF_list[[i]][cont]<-svr_test[1]
} 

#Constroi o Grafico do MSE do SVR sobre a RF
pdf(file="Grafico_MSE_SVR_bessel_4.pdf",height=5,width=9)
plot(seq(passo,limite*(passo)+10,by=passo),
     mean_svr_RF_list[[i]],
     type="o",
     lwd=2,
     xlab="Grupo de Marcadores",
     ylab="MSE do SVR")

#ESCOLHE O PRIMEIRO PONTO DE CORTE
minimo[[i]]<-which.min(mean_svr_RF_list[[i]])
corte[[i]]<-(minimo[[i]]+1)*passo #Adotou-se o segundo menor MSE do SVR, pois o menor foi o referente ao primeiro marcador SNP2
abline(v=corte[[i]],col="black",lty=2)
dev.off()


snps_selec_corte[[i]]<-names(rank_RF[1:corte[[i]]])
snps_selec_corte[[i]]

#Construindo o dataframe com rank_RF juntamente com o p_valor dos SNPs selecionados na etapa de corte
rank_global_selec<-list()
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]][c(snps_selec_corte[[i]],'fenotipo')])
rank_global_selec[[i]]<-cbind(rank_RF[1:corte[[i]]],valor_p)
#View(rank_global_selec[[i]])

#Construindo o dataframe com rank_RF juntamente com o p_valor de todos os SNPs da base de dados inicial.
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]])
rank_global<-cbind(importance(RF)[,1],valor_p)
#View(rank_global)

#Selecionando somente o genotipo apos o corte
genotipo[[2]]<-dados[[1]][,names(rank_RF[1:corte[[i]]])]

#Selecionando o genotipo e o fenotipo apos o corte
dados[[2]]<-cbind(genotipo[[2]],dados[[1]]$fenotipo)
names(dados[[2]])[ncol(dados[[2]])]<-"fenotipo"


###AG para a SEGUNDA selecao dos marcadores na base de dados inicial (REFINAMENTO)

#Declara o objeto GA como uma lista

f<-function(x)
{
  inc<-which(x==1)
  
  cat("inc = ",inc, "\n")
  
  dados_validacao<- dados[[2]][,c(inc,ncol(dados[[2]]))]
  model<-validacao_cruzada_bessel(dados_validacao,
                                  folds,
                                  sigma,
                                  order,
                                  degree,
                                  cost,
                                  epsilon);
  media<-model[5];
  return(media);       
}
fitness<-function(x) {f(x)}


#Parametros do GA
run = 30            #Numero maximo de geracoes do GA sem melhoria.
maxiter = 10        #Numero maximo de geracoeso GA.
pcross = 0.8        #Probabilidade de crossing over
pmut = 0.1          #Probabilidade de mutacao.
elitism = 5         #Numero de melhores individuos do GA que irao para proxima geracaos em alteracao alguma.
popSize=100         #Tamanho da populacao em cada geracao do GA.

GA[[i]]<-ga(type="binary",
            fitness=fitness, 
            nBits=ncol(genotipo[[2]]),
            popSize=popSize,
            names=colnames(genotipo[[2]]),
            maxiter=maxiter,
            seed=i,
            parallel=TRUE,
            run=run, 
            suggestions=matrix(rep(1,ncol(genotipo[[2]])),
                               ncol=ncol(genotipo[[2]])))

snps_selec_ref[[i]]<-names(GA[[i]]@solution[,which(GA[[i]]@solution==1)])

snps_selec_SVR_bessel_4_GA_corr<-snps_selec_ref[[i]]
snps_selec_SVR_bessel_4_GA_corr

#Relatorio do AG
summary(GA[[i]])

#Grafico do AG
plot(GA[[i]])

#Salvando o grafico do AG com ggplot2
pdf(file="Grafico_GA_SVR_bessel_4.pdf",height=5,width=9)
geracao<-seq(1,GA[[1]]@iter,by=1)
mean_fitness   <-  GA[[i]]@summary[,2]
median_fitness <- GA[[i]]@summary[,4]
best_fitness   <- GA[[i]]@summary[,1]
Estatisticas <- c(rep("Mediana", length(geracao)), rep("Media", length(geracao)),rep("Melhor", length(geracao)) )
data_grafico_1 <- data.frame(
  Geracao = rep(geracao, 3),
  Aptidao = c(mean_fitness, median_fitness, best_fitness),
  Estatisticas = Estatisticas
)

ggplot(data_grafico_1, aes(x = Geracao, y = Aptidao, group = Estatisticas)) +
  geom_line(aes(colour = Estatisticas, linetype = Estatisticas),size=2) + 
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Geracao), max(data_grafico_1$Geracao), by = 1))
dev.off()




############ SMS com SVR com kernel bessel 5 #############

i=6 #Contador do kernel "bessel_5" do SVR utilizado



#Erro MSE do SVR em relacao ao rank da RF

#Inicializacao dos parametros do SVR com kernel Bessel 5
folds = 10
sigma = 0.50
order = 1.80
degree = 1
cost = 1.0
epsilon = 0.1

#Etapa de Corte com o SVR sobre o Rank da RF
mean_svr_RF_list[[i]]<-vector()
passo<-10;

limite<-floor((length(names(dados[[1]][,-dim(dados[[1]])[2]]))/passo)*percentual_snps);


#Constroi a sequencia crescente de SNPs mais importantes a partir do rank da RF
for (cont in 1:(limite+1))
{
  
  print(cont)
  j<-cont*passo
  #Ajuste quando j==1, pois o R demonstra um erro quando j==1.
  if (j==1) {
    
    var_sel<-names(rank_RF)[1]
    
  }else{
    
    var_sel<-names(rank_RF)[1:j]
    
  }
  
  #Avalia os 10-fold para o SVR na sequencia crescente de SNPs a partir do rank da RF
  svr_test<- validacao_cruzada_bessel(data = dados[[1]][c(var_sel,"fenotipo")],
                                      folds = folds, 
                                      sigma = sigma,
                                      order = order,
                                      degree = degree,
                                      cost = cost,
                                      epsilon = epsilon)
  mean_svr_RF_list[[i]][cont]<-svr_test[1]
} 

#Constroi o Grafico do MSE do SVR sobre a RF
pdf(file="Grafico_MSE_SVR_bessel_5.pdf",height=5,width=9)
plot(seq(passo,limite*(passo)+10,by=passo),
     mean_svr_RF_list[[i]],
     type="o",
     lwd=2,
     xlab="Grupo de Marcadores",
     ylab="MSE do SVR")

#ESCOLHE O PRIMEIRO PONTO DE CORTE
minimo[[i]]<-minimo[[i]]<-which.min(mean_svr_RF_list[[i]][1]) ######## Tive que escolher manualmente o segundo menor MSE, pois o menor valor foi a base de dados com os 100 SNPs.
corte[[i]]<-(minimo[[i]]+1)*passo #Adotou-se o segundo menor MSE do SVR, pois o menor foi o referente ao primeiro marcador SNP2
abline(v=corte[[i]],col="black",lty=2)
dev.off()


snps_selec_corte[[i]]<-names(rank_RF[1:corte[[i]]])
snps_selec_corte[[i]]

#Construindo o dataframe com rank_RF juntamente com o p_valor dos SNPs selecionados na etapa de corte
rank_global_selec<-list()
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]][c(snps_selec_corte[[i]],'fenotipo')])
rank_global_selec[[i]]<-cbind(rank_RF[1:corte[[i]]],valor_p)
#View(rank_global_selec[[i]])

#Construindo o dataframe com rank_RF juntamente com o p_valor de todos os SNPs da base de dados inicial.
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]])
rank_global<-cbind(importance(RF)[,1],valor_p)
#View(rank_global)

#Selecionando somente o genotipo apos o corte
genotipo[[2]]<-dados[[1]][,names(rank_RF[1:corte[[i]]])]

#Selecionando o genotipo e o fenotipo apos o corte
dados[[2]]<-cbind(genotipo[[2]],dados[[1]]$fenotipo)
names(dados[[2]])[ncol(dados[[2]])]<-"fenotipo"


###AG para a SEGUNDA selecao dos marcadores na base de dados inicial (REFINAMENTO)

#Declara o objeto GA como uma lista

f<-function(x)
{
  inc<-which(x==1)
  
  cat("inc = ",inc, "\n")
  
  dados_validacao<- dados[[2]][,c(inc,ncol(dados[[2]]))]
  model<-validacao_cruzada_bessel(dados_validacao,
                                  folds,
                                  sigma,
                                  order,
                                  degree,
                                  cost,
                                  epsilon);
  media<-model[5];
  return(media);       
}
fitness<-function(x) {f(x)}


#Parametros do GA
run = 30            #Numero maximo de geracoes do GA sem melhoria.
maxiter = 10        #Numero maximo de geracoeso GA.
pcross = 0.8        #Probabilidade de crossing over
pmut = 0.1          #Probabilidade de mutacao.
elitism = 5         #Numero de melhores individuos do GA que irao para proxima geracaos em alteracao alguma.
popSize=100         #Tamanho da populacao em cada geracao do GA.

GA[[i]]<-ga(type="binary",
            fitness=fitness, 
            nBits=ncol(genotipo[[2]]),
            popSize=popSize,
            names=colnames(genotipo[[2]]),
            maxiter=maxiter,
            seed=i,
            parallel=TRUE,
            run=run, 
            suggestions=matrix(rep(1,ncol(genotipo[[2]])),
                               ncol=ncol(genotipo[[2]])))

snps_selec_ref[[i]]<-names(GA[[i]]@solution[,which(GA[[i]]@solution==1)])

snps_selec_SVR_bessel_1_GA_corr<-snps_selec_ref[[i]]
snps_selec_SVR_bessel_1_GA_corr

#Relatorio do AG
summary(GA[[i]])

#Grafico do AG
plot(GA[[i]])

#Salvando o grafico do AG com ggplot2
pdf(file="Grafico_GA_SVR_bessel_5.pdf",height=5,width=9)
geracao<-seq(1,GA[[1]]@iter,by=1)
mean_fitness   <-  GA[[i]]@summary[,2]
median_fitness <- GA[[i]]@summary[,4]
best_fitness   <- GA[[i]]@summary[,1]
Estatisticas <- c(rep("Mediana", length(geracao)), rep("Media", length(geracao)),rep("Melhor", length(geracao)) )
data_grafico_1 <- data.frame(
  Geracao = rep(geracao, 3),
  Aptidao = c(mean_fitness, median_fitness, best_fitness),
  Estatisticas = Estatisticas
)

ggplot(data_grafico_1, aes(x = Geracao, y = Aptidao, group = Estatisticas)) +
  geom_line(aes(colour = Estatisticas, linetype = Estatisticas),size=2) + 
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Geracao), max(data_grafico_1$Geracao), by = 1))
dev.off()


############ Termino do SMS #############


####Uniao dos SNPs selecionados####
#O argumento lista_SNPs tem que estar no formato lista com pelo menos 2 elementos.
uniao_snps<-function(lista_SNPs){
  uniao <- list()
  uniao[[1]]<-union(lista_SNPs[[1]],lista_SNPs[[2]])
  for (a in 1:(length(lista_SNPs)-2)){
    uniao[[a+1]]<-union(uniao[[a]],lista_SNPs[[a+2]])
  }
  uniao_final <- uniao[[a+1]]
  return(uniao_final)
}

uniao_snps(snps_selec_ref)

uniao_final<-uniao_snps(snps_selec_ref)

####Intersecao dos SNPs selecionados####
#O argumento lista_SNPs tem que estar no formato lista com pelo menos 2 elementos.
intersecao_snps<-function(lista_SNPs){
  intersecao <- list()
  intersecao[[1]]<-intersect(lista_SNPs[[1]],lista_SNPs[[2]])
  for (a in 1:(length(lista_SNPs)-2)){
    intersecao[[a+1]]<-intersect(intersecao[[a]],lista_SNPs[[a+2]])
  }
  intersecao_final <- intersecao[[a+1]]
  if (length(intersecao_final) == 0L) {
    return("Conjunto vazio")
  } else
    return(intersecao_final)
}

intersecao_snps(snps_selec_ref)

intersecao_final<-intersecao_snps(snps_selec_ref)


#####Valor-p bruto e ajustado - Selecao#####
#Calculo do valor p bruto e ajustado pela correcao de Bonferroni
valor.p<-function(genotipo_fenotipo)
{
  valor_p_bruto<-vector();
  valor_p_ajustado<-vector();
  y<-genotipo_fenotipo[,ncol(genotipo_fenotipo)];
  m<-ncol(genotipo_fenotipo)-1;
  model_regression<-list();
  saida<-data.frame();
  
  for (i in 1:(ncol(genotipo_fenotipo)-1))
  {
    x<-genotipo_fenotipo[,i];
    model_regression[[i]]<-lm(y~x,data=genotipo_fenotipo);
    valor_p_bruto[i]<-ifelse(all(x==3)|all(x==2)|all(x==1),1,summary(model_regression[[i]])[[4]][2,4]);
    
  }
  valor_p_ajustado<- m*valor_p_bruto;  
  valor_p<-data.frame();
  valor_p<-cbind(valor_p_bruto,valor_p_ajustado);
  colnames(valor_p)<-c("Valor p bruto","Valor p ajustado")
  rownames(valor_p)<-names(genotipo_fenotipo)[1:(ncol(genotipo_fenotipo)-1)];
  return(valor_p)
}

valor_p<-as.data.frame(valor.p(dados[[1]]))
valor_p_order<-valor_p[order(valor_p[,1], decreasing = FALSE),]
#View(valor_p_order)
valor_p_bruto_order <- valor_p_order[valor_p_order$`Valor p bruto`<= 0.05,] # Filtrando somente os SNPs com valor-p bruto <= 0.05
valor_p_corrigido_order <- valor_p_order[valor_p_order$`Valor p ajustado`<= 0.05,] # Filtrando somente os SNPs com valor-p corrigido <= 0.05
valor_p_bruto_selecao<-rownames(valor_p_bruto_order) # SNPs selecionados pelo valor-p bruto <= 0.05
valor_p_corrigido_selecao<-rownames(valor_p_corrigido_order) # SNPs selecionados pelo valor-p bruto <= 0.05

sink("my_output_SMS_Rank_Corr_6.txt")

cat('Seleção valor-p bruto =', valor_p_bruto_selecao, '\n') # Exibe os SNPs selecionados pelo valor-p bruto <= 0.05
cat('Seleção valor-p corrigido =', valor_p_corrigido_selecao, '\n')
rank_RF_df = as.data.frame(rank_RF)
cat('Seleção rank da Random Forest', View(rank_RF_df))

sink()


rank_RF_df = as.data.frame(rank_RF)

write.csv(rank_RF_df, "rank_Random_Forest_6.txt", row.names = TRUE)

write.csv(as.data.frame(rank_global_selec[[i]]),"rank_global_6.txt", row.names = TRUE )

#write this string to file
#define file name
sink("my_output_SMS_Bessel_Corr_6.txt")

#Imprimindo as informacoes finais para o texto

cat('SMS Linear =',snps_selec_ref[[1]], "\n")

cat('SMS Bessel 1 =',snps_selec_ref[[2]], "\n")

cat('SMS Bessel 2 =',snps_selec_ref[[3]], "\n")

cat('SMS Bessel 3 =',snps_selec_ref[[4]], "\n")

cat('SMS Bessel 4 =',snps_selec_ref[[5]], "\n")

cat('SMS Bessel 5 =',snps_selec_ref[[6]], "\n")

cat('Uniao =',uniao_final, "\n")

cat('Intersecao =',intersecao_final, "\n")

cat('Valor-p bruto =',valor_p_bruto_selecao, "\n")

cat('Valor-p corrigido =',valor_p_corrigido_selecao, "\n")

################# Medidas de Otimalidade dos SNPs causais###################

var_sel_causais<-c("SNP1","SNP10","SNP20","SNP30","SNP40","SNP50","SNP60")

# Parâmetros do SVR Kernel Linear
folds = 10
gamma = 0.001
cost = 1.0
epsilon = 0.1
kernel = 'linear'
resultado_causais_linear<-validacao_cruzada(data = dados[[1]][c(var_sel_causais,"fenotipo")],
                                            folds = folds, 
                                            gamma = gamma,
                                            cost = cost,
                                            epsilon = epsilon,
                                            kernel = kernel)


# Parâmetros do SVR Kernel Bessel 1
folds = 10
sigma = 0.10
order = 1.00
degree = 1
cost = 1.0
epsilon = 0.1

resultado_causais_bessel_1<-validacao_cruzada_bessel(data = dados[[1]][c(var_sel_causais,"fenotipo")],
                                                     folds = folds, 
                                                     sigma = sigma,
                                                     order = order,
                                                     degree = degree,
                                                     cost = cost,
                                                     epsilon = epsilon)

# Parâmetros do SVR Kernel Bessel 2
folds = 10
sigma = 0.10
order = 1.50
degree = 1
cost = 1.0
epsilon = 0.1

resultado_causais_bessel_2<-validacao_cruzada_bessel(data = dados[[1]][c(var_sel_causais,"fenotipo")],
                                                     folds = folds, 
                                                     sigma = sigma,
                                                     order = order,
                                                     degree = degree,
                                                     cost = cost,
                                                     epsilon = epsilon)


# Parâmetros do SVR Kernel Bessel 3
folds = 10
sigma = 0.30
order = 1.00
degree = 1
cost = 1.0
epsilon = 0.1

resultado_causais_bessel_3<-validacao_cruzada_bessel(data = dados[[1]][c(var_sel_causais,"fenotipo")],
                                                     folds = folds, 
                                                     sigma = sigma,
                                                     order = order,
                                                     degree = degree,
                                                     cost = cost,
                                                     epsilon = epsilon)


# Parâmetros do SVR Kernel Bessel 4 
folds = 10
sigma = 0.30
order = 1.50
degree = 1
cost = 1.0
epsilon = 0.1

resultado_causais_bessel_4<-validacao_cruzada_bessel(data = dados[[1]][c(var_sel_causais,"fenotipo")],
                                                     folds = folds, 
                                                     sigma = sigma,
                                                     order = order,
                                                     degree = degree,
                                                     cost = cost,
                                                     epsilon = epsilon)

# Parâmetros do SVR Kernel Bessel 5 
folds = 10
sigma = 0.30
order = 1.50
degree = 1
cost = 1.0
epsilon = 0.1

resultado_causais_bessel_5<-validacao_cruzada_bessel(data = dados[[1]][c(var_sel_causais,"fenotipo")],
                                                     folds = folds, 
                                                     sigma = sigma,
                                                     order = order,
                                                     degree = degree,
                                                     cost = cost,
                                                     epsilon = epsilon)
###########################################################################

#Saida do SVR kernel linear com somente SNPs causais
cat('Correlacao do SVR Linear = ',resultado_causais_linear[5], "\n")
cat('R2 ajustado do SVR Linear = ',resultado_causais_linear[9], "\n")
cat('MSE do SVR Linear = ',resultado_causais_linear[1], "\n")

#Saida do SVR kernel radia gamma = 0.001 com somente SNPs causais
cat('Correlacao do SVR Bessel 1 = ',resultado_causais_bessel_1[5], "\n")
cat('R2 ajustado do SVR Bessel 1 = ',resultado_causais_bessel_1[9], "\n")
cat('MSE do SVR do SVR Bessel 1 = ',resultado_causais_bessel_1[1], "\n")

#Saida do SVR kernel radia gamma = 0.01 com somente SNPs causais
cat('Correlacao do SVR Bessel 2 = ',resultado_causais_bessel_2[5], "\n")
cat('R2 ajustado do SVR Bessel 2 = ',resultado_causais_bessel_2[9], "\n")
cat('MSE do SVR Bessel 2 = ',resultado_causais_bessel_2[1], "\n")

#Saida do SVR kernel radia gamma = 0.1 com somente SNPs causais
cat('Correlacao do SVR Bessel 3 = ',resultado_causais_bessel_3[5], "\n")
cat('R2 ajustado do SVR Bessel 3 = ',resultado_causais_bessel_3[9], "\n")
cat('MSE do SVR Bessel 3 = ',resultado_causais_bessel_3[1], "\n")

#Saida do SVR kernel radia gamma = 1 com somente SNPs causais
cat('Correlacao do SVR Bessel 4 = ',resultado_causais_bessel_4[5], "\n")
cat('R2 ajustado do SVR Bessel 4 = ',resultado_causais_bessel_4[9], "\n")
cat('MSE do SVR Bessel 4 = ',resultado_causais_bessel_4[1], "\n")

#Saida do SVR kernel radia gamma = 1 com somente SNPs causais
cat('Correlacao do SVR Bessel 5 = ',resultado_causais_bessel_5[5], "\n")
cat('R2 ajustado do SVR Bessel 5 = ',resultado_causais_bessel_5[9], "\n")
cat('MSE do SVR Bessel 5 = ',resultado_causais_bessel_5[1], "\n")

#close the external connection
sink()

save.image(file="SMS_corr_6.RData")


################# Parâmetros escolhidos para o kernel Bessel ###############













#bessel <- besseldot(sigma = 1, order = 1, degree = 2 ) # Estabelecer os parâmetros do kernel Bessel.
#bessel_1 <- ksvm(fenotipo~., data = data_temp) # Treinar o modelo SVM com Kernel Bessel na base de dados de SNPs e Fenótipo.
#bessel_1 # Mostrar o modelo SVM treinado.

folds = 10
sigma = 0.10
order = 1.00
degree = 1
cost = 1.0
epsilon = 0.1


svr_test_bessel_1 <- validacao_cruzada_bessel(data = dados[[1]],
                                              folds = folds, 
                                              sigma = sigma,
                                              order = order,
                                              degree = degree,
                                              cost = cost,
                                              epsilon = epsilon)

svr_test_bessel_1[5]

folds = 10
sigma = 0.10
order = 1.50
degree = 1
cost = 1.0
epsilon = 0.1

svr_test_bessel_2 <- validacao_cruzada_bessel(data = dados[[1]],
                                              folds = folds, 
                                              sigma = sigma,
                                              order = order,
                                              degree = degree,
                                              cost = cost,
                                              epsilon = epsilon)

svr_test_bessel_2[5]


folds = 10
sigma = 0.30
order = 1.00
degree = 1
cost = 1.0
epsilon = 0.1

svr_test_bessel_3 <- validacao_cruzada_bessel(data = dados[[1]],
                                              folds = folds, 
                                              sigma = sigma,
                                              order = order,
                                              degree = degree,
                                              cost = cost,
                                              epsilon = epsilon)

svr_test_bessel_3[5]



folds = 10
sigma = 0.30
order = 1.50
degree = 1
cost = 1.0
epsilon = 0.1

svr_test_bessel_4 <- validacao_cruzada_bessel(data = dados[[1]],
                                              folds = folds, 
                                              sigma = sigma,
                                              order = order,
                                              degree = degree,
                                              cost = cost,
                                              epsilon = epsilon)

svr_test_bessel_4[5]


folds = 10
sigma = 0.50
order = 1.80
degree = 1
cost = 1.0
epsilon = 0.1

svr_test_bessel_5 <- validacao_cruzada_bessel(data = dados[[1]],
                                              folds = folds, 
                                              sigma = sigma,
                                              order = order,
                                              degree = degree,
                                              cost = cost,
                                              epsilon = epsilon)

svr_test_bessel_5[5]
