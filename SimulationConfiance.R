#######################################################################################################
#                           ISSEA 2020-2021                                                           #               
#                                                                                                     #
# PROJET DE LOGICIELS STATISTQIUES: SIMULATION D'UN INTERVALLE DE CONFIANCE PAR REECHANTILLONAGE      # 
#                                                                                                     #
#                                                                                                     #
#                                                                                                     #
#  Redigé par : KETEGAZA PRINCE HERITIER GODBLESS               Examinateur: LUC B. DIMAI             #
#               KEUBOU KENFACK DILANE DESIRE                                                          #
#                                                                                                     #
#######################################################################################################

#fonction de calcul des mesures de pauvrétés

povAF1 = function(X,z,w,k){
  #verification de la conformité des paramètres
  #idee : faire celà dans la premiere fonction ou pas
  if(!is.data.frame(X)&&!is.matrix(X)){return(-1)}
  if(!(length(z)==ncol(X))&&!(length(w)==ncol(X))&&!(length(z)==length(w))){return(-1)}
  if(!(k>=1&& k<=d)) {return(-1)}
  d = ncol(X)
  n = nrow(X)
  g = matrix(c(rep(0,n*d)),nrow = n,byrow = FALSE)
  c1 = vector()
  pk = vector()
  ck = vector()
  #evaluation de la matrice de privation g
  temp=c(t(X))<z
  g=rep(w,n)
  g=matrix(g*temp,nrow = n,byrow = TRUE)
  #}"
  #g=matrix(if(X >=z) 0 else w,nrow = n,ncol = d,byrow = FALSE)
  #evaluation du vecteur d'intensite"
  c1 = apply(g, 1,sum)
  #definition de la fonction d'identification 
  pk=as.numeric(c1>=k)
  #evaluation de l'indice de pauvrete H
  q=length(pk[pk == 1])
  H=q/n
  #calcule du vecteur d'intensite de privations censure
  ck = c1*pk
  #calcul de la proportion d'intensite parmis les pauvres
  A=sum(ck)/(q*d)
  #calcul de la mesure de pauvrete ajustee
  M0= A*H
  #vecteur en sortie
  Res = c("A"=A,"H"=H,"M0"=M0 )
  return(Res)
}
#####################################################################################################################
####################################################################################################################???########
                   ########fonction simulation intervalle de confiance#######

simulerC = function(database,typeSondage,strate,FUN,methoSim,times=1000,er=0.05,z,w,k){
  databaseTemp=database
  X=database[,-strate]
  n=nrow(database)
  M01= vector() #vecteur resultat des variables 
  #echantillonnage
  if(methoSim  == "bootstrap" ){
      if(typeSondage=="SAT"){
    #conversion de la variable de stratification en facteur
    databaseTemp[,strate] = as.factor(databaseTemp[,strate])
    nstrate = length(levels(databaseTemp[,strate]))
    #division en sous ensemble de la base de donnees chaque element de la liste est un dataframe, associe à un ses
    SousensDataFrame=split(database,database[,strate])
    #fonction pour l'echantillonnage de niveau tour 
    echanttillonnnage=function(tour,SousensDataFrame){
      # on va parcourir chaque couche et exttraire le nombre d'individus necessaires 
      echantillonI=function(i,SousensDataFrame){
        #on va parcourir le sous enssemble i et extraire avec remise  le nombre d'echantillon
        return(sample(row.names(SousensDataFrame[[i]]),nrow(SousensDataFrame[[i]]),replace = TRUE))
      }
      # on applique maintenant l'echantillonnage à chaque sous groupe
      echantillon_Ind=unlist(lapply(X=1:nstrate,FUN = echantillonI,SousensDataFrame=SousensDataFrame))
      #fabrication de l'echantillon complet au tour tour
      echantillon= X[as.numeric(echantillon_Ind),]
      #application de la fonction mesurePauvrete
      res=FUN(echantillon,z=z,w=w,k=k)
      return(res[3])
    }
    #calcul de la grandeur de pauvrete pour chaque tour 
    M01=vapply(1:times,FUN = echanttillonnnage,SousensDataFrame=SousensDataFrame,FUN.VALUE = numeric(1))
      }
    if(typeSondage=="SAS"){
      echanttillonnnage3=function(tour){
       individu=sample(1:n,n,replace = TRUE)
       echantillon=database[individu,-10]
       return(povAF1(X=echantillon,z,w,k)[3])
      }
      M01=vapply(1:times,FUN = echanttillonnnage3,FUN.VALUE = numeric(1))
    }
    }
  # methode jacknifff
  if(methoSim == "jackknife"){
    echanttillonnnage2=function(individu){
      X1=X[-individu,]
      return(FUN(X1,z,w,k)[3])
    }
    #application de cette fonction à l'ensembe de tous les donnees 
    M01=vapply(1:nrow(X), FUN = echanttillonnnage2, FUN.VALUE = numeric(1))
  }
  #calcule de l'intervalle interquartille 
  intervalle=quantile(M01,probs=c(er/2,1-er/2))
  #valeurs retournees
  Res = list()
  Res[[1]] = database
  Res[[2]] = methoSim
  Res[[3]] = "nom de la"
  Res[[4]] = M01
  Res[[5]] = intervalle
  return(Res)
}