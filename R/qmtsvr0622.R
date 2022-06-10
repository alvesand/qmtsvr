#--------------------------------------------------------------------#
# (Quasi) Multi-task Support Vector Regression with GA optimization
#--------------------------------------------------------------------#

#--------------------------------------------------------------------#
#----------------- Weighted Euclidean Distances ---------------------#
#--------------------------------------------------------------------#
qmtsvr.dist<-function (x, w, verbose = verbose, u, scale, vardiag){ 
  start.time <- Sys.time()
  if(missing(verbose)==T){verbose = F}
  if(missing(scale)==T){scale = F}
  if(missing(u)==T){u = 1}
  if(missing(vardiag)==T){vardiag = F}
  if(typeof(vardiag)!="logical"){stop("Error ... vardiag must be T or F")}
  if(typeof(scale)!="logical"){stop("Error ... scale must be T or F")}
  if(typeof(u)!="double"){stop("Error ... u parameter must be numeric")}
  if(verbose == T){welcome()}
  
  #------------- CREATING MML ---------------------------------------#
  if (missing(w)==T){
    MML<-tcrossprod(x)
    len = dim(x)[2]
  }else{
    n = dim(x)[1]
    p = dim(x)[2]
    len = sum(w)
    Xw = matrix(NA, n,p)
    for (k in 1:p){
      Xw[,k] <-  x[,k]*w[k]
      
    }
    
    MML<-tcrossprod(Xw, x)}
  
  ## DIAGONAL OF MML ###
  S<-diag(MML)
  
  ### CREATING THE EUCLIDEAN DISTANCE MATRIX (Winkelman et al. 2015) ###
  EDM<-matrix(NA,nrow=nrow(MML),ncol=ncol(MML))
  dimnames(EDM)[[1]] <- rownames(MML)
  dimnames(EDM)[[2]] <- colnames(MML)
  for(i in 1:dim(MML)[1]){
    for(j in 1:dim(MML)[2]){
      EDM[i,j] <- sqrt(S[i]+S[j]-2*MML[i,j])
      #print(i)
    }
  }
  
  
  if(scale==T){EDM=(EDM/len^(1/u))^u}else{EDM = EDM^u}
  if(vardiag==T){diag(EDM) = 1-(S/(sqrt(S)*sqrt(median(S))))}
  
  
  end.time <- Sys.time()
  if(verbose == T){
    mytime = as.character(round(difftime(end.time, start.time, unit = "mins"),2))
    mytime = paste(mytime, "min", sep = " ")
    star = paste(rep("*",56),collapse = "")
    pad = paste(rep(" ",9), collapse="")
    msg1 = "Euclidean distance matrix ready"
    msg2 = paste(pad, 'Elapsed time: ', collapse="")
    msg2 = paste(msg2, mytime, collapse="")
    n_elem = nchar(msg2)
    npad2 = 56 - (n_elem) - 5
    pad2 = paste(rep(" ",npad2), collapse="")
    cat(star, '\n')
    cat('*',pad,msg1,pad,' *', '\n')
    cat('*',msg2,pad2,'*', '\n')
    cat(star, '\n')
    }                                       
  return(EDM)
}

#--------------------- End of subroutine ----------------------------------------------------#
#--------------------------------------------------------------------------------------------#
#qmtsvr.fit: This function fits a (Quasi) Multitask SVR with user pre-defined hyper parameters
#requires the package kernlab
#Y: is an n x t matrix of correlated observations with the target trait on the 1st column
#X: is an n x p matrix of predictor variables trait-specific and common to all observations
#set_hyper: a list with all required hyper parameters 
#D: a list with all pre-computed EDM  can be passed as an argument alternatively to the matrix X 
#-------------------------------------------------------------------------------------------#
qmtsvr.fit = function(Y, X, w, set_hyper, D, verbose, vardiag){
  
  #--- Requirement checking! ---#
  caution()
  #-----------------------------#
  
  #--------- Normalizes y to lie between 0 and 1 ----------------------#
  normalize <- function(x) {
    num <- x - min(x, na.rm = T)
    denom <- max(x, na.rm = T) - min(x, na.rm = T)
    return (num/denom)
  }
  
  y = NULL
  ntrait = ncol(Y)
  N = nrow(Y)
  
  if (missing(verbose)==T){verbose=F}
  if(missing(vardiag)==T){vardiag = F}
  if(typeof(vardiag)!="logical"){stop("Error ... vardiag must be T or F")}
  if(verbose == T){welcome()}
  
  if(missing(X)==F){p = ncol(X)}
  if(ntrait>1){
    nb = ntrait*(ntrait+1)/2
    nr = ntrait*(ntrait+1)/2 - ntrait
    band_id = paste("b",seq(1,nb,1), sep = "")
    r_id = paste("r",seq(1,nr,1), sep = "")
    hyper_dic = c("C","eps",band_id,r_id)}else{hyper_dic = c("C","eps","b1");nb=1}
  
  #Checking the sequence of provided hyperparameters
  for (i in 1:length(hyper_dic)){
    if(names(set_hyper)[i]!=hyper_dic[i]){
      stop("Invalid hyperparameters list provided")}}
  
  for (i in 1:ntrait){y = c(y, normalize(Y[,i]))}
  bandwidth = set_hyper[3:(2+nb)]
  
  if (missing(D)==T){
    if(verbose==T){
      cat('#--| Euclidean Distance Matrices (EDM) not provided |--#', '\n')
      cat('Computing EDMs...', '\n')}
    if(missing(w)==T){
      D = list();
      Dtmp = qmtsvr.dist(x=X, verbose = F, u =2, scale = T, vardiag = vardiag)
      D[[1]] = Dtmp; rm(Dtmp)
      if(verbose==T){
        cat('Done...', '\n')}}
    else{
      D = list();
      for (i in 1:length(w)){Dtmp = qmtsvr.dist(x=X, w = w[[i]], verbose = F, u =2, scale = T, vardiag=vardiag)
      D[[i]] = Dtmp; rm(Dtmp)}
      if(verbose == T){
        cat('Done...', '\n')}
    }
  }
  
  #Building the kernel blocks
  if (verbose==T){
    cat('Building the kernel blocks...', '\n')}
  nD = length(D)
  if (nD>1&nD!=nb){stop("Error! Number of EMDs do not match number of bandwidth parameters")}
  K = list()
  if (nD == 1){
    for (i in 1:nb){
      K[[i]] = exp(-1*bandwidth[[i]]*D[[1]])}}
  else{
    for (i in 1:nb){
      K[[i]] = exp(-1*bandwidth[[i]]*D[[i]])}}
  
  #Getting the indexes for building Q
  ind = matrix(0, nrow = ntrait,ntrait)
  j=1
  first = 1
  last = ntrait
  for (i in 1:nrow(ind)){
    doseq = seq(first,last,1)
    ind[i,j:ncol(ind)] = doseq
    first = seq(first,last,1)[length(seq(first,last,1))]+1
    j=j+1
    last = first + ncol(ind) - j}
  
  index = c(ind[upper.tri(ind, diag = F)])
  ind2 = t(ind)
  diag(ind2) = 0
  ind = ind + ind2
  rm(ind2)
  
  if(ntrait>1){
    phi = set_hyper[(nb+3):(nb+2+nr)]
    for (i in 1:length(index)){
      K[[index[i]]] = phi[[i]]*K[[index[i]]]}
  }
  
  L = list()
  Krow = NULL
  for (i in 1:nrow(ind)){
    ki = ind[i,]
    for (o in 1:length(ki)){
      v = ki[o]
      Krow = cbind(Krow,K[[v]])}
    L[[i]] = Krow
    Krow = NULL
  }
  
  Q = NULL
  for (i in 1:length(L)){Q = rbind(Q,L[[i]])}
  if (verbose == T){cat('Done...', '\n')}
  mod = ksvm(Q, y, kernel = 'matrix', 
             type = "eps-svr", C = set_hyper[[1]], e = set_hyper[[2]])
  
  yhat = predict(mod)
  YHAT = matrix(0,nrow = N, ncol = ntrait)
  ni = 0
  nf = 0
  for (i in 1:ntrait){
    max_train = max(Y[,i], na.rm = T)
    min_train = min(Y[,i], na.rm = T)
    yhat1 = (max_train - min_train)*yhat[(1+ni):(nf+N)] + min_train
    ni = ni+N
    nf = nf+N
    YHAT[,i] = yhat1
  }
  return(YHAT)
}

#--------------------------------------------------------------------#
#----------------- GA for the QMTSVR optimization -------------------#
#--------------------------------------------------------------------#
qmtsvr.GA = function(Y, X, w = w, hyper, ngen, popsize, mut_rate, cross_rate, elitism, vartype = vartype,
                     verbose, cost, nfolds, val_pop, tsize, custom_val, k, MRCR, vardiag, lambda, dopar){
  if(missing(dopar)==T){dopar=F}
  if(typeof(dopar)!="logical"){stop("Error ... dopar parameter must be T or F")}
  if(dopar == T){
    require("foreach")
    require("parallel")
    require("doParallel")
  }
  require ("kernlab")
  
  Y = Y
  X = X
  N = nrow(X)
  nsnp = ncol(X)
  if(missing(w)==T){w=NULL}
  if(missing(verbose)==T){verbose=NULL}
  if(missing(nfolds)==T){nfolds=NULL}
  if(missing(tsize)==T){tsize = 5}
  if(missing(lambda)==T){lambda = 0}
  if(missing(elitism)==T){elitism = 1}
  if(missing(vartype)==T){vartype = "continuous"}
  if(vartype!="continuous"&vartype!="binary"){stop("Error: vartype must be 'continuous' or 'binary'...")}
  if(missing(MRCR)==T){MRCR="fixed"}
  if(MRCR!="fixed"&MRCR!="dynamic"){stop("Error MCRC must be 'fixed' or 'dynamic'...")}
  if(missing(vardiag)==T){vardiag = F}
  if(typeof(vardiag)!="logical"){stop("Error ... vardiag must be T or F")}
  
  #*******************************************************************#
  #*                  NESTED FUNCTIONS                               *#
  #*******************************************************************#
  #1.cost metrics 
  #Classification
  class_metrics = function(df, t){
    if(missing(t)==T){t=0.5}
    db = data.frame(df)
    db[,2] = ifelse(db[,2]>=t,1,0)
    colnames(db) = c("Target", "Classification")
    target = factor(db$Target,levels = c(0,1))
    classification = factor(db$Classification, levels = c(0,1))
    N = nrow(db)
    CM = table(target, classification)
    accuracy = (CM[2,2] + CM[1,1]) / N
    precision = (CM[2,2] / (CM[2,2] + CM[1,2]))
    sensitivity = CM[2,2] / (CM[2,2] + CM[2,1])
    specificity = CM[1,1] / (CM[1,1] + CM[1,2])
    f1 = (2*precision*sensitivity) / (precision + sensitivity)
    if(is.nan(precision)==T){precision = 0}
    if(is.nan(f1)==T){f1 = 0}
    
    for(l in 1:N){
      df[l,2] = max(0,df[l,2])
    }
    cross_entropy = 1/N*sum(-1*(df[,1]*log(df[,2])+(1-df[,1])*log(1-df[,2])))
    
    return(list(CM = CM, acc = accuracy, cross_entropy = cross_entropy, precision = precision, 
                sensitivity = sensitivity, specificity = specificity, f1 = f1))
    
  }
  #Regression
  reg_metrics = function(y, yhat){
    acc = cor(yhat, y)
    rmse = sqrt(mean((yhat - y)^2))
    mae = mean(abs(yhat - y))
    return(list(cor = acc, rmse = rmse, mae = mae))
  }
  
  #Convert a binary array to an integer value
  BinToDec <- function(x) 
    sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))+1
  #*************************************************************************#
  
  popsize = popsize
  pop_average = NULL
  best_average = NULL
  pop = NULL
  gen = 0
  hyper = hyper
  cost = cost
  if(missing(nfolds)==T){nfolds=3}
  if(missing(custom_val)==T&val_pop=="custom"){stop("Error... Please provide a valid array for custom validation pop!")}
  if(vartype=="continuous"&cost!="rmse"&cost!="cor"&cost!="mae"){stop("Error... Please provide a valid cost character for continuous target variable")}
  if(vartype=="binary"&cost!="accuracy"&cost!="cross_entropy"&cost!="sen"&cost!="spec"
     &cost!="f1"&cost!="precision"){stop("Error... Please provide a valid cost character for binary target variable")}
  if(cost =="cor"|cost == "precision"| cost == "sen" | cost == "spec" | cost == "f1"| cost == "accuracy"){argmax = T}else{argmax=F} 
  if (verbose == T){welcome()}
  
  #----------------------------------------------------------------------------------------#
  #-------------------------   Data pre-processing   --------------------------------------#
  #----------------------------------------------------------------------------------------#
  ntrait = ncol(Y)
  
  
  #--------------------------------------------------------
  #Building the hyperparameter space
  if(ntrait>1){
    nb = ntrait*(ntrait+1)/2
    nr = ntrait*(ntrait+1)/2 - ntrait
    band_id = paste("b",seq(1,nb,1), sep = "")
    r_id = paste("r",seq(1,nr,1), sep = "")
    hyper_dic = c("C","eps",band_id,r_id)}else{hyper_dic = c("C","eps","b1");nb=1}
  
  hyperspace = list()
  for (i in 1:length(hyper_dic)){
    par = hyper[[i]]
    if(sum(as.numeric(par[1] == hyper_dic))==1){
      start = as.numeric(par[2])
      end = as.numeric(par[3])
      interval = (end - start)/(as.numeric(par[4])-1)
      hyperspace[[i]] = c(seq(start,end,interval))
    }
  }
  bin_size = 0
  for (i in 1:length(hyperspace)){bin_size = bin_size+log(length(hyperspace[[i]]),2)}
  n_par = length(hyper_dic)
  if(verbose==T){
    cat('Computing EDMs...', '\n')
    ini = Sys.time()}
    
  #Building the EDMs
  if(length(w)==0){
    D = list()
    D[[1]] = qmtsvr.dist(X, verbose = F, scale = T, u = 2, vardiag = vardiag); rm(X)
  }else{
    D = list()
    for(j in 1:nb){
      D[[j]] = qmtsvr.dist(X,w[[j]],verbose = F, scale = T, u = 2, vardiag=vardiag)};rm(X)}
  if(verbose==T){
    fim = Sys.time()
    cat('Done...', '\n')
    a=utils::capture.output(fim-ini)
    cat('Time elapsed:', substr(a, 19, 100), '\n')}
  #---------------------------- End of subroutine ----------------------------------------#
  
  #----------------------------------------------------------------------------------------#
  #-This function computes the cost function for the GA -----------------------------------#
  #----------------------------------------------------------------------------------------#
  
  fun_fit = function(x){
    xp = x
    genes = list()
    for(i in 1:length(hyper_dic)){
      assign(hyper_dic[i],hyperspace[[i]])
      slice = log(length(get(hyper_dic[i])),2)
      genes[[i]] = assign(paste("x",i,sep=""),xp[1:slice])
      xp = xp[slice+1:length(xp)]}
    rm(xp)
    
    set_hyper = list()
    for (i in 1:length(hyper_dic)){
      set_hyper[[i]] = get(hyper_dic[i])[BinToDec(genes[[i]])]
    }
    names(set_hyper) = hyper_dic
    if(gen>0){penalty = sum(x!=pop[best,])}else{penalty=0}
    if(val_pop!="cross"){
      y = Y
      holdout=1
      y[which(folds==holdout),1] = NA
      YHAT = qmtsvr.fit(Y=y, D = D, set_hyper = set_hyper, verbose = F)
      if(vartype=="continuous"){
        metrics = reg_metrics(y=Y[which(folds==holdout),1], yhat=YHAT[which(folds==holdout),1])
        if(cost=="cor"){return(list(metrics$cor-lambda*exp(-penalty),set_hyper))}
        else if(cost=="rmse"){return(list(metrics$rmse+lambda*exp(-penalty),set_hyper))}
        else{return(list(metrics$mae+lambda*exp(-penalty),set_hyper))}
      }else {
        metrics = class_metrics(df = cbind(Y[which(folds==holdout),1], YHAT[which(folds==holdout),1]))
        if (cost == "accuracy"){return(list(metrics$acc+lambda*exp(-penalty),set_hyper))}
        else if(cost=="cross_entropy"){return(list(metrics$cross_entropy+lambda*exp(-penalty),set_hyper))}
        else if (cost == "sen"){return(list(metrics$sensitivity+lambda*exp(-penalty),set_hyper))}
        else if (cost == "spec"){return(list(metrics$specificity+lambda*exp(-penalty),set_hyper))}
        else if (cost == "precision"){return(list(metrics$precision+lambda*exp(-penalty),set_hyper))}
        else if (cost == "f1"){return(list(metrics$f1+lambda*exp(-penalty),set_hyper))}
      }
    }
    else{
      RESUL = matrix(NA,nfolds, 3)
      for (m in 1:nfolds){
        y = Y
        holdout = m
        y[which(folds==holdout),1] = NA
        YHAT = qmtsvr.fit(Y=y, D = D, set_hyper = set_hyper, verbose = F)
        if(vartype=="continuous"){
          metrics = reg_metrics(y=Y[which(folds==holdout),1], yhat=YHAT[which(folds==holdout),1])
          if(cost=="cor"){RESUL[m,1] = metrics$cor-lambda*exp(-penalty)}
          else if(cost=="rmse"){RESUL[m,1]=metrics$rmse+lambda*exp(-penalty)}
          else {RESUL[m,1]=metrics$mae+lambda*exp(-penalty)}
        }
        else{
          metrics = class_metrics(df = cbind(Y[which(folds==holdout),1],
                                             YHAT[which(folds==holdout),1]))
          if (cost == "accuracy"){RESUL[m,1]=metrics$acc+lambda*exp(-penalty)}
          else if(cost=="cross_entropy"){RESUL[m,1]=metrics$cross_entropy+lambda*exp(-penalty)}
          else if (cost == "sen"){RESUL[m,1]=metrics$sensitivity+lambda*exp(-penalty)}
          else if (cost == "spec"){RESUL[m,1]=metrics$specificity+lambda*exp(-penalty)}
          else if (cost == "precision"){RESUL[m,1]=metrics$precision+lambda*exp(-penalty)}
          else if (cost == "f1"){RESUL[m,1]=metrics$f1+lambda*exp(-penalty)}
        }
      }  
    return(list(mean(RESUL[,1]),set_hyper))}
  }
  
  #--------------------- End of subroutine ------------------------------------------------------#
  
  #-------------- The crossing-over function --------------------------#
  cross_over = function(population, mut_rate, cross_rate, tsize, elitism, score){
    popsize =  nrow(population)
    nchildren = popsize - elitism
    chr_len = dim(population)[2]
    children = NULL
    for (i in 1:nchildren){
      #Tournament selection 
      index1 = sample(popsize, tsize, replace = F)
      index2 = sample(popsize, tsize, replace = F)
      while(sum(index1%in%index2)>0){index2 = sample(popsize, tsize, replace = F)}
      candidates1 = population[index1,]
      candidates2 = population[index2,]
      score1 = c(score[index1])
      score2 = c(score[index2])
      best1 =  whichpart(x=score1, n = 1, argmax = argmax)
      best2 =  whichpart(x=score2, n = 1, argmax = argmax)
      parent1 = candidates1[best1,]
      parent2 = candidates2[best2,]
      #one-point crossing over
      point = sample(2:(chr_len-1),1)
      do_cross = rbinom(1,1,cross_rate)
      if(do_cross==1){
        child = c(rep(0,chr_len))
        child[1:point] = parent1[1:point]
        child[(point+1):chr_len] = parent2[(point+1):chr_len]
      }else{child = parent1}
      #Mutations
      z = 1
      for (z in 1:chr_len){
        child[z] = abs(rbinom(1, prob = mut_rate, size = 1) - child[z])}
      children = rbind(children,child)}
    elit =  whichpart(score, n = elitism, argmax = argmax)
    children = rbind(children,population[elit,])
    return (children)}
  #--------------------- End of subroutine ----------------------------#
  
  #-------------------- Selects the best n individuals -------------------------------------#
  whichpart <- function(x, n, argmax) {
    ind = order(x, decreasing = argmax)[1:n]
  }
  
  #------------ Create the population for generation 0 ----------------#
  for (i in 1:popsize){
    pop = rbind(pop,matrix(rbinom(n=bin_size, prob = 0.5, size = 1),1,bin_size))}
  
  if(val_pop=="cross"){
    folds = sample(1:nfolds,N,replace = T)
    folds[which(is.na(Y[,1]))] = nfolds+1
  }else if(val_pop=="custom"){
    folds = custom_val; holdout = 1}
  else{
    target = which(is.na(Y[,1]))
    train_ind = NULL
    for (z in 1:length(target)){
      ind = target[z]
      notuse = unique(c(-ind,-target))
      nearest = order(D[[1]][ind,notuse])[1:k]
      train_ind = c(train_ind,nearest)
    }
    train_ind = unique(train_ind)
    cp = length(train_ind)
    if(verbose==T){cat("\nUsing", cp, "unique nearest points as target training...\n")}
    
    folds = rep(2,N)
    folds[train_ind] = 1
    holdout = 1
  }
  
  #------------ Begin looping for the n generations -------------------#
  if (verbose == T){
    ini = Sys.time()
    cat('#------------ Starting GA generations!-----------------#')
    cat('\n')} 
  
  if (dopar == T){
    n.cores <- parallel::detectCores() - 1
    my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "FORK")
    #Register the cluster
    doParallel::registerDoParallel(cl = my.cluster)
    if (verbose == T){cat("\n", "Number of registered cores:", foreach::getDoParWorkers(), "\n" )}
  }
  
  MR=mut_rate
  CR=cross_rate
  while (gen <= ngen){
    tryCatch({
      score = NULL
      hypercomb = list()
      
      
      if(dopar == T){
        resul=NULL
        resul = foreach(z=1:popsize)%dopar%{
        fun_fit(pop[z,])}
      }else{
        resul=list()
        for (z in 1:popsize){
          resul[[z]] = fun_fit(pop[z,])}
      }
      
      for (i in 1:popsize){
        score =c(score, unlist(resul[[i]][1]))
        hypercomb[[i]] = unlist(resul[[i]][2])}
      
      best =  whichpart(x=score, n = 1, argmax = argmax)
      pop_average =c(pop_average,mean(score, na.rm = T))
      best_average = c(best_average,score[best])
      best_hyper = hypercomb[[best]]
      
      
      e = sum(best_average[length(best_average)] == best_average)
      if(gen==0){e=0}
        
      if(MRCR=="dynamic"){
        mut_rate = MR+(1.01^e-1)
        cross_rate = CR-(1.01^e-1)
        if(mut_rate>0.50){mut_rate=0.50}
        if(cross_rate<0.50){cross_rate=0.50}
      }
      
      pop=cross_over(population=pop, mut_rate=mut_rate, cross_rate=cross_rate, tsize=tsize, elitism=elitism, score=score)
      
      if (verbose == T){
        cat('\n')
        cat('Generation:', gen,  "\n")
        cat('Population average:', pop_average[gen+1], "\n")
        cat('Best performance:', best_average[gen+1], "\n")
        cat('MR:', mut_rate, " CR: ", cross_rate, "\n")
        cat('stuck:', e, " generations", "\n")
        cat('Hyperparameters:', "\n")
        print(unlist(best_hyper))}
      
      gen=gen+1
    }, interrupt = function(x) {
      message("\nGA interromped by user, closing cluster...\n")
      stopCluster(cl=my.cluster)
      stop()
    })}
  if (verbose == T){
    fim = Sys.time()
    a=utils::capture.output(fim-ini)
    cat('Time elapsed:', substr(a, 19, 100), '\n')}
  
  stopCluster(cl=my.cluster)
  return (list(set_hyper = best_hyper, actual_population = pop, population_scores = score, 
               best_index = best[1], pop_average = pop_average, best_performance = best_average))
  
}

#--------------------- End of subroutine ----------------------------------------------------#
#

#--------------------------------------------------------------------#
#-------------------- Nested functions ------------------------------#
#--------------------------------------------------------------------#

#-------------------- Welcome panel ---------------------------------#
welcome = function(){
 cat('\n
 #------- The University of Wisconsin - Madison --------#
 #------ Department of Animal and Dairy Sciences -------#
 #                     /)  (\                            # 
 #                .-._((.~~.))_.-.                      #
 #                `-.   @@   .-.-                       #
 #                  / .o--o. \                           #
 #                 ( ( .__. ) )                         #
 #------------------------------------------------------#
 #       QMTSVR v.0.1.4 (beta) - June 2022              #
 #   (Quasi)Multitask support vector regression         #
 #  Anderson A.C. Alves (alves.zootecnista@gmail.com)   #
 #------------------------------------------------------#
      \n')}


#--------- Requirement checking!!! ---------------------------------#
caution = function(){
  suppressWarnings(if (!require("kernlab"))
  {use_rsp =  readline("kernlab package not installed, would you like to install it? (0 = 'No', 1 = 'Yes': ");
  while (use_rsp!="0"&use_rsp!="1"){use_rsp =  readline("Character invalid, please inform 0 (No) or 1 (Yes):")}
  if(use_rsp == "0"){stop("Error!!! kernalb must be installed before using qmtsvr.fit()")}
  else{install.packages("kernlab", dependencies = TRUE);require(kernlab)}})}

#--------------------- End of subroutine ----------------------------#

#-------------------#
#  END of program   #
#-------------------#
