
BinBayes <- function(m_data,model,link)
{
  
  library(coda)
  library(lme4)
  library(rjags)
  
  n.obs<-nrow(m_data)
  n.subject<-length(unique(m_data$subj))
  n.item<-length(unique(m_data$itemID))
  n.cond<-length(unique(m_data$cond))
  
  #store the labels used for subjects, conditions, items
  labels.subject<-unique(m_data$subj)
  labels.cond<-unique(m_data$cond)
  labels.item<-unique(m_data$itemID)
  #setup m_data vectors as needed for jags
  subject<-match(m_data$subj,labels.subject)
  item<-match(m_data$itemID,labels.item)
  condition<-match(m_data$cond,labels.cond)
  y<-as.numeric(m_data$Acc)
  
  
  n.iter.marginal<-250000
  n.thin.marginal<-5
  n.iter.posterior<-20000
  n.iter.WAIC<-10000
  
  
  M  <- paste( model,"_",link,sep="")
  
  if( M == "M1_Logit")	
  {
    
    #### Compute the BIC first
    L1<-glmer(Acc~(1|subj)+(1|itemID),data=m_data,family=binomial(logit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
    L1.bic<-BIC(L1)
    
    
    ######## Compute the WAIC
    L1<-'model {
    for (l in 1:n_obs) {
    y[l] ~ dbern(p[l])
    logit(p[l])<-beta0+a[item[l]]+b[subject[l]]
    loglik[l]<-y[l]*log(p[l]) + (1-y[l])*log(1-p[l])
    }
    loglik_full<-sum(loglik)
    #priors
    for (k in 1:n_subject) {
    b[k] ~ dnorm(0,tau_b) #random subject effects
    }
    for (j in 1:n_item) {
    a[j] ~ dnorm(0,tau_a) #random item effects
    }
    beta0 ~ dt(0,0.01,1) #overall intercept
    #hyper-priors
    tau_a ~ dgamma(0.01,0.01)
    tau_b ~ dgamma(0.01,0.01)
    sigma_a <-sqrt(1/tau_a)
    sigma_b <-sqrt(1/tau_b)
  }'

    
    jags <- jags.model(textConnection(L1),data = list('n_item'=n.item,'n_obs'=n.obs,'n_subject'=n.subject, 'y'=y,'item'=item,'subject'=subject),n.chains = 1,n.adapt = 1000)
    #run markov chain to burnin
    update(jags, 1000)
    
    ### Compute WAIC
    fit<-coda.samples(model=jags, variable.names=c('loglik'), n.iter= n.iter.WAIC)
    log_lik<-fit[[1]]
    S.waic<-nrow(log_lik)
    n.waic<-ncol(log_lik)
    lpd<-sum(log(colMeans(exp(log_lik))))
    p_waic<-sum(apply(log_lik, 2, var))
    L1.waic<--2*(lpd-p_waic)
    rm(log_lik)
    rm(fit)
    post<-coda.samples(model=jags,variable.names=c('beta0','sigma_a','sigma_b','a','b'), n.iter= n.iter.posterior)
    
    return(list("bic" =L1.bic,"waic"=L1.waic,"poster_distribution" = post))
    }
  
  else if (M == "M2_Logit")
  {
    ##### Compute the BIC
    L2<-glmer(Acc~cond+(1|subj)+(1|itemID),data=accuracy,family=binomial(logit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
    L2.bic<-BIC(L2)
    
    
    #####Compute the WAIC
    #Model 2, LOGISITC; 
    #fixed effect for condition and random subject and item effects
    L2<-'model {
    for (l in 1:n_obs) {
    y[l] ~ dbern(p[l])
    logit(p[l])<-beta0+alpha[condition[l]]+a[item[l]]+b[subject[l]]
    loglik[l]<-y[l]*log(p[l]) + (1-y[l])*log(1-p[l])
    }
    loglik_full<-sum(loglik)
    #priors
    for (k in 1:n_subject) {
    b[k] ~ dnorm(0,tau_b) #random subject effects
    }
    for (j in 1:n_item) {
    a[j] ~ dnorm(0,tau_a) #random item effects
    }
    alpha[1]<-0 #identification constraint
    for (i in 2:n_cond) {
    alpha[i] ~dt(0,0.16,1) #prior for fixed effects
    }
    beta0 ~ dt(0,0.01,1) #overall intercept
    #hyper-priors
    tau_a ~ dgamma(0.01,0.01)
    tau_b ~ dgamma(0.01,0.01)
    sigma_a <-sqrt(1/tau_a)
    sigma_b <-sqrt(1/tau_b)
  }'

    
    ###################FIT L2 - EXAMINE POSTERIOR
    jags <- jags.model(textConnection(L2),data = list('n_item'=n.item,'n_obs'=n.obs,'n_subject'=n.subject, 'y'=y,'item'=item,'subject'=subject),n.chains = 1,n.adapt = 1000)
    #run markov chain to burnin
    update(jags, 1000)
    
    ### Compute WAIC
    fit<-coda.samples(model=jags, variable.names=c('loglik'), n.iter= n.iter.WAIC)
    log_lik<-fit[[1]]
    S.waic<-nrow(log_lik)
    n.waic<-ncol(log_lik)
    lpd<-sum(log(colMeans(exp(log_lik))))
    p_waic<-sum(apply(log_lik, 2, var))
    L2.waic<--2*(lpd-p_waic)
    rm(log_lik)
    rm(fit)
    
    ########Posterior distribution
    jags <- jags.model(textConnection(L2),data = list('n_item'=n.item,'n_obs'=n.obs,'n_subject'=n.subject, 'y'=y,'item'=item,'subject'=subject),n.chains = 1,n.adapt = 1000)
    #run markov chain to burnin
    update(jags, 1000)
    #sample from stationary distribution
    fit<-coda.samples(model=jags,variable.names=c('beta0','sigma_a','sigma_b',,'a','b'), n.iter= n.iter.posterior)
    post <- fit
    
    return(list("bic" =L2.bic,"waic"=L2.waic,"post_summary" = post))	
    }
  
  
  else if (M == "M3_Logit")
  { 
    #Model 3, LOGISITC; 
    #random subject and item effects and random effect for condition depending on subjects
    
    # Compute BIC
    L3<-glmer(Acc~cond+(1+cond|subj)+(1|itemID),data=accuracy,family=binomial(logit),nAGQ=1)
    L3.bic<-BIC(L3)
    
    # Compute WAIC
    L3<-'model {
    for (l in 1:n_obs) {
    y[l] ~ dbern(p[l])
    logit(p[l])<-beta0+alpha[condition[l]]+a[item[l]]+b[subject[l]] + alpha_b[condition[l],subject[l]]
    loglik[l]<-y[l]*log(p[l]) + (1-y[l])*log(1-p[l])
    }
    loglik_full<-sum(loglik)
    #priors
    for (k in 1:n_subject) {
    b[k] ~ dnorm(0,tau_b) #random subject effects
    }
    for (j in 1:n_item) {
    a[j] ~ dnorm(0,tau_a) #random item effects
    }
    alpha[1]<-0 #identification constraint
    for (i in 2:n_cond) {
    alpha[i] ~dt(0,0.16,1) #prior for fixed effects
    }
    beta0 ~ dt(0,0.01,1) #overall intercept
    
    #condition by subject interaction
    for (k in 1:n_subject){
    alpha_b[1,k] <-0
    }
    for (i in 2:n_cond){
    for (k in 1:n_subject){
    alpha_b[i,k] ~ dnorm(0,tau_alpha_b)
    }
    }
    
    #hyper-priors
    tau_a ~ dgamma(0.01,0.01)
    tau_b ~ dgamma(0.01,0.01)
    tau_alpha_b ~ dgamma(0.01,0.01)
    sigma_a <-sqrt(1/tau_a)
    sigma_b <-sqrt(1/tau_b)
    sigma_alpha_b <-sqrt(1/tau_alpha_b)

  }'
    
   
    
    
    ###################FIT L3 - EXAMINE POSTERIOR AND COMPUTE WAIC
    #initialize the model; set tuning parameters based on n.adapt runs
    jags <- jags.model(textConnection(L3),data = list('n_item'=n.item,'n_cond'=n.cond,'n_obs'=n.obs,'n_subject'=n.subject,'y'=y,'condition'=condition,'item'=item,'subject'=subject),
    n.chains = 1,n.adapt = 1000)
    #run markov chain to burnin
    update(jags, 1000)
    ## compute WAIC
    fit<-coda.samples(model=jags, variable.names=c('loglik'), n.iter= n.iter.WAIC)
    log_lik<-fit[[1]]
    S.waic<-nrow(log_lik)
    n.waic<-ncol(log_lik)
    lpd<-sum(log(colMeans(exp(log_lik))))
    p_waic<-sum(apply(log_lik, 2, var))
    L3.waic<--2*(lpd-p_waic)
    rm(log_lik)
    rm(fit)
    post <- coda.samples(model=jags, variable.names=c('beta0','sigma_a','sigma_b','alpha','a','b','sigma_alpha_b','alpha_b'), n.iter= n.iter.posterior)
    return(list("bic" =L3.bic,"waic"=L3.waic,"post_summary" = post))
  }
  
  else if (M == "M4_Logit")
  { 
  		L4<- glmer(Acc~cond+(1|subj)+(1+cond|itemID),data=accuracy,family=binomial(logit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))

L4.bic<- BIC(L4)
		
#Model 4, LOGISITC; 
#random subject and item effects and random effect for condition depending on items
L4<-'model {
for (l in 1:n_obs) {
y[l] ~ dbern(p[l])
logit(p[l])<-beta0+alpha[condition[l]]+a[item[l]]+b[subject[l]] + alpha_a[condition[l],item[l]]
loglik[l]<-y[l]*log(p[l]) + (1-y[l])*log(1-p[l])
}
  loglik_full<-sum(loglik)
#priors
for (k in 1:n_subject) {
b[k] ~ dnorm(0,tau_b) #random subject effects
}
for (j in 1:n_item) {
a[j] ~ dnorm(0,tau_a) #random item effects
}
alpha[1]<-0 #identification constraint
for (i in 2:n_cond) {
alpha[i] ~dt(0,0.16,1) #prior for fixed effects
}
beta0 ~ dt(0,0.01,1) #overall intercept

#condition by item interaction
for (k in 1:n_item){
alpha_a[1,k] <-0
}
for (i in 2:n_cond){
for (k in 1:n_item){
alpha_a[i,k] ~ dnorm(0,tau_alpha_a)
}
}

#hyper-priors
tau_a ~ dgamma(0.01,0.01)
tau_b ~ dgamma(0.01,0.01)
tau_alpha_a ~ dgamma(0.01,0.01)
sigma_a <-sqrt(1/tau_a)
sigma_b <-sqrt(1/tau_b)
sigma_alpha_a <-sqrt(1/tau_alpha_a)
}'


###################FIT L4 - EXAMINE POSTERIOR
#initialize the model; set tuning parameters based on n.adapt runs
jags <- jags.model(textConnection(L4),data = list('n_item'=n.item,'n_cond'=n.cond,'n_obs'=n.obs,'n_subject'=n.subject, 'y'=y,'condition'=condition,'item'=item,'subject'=subject),n.chains = 1,n.adapt = 1000)
#run markov chain to burnin
update(jags, 1000)

	## compute WAIC
	fit<-coda.samples(model=jags, variable.names=c('loglik'), n.iter= n.iter.WAIC)
	log_lik<-fit[[1]]
	S.waic<-nrow(log_lik)
	n.waic<-ncol(log_lik)
	lpd<-sum(log(colMeans(exp(log_lik))))
	p_waic<-sum(apply(log_lik, 2, var))
	L4.waic<--2*(lpd-p_waic)
	rm(log_lik)
	rm(fit)
	post <- coda.samples(model=jags, variable.names=c('beta0','sigma_a','sigma_b','alpha','a','b','sigma_alpha_a','alpha_a'), n.iter= n.iter.posterior)
  	return(list("bic" =L4.bic,"waic"=L4.waic,"post_summary" = post))	
  }

	else if (M == "M5_Logit")
   { 
	#Model 5, LOGISITC; 
	#random subject and item effects and random interaction effect for condition depending on items as well as condition depending on subjects
	L5<-glmer(Acc~cond+(1+cond|subj)+(1+cond|itemID),data=accuracy,family=binomial(logit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
	L5.bic<-BIC(L5)
	
	#Model 5, LOGISITC; 
#random subject and item effects and random interaction effect for condition depending on items as well as condition depending on subjects
L5<-'model {
for (l in 1:n_obs) {
y[l] ~ dbern(p[l])
logit(p[l])<-beta0+alpha[condition[l]]+a[item[l]]+b[subject[l]] + alpha_a[condition[l],item[l]] + alpha_b[condition[l],subject[l]]
loglik[l]<-y[l]*log(p[l]) + (1-y[l])*log(1-p[l])
}
  loglik_full<-sum(loglik)
#priors
for (k in 1:n_subject) {
b[k] ~ dnorm(0,tau_b) #random subject effects
}
for (j in 1:n_item) {
a[j] ~ dnorm(0,tau_a) #random item effects
}
alpha[1]<-0 #identification constraint
for (i in 2:n_cond) {
alpha[i] ~dt(0,0.16,1) #prior for fixed effects
}
beta0 ~ dt(0,0.01,1) #overall intercept

#condition by item interaction
for (k in 1:n_item){
alpha_a[1,k] <-0
}
for (i in 2:n_cond){
for (k in 1:n_item){
alpha_a[i,k] ~ dnorm(0,tau_alpha_a)
}
}

#condition by subject interaction
for (k in 1:n_subject){
alpha_b[1,k] <-0
}
for (i in 2:n_cond){
for (k in 1:n_subject){
alpha_b[i,k] ~ dnorm(0,tau_alpha_b)
}
}

#hyper-priors
tau_a ~ dgamma(0.01,0.01)
tau_b ~ dgamma(0.01,0.01)
tau_alpha_a ~ dgamma(0.01,0.01)
tau_alpha_b ~ dgamma(0.01,0.01)
sigma_a <-sqrt(1/tau_a)
sigma_b <-sqrt(1/tau_b)
sigma_alpha_a <-sqrt(1/tau_alpha_a)
sigma_alpha_b <-sqrt(1/tau_alpha_b)
}'

###################FIT L5 - EXAMINE POSTERIOR
#initialize the model; set tuning parameters based on n.adapt runs
jags <- jags.model(textConnection(L5),data = list('n_item'=n.item,'n_cond'=n.cond,'n_obs'=n.obs,'n_subject'=n.subject,'y'=y,'condition'=condition,'item'=item,'subject'=subject), n.chains = 1,n.adapt = 1000)
#run markov chain to burnin
update(jags, 1000)

## compute WAIC
fit<-coda.samples(model=jags, variable.names=c('loglik'), n.iter= n.iter.WAIC)
log_lik<-fit[[1]]
S.waic<-nrow(log_lik)
n.waic<-ncol(log_lik)
lpd<-sum(log(colMeans(exp(log_lik))))
p_waic<-sum(apply(log_lik, 2, var))
L5.waic<--2*(lpd-p_waic)
rm(log_lik)
rm(fit)

post <- coda.samples(model=jags, variable.names=c('beta0','sigma_a','sigma_b','alpha','a','b','sigma_alpha_a','sigma_alpha_b','alpha_a','alpha_b'), n.iter= n.iter.posterior)
  return(list("bic" =L5.bic,"waic"=L5.waic,"post_summary" = post))	
	}
	
	else if (M == "M1_Probit")
   { 
 
 #Model 1, PROBIT; 
#no effect for condition and random subject and item effects
P1<-glmer(Acc~(1|subj)+(1|itemID),data=accuracy,family=binomial(probit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
P1.bic<-BIC(P1)

#Model 1, PROBIT; 
#no effect for condition and random subject and item effects
P1<-'model {
for (l in 1:n_obs) {
y[l] ~ dbern(p[l])
probit(p[l])<-beta0+a[item[l]]+b[subject[l]]
loglik[l]<-y[l]*log(p[l]) + (1-y[l])*log(1-p[l])
}
  loglik_full<-sum(loglik)
#priors
for (k in 1:n_subject) {
b[k] ~ dnorm(0,tau_b) #random subject effects
}
for (j in 1:n_item) {
a[j] ~ dnorm(0,tau_a) #random item effects
}
beta0 ~ dt(0,0.01,1) #overall intercept
#hyper-priors
tau_a ~ dgamma(0.01,0.01)
tau_b ~ dgamma(0.01,0.01)
sigma_a <-sqrt(1/tau_a)
sigma_b <-sqrt(1/tau_b)
}'

###################FIT P1 - EXAMINE POSTERIOR
#initialize the model; set tuning parameters based on n.adapt runs
jags <- jags.model(textConnection(P1),data = list('n_item'=n.item,'n_obs'=n.obs,'n_subject'=n.subject,'y'=y,'item'=item,'subject'=subject),n.chains = 1,n.adapt = 1000)
#run markov chain to burnin
update(jags, 1000)

## compute WAIC
fit<-coda.samples(model=jags, variable.names=c('loglik'), n.iter= n.iter.WAIC)
log_lik<-fit[[1]]
S.waic<-nrow(log_lik)
n.waic<-ncol(log_lik)
lpd<-sum(log(colMeans(exp(log_lik))))
p_waic<-sum(apply(log_lik, 2, var))
P1.waic<--2*(lpd-p_waic)
rm(log_lik)
rm(fit)

  post<-coda.samples(model=jags,variable.names=c('beta0','sigma_a','sigma_b','a','b'), n.iter= n.iter.posterior)

  return(list("bic" =P1.bic,"waic"=P1.waic,"post_summary" = post))	
	}
	
	else if (M == "M2_Probit")
   { 
#Model 2, PROBIT; 
#fixed effect for condition and random subject and item effects
P2<-glmer(Acc~cond+(1|subj)+(1|itemID),data=accuracy,family=binomial(probit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
P2.bic<-BIC(P2)

#Model 2, PROBIT; 
#fixed effect for condition and random subject and item effects
P2<-'model {
  for (l in 1:n_obs) {
    y[l] ~ dbern(p[l])
    probit(p[l])<-beta0+alpha[condition[l]]+a[item[l]]+b[subject[l]]
    loglik[l]<-y[l]*log(p[l]) + (1-y[l])*log(1-p[l])
  }
  loglik_full<-sum(loglik)
  #priors
  for (k in 1:n_subject) {
  b[k] ~ dnorm(0,tau_b) #random subject effects
  }
  for (j in 1:n_item) {
  a[j] ~ dnorm(0,tau_a) #random item effects
  }
  alpha[1]<-0 #identification constraint
  for (i in 2:n_cond) {
  alpha[i] ~dt(0,0.16,1) #prior for fixed effects
  }
  beta0 ~ dt(0,0.01,1) #overall intercept
  #hyper-priors
  tau_a ~ dgamma(0.01,0.01)
  tau_b ~ dgamma(0.01,0.01)
  sigma_a <-sqrt(1/tau_a)
  sigma_b <-sqrt(1/tau_b)
}'

###################FIT P2 - EXAMINE POSTERIOR
#initialize the model; set tuning parameters based on n.adapt runs
jags <- jags.model(textConnection(P2),data = list('n_item'=n.item,'n_cond'=n.cond,'n_obs'=n.obs,'n_subject'=n.subject,'y'=y,'condition'=condition,'item'=item,'subject'=subject),n.chains = 1,n.adapt = 1000)
#run markov chain to burnin
update(jags, 1000)

## Compute WAIC
fit<-coda.samples(model=jags, variable.names=c('loglik'), n.iter= n.iter.WAIC)
log_lik<-fit[[1]]
S.waic<-nrow(log_lik)
n.waic<-ncol(log_lik)
lpd<-sum(log(colMeans(exp(log_lik))))
p_waic<-sum(apply(log_lik, 2, var))
P2.waic<--2*(lpd-p_waic)
rm(log_lik)
rm(fit)

 post<-coda.samples(model=jags,variable.names=c('beta0','sigma_a','sigma_b','a','b','alpha'), n.iter= n.iter.posterior)
 
 return(list("bic" =P2.bic,"waic"=P2.waic,"post_summary" = post))	
   	}
   	
	else if (M == "M3_Probit")
   { 	
   	#Model 3, PROBIT; 
#random subject and item effects and random effect for condition depending on subjects
P3<-glmer(Acc~cond+(1+cond|subj)+(1|itemID),data=accuracy,family=binomial(probit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
P3.bic<-BIC(P3)

#Model 3, PROBIT; 
#random subject and item effects and random effect for condition depending on subjects
P3<-'model {
for (l in 1:n_obs) {
y[l] ~ dbern(p[l])
probit(p[l])<-beta0+alpha[condition[l]]+a[item[l]]+b[subject[l]] + alpha_b[condition[l],subject[l]]
loglik[l]<-y[l]*log(p[l]) + (1-y[l])*log(1-p[l])
}
  loglik_full<-sum(loglik)
#priors
for (k in 1:n_subject) {
b[k] ~ dnorm(0,tau_b) #random subject effects
}
for (j in 1:n_item) {
a[j] ~ dnorm(0,tau_a) #random item effects
}
alpha[1]<-0 #identification constraint
for (i in 2:n_cond) {
alpha[i] ~dt(0,0.16,1) #prior for fixed effects
}
beta0 ~ dt(0,0.01,1) #overall intercept

#condition by subject interaction
for (k in 1:n_subject){
alpha_b[1,k] <-0
}
for (i in 2:n_cond){
for (k in 1:n_subject){
alpha_b[i,k] ~ dnorm(0,tau_alpha_b)
}
}

#hyper-priors
tau_a ~ dgamma(0.01,0.01)
tau_b ~ dgamma(0.01,0.01)
tau_alpha_b ~ dgamma(0.01,0.01)
sigma_a <-sqrt(1/tau_a)
sigma_b <-sqrt(1/tau_b)
sigma_alpha_b <-sqrt(1/tau_alpha_b)
}'

	###################FIT P3 - EXAMINE POSTERIOR
#initialize the model; set tuning parameters based on n.adapt runs
jags <- jags.model(textConnection(P3),data = list('n_item'=n.item,'n_cond'=n.cond,'n_obs'=n.obs,'n_subject'=n.subject,'y'=y,'condition'=condition,'item'=item,'subject'=subject), n.chains = 1,n.adapt = 1000)
#run markov chain to burnin
update(jags, 1000)

## compute WAIC
fit<-coda.samples(model=jags, variable.names=c('loglik'), n.iter= n.iter.WAIC)
log_lik<-fit[[1]]
S.waic<-nrow(log_lik)
n.waic<-ncol(log_lik)
lpd<-sum(log(colMeans(exp(log_lik))))
p_waic<-sum(apply(log_lik, 2, var))
P3.waic<--2*(lpd-p_waic)
rm(log_lik)
rm(fit)
   	
 post<-coda.samples(model=jags,variable.names=c('beta0','sigma_a','sigma_b','a','b','alpha_b','sigma_alpha_b'), n.iter= n.iter.posterior)
 
 return(list("bic" =P3.bic,"waic"=P3.waic,"post_summary" = post))	
   	}
   	
   	else if (M == "M4_Probit")
   { 
   	#Model 4, PROBIT; 
#random subject and item effects and random effect for condition depending on items
P4<-glmer(Acc~cond+(1|subj)+(1+cond|itemID),data=accuracy,family=binomial(probit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
P4.bic<-BIC(P4)

#Model 4, PROBIT; 
#random subject and item effects and random effect for condition depending on items
P4<-'model {
for (l in 1:n_obs) {
y[l] ~ dbern(p[l])
probit(p[l])<-beta0+alpha[condition[l]]+a[item[l]]+b[subject[l]] + alpha_a[condition[l],item[l]]
loglik[l]<-y[l]*log(p[l]) + (1-y[l])*log(1-p[l])
}
  loglik_full<-sum(loglik)
#priors
for (k in 1:n_subject) {
b[k] ~ dnorm(0,tau_b) #random subject effects
}
for (j in 1:n_item) {
a[j] ~ dnorm(0,tau_a) #random item effects
}
alpha[1]<-0 #identification constraint
for (i in 2:n_cond) {
alpha[i] ~dt(0,0.16,1) #prior for fixed effects
}
beta0 ~ dt(0,0.01,1) #overall intercept

#condition by item interaction
for (k in 1:n_item){
alpha_a[1,k] <-0
}
for (i in 2:n_cond){
for (k in 1:n_item){
alpha_a[i,k] ~ dnorm(0,tau_alpha_a)
}
}

#hyper-priors
tau_a ~ dgamma(0.01,0.01)
tau_b ~ dgamma(0.01,0.01)
tau_alpha_a ~ dgamma(0.01,0.01)
sigma_a <-sqrt(1/tau_a)
sigma_b <-sqrt(1/tau_b)
sigma_alpha_a <-sqrt(1/tau_alpha_a)
}'

###################FIT P4 - EXAMINE POSTERIOR
#initialize the model; set tuning parameters based on n.adapt runs
jags <- jags.model(textConnection(P4),data = list('n_item'=n.item,'n_cond'=n.cond,'n_obs'=n.obs,'n_subject'=n.subject,'y'=y,'condition'=condition,'item'=item,'subject'=subject), n.chains = 1,n.adapt = 1000)
#run markov chain to burnin
update(jags, 1000)

## compute WAIC
fit<-coda.samples(model=jags, variable.names=c('loglik'), n.iter= n.iter.WAIC)
log_lik<-fit[[1]]
S.waic<-nrow(log_lik)
n.waic<-ncol(log_lik)
lpd<-sum(log(colMeans(exp(log_lik))))
p_waic<-sum(apply(log_lik, 2, var))
P4.waic<--2*(lpd-p_waic)
rm(log_lik)
rm(fit)

post<-coda.samples(model=jags,variable.names=c('beta0','sigma_a','sigma_b','a','b','alpha_a','sigma_alpha_a'), n.iter= n.iter.posterior)
 
 return(list("bic" =P4.bic,"waic"=P4.waic,"post_summary" = post))	
   	}
   	
   	 	else if (M == "M5_Probit")
   { 
   	
 #Model 5, PROBIT; 
#random subject and item effects and random interaction effect for condition depending on items as well as condition depending on subjects
P5<-glmer(Acc~cond+(1+cond|subj)+(1+cond|itemID),data=accuracy,family=binomial(probit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
P5.bic<-BIC(P5)

#Model 5, PROBIT; 
#random subject and item effects and random interaction effect for condition depending on items as well as condition depending on subjects
P5<-'model {
for (l in 1:n_obs) {
y[l] ~ dbern(p[l])
probit(p[l])<-beta0+alpha[condition[l]]+a[item[l]]+b[subject[l]] + alpha_a[condition[l],item[l]] + alpha_b[condition[l],subject[l]]
loglik[l]<-y[l]*log(p[l]) + (1-y[l])*log(1-p[l])
}
  loglik_full<-sum(loglik)
#priors
for (k in 1:n_subject) {
b[k] ~ dnorm(0,tau_b) #random subject effects
}
for (j in 1:n_item) {
a[j] ~ dnorm(0,tau_a) #random item effects
}
alpha[1]<-0 #identification constraint
for (i in 2:n_cond) {
alpha[i] ~dt(0,0.16,1) #prior for fixed effects
}
beta0 ~ dt(0,0.01,1) #overall intercept

#condition by item interaction
for (k in 1:n_item){
alpha_a[1,k] <-0
}
for (i in 2:n_cond){
for (k in 1:n_item){
alpha_a[i,k] ~ dnorm(0,tau_alpha_a)
}
}

#condition by subject interaction
for (k in 1:n_subject){
alpha_b[1,k] <-0
}
for (i in 2:n_cond){
for (k in 1:n_subject){
alpha_b[i,k] ~ dnorm(0,tau_alpha_b)
}
}

#hyper-priors
tau_a ~ dgamma(0.01,0.01)
tau_b ~ dgamma(0.01,0.01)
tau_alpha_a ~ dgamma(0.01,0.01)
tau_alpha_b ~ dgamma(0.01,0.01)
sigma_a <-sqrt(1/tau_a)
sigma_b <-sqrt(1/tau_b)
sigma_alpha_a <-sqrt(1/tau_alpha_a)
sigma_alpha_b <-sqrt(1/tau_alpha_b)
}'

###################FIT P5 - EXAMINE POSTERIOR
#initialize the model; set tuning parameters based on n.adapt runs
jags <- jags.model(textConnection(P5),data = list('n_item'=n.item,'n_cond'=n.cond,'n_obs'=n.obs,'n_subject'=n.subject,'y'=y,'condition'=condition,'item'=item,'subject'=subject), n.chains = 1,n.adapt = 1000)
#run markov chain to burnin
update(jags, 1000)

## compute WAIC
fit<-coda.samples(model=jags, variable.names=c('loglik'), n.iter= n.iter.WAIC)
log_lik<-fit[[1]]
S.waic<-nrow(log_lik)
n.waic<-ncol(log_lik) 
lpd<-sum(log(colMeans(exp(log_lik))))
p_waic<-sum(apply(log_lik, 2, var))
WAIC<--2*(lpd-p_waic)
rm(log_lik)
rm(fit)

post<-coda.samples(model=jags,variable.names=c('beta0','sigma_a','sigma_b','a','b','alpha_a','alpha_b','sigma_alpha_a','sigma_alpha_b'), n.iter= n.iter.posterior)

 return(list("bic" = P5.bic,"waic" = P5.waic,"post_summary" = post))	

   	}
   	
}
