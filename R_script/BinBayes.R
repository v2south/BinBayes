BinBayes <- function(factor, m_data, model_struct, link,baseline=NULL)
{
	library(coda)
	library(lme4)
 	library(rjags)

	if  (factor ==1)
	{
			# Rename the columns in the data. 
			colnames(m_data) <- c("subj","itemID", "cond","Acc")
			n.obs<-nrow(m_data)
			n.subject<-length(unique(m_data$subj))
			n.item<-length(unique(m_data$itemID))
			n.cond<-length(unique(m_data$cond))
  			
  			#store the labels used for subjects, conditions, items
  			labels.subject<-unique(m_data$subj)
  			labels.item<-unique(m_data$itemID)
  			labels.cond<-unique(m_data$cond)
			
			# Determine the baseline condition 
			if (is.null(baseline))
					{
	 					labels.cond<-unique(m_data$cond)
						baseline  <- labels.cond[1]
					} 
					else
					{	
	  					labels.cond<-unique(m_data$cond)
	  					index <- which(baseline == labels.cond)
	  					temp <- labels.cond[1]
	 					labels.cond[index] <- temp
	  					labels.cond[1]  <- baseline
					}
					
			#setup m_data vectors as needed for jags
  			subject<-match(m_data$subj,labels.subject)
  			item<-match(m_data$itemID,labels.item)
  			condition<-match(m_data$cond,labels.cond)
  			y<-as.numeric(m_data$Acc)		
			
			n.iter.marginal<-250000
  			n.thin.marginal<-5
  			n.iter.posterior<-20000
  			n.iter.WAIC<-10000
			
		M  <- paste( "M",model_struct,"_",link,sep="")
		if( M == "M1_Logit")	
  {
    
    #### Compute the BIC first
    L1<-glmer(Acc~(1|subj)+(1|itemID),data=m_data,family=binomial(logit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
    L1.bic<-BIC(L1)
    
    
    #### Compute the WAIC
    L1<-'model 
    {
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
  }
  '

    
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
    
    return(list("bic" =L1.bic,"waic"=L1.waic,"post_summary" = post, "condition_level" = labels.cond, "baseline"= baseline))
    }

		else if (M == "M2_Logit")
  {
    ##### Compute the BIC
    L2<-glmer(Acc~cond+(1|subj)+(1|itemID),data=m_data,family=binomial(logit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
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
  jags <- jags.model(textConnection(L2),data = list('n_item'=n.item,'n_cond'=n.cond,'n_obs'=n.obs,'n_subject'=n.subject, 
                                                  'y'=y,'condition'=condition,'item'=item,'subject'=subject),
                   n.chains = 1,n.adapt = 1000)
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

    #sample from stationary distribution
    post<-coda.samples(model=jags,variable.names=c('beta0','sigma_a','sigma_b','alpha','a','b'), n.iter= n.iter.posterior)
        
    return(list("bic" =L2.bic,"waic"=L2.waic,"post_summary" = post, "condition_level" = labels.cond, "baseline"= baseline))
    }
  
  
  else if (M == "M3_Logit")
  { 
    #Model 3, LOGISITC; 
    #random subject and item effects and random effect for condition depending on subjects
    
    # Compute BIC
    L3<-glmer(Acc~cond+(1+cond|subj)+(1|itemID),data=m_data,family=binomial(logit),nAGQ=1)
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
    
    return(list("bic" =L3.bic,"waic"=L3.waic,"post_summary" = post, "condition_level" = labels.cond, "baseline"= baseline))
  }
  
  else if (M == "M4_Logit")
  { 
  		L4<- glmer(Acc~cond+(1|subj)+(1+cond|itemID),data=m_data,family=binomial(logit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))

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
    return(list("bic" =L4.bic,"waic"=L4.waic,"post_summary" = post, "condition_level" = labels.cond, "baseline"= baseline))
  }

	else if (M == "M5_Logit")
   { 
	#Model 5, LOGISITC; 
	#random subject and item effects and random interaction effect for condition depending on items as well as condition depending on subjects
	L5<-glmer(Acc~cond+(1+cond|subj)+(1+cond|itemID),data=m_data,family=binomial(logit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
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
    return(list("bic" =L5.bic,"waic"=L5.waic,"post_summary" = post, "condition_level" = labels.cond, "baseline"= baseline))
	}
	
	else if (M == "M1_Probit")
   { 
 
 #Model 1, PROBIT; 
#no effect for condition and random subject and item effects
P1<-glmer(Acc~(1|subj)+(1|itemID),data=m_data,family=binomial(probit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
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

    return(list("bic" =P1.bic,"waic"=P1.waic,"post_summary" = post, "condition_level" = labels.cond, "baseline"= baseline))
	}
	
	else if (M == "M2_Probit")
   { 
#Model 2, PROBIT; 
#fixed effect for condition and random subject and item effects
P2<-glmer(Acc~cond+(1|subj)+(1|itemID),data=m_data,family=binomial(probit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
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
 
    return(list("bic" =P2.bic,"waic"=P2.waic,"post_summary" = post, "condition_level" = labels.cond, "baseline"= baseline))
   	}
   	
	else if (M == "M3_Probit")
   { 	
   	#Model 3, PROBIT; 
#random subject and item effects and random effect for condition depending on subjects
P3<-glmer(Acc~cond+(1+cond|subj)+(1|itemID),data=m_data,family=binomial(probit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
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
 
    return(list("bic" =P3.bic,"waic"=P3.waic,"post_summary" = post, "condition_level" = labels.cond, "baseline"= baseline))
   	}
   	
   	else if (M == "M4_Probit")
   { 
   	#Model 4, PROBIT; 
#random subject and item effects and random effect for condition depending on items
P4<-glmer(Acc~cond+(1|subj)+(1+cond|itemID),data=m_data,family=binomial(probit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
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
 
    return(list("bic" =P4.bic,"waic"=P4.waic,"post_summary" = post, "condition_level" = labels.cond, "baseline"= baseline))
   	}
   	
   	 	else if (M == "M5_Probit")
   { 
   	
 #Model 5, PROBIT; 
#random subject and item effects and random interaction effect for condition depending on items as well as condition depending on subjects
P5<-glmer(Acc~cond+(1+cond|subj)+(1+cond|itemID),data=m_data,family=binomial(probit),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
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
P5.waic<--2*(lpd-p_waic)
rm(log_lik)
rm(fit)

post<-coda.samples(model=jags,variable.names=c('beta0','sigma_a','sigma_b','a','b','alpha_a','alpha_b','sigma_alpha_a','sigma_alpha_b'), n.iter= n.iter.posterior)

    return(list("bic" =P5.bic,"waic"=P5.waic,"post_summary" = post, "condition_level" = labels.cond, "baseline"= baseline))

   	}
	
	
	}
	
	else if(factor ==2)
	{
			# Rename the columns 
			colnames(m_data) <- c("subj","itemID", "cond1","cond2","Acc")
			# number of observation
			n.obs<-nrow(m_data)
			n.subject<-length(unique(m_data$sub))
			n.item<-length(unique(m_data$itemID))
			n.cond.factor1<-length(unique(m_data$cond1))
			n.cond.factor2<-length(unique(m_data$cond2))


			#store the labels used for subjects, conditions, items
	  		labels.subject<-unique(m_data$subj)
  			labels.item<-unique(m_data$itemID)
  			labels.cond.factor1<-unique(m_data$cond1)
  			labels.cond.factor2<-unique(m_data$cond2)
			condition  <- list(labels.cond.factor1, labels.cond.factor2)
			# Determine the baseline condition 
			if (is.null(baseline))
					{
	 					labels.cond.factor1<-unique(m_data$cond1)
	 					labels.cond.factor2<-unique(m_data$cond2)
					    baseline.factor1<- labels.cond.factor1[1]
						baseline.factor2 <- labels.cond.factor2[1]
						baseline  <- list(baseline.factor1, baseline.factor2)
					} 
				else if (is.na(baseline[1]))
					{
						# Factor 1 is null, will use the default baseline
						labels.cond.factor1<-unique(m_data$cond1)
					    baseline.factor1<- labels.cond.factor1[1]
					    
					    # Factor 2 is defined. 
					    labels.cond.factor2<-unique(m_data$cond2)
	  					index.factor2 <- which(baseline[2] == 	labels.cond.factor2)
	  					temp2 <- labels.cond.factor2[1]
	 					labels.cond.factor2[index.factor2] <- temp2
	  					labels.cond.factor2[1]  <- baseline[2]
						
						baseline  <- list(baseline.factor1, baseline[2])
					}
						else if (is.na(baseline[2]))
					{
						# Factor 1 is defined baseline
						labels.cond.factor1<-unique(m_data$cond1)
	 					index.factor1 <- which(baseline[1] == 	labels.cond.factor1)
	  					temp1 <- labels.cond.factor1[1]
	 					labels.cond.factor1[index.factor1] <- temp1
	  					labels.cond.factor1[1]  <- baseline[1]
	  					
	  					# Factor 2 is null, will use the default baseline
	  					labels.cond.factor2<-unique(m_data$cond2)
						baseline.factor2 <- labels.cond.factor2[1]

						baseline  <- list(baseline[1], baseline.factor2)
					}
					
				else 
					{		
						# Factor 1
						labels.cond.factor1<-unique(m_data$cond1)
	 					index.factor1 <- which(baseline[1] == 	labels.cond.factor1)
	  					temp1 <- labels.cond.factor1[1]
	 					labels.cond.factor1[index.factor1] <- temp1
	  					labels.cond.factor1[1]  <- baseline[1]
	  					#Factor 2
	  					labels.cond.factor2<-unique(m_data$cond2)
	  					index.factor2 <- which(baseline[2] == 	labels.cond.factor2)
	  					temp2 <- labels.cond.factor2[1]
	 					labels.cond.factor2[index.factor2] <- temp2
	  					labels.cond.factor2[1]  <- baseline[2]
	  					baseline <- list(baseline[1], baseline[2])
					}

subject<-match(m_data[,1],labels.subject)
item<-match(m_data[,2],labels.item)
condition.factor1<-match(m_data[,3],labels.cond.factor1)
condition.factor2<-match(m_data[,4],labels.cond.factor2)
y<-m_data[,5]

# Two factor model
#input: link - either 'P' for probit link or 'L' for logit link
#input: M - value in {1,2,3,4,5} indicating model structure for first factor as described in paper
#input: N - value in {1,2,3,4,5} indicating model structure for the second factor as described in paper
#input: I - value in {0,1} indicating presense of interaction between factors, I=1 gives the interaction
#input: m_data - dataframe with columns: "subj","itemID", "cond1","cond2","Acc"
#output: a list, the first component is the model string for jags, 
#output: a list, the second component is a vector of strings giving that can be used for the variable.names input for coda.samples 
#output a list, the third component is a list that can be used for the data argument of jags.model
#output a list, the fourth component is a list to be used when calling glmer to compute the BIC

two.factor.model<-function(link,M,N,I,m_data)
			{ 
  final.list<-list()
  model.string <- 'model {\n for (l in 1:n_obs) {\ny[l] ~ dbern(p[l])\n'
  #set the link function
  if (link=='Logit')
  {model.string<-paste(model.string,'logit(p[l])<-beta0 + a[item[l]] + b[subject[l]]')}
  else
  {model.string<-paste(model.string,'probit(p[l])<-beta0 + a[item[l]] + b[subject[l]]')}
  
  #add the rest of the regression structure for factor 1
  if (M==2)
  {model.string<-paste(model.string,'+ alpha[condition.factor1[l]]')}
  else if (M==3)
  {model.string<-paste(model.string,'+ alpha[condition.factor1[l]] + alpha_b[condition.factor1[l],subject[l]]')}
  else if (M==4)
  {model.string<-paste(model.string,'+ alpha[condition.factor1[l]] + alpha_a[condition.factor1[l],item[l]]')}
  else if (M==5)
  {model.string<-paste(model.string,'+ alpha[condition.factor1[l]] + alpha_b[condition.factor1[l],subject[l]]+ alpha_a[condition.factor1[l],item[l]]')}
  
  #add the rest of the regression structure for factor 1
  if (N==2)
  {model.string<-paste(model.string,'+ gamma[condition.factor2[l]]')}
  else if (N==3)
  {model.string<-paste(model.string,'+ gamma[condition.factor2[l]] + gamma_b[condition.factor2[l],subject[l]]')}
  else if (N==4)
  {model.string<-paste(model.string,'+ gamma[condition.factor2[l]] + gamma_a[condition.factor2[l],item[l]]')}
  else if (N==5)
  {model.string<-paste(model.string,'+ gamma[condition.factor2[l]] + gamma_b[condition.factor2[l],subject[l]] + gamma_a[condition.factor2[l],item[l]]')}
  
  #add the interaction 
  if (I==1)
  {model.string<-paste(model.string,'+ gamma_alpha[condition.factor2[l],condition.factor1[l]]')}
  
  #add terms required for all models
  model.string<-paste(model.string,'\nloglik[l]<-y[l]*log(p[l]) + (1-y[l])*log(1-p[l])\n}\nloglik_full<-sum(loglik)\n\n#priors\n\n#subject effects\nfor (k in 1:n_subject) {\nb[k] ~ dnorm(0,tau_b) #random subject effects\n}\n\n#item effects\nfor (j in 1:n_item) {\na[j] ~ dnorm(0,tau_a) #random item effects\n}\n\n#intercept\nbeta0 ~ dt(0,0.01,1) #overall intercept\n\n')
  
  #add terms for prior for first factor
  if (M==2)
  {model.string<-paste(model.string,'\n\n alpha[1]<-0 \n for (i in 2:n_cond_factor1) {\nalpha[i] ~dt(0,0.16,1) \n}\n\n')}
  else if (M==3)
  {model.string<-paste(model.string,'\n\n alpha[1]<-0 \n for (i in 2:n_cond_factor1) {\nalpha[i] ~dt(0,0.16,1) \n}\n\n \n\n#factor1 by subject interaction\nfor (k in 1:n_subject){\nalpha_b[1,k] <-0\n}\nfor (i in 2:n_cond_factor1)\n  {\n  for (k in 1:n_subject){\n  alpha_b[i,k] ~ dnorm(0,tau_alpha_b)\n  }\n}\n\n')}
  else if (M==4)
  {model.string<-paste(model.string,'\n\n alpha[1]<-0 \n for (i in 2:n_cond_factor1) {\nalpha[i] ~dt(0,0.16,1) \n}\n\n \n\n#factor1 by item interaction\nfor (k in 1:n_item){\nalpha_a[1,k] <-0\n}\nfor (i in 2:n_cond_factor1){\n  for (k in 1:n_item){\n  alpha_a[i,k] ~ dnorm(0,tau_alpha_a)\n  }\n}\n')}
  else if (M==5)
  {model.string<-paste(model.string,'\n\n alpha[1]<-0 \n for (i in 2:n_cond_factor1) {\nalpha[i] ~dt(0,0.16,1) \n}\n\n \n\n#factor1 by subject interaction\nfor (k in 1:n_subject){\nalpha_b[1,k] <-0\n}\nfor (i in 2:n_cond_factor1)\n  {\n  for (k in 1:n_subject){\n  alpha_b[i,k] ~ dnorm(0,tau_alpha_b)\n  }\n}\n\n \n\n#factor1 by item interaction\nfor (k in 1:n_item){\nalpha_a[1,k] <-0\n}\nfor (i in 2:n_cond_factor1){\n  for (k in 1:n_item){\n  alpha_a[i,k] ~ dnorm(0,tau_alpha_a)\n  }\n}\n')}  

  #add terms for prior for second factor
  if (N==2)
  {model.string<-paste(model.string,'\n#fixed effect factor2 \n gamma[1]<-0 #identification constraint \n for (i in 2:n_cond_factor2) {\n gamma[i] ~dt(0,0.16,1) #prior for fixed effects\n}\n')}
  else if (N==3)
  {model.string<-paste(model.string,'\n#fixed effect factor2 \n gamma[1]<-0 #identification constraint \n for (i in 2:n_cond_factor2) {\n gamma[i] ~dt(0,0.16,1) #prior for fixed effects\n}\n \n\n#factor2 by subject interaction\nfor (k in 1:n_subject){\ngamma_b[1,k] <-0\n}\nfor (i in 2:n_cond_factor2)\n  {\n  for (k in 1:n_subject){\n  gamma_b[i,k] ~ dnorm(0,tau_gamma_b)\n  }\n}\n')}
  else if (N==4)
  {model.string<-paste(model.string,'\n#fixed effect factor2 \n gamma[1]<-0 #identification constraint \n for (i in 2:n_cond_factor2) {\n gamma[i] ~dt(0,0.16,1) #prior for fixed effects\n}\n \n#factor2 by item interaction\nfor (k in 1:n_item){\ngamma_a[1,k] <-0\n}\nfor (i in 2:n_cond_factor2){\nfor (k in 1:n_item){\ngamma_a[i,k] ~ dnorm(0,tau_gamma_a)\n}\n}\n')}
  else if (N==5)
  {model.string<-paste(model.string,'\n#fixed effect factor2 \n gamma[1]<-0 #identification constraint \n for (i in 2:n_cond_factor2) {\n gamma[i] ~dt(0,0.16,1) #prior for fixed effects\n}\n \n\n#factor2 by subject interaction\nfor (k in 1:n_subject){\ngamma_b[1,k] <-0\n}\nfor (i in 2:n_cond_factor2)\n  {\n  for (k in 1:n_subject){\n  gamma_b[i,k] ~ dnorm(0,tau_gamma_b)\n  }\n}\n \n#factor2 by item interaction\nfor (k in 1:n_item){\ngamma_a[1,k] <-0\n}\nfor (i in 2:n_cond_factor2){\nfor (k in 1:n_item){\ngamma_a[i,k] ~ dnorm(0,tau_gamma_a)\n}\n}\n')}
  
  #add terms for prior for interaction term
  if (I==1)
  {model.string<-paste(model.string,'\n#fixed interaction for two factors\nfor (i in 1:n_cond_factor1) {\ngamma_alpha[1,i] <-0 #prior for fixed effects\n}\nfor (i in 2:n_cond_factor2) {\ngamma_alpha[i,1] <-0 #prior for fixed effects\n}\nfor (i in 2:n_cond_factor1){\n  for (k in 2:n_cond_factor2){\n  gamma_alpha[k,i] ~dt(0,0.16,1) #prior for fixed effects\n  }\n}\n\n')}
  
  #add terms for hyperpriors
  model.string<-paste(model.string,'\n\n#hyper-priors\ntau_a ~ dgamma(0.01,0.01)\ntau_b ~ dgamma(0.01,0.01)\nsigma_a <-sqrt(1/tau_a)\nsigma_b <-sqrt(1/tau_b)\n')
  
  if (M==3)
  {model.string<-paste(model.string,'\ntau_alpha_b ~ dgamma(0.01,0.01)\nsigma_alpha_b <-sqrt(1/tau_alpha_b)\n')}
  else if (M==4)
  {model.string<-paste(model.string,'\ntau_alpha_a ~ dgamma(0.01,0.01)\nsigma_alpha_a <-sqrt(1/tau_alpha_a)\n')}
  else if (M==5)
  {model.string<-paste(model.string,'\ntau_alpha_b ~ dgamma(0.01,0.01)\nsigma_alpha_b <-sqrt(1/tau_alpha_b)\n \ntau_alpha_a ~ dgamma(0.01,0.01)\nsigma_alpha_a <-sqrt(1/tau_alpha_a)\n')}
  
  if (N==3)
  {model.string<-paste(model.string,'\ntau_gamma_b ~ dgamma(0.01,0.01)\nsigma_gamma_b <-sqrt(1/tau_gamma_b)\n')}
  else if (N==4)
  {model.string<-paste(model.string,'\ntau_gamma_a ~ dgamma(0.01,0.01)\nsigma_gamma_a <-sqrt(1/tau_gamma_a)\n')}
  else if (N==5)
  {model.string<-paste(model.string,'\ntau_gamma_b ~ dgamma(0.01,0.01)\nsigma_gamma_b <-sqrt(1/tau_gamma_b)\n \ntau_gamma_a ~ dgamma(0.01,0.01)\nsigma_gamma_a <-sqrt(1/tau_gamma_a)\n')}
  
  #add final closing bracket
  model.string<-paste(model.string,'\n}')
  
  final.list[[1]]<-model.string
  
  #second component
  #all models have these parameters
  variable.names<-c('beta0','a','b','sigma_a','sigma_b')
  
  if (M==2)
  {variable.names<-c(variable.names,'alpha')}
  else if (M==3)
  {variable.names<-c(variable.names,'alpha','alpha_b','sigma_alpha_b')}
  else if (M==4)
  {variable.names<-c(variable.names,'alpha','alpha_a','sigma_alpha_a')}
  else if (M==5)
  {variable.names<-c(variable.names,'alpha','alpha_b','sigma_alpha_b','alpha_a','sigma_alpha_a')}
  
  #add the rest of the regression structure for factor 1
  if (N==2)
  {variable.names<-c(variable.names,'gamma')}
  else if (N==3)
  {variable.names<-c(variable.names,'gamma','gamma_b','sigma_gamma_b')}
  else if (N==4)
  {variable.names<-c(variable.names,'gamma','gamma_a','sigma_gamma_a')}
  else if (N==5)
  {variable.names<-c(variable.names,'gamma','gamma_b','sigma_gamma_b','gamma_a','sigma_gamma_a')}
  
  #add the interaction 
  if (I==1)
  {variable.names<-c(variable.names,'gamma_alpha')}
  
 #add to final output
  final.list[[2]]<-variable.names  
 
 # create the data list
 data.list<-list('n_item'=n.item,'n_obs'=n.obs,'n_subject'=n.subject,'y'=y,'item'=item,'subject'=subject)
 if (M>=2)
 {data.list$n_cond_factor1 = n.cond.factor1
  data.list$condition.factor1=condition.factor1
 }
 if (N>=2)
 {
   data.list$n_cond_factor2 = n.cond.factor2
   data.list$condition.factor2=condition.factor2
 }
 
 final.list[[3]]<-data.list
 
 
 model.string.bic<-'Acc~'
 
 #get string for fixed effects structure
 if (I==1)
  {model.string.bic<-paste(model.string.bic,'cond1*cond2')}
 else
 {
  if ((M!=1)&&(N==1))  
  {model.string.bic<-paste(model.string.bic,'cond1')}
  else if((M==1)&&(N!=1))
  {model.string.bic<-paste(model.string.bic,'cond2')}
  else if((M!=1)&&(N!=1))
  {model.string.bic<-paste(model.string.bic,'cond1+cond2')}  
 }
 
 #add random effect structure for subjects
 if ((M!=3)&&(M!=5)&&(N!=3)&&(N!=5))
 {model.string.bic<-paste(model.string.bic,'+ (1|subj)')}
 else if (((M==5)||(M==3))&&(N!=3)&&(N!=5))
 {model.string.bic<-paste(model.string.bic,'+ (1+cond1|subj)')}
 else if((M!=3)&&(M!=5)&&((N==3)||(N==5)))
 {model.string.bic<-paste(model.string.bic,'+ (1+cond2|subj)')}
 else if(((M==5)||(M==3))&&((N==3)||(N==5)))
 {model.string.bic<-paste(model.string.bic,'+ (1+cond1+cond2|subj)')}
   
 #add random effect structure for items
 if ((M!=4)&&(M!=5)&&(N!=4)&&(N!=5))
 {model.string.bic<-paste(model.string.bic,'+ (1|itemID)')}
 else if (((M==5)||(M==4))&&(N!=4)&&(N!=5))
 {model.string.bic<-paste(model.string.bic,'+ (1+cond1|itemID)')}
 else if((M!=4)&&(M!=5)&&((N==4)||(N==5)))
 {model.string.bic<-paste(model.string.bic,'+ (1+cond2|itemID)')}
 else if(((M==5)||(M==4))&&((N==4)||(N==5)))
 {model.string.bic<-paste(model.string.bic,'+ (1+cond1+cond2|itemID)')}
 

 call.bic<-list()
 call.bic$formula<-as.formula(model.string.bic)
 call.bic$data<-m_data
 if(link=='Logit')
  {
  call.bic$family<-binomial(logit)
  }
 else{
   call.bic$family<-binomial(probit)
 }
 call.bic$nAGQ<-1
 call.bic$control<-glmerControl(optimizer="bobyqa") 
 final.list[[4]]<-call.bic
 
 final.list
}
		    M <- model_struct[1]
		    N <- model_struct[2]
		    I <- model_struct[3]
		    
		     n.iter.posterior<-20000
			 n.iter.WAIC<-10000
		    
#### Fit the two factor model		    
		    M.fit<-two.factor.model(link,M,N,I ,m_data)
###################FIT Model: setup, adapt, burnin, WAIC, posterior sampling of parameters
#initialize the model; set tuning parameters based on n.adapt runs
#M.fit[[1]] contains the jags code for the model
#M.fit[[3]] holds the data required
jags <- jags.model(textConnection(M.fit[[1]]),data = M.fit[[3]], n.chains = 1,n.adapt = 1000)
#run markov chain to burnin
update(jags, 1000)
### Compute WAIC
fit<-coda.samples(model=jags, variable.names=c('loglik'), n.iter= n.iter.WAIC)
log_lik<-fit[[1]]
S.waic<-nrow(log_lik)
n.waic<-ncol(log_lik)
lpd<-sum(log(colMeans(exp(log_lik))))
p_waic<-sum(apply(log_lik, 2, var))


model.waic<--2*(lpd-p_waic) ##THIS HOLDS WAIC FOR THE MODEL

rm(log_lik)
rm(fit)
###sample posterior: M.fit[[2]] is a list of all the parameters that are sampled in the mcmc

post<-coda.samples(model=jags,variable.names=M.fit[[2]], n.iter= n.iter.posterior)##THIS HOLDS POSTERIOR SAMPLES FOR MODEL

###############################################################################################Compute BIC

fit<-do.call(glmer,M.fit[[4]])
model.bic<-BIC(fit) ## THIS HOLDS THE BIC FOR THE MODEL


	return(list("bic" = model.bic, "waic" = model.waic, "post_summary" = post, "condition_level" = condition, "baseline" = baseline))
	}
	

}
