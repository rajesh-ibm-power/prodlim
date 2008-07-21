#library(dental,lib.loc='/home/bs/tag/R/library')
#library(prodlim,lib.loc='/home/bs/tag/R/library')
## source('/home/bs/ragr/ThomasGerds/ExtendedModel/PetoInt.R')
## source('/home/bs/ragr/ThomasGerds/ExtendedModel/IntIndextirsdag.R')
#source('/home/bs/ragr/ThomasGerds/IntIndex.R')
## compGMLE arguments:
##
## L and R are left and right end-points of the interval censored observations
## Exact observations L=R.
##
## Status: denotes right censored by 0, 
## exact time by 1 and interval censored by 2. 
##
## Cause is a vector of 0,1,2 denoting the type of event. 0 is right censoring,
## 1 is fracture failur, 2 is exfoliation.
##
## maxstep: is the maximum number of iterations of the estimating algorithm.
##
## tol: denotes the tolerance level 10^{-tol}. If there is less than 10^{-tol} in difference between two succesive estimates then the algorithm is stoppped.
## 
## An initial value of the GMLE (Z) are choosen to be a vector of length as the number of peto-intervals (M)  and with equal mass (1/M) in each interval.

## compGMLE returns:
##
## Int: Peto-intervals found by using the PetoInt function.
## Z: Generalised Maximum Likelihood Estimator.
## steps: Number of iterations.

CompGMLEdental <- function(L,R,cause, maxstep=1000,tol=4, competing=FALSE){
  ## creating a status variable indicating
  ## if the obs are excact,interval censored ore right censored.
  status <- rep(NA,length(L))
  status[(R-L)>0] <- 2
  status[(R-L)==0] <- 1
  status[cause==0] <- 0
  
  oL <- order(L)
  L <- L[oL]
  R <- R[oL]
  status <- status[oL]
  cause <- cause[oL]
  R[status==0] <- max(L,R[status!=0],na.rm=T)+1 #to ensure that the right censored observation behaves as [s_c;inf)

  dyn.load('/home/bs/ragr/ThomasGerds/ExtendedModel/indicator.so')
  dyn.load('/home/bs/ragr/ThomasGerds/ExtendedModel/sumprod.so')


  N <- length(L) #number of observations
  ntol <- 10^{-tol}
  s_max <- max(L[status==0])-1  # the int Q_{I'+1} is defined as [s_max;inf]

  ##############################
  ## Possible support points ###
  ##############################
  
  #Construction of the O and A set here
  Lm <- L[cause==1]
  Rm <- R[cause==1]

  A <- cbind(Lm,Rm) 
  A <- unique(unlist(as.vector(apply(A,1,function(x){min(x):max(x)})))) # takes Lm:Rm for all m and then contruct a set by using unique funtion

  Bk <- L[cause==2]
  Ek <- R[cause==2]

  O <- cbind(Bk,Ek-1)
  O <- unique(unlist(as.vector(apply(O,1,function(x){min(x):max(x)}))))

  AuO <- sort(as.vector(as.numeric(union(A,O))))
  print(class(AuO))
  print(AuO)
  #Constructs the L01 and R01 sets and returns the possible support points of transition 01.
  if(s_max>max(c(Rm,Ek-1))){Support01 <- as.matrix(rbind(c(AuO,s_max),c(AuO,max(R))))}
  else{Support01 <- as.matrix(rbind(AuO,AuO))}

    #Construction of the D set here:
  D <- cbind(Bk+1,Ek) #the Dk sets are open on the left 
  D <- unique(unlist(as.vector(apply(D,1,function(x){min(x):max(x)})))) # takes Bk:Ek for all k and then contruct a set by using unique funtion
  

  Sc <- L[status==0]
  unionSc <- Sc[Sc%in%D]-1 #right censored observations which is included in at least one of the (B_k,E_k} intervals

  #constructs the two sets L_02 and R_02. 
  if( s_max>max(R[cause==2])){ L02 <- c(L[cause==2]-1,unionSc,s_max)
                               R02 <- c(R[cause==2],rep(max(L)+1,length(unionSc)+1))
                               oL02 <- order(L02)
                               status02 <- c(status[cause==2],rep(0,length(unionSc)+1)) } 
  else{L02 <- c(L[cause==2]-1,unionSc)
       R02 <- c(R[cause==2],rep(max(R)+1,length(unionSc)))
       oL02 <- order(L02)
       status02 <- c(status[cause==2],rep(0,length(unionSc)))}

  L02 <- L02[oL02]
  R02 <- R02[oL02]
  status02 <-status02[oL02]
  Support02 <- as.matrix(PetoInt(L02,R02,status02)) #This returns the Q_i, I<i<=I'


  NoInt <- ncol(Support01)+ncol(Support02) #number of intervals with potential increase in the subdistributionfunctions

  #Indicator functions, make c function that returns relevatn indicator functions
  I <- ncol(Support01)
  Imark <- NoInt-I

  M <- length(L[cause==1])#Number of filling fracture
  K <- length(L[cause==2])#Number of exfoliations
  C <- length(L[cause==0])#Number of right censored observations

  #vectors to contain indicators for the intervals.
  ind1 <- rep(0,I*M) 
  ind2 <- rep(0,I*K)
  ind3 <- rep(0,Imark*K)
  ind4 <- rep(0,NoInt*C)

  beta  <- .C('Indicator',
              as.integer(I),
              as.integer(M),
              as.integer(Lm),
              as.integer(Rm),
              as.integer(Support01[1,]),
              as.integer(Support01[1,]),
              ind=as.integer(ind1))
              
  beta <- beta$ind


  
  gamma  <- .C('Indicator',
              as.integer(I),
              as.integer(K),
              as.integer(Bk),
              as.integer(Ek-1),
              as.integer(Support01[1,]),
              as.integer(Support01[1,]),
              ind=as.integer(ind2))
              
  gamma <- gamma$ind

  
  eta  <- .C('Indicator',
             as.integer(Imark),
             as.integer(K),
             as.integer(Bk),
             as.integer(Ek),
             as.integer(Support02[1,]),
             as.integer(Support02[2,]),
             ind=as.integer(ind3))
              
  eta <- eta$ind

  Sc <- L[cause==0]
  ScEnd <- R[cause==0]
  
  alpha  <- .C('Indicator',
             as.integer(NoInt),
             as.integer(C),
             as.integer(Sc),
             as.integer(ScEnd),
             as.integer(c(Support01[1,],Support02[1,])),
             as.integer(c(Support01[1,],Support02[2,])),
             ind=as.integer(ind4))
              
  alpha <- alpha$ind


################################
#### The estimation algorithm.##
################################

  maxD <- max(Ek,D)     #to ensure that the the vector lambda is long enough in the sumprod function.
  Dvector <- rep(0,maxD) #Bliver denne brugt??

  
  lambda <- rep(0,maxD)
    
  if(competing==FALSE){Ldenom <- rep(0,K)#vektor for the L_k values that are calculated in eache step in the algorithm. THIS IS NOT LEFT ENDPOINTS
                       lambda[D] <- 0.5  #initial value for the lambda is set to 0.5 only for values in the D set.
                       comp <- 0}        #variable used in c function to distinguish the two models
  else{Ldenom <- rep(0,K)                
       comp <- 1}

  Zold <- rep(1/NoInt,NoInt) #Initial value of z  
  Z <- rep(1/NoInt,NoInt)    #vector to contain the resulting z vector.

  print(Ldenom)
  print(maxD)
  

  dental <- .C('DentalGMLE',
               as.integer(comp),
               as.double(Zold),
               as.double(lambda),
               Lambda=as.double(lambda),
               Ldenom=as.double(Ldenom),
               as.integer(K),
               as.integer(M),
               as.integer(C),
               as.integer(NoInt),
               as.integer(I),
               as.integer(maxD),
               as.integer(Ek),
               as.integer(Rm),
               as.integer(Support01[1,]),#There is not taken care of the situation where the last r_i is an interval and not a singelton
               as.integer(eta),
               as.integer(gamma),
               as.integer(beta),
               as.integer(alpha),
               as.integer(maxstep),
               iter=integer(1),
               as.double(ntol),
               z=as.double(Z))

  z01 <- dental$z[1:I]
  z02 <- dental$z[(I+1):NoInt]
    
     
  out <- list(z01=z01,z02=z02,Iterations=dental$iter,res=dental$Ldenom,Support01=Support01,Support02=Support02, beta=beta, gamma=gamma, eta=eta,
              alpha=alpha,D=D, lambda=dental$Lambda, lambdastart=lambda, obsordered=rbind(L,R,cause))
  out
}


## L <- c(2,3,5,4,4,5)
## R <- c(5,6,8,7,20,20)
## status <- c(2,2,2,2,0,0)
## cause <- c(1,1,2,2,0,0)
## out <- CompGMLEdental(L,R,cause, maxstep=1)
## out
## sum(out$z01)+sum(out$z02)


## out2 <- CompGMLEdental(L,R,status,cause, competing=TRUE)
## sum(out2$z01)+sum(out2$z02)



## L <- c(2,3,1,4,4,5)
## R <- c(5,6,8,7,20,20)
## status <- c(2,2,2,2,0,0)
## cause <- c(1,1,2,2,0,0)
## out <- CompGMLEdental(L,R,cause, maxstep=1)
## out
## sum(out$z01)+sum(out$z02)

## L <- c(4,3,8,4,4,5)
## R <- c(5,6,8,6,20,20)
## status <- c(2,2,2,2,0,0)
## cause <- c(1,1,2,2,0,0)
## out <- CompGMLEdental(L,R,cause, maxstep=1)
## out
## sum(out$z01)+sum(out$z02) #does not sum to one. there is missing 2/6 in the sum.
## #can this be because the contribution from the Bk;Ek intervals are 'removed'.









## L <- c(2,3,4,5,4,9)
## R <- c(5,6,7,8,20,20)
## status <- c(2,2,2,2,0,0)
## cause <- c(1,1,2,2,0,0)
## out <- CompGMLEdental(L,R,status,cause, maxstep=1)
## out
## sum(out$z01)+sum(out$z02)

## ############## Problems!!!!!!!!!!!!!!!!

## L<-c(1,4,6,8,4,7,6) 
## R<-c(6,5,7,9,20,20,9)
## status<-c(2,2,2,2,0,0,2)
## cause <- c(1,1 ,2 ,2 ,0 ,0 ,1)
## out <- CompGMLEdental(L,R,status,cause,maxstep=1)
## out
## sum(out$z01)+sum(out$z02)


## L<-c(1,4,3,8,6,5,4,10,6) 
## R<-c(6,5,3,9,7,5,20,20,9)
## status<-c(2,2,2,2,2,2,0,0,2)
## PetoInt(L,R,status)
## out <- CompGMLEdental(L,R,status,c(1,2 ,1 ,1 ,2 ,2 ,0 ,0 ,1))
## out

