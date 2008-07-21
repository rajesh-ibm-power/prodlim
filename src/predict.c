void findex(int *findex,
	    int *type,
	    int *S,
	    int *freq_strata,
	    double *Z,
	    double *NN,
	    int *NR,
	    int *NT){
  
  int i,x,last,f;
  
  for (i=0;i<*NR;i++){
    
    /* goto strata of i */
    if (S[i]==1)
      x=0;
    else
      x = freq_strata[S[i]-2];
    last = freq_strata[S[i]-1] -1;
    
    /* find the closest neighbor */
    if (*type==0)
      findex[i]=x;
    else{
      if (Z[i] <= NN[x]) /* <= first */
	findex[i] = x; 
      else{
	if (Z[i] >= NN[last]){/* >= last */
	  findex[i] = last;
	}
	else { /* sitting between two neighbors*/
	  while (Z[i] >= NN[x]) x++;
	  if ((NN[x] - Z[i]) < (Z[i] - NN[x-1]))
	    findex[i] = x;
	  else
	    findex[i] = x-1;
	}
      }
    }
    findex[i]+=1; /* in `R' counting starts at 1 */
  }
}

void pred_index(int *pindex,
		double *Y,
		double *time,
		int *first,
		int *size,
		int *NR,
		int *NT){
  
  int i,t,f;
    
  for (i=0;i<*NR;i++){    
    f=0;
    for (t=0;t<(*NT);t++){
      
      if (Y[t] < time[first[i]-1]){ /* < first */
	pindex[t + i * (*NT)] = 0;
      }
      else{
	if (Y[t] > time[first[i]-1 + size[i]-1]){ /* > last */
	  while(t<(*NT)){ 
	    pindex[t + i * (*NT)] = -1;
	    t++; 
	  } 
	}
	else{ /* sitting between to jump times */
	  
	  while (Y[t] >= time[first[i]-1 + f]
		 && f <= size[i]-1)
	    f++;
	  pindex[t + i * (*NT)] = first[i] -1 + f;
	  /* do NOT reset f because the next requested time
	     is greater or equal to the current time */
	}
      }
    }
  }
}


void predict_individual_survival(double *pred,
				 double *surv,
				 double *jumptime,
				 double *Y,
				 int *first,
				 int *size,
				 int *n,
				 int *lag){
  int j,i; /* start at index 0 */

  /* predicted survival probabilities at or just before the
     individual event times Y[i] */  
  for (i=0;i<(*n);i++){
    j=0;
    /* index j is in stratum i if j < size[i] */
    while(j < size[i] - 1 &&
	  jumptime[first[i] - 1 + j] != Y[i])
      j++;
    if (j - *lag < 0)
      pred[i]=1;
    else
      pred[i] = surv[first[i] - 1 + j - *lag];
  }

}
void life_table(int *pred_nrisk,
		int *pred_nevent,
		int *pred_nlost,
		int *nrisk,
		int *nevent,
		int *nlost,
		double *Y,
		double *jump,
		int *first,
		int *size,
		int *NR,
		int *NT,
		int *intervals){
  
  int i,t,f,count_e,count_l;
  double min_jump, max_jump;

  /*  count events, right censored (lost) and numbers at risk
      at the given time points or when intervals=1 in between
      the given time points. the requested time points are given
      by Y. the event times are given by jump.
  */
  for (i=0;i<*NR;i++){    
    min_jump = jump[first[i]-1];
    max_jump = jump[first[i]-1 + size[i]-1];
    f=0;
    for (t=0;t<(*NT);t++){
      count_e =0;
      count_l =0;
      if (Y[t] < min_jump){ /*
			      first set everything to zero
			       at all requested times that are 
			       before the smallest event time
			    */
	pred_nrisk[t + i *(*NT)] = nrisk[first[i]-1];
	pred_nevent[t + i *(*NT)] = 0;
	pred_nlost[t + i *(*NT)] = 0;
      }
      else{
	if (Y[t] > max_jump && t==0 || Y[t] > max_jump && Y[t-1] >= max_jump){
	  /*
	    once the current time and the one before are after the last jump
	     nothing may happen at later times and hence all is zero. do the same
	     if all requested times are after the last jump.
	  */
	  while(t<(*NT)){ 
	    pred_nrisk[t + i *(*NT)] = 0;
	    pred_nevent[t + i *(*NT)] = 0;
	    pred_nlost[t + i *(*NT)] = 0;
	    t++; 
	  } 
	}
	else{ /*
		for all jumps that are smaller or equal to
		the currently requested time count
		events and losts.
	      */
	  while (Y[t] >= jump[first[i]-1 + f] 
		 && f <= size[i]-1){
	    count_e +=nevent[first[i] -1 +f];
	    count_l +=nlost[first[i] -1 +f];
	    f++;
	  }
	  
	  pred_nevent[t + i *(*NT)] = count_e;
	  pred_nlost[t + i *(*NT)] = count_l;
	  
	  if (Y[t] > max_jump) {
	    pred_nrisk[t + i *(*NT)] = 0;
	  }
	  else{
	    if (Y[t]==jump[first[i] -1 -1 +f]){ /* looking back */
	      if (*intervals!=1)
		pred_nrisk[t + i *(*NT)] = nrisk[first[i] -1 +f -1];
	      else{
		if (jump[first[i] -1 +f -1]>=max_jump)
		  pred_nrisk[t + i *(*NT)] = 0;
		else
		  pred_nrisk[t + i *(*NT)] = nrisk[first[i] -1 +f -1];
	      }
	    }
	    else{ /* really inbetween */
	      pred_nrisk[t + i *(*NT)] = nrisk[first[i] -1 +f];
	    }
	  }
	  /* do NOT reset f, count_l and count_e because the next requested time
	     is greater or equal to the current time. */
	}
      }
    }
  }
}


