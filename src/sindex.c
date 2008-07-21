/* compute the values of a step function, 
   ie how many of the jumps are smaller or
   equal to the eval points  */

void sindex(int *index,
	    double *jump,
	    double *eval,
	    int *njump,
	    int *neval,
	    int *strict){
  int n,nt,i,t;
  
  n = *njump;
  
  nt = *neval;

  index[0] = 0;

  i = 0;
  if (*strict==0)
    for (t=0;t<nt;t++){
    
    while(jump[i]<=eval[t]
	  && i<n)
      i++;
    
    index[t] = i;
    }
  else
    for (t=0;t<nt;t++){
      
      while(jump[i] < eval[t]
	    && i<n)
	i++;
      
      index[t] = i;
  }
}

