#include <math.h>
#include <R.h>

#define max(A,B) ((A) > (B) ? (A):(B))
#define min(A,B) ((A) < (B) ? (A):(B))

void icens_prodlim(double *L,
		   double *R,
		   double *grid,
		   int *indexL,
		   int *indexR,
		   int *iindex,
		   int *imax,
		   int *status,
		   double *N,
		   double *NS,
		   double *nrisk,
		   double *nevent,
		   double *ncens,
		   double *hazard,
		   double *var_hazard,
		   double *surv,
		   double *oldsurv,
		   double *tol,
		   int *maxstep,
		   int *niter) {
  
  int i, j, s, done=0, step=0, n, ns, start, stop;
  double atrisk, pl, haz, varhaz, diff, survL, survR, lenOBS;

  
  n = (int) *N; /* number of interval censored observations */
  ns = (int) *NS;  /* number of grid points + 1 */

  /*
    loop until convergence or the maximum stepsize reached
  */
  while (done==0 && step < *maxstep){
    /*
      initialize some parameters before each step
    */
    diff=0;
    atrisk = *N;
    pl=1; 
    haz=0;
    varhaz=0;
    nevent[0] = 0;
    ncens[0] = 0;
    start=0;
    stop=max(0,imax[0]);
    /*
      loop over all grid points
      starting with the interval
      [grid[0] ; grid[1]]
    */ 
    for (s=0; s < (ns-1); s++){

      /* initialize the survival probability at the */
      /* first grid point and the number at risk */
      surv[s]=1;
      nrisk[s]=atrisk;
      
      for (j=start; j < stop; j++){ 
	/*
	  use only those intervals that overlap the current
	  grid interval: [grid[s] ; grid[s+1]]
	  the overlap is determined by iindex
	*/
	i=iindex[j]-1;
	/*
	  right censored observations are handled
	  as for the usual Kaplan-Meier method
	  here they are counted at the right end grid[s+1] 
	*/
	if (status[i]==0 && L[i] == grid[s+1]) ncens[s]++;
	/*
	  interval censored observations
	*/
	if (status[i]>0){
	  lenOBS = R[i] - L[i];
	  /*
	    are either exact observations and  
	    counted at the right end grid[s+1] of
	    the current interval
	  */
	  if (lenOBS==0 && L[i] == grid[s+1]) nevent[s] ++; 
	  /*
	    or really interval censored and contribute
	    relative to the overlap with the current
	    grid-interval
	  */
	  if (lenOBS > 0){
	    /*
	      at the very first step
	      assume a uniform distribution
	    */
	    if (step==0){
	      nevent[s] += max(0,min(R[i],grid[s+1]) - max(grid[s],L[i]))/lenOBS;
	    }
	    else{
	      /*
		at subsequent steps
		use the survival probability obtained
		in the previous step
	      */
	      if (indexL[i]<=1)
		survL=1; else survL=surv[indexL[i]-2];
	      
	       if (indexR[i]<=1)
		 survR=1; else survR=surv[indexR[i]-2];
	       
	       /* if (s==0) survG=1; else survG=surv[s-1]; */
	       nevent[s] += max(0,min(survL,surv[s]) - max(surv[s+1],survR))/(survL-survR);
	    }
	  }
	}
      }
      start=max(0,imax[s]);
      stop=max(imax[s+1],0);
      if (nevent[s]>0){
	haz = nevent[s] / (double) atrisk;
	pl*=(1 - (nevent[s] / (double) atrisk)); 
	varhaz += nevent[s] / (double) (atrisk * (atrisk - nevent[s])); 
      }
      /*  store the current estimate in oldsurv */
      if (step>0) oldsurv[s]= surv[s];
      /*
	update the survival probability
	and the hazard at grid s
      */
      surv[s]=pl;
      hazard[s] = haz;
      var_hazard[s] = varhaz;
      /*
	update the number at risk
      */
      atrisk-=(nevent[s]+ncens[s]);
      /*
	initialize the event and censored counts
      */
      nevent[s+1] = 0;
      ncens[s+1] = 0;
    }
    /*
      check if the algorithm converged
    */
    for (s=0;s<(ns-1);s++){
      diff=max(max(surv[s]-oldsurv[s],oldsurv[s]-surv[s]),diff);
    }
    if (diff < *tol) done=1;
    step++;
  }
  /*     Rprintf("Step %d\n",step); */
  niter[0]=step; 
}



