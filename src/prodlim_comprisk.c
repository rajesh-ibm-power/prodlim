#include <math.h>
#include <R.h>
#include "prodlim.h"

void prodlim_comprisk(double* y, 
		      int* status,
		      int* cause,
		      int* NS, /* number of causes (states) */
		      double* time,
		      double* nrisk,
		      int* event, 
		      int* loss, 
		      double* surv,
		      double* cuminc,
		      double* cause_hazard,
		      double* varcuminc,
		      double* I, /* current cumulative incidence */ 
		      double*I_lag, /* time lagged cumulative incidence */ 
		      double* v1,
		      double* v2,
		      int *t,
		      int start,
		      int stop) {
  
  int i,j,s,d,d1,d2;
  double S,S_lag,H,varH,n;
  
  /*   initialize */
  s=(*t);
  S=1;
  H=0;
  for(j=0; j < (*NS); ++j) {
    I[j]=0;
    I_lag[j]=0;
    v1[j]=0;
    v2[j]=0;
  }
  varH=0;
  n=(double) stop-start; /* (sub-)sample size */ 
  
  if (status[start]>0)
    event[s *(*NS) + cause[start]]=1;
  else
    loss[s]=1;
  
  for (i=(1+start);i<=stop;i++){
    if (y[i]==y[i-1] && i<stop){
      if (status[i]>0)
	event[s * (*NS) + cause[i]] +=1;
      else
	loss[s]+=1;
    }
    else {
      time[s]=y[i-1];
      nrisk[s]=n;
      d = 0;
      for(j=0; j < (*NS); ++j) {
	/* 	Rprintf("i=%d\ts=%d\tcause=%d\teventj=%d\teventi=%d\n",i,s,cause[i],event[s * (*NS) + j],event[s * (*NS) + cause[i]]);   */
	cause_hazard[s * (*NS) + j] = (double) (event[s * (*NS) + j] / n);
	I_lag[j] = I[j];
	I[j] += S * cause_hazard[s * (*NS) + j];
	cuminc[s * (*NS) + j] = I[j];
	d += (double) event[s * (*NS) + j];
      }
      S_lag = S;
      pl_step(&S, &H, &varH, n, d,0);
      
      /* variance estimate Marubini & Valsecchi (1995), Wiley, chapter 10, page 341 */      
      for (j=0; j < (*NS); ++j){

	d1 = event[s * (*NS) + j];
	d2 = d - d1;
	
	v1[j] += I[j] * (d / (n * (n - d))) + (S_lag * d1) / (n * n);
	v2[j] += (I[j] * I[j]) * (d / (n * (n - d)))
	  + ((S_lag * S_lag) * (n - d1) * d1) / (n * n * n)
	  + (2 * I[j] * S_lag * d1) / (n * n);
	varcuminc[s * (*NS) + j] = (I[j] * I[j]) *  varH - 2 * I[j] * v1[j] + v2[j];
	/* varH is greenwood's formula */
	
	/* Rprintf("d=%d\tn=%1.0f\tt=%1.0f\tI_t=%1.4f\tv1=%1.2f\tv2=%1.2f\tS=%1.3f\tvar=%1.5f\n",d,n,time[s],I[j],v1[j],v2[j],S_lag,varcuminc[s * (*NS) + j]);    */
	/* variance estimate Korn & Dorey (1992), Stat in Med, Vol 11, page 815 */
	/* I1 = (I[j] - I_lag[j]) / 2; */

      }
      surv[s] = S;
      
      if (i<stop){
	n -= (d + loss[s]);
	s++;
	if (status[i]>0){
	  event[s *(*NS) + cause[i]]=1;
	  /* 	  Rprintf("i=%d\ts=%d\tcause=%d\tevent=%d\n",i,cause[i],event[s * (*NS) + cause[i]]);  */
	}
	else
	  loss[s]=1;
      }
    }
  }
  *t=(s+1); /* for the next strata  */
}


