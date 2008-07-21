#include <math.h>
#include <R.h>
#include "prodlim.h"

void prodlim_clustersurv(double *y,
			 int *status,
			 int *cluster,
			 int *NC,
			 double *time,
			 double *nrisk,
			 double *cluster_nrisk,
			 int *nevent,
			 int *loss,
			 int *ncluster_with_event,
			 int *sizeof_cluster,
			 int *nevent_in_cluster,
			 double *surv,
			 double *hazard,
			 double *varhazard,
			 double *adj1,
			 double *adj2,
			 double *adjvarhazard,
			 int *t,
			 int start,
			 int stop){
  
  int s,i,l,k;
  double surv_step, hazard_step, V1, V2, atrisk, cluster_atrisk;
  
  s = (*t);			/* set the counter of time points */
  

  /*  fill the clusters */
  for (i=start;i<stop;i++) {
    sizeof_cluster[cluster[i]-1]++;
  }
  
  /* initialize  */
  surv_step=1; hazard_step=0; V1=0; V2=0;

  atrisk=(double) stop-start;
  
  cluster_atrisk= (double) *NC;
  nevent[s] = status[start];
  ncluster_with_event[s] = status[start];
  nevent_in_cluster[cluster[start]-1] = status[start];
  loss[s] = (1-status[start]);
  
  for(i=(1+start);i <=stop;i++){ /* start is i=1 for the case with ties */
    
    if (y[i]==y[i-1] && i<stop){
      
      nevent[s] += status[i];
      loss[s] += (1 - status[i]);
      nevent_in_cluster[cluster[i]-1] += status[i];
      if (cluster[i]!=cluster[i-1])
	ncluster_with_event[s]+= status[i];
    }
    
    else {
      time[s] = y[i-1];
      nrisk[s] = atrisk;
      cluster_nrisk[s] = cluster_atrisk;

      /* marginal Kaplan-Meier and naive variance estimator */
      pl_step(&surv_step, &hazard_step, &V1, atrisk, nevent[s],0);
	
      surv[s]=surv_step;
      hazard[s]=hazard_step;
      varhazard[s] = V1;
	
      /* adjusted variance estimator of Ying and Wei (1994)  */
	
      V2=0;
      for (k=0;k<*NC;k++) {
	adj1[k] += nevent_in_cluster[k] / (double) atrisk;
	adj2[k] += sizeof_cluster[k] * nevent[s] / (double) (atrisk * atrisk);	
	V2 += (adj1[k]-adj2[k]) * (adj1[k]-adj2[k]);
      }	

      /* collect the results for unique time points */
      surv[s] = surv_step; 
      varhazard[s]=V1;
      adjvarhazard[s]=V2; 
      
      /* initialize the next time point */
      if (i < stop) {
	atrisk-=(nevent[s]+loss[s]);
	
	for (l=1;l<=(nevent[s]+loss[s]);l++) {
	  sizeof_cluster[cluster[i-l]-1]--;
	  if (sizeof_cluster[(cluster[i-l]-1)]==0) cluster_atrisk--; /* if the last obs in a cluster is gone ...  */
	  nevent_in_cluster[cluster[i-l]-1]=0; /* reset for next time point  */
	}
	  
	s++;
	nevent_in_cluster[cluster[i]-1] = status[i];
	nevent[s] = status[i];
	ncluster_with_event[s] = status[i];
	loss[s] = (1-status[i]);
      }
    }
  }
  *t=(s+1); /* for the next strata and finally for R */
}


/* void prodlim_cluster(double *y,int *status,int *cluster,int *nc,int *NU,int *size,double *time,double *nrisk,double *cluster_nrisk,int *nevent,int *loss,int *ncluster_with_event,int *sizeof_cluster,int *nevent_in_cluster,double *surv,double *hazard,double *varhazard,double *adj1,double *adj2,double *adjvarhazard,int *len,int *size_strata,int *first_strata){ */
/*   int t,u, start, stop, size_temp; */
/*   t=0; */
/*   start=0; */
/*   size_temp=0; */
/*   for (u=0;u<*NU;u++){ */
/*     stop=start+size[u]; */
/*     clustersurv_intern(y,status,cluster,nc,time,nrisk,cluster_nrisk,nevent,ncluster_with_event,loss,sizeof_cluster,nevent_in_cluster,surv,hazard,varhazard,adj1,adj2,adjvarhazard,&t,start,stop); */
/*     start+=size[u]; */
/*     size_strata[u] = t - size_temp; */
/*     first_strata[u] = t + 1 - size_strata[u]; */
/*     size_temp += size_strata[u]; */
/*   } */
/*   *len=t; */
/* } */
