data{
  int N;
  int N1;
  int N2;
  int n;
  int nbetas;
  vector[N] y;
  matrix[n,nbetas] X;
  int<lower=1,upper=n> ID1[N1];
  int<lower=1,upper=n> ID2[N2];
  vector[N1] times1;
  vector[N2] times2;
  vector[N2] treat_times;
  vector[n] time;
  vector[n] status;
  vector[3] theta;
  matrix[n,3] bi;
}


parameters{
  vector[nbetas] beta;
  real alpha;
  real<lower=0> gamma;
}

model{
// ------------------------------------------------------
//        LOG-LIKELIHOOD FOR SURVIVAL SUBMODEL                
// ------------------------------------------------------
{
  vector[n] logHaz;
  vector[n] cumHaz;

  for(i in 1:n){
       // log-Hazard
       logHaz[i] = log(gamma) + (gamma-1)*log(time[i]) + X[i,] * beta + alpha * exp(theta[2] + bi[i,2]);
       // Cumulative hazard H[t] = int_0^t h[u] du
       cumHaz[i] = (time[i]^gamma) * exp( X[i,] * beta + alpha * exp(theta[2] + bi[i,2]) );

       target += status[i] * logHaz[i] - cumHaz[i];
  }
}   
// ------------------------------------------------------
//                       LOG-PRIORS                       
// ------------------------------------------------------
   // Survival fixed effects
   target += normal_lpdf(beta | 0, 100);

   // Association parameter
   target += normal_lpdf(alpha | 0, 100);

   // Shape parameter (Weibull hazard)
   target += cauchy_lpdf(gamma | 0, 25);

}
