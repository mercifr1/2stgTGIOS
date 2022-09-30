functions{
// ------------------------------------------------------
//     NONLINEAR PREDICTOR FOR THE LONGITUDINAL SUBMODEL                
// ------------------------------------------------------ 
    vector nonlinear_predictor(int[] ID1, int[] ID2, vector times1, vector times2, vector treat_times, vector theta, matrix bi){
         int N1 = num_elements(times1);
         int N2 = num_elements(times2);
         int N = N1 + N2;
         vector[N1] omega1 = exp(theta[1] + bi[ID1,1]);
         vector[N2] omega2 = exp(theta[1] + bi[ID2,1]);
         vector[N1] g1 = exp(theta[2] + bi[ID1,2]);
         vector[N2] g2 = exp(theta[2] + bi[ID2,2]);
         vector[N2] d2 = exp(theta[3] + bi[ID2,3]);
         vector[N2] part1 = exp(-rows_dot_product( d2, times2-treat_times ));
         vector[N2] part2 = exp(rows_dot_product( g2, times2-treat_times ));
         vector[N1] out1;
         vector[N2] out2;
         vector[N] out;

         out1 = rows_dot_product( omega1, exp( rows_dot_product(times1, g1) ) );
         out2 = rows_dot_product( rows_dot_product( omega2, exp( rows_dot_product(treat_times, g2) ) ), part1 + part2 - 1 );

         out = append_row(out1, out2);
         return out;
    }
// ------------------------------------------------------ 
}


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
  real<lower=0> Var_b[3];
  real<lower=0> Var_e;
}


parameters{  
  matrix[n,3] bi;
  vector[nbetas] beta;
  real alpha;
  real<lower=0> gamma;
}

model{
// ------------------------------------------------------
//        LOG-LIKELIHOOD FOR LONGITUDINAL SUBMODEL                
// ------------------------------------------------------
{
   vector[N] nonlinpred; 

   // Nonlinear predictor
   nonlinpred = nonlinear_predictor(ID1, ID2, times1, times2, treat_times, theta, bi);

   // Longitudinal Normal log-likelihood
   target += normal_lpdf(y | nonlinpred, nonlinpred * sqrt(Var_e));
}
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
   // Random-effects
   target += normal_lpdf(bi[1:n,1] | 0, sqrt(Var_b[1])); 
   target += normal_lpdf(bi[1:n,2] | 0, sqrt(Var_b[2]));
   target += normal_lpdf(bi[1:n,3] | 0, sqrt(Var_b[3]));

   // Survival fixed effects
   target += normal_lpdf(beta | 0, 100);

   // Association parameter
   target += normal_lpdf(alpha | 0, 100);

   // Shape parameter (Weibull hazard)
   target += cauchy_lpdf(gamma | 0, 25);

}
