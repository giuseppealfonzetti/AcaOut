#ifndef cr_H
#define cr_H
#include "extractParams.h"


namespace cr
{

  //' Evaluate hazard function based on outcome and year
   //'
   //' @param OUTCOME 1 for graduation, 2 for dropout, 3 for transfer
   //' @param YEAR Possible values 1:YB in case of dropout/transfer;
   //' @param THETA Parameter vector
   //' @param COVARIATES The last 2 values refers to ability and speed respectively. Remaining values are external covariates
   //' @param YB Maximum number of years allowed before graduation.
   //' @param LOGFLAG Set TRUE to return log value.
   //' @returns It returns the hazard probability of the specific outcome and year.
  double hazard(
    const int OUTCOME,
    const int YEAR,
    const Eigen::Ref<const Eigen::VectorXd> THETA,
    const Eigen::Ref<const Eigen::VectorXd> COVARIATES,
    const int YB,
    const bool LOGFLAG = false
  ){

    const int q = COVARIATES.size();
    Eigen::VectorXd theta_cr = THETA.tail(2*(YB+q)+1);
    double out;
    if((OUTCOME == 2 | OUTCOME == 3) && YEAR > YB) Rcpp::stop("`YEAR` larger than `YB`");

    if(OUTCOME == 2 | OUTCOME == 3){
      const double int_d = theta_cr(YEAR);
      const double int_t = theta_cr(YEAR+YB+q);
      const Eigen::VectorXd beta_d = theta_cr.segment(YB+1, q);
      const Eigen::VectorXd beta_t = theta_cr.segment(2*YB+q+1, q);
      const double eta_d = int_d + beta_d.dot(COVARIATES);
      const double eta_t = int_t + beta_t.dot(COVARIATES);

      if(OUTCOME == 2){
        out = eta_d-log1pexp(R::logspace_add(eta_d, eta_t));
      }else if(OUTCOME == 3){
        out = eta_t-log1pexp(R::logspace_add(eta_d, eta_t));
      }
    }

    if(OUTCOME == 1){
      const double int_g = theta_cr(0);
      out = int_g - log1pexp(int_g);
    }

    if(LOGFLAG){
      return(out);
    }else{
      return(exp(out));
    }
  }


  //' Evaluate survival function given a the range of years of interest
  //'
  //' @param YEAR_LAST Last year to evaluate.
  //' @param THETA Parameter vector
  //' @param COVARIATES The last 2 values refers to ability and speed respectively. Remaining values are external predictors.
  //' @param YB Total number of years in the non-graduatable regime. Needed for determining how many time-related intercepts.
  //' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time
  //' @param LOGFLAG Set TRUE to return log value.
  //'
  //' @returns It returns the probability of survival from `YEAR FIRST` to `YEAR_LAST` included.
  double survival(
    const int YEAR_FIRST,
    const int YEAR_LAST,
    const Eigen::Ref<const Eigen::VectorXd> THETA,
    const Eigen::Ref<const Eigen::VectorXd> COVARIATES,
    const int YB,
    const int YEAR_LAST_EXAM = 100,
    const bool LOGFLAG = false
  ){

  double logout = 0;
  // double out = 1;
  if(YEAR_LAST_EXAM > YEAR_LAST){

    // Regime where graduation is not possible
    if(YEAR_LAST > YB) Rcpp::stop("`YEAR_LAST` > `YB`");

    // Remain enrolled until YEAR_LAST (conditioned on not having all exams)
    for(unsigned int year = YEAR_FIRST; year<=YEAR_LAST; year++){
      logout += log1mexp(-R::logspace_add(hazard(2, year, THETA, COVARIATES, YB, true),
                                          hazard(3, year, THETA, COVARIATES, YB, true)));
    }

  }else if(YEAR_LAST_EXAM <= YEAR_LAST){

    // Regime where graduation is possible from year `YEAR_LAST_EXAM`

    // Remain enrolled until YEAR_LAST_EXAM
    for(unsigned int year = YEAR_FIRST; year < YEAR_LAST_EXAM; year++){
      logout += log1mexp(-R::logspace_add(hazard(2, year, THETA, COVARIATES, YB, true),
                                          hazard(3, year, THETA, COVARIATES, YB, true)));
    }

    // Remain enrolled from YEAR_LAST_EXAM to YEAR_LAST
    for(unsigned int year = YEAR_LAST_EXAM; year <= YEAR_LAST; year++){
      logout += log1mexp(-hazard(1, year - YEAR_LAST_EXAM + 1, THETA, COVARIATES, YB, true));
    }

  }

  logout=std::max(-10000.0, logout);
  if(LOGFLAG){
    return(logout);
  }else{
    return(exp(logout));
  }
 }

  //' Evaluate Outcome Likelihood
  //'
  //' @param OUTCOME  `1` for graduation, `2` for dropout, `3` for transfer. `0` if no outcome is observed.
  //' @param YEAR_LAST Last year to evaluate.
  //' @param THETA PParameter vector
  //' @param COVARIATES The first 2 values refers to ability and speed respectively. Remaining values are external predictors.
  //' @param YB Total number of years in the non-graduatable regime. Needed for determining how many time-related intercepts.
  //' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time.
  //' @param LOGFLAG Set TRUE to return log value.
  //'
  double outcomeLik(
      const unsigned int OUTCOME,
      const unsigned int YEAR_FIRST,
      const unsigned int YEAR_LAST,
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const Eigen::Ref<const Eigen::VectorXd> COVARIATES,
      const unsigned int YB,
      const unsigned int YEAR_LAST_EXAM = 100,
      const bool LOGFLAG = false
  ){
    double logout;

    if(OUTCOME==0){
      // The student is still enrolled
      logout = survival(YEAR_FIRST, YEAR_LAST, THETA, COVARIATES, YB, YEAR_LAST_EXAM, true);
    }else if(OUTCOME==2|OUTCOME==3){
      if((YEAR_LAST_EXAM<=YEAR_LAST) & !LOGFLAG) return 0;
      if((YEAR_LAST_EXAM<=YEAR_LAST) & LOGFLAG) return -10000;//return R_NegInf;

      //The student remain enrolled until YEAR_LAST-1 and the experience OUTCOME
      logout  = survival(YEAR_FIRST, YEAR_LAST-1, THETA, COVARIATES, YB, YEAR_LAST_EXAM, true);
      logout += hazard(OUTCOME, YEAR_LAST, THETA, COVARIATES, YB, true);
    }else if(OUTCOME==1){
      if((YEAR_LAST_EXAM>YEAR_LAST) & !LOGFLAG) return 0;
      if((YEAR_LAST_EXAM>YEAR_LAST) & LOGFLAG) return -10000;//return R_NegInf;

      //The student remain enrolled until YEAR_LAST-1 and the experience OUTCOME
      logout  = survival(YEAR_FIRST, YEAR_LAST-1, THETA, COVARIATES, YB, YEAR_LAST_EXAM, true);
      logout += hazard(1, YEAR_LAST-YEAR_LAST_EXAM+1, THETA, COVARIATES, YB, true);
    }

    if(LOGFLAG){
      return(logout);
    }else{
      return(exp(logout));
    }
  }
}

namespace cr::grad{
//' Evaluate hazard function based on outcome and year
 //'
 //' @param OUTCOME 1 for graduation, 2 for dropout, 3 for transfer
 //' @param YEAR Possible values 1:YB in case of dropout/transfer;
 //' @param THETA Parameter vector
 //' @param COVARIATES The last 2 values refers to ability and speed respectively. Remaining values are external covariates
 //' @param YB Maximum number of years allowed before graduation.
 //' @param LOGFLAG Set TRUE to return log value.
 //' @returns It returns the hazard probability of the specific outcome and year.
 //'
 Eigen::VectorXd gr_hazard(
     const int OUTCOME,
     const int YEAR,
     const Eigen::Ref<const Eigen::VectorXd> THETA,
     const Eigen::Ref<const Eigen::VectorXd> COVARIATES,
     const int YB,
     const bool LATPARFLAG
 ){

   const int q = COVARIATES.size();
   const int dim_cr = 2*(YB+q)+1;
   const int dim_irtlat = THETA.size()-dim_cr;
   Eigen::VectorXd theta_cr = THETA.tail(dim_cr);
   Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA.size());
   double lpr;
   if((OUTCOME == 2 | OUTCOME == 3) && YEAR > YB) Rcpp::stop("`YEAR` larger than `YB`");

   if(OUTCOME == 2 | OUTCOME == 3){
     const double int_d = theta_cr(YEAR);
     const double int_t = theta_cr(YEAR+YB+q);
     const Eigen::VectorXd beta_d = theta_cr.segment(YB+1, q);
     const Eigen::VectorXd beta_t = theta_cr.segment(2*YB+q+1, q);
     const double eta_d = int_d + beta_d.dot(COVARIATES);
     const double eta_t = int_t + beta_t.dot(COVARIATES);

     double logConstNorm = log1pexp(R::logspace_add(eta_d, eta_t));
     double lpr_d = eta_d-logConstNorm;
     double lpr_t = eta_t-logConstNorm;
     double nprod = -exp(lpr_d+lpr_t);

     if(OUTCOME == 2){
       const double tmp = exp(lpr_d) - pow(exp(lpr_d), 2);
       gr(dim_irtlat+YEAR) = tmp;
       gr.segment(dim_irtlat+YB+1, q) = tmp*COVARIATES;
       gr.segment(dim_irtlat+2*YB+q+1, q) = nprod*COVARIATES;
       gr(dim_irtlat+YEAR+YB+q) = nprod;

       if(LATPARFLAG){
         const double d1 = COVARIATES(q-2); //quadrature point for ability does not rotate
         const double d2 = (COVARIATES(q-1)-THETA(dim_irtlat-2)*d1)/THETA(dim_irtlat-1); // quadrature grid point for speed

         // derivative wrt l1 (L(1,0) where L is the lower Cholesky decomp of lat covariance matrix)
         gr(dim_irtlat-2) = tmp*beta_d(q-1)*d1+nprod*beta_t(q-1)*d1;
         // derivative wrt l1 (L(1,1) where L is the lower Cholesky decomp of lat covariance matrix)
         gr(dim_irtlat-1) = tmp*beta_d(q-1)*d2+nprod*beta_t(q-1)*d2;
       }


     }else if(OUTCOME == 3){
       const double tmp = exp(lpr_t) - pow(exp(lpr_t), 2);
       gr(dim_irtlat+YEAR+YB+q) = tmp;
       gr.segment(dim_irtlat+2*YB+q+1, q) = tmp*COVARIATES;
       gr.segment(dim_irtlat+YB+1, q) = nprod*COVARIATES;
       gr(dim_irtlat+YEAR) = nprod;

       if(LATPARFLAG){
         const double d1 = COVARIATES(q-2); //quadrature point for ability does not rotate
         const double d2 = (COVARIATES(q-1)-THETA(dim_irtlat-2)*d1)/THETA(dim_irtlat-1); // quadrature grid point for speed

         // derivative wrt l1 (L(1,0) where L is the lower Cholesky decomp of lat covariance matrix)
         gr(dim_irtlat-2) = tmp*beta_t(q-1)*d1+nprod*beta_d(q-1)*d1;
         // derivative wrt l1 (L(1,1) where L is the lower Cholesky decomp of lat covariance matrix)
         gr(dim_irtlat-1) = tmp*beta_t(q-1)*d2+nprod*beta_d(q-1)*d2;
       }
     }


   }

   if(OUTCOME == 1){
     const double tmp = exp(theta_cr(0) - log1pexp(theta_cr(0)));
     gr(dim_irtlat) = tmp-pow(tmp, 2); // Probability of graduation when all exams are done
   }


   return gr;
 }

//' Evaluate gradient of survival function given  the range of years of interest
 //'
 //' @param YEAR_LAST Last year to evaluate.
 //' @param THETA Parameter vector
 //' @param COVARIATES The first 2 values refers to ability and speed respectively. Remaining values are external predictors.
 //' @param YB Total number of years in the non-graduatable regime. Needed for determining how many time-related intercepts.
 //' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time
 //' @param LOGFLAG Set TRUE to return log value.
 //'
 //' @returns It returns the probability of survival from `YEAR FIRST` to `YEAR_LAST` included.
 //'
 Eigen::VectorXd gr_survival(
     const unsigned int YEAR_FIRST,
     const unsigned int YEAR_LAST,
     const Eigen::Ref<const Eigen::VectorXd> THETA,
     const Eigen::Ref<const Eigen::VectorXd> COVARIATES,
     const unsigned int YB,
     const unsigned int YEAR_LAST_EXAM,
     const bool ROTATED
 ){

   const int q = COVARIATES.size();
   const int dim_cr = 2*(YB+q)+1;
   const int dim_irtlat = THETA.size()-dim_cr;
   Eigen::VectorXd grl = Eigen::VectorXd::Zero(THETA.size());
   double logout = 0;


   if(YEAR_LAST_EXAM > YEAR_LAST){
     // Regime where graduation is not possible
     if(YEAR_LAST > YB) Rcpp::stop("Regime 1 mismatch: `YEAR_LAST` > `NYB`");

     // Remain enrolled until YEAR_LAST (conditioned on not having all exams)
     for(unsigned int year = YEAR_FIRST; year<=YEAR_LAST; year++){
       double logouty = log1mexp(-R::logspace_add(hazard(2, year, THETA, COVARIATES, YB, true), hazard(3, year, THETA, COVARIATES, YB, true)));
       logout += logouty;
       Eigen::VectorXd gry = - gr_hazard(2, year, THETA, COVARIATES, YB, ROTATED) - gr_hazard(3, year, THETA, COVARIATES, YB, ROTATED);
       grl += gry/exp(logouty);
     }
   }else if(YEAR_LAST_EXAM <= YEAR_LAST){

     // Regime where graduation is possible from year `YEAR_LAST_EXAM`

     // Remain enrolled until YEAR_LAST_EXAM
     for(unsigned int year = YEAR_FIRST; year < YEAR_LAST_EXAM; year++){
       double logouty = log1mexp(-R::logspace_add(hazard(2, year, THETA, COVARIATES, YB, true), hazard(3, year, THETA, COVARIATES, YB, true)));
       logout += logouty;
       Eigen::VectorXd gry = - gr_hazard(2, year, THETA, COVARIATES, YB, ROTATED) - gr_hazard(3, year, THETA, COVARIATES, YB, ROTATED);
       grl += gry/exp(logouty);
     }
     // Remain enrolled from YEAR_LAST_EXAM to YEAR_LAST
     for(unsigned int year = YEAR_LAST_EXAM; year <= YEAR_LAST; year++){
       double logouty = log1mexp(-hazard(1, year - YEAR_LAST_EXAM + 1, THETA, COVARIATES, YB, true));
       logout += logouty;
       grl(dim_irtlat)+=exp(logouty)-1;
     }


   }

   return(grl*exp(logout));
 }



//' Evaluate the gradient of the outcome log-likelihood
 //'
 //' @param OUTCOME  `1` for graduation, `2` for dropout, `3` for transfer. `0` if no outcome is observed.
 //' @param YEAR_LAST Last year to evaluate.
 //' @param THETA Parameter vector.
 //' @param COVARIATES The first 2 values refers to ability and speed respectively. Remaining values are external predictors.
 //' @param YB Total number of years in the non-graduatable regime. Needed for determining how many time-related intercepts.
 //' @param YEAR_LAST_EXAM Year at which the all exams are completed for the first time.
 //'
 Eigen::VectorXd grl_outcomeLik(
     const unsigned int OUTCOME,
     const unsigned int YEAR_FIRST,
     const unsigned int YEAR_LAST,
     const Eigen::Ref<const Eigen::VectorXd> THETA,
     const Eigen::Ref<const Eigen::VectorXd> COVARIATES,
     const unsigned int YB,
     const unsigned int YEAR_LAST_EXAM,
     const bool ROTATED
 ){
   Eigen::VectorXd grl = Eigen::VectorXd::Zero(THETA.size());
   double logout;

   if(OUTCOME==0){

     logout = survival(YEAR_FIRST, YEAR_LAST, THETA, COVARIATES, YB, YEAR_LAST_EXAM, true);
     grl = gr_survival(YEAR_FIRST, YEAR_LAST, THETA, COVARIATES, YB, YEAR_LAST_EXAM, ROTATED)/exp(logout);

   }else if(OUTCOME==2 | OUTCOME==3){
     if(YEAR_LAST_EXAM<=YEAR_LAST){
       logout = -10000; //return R_NegInf;
     } else {
       logout  = survival(YEAR_FIRST, YEAR_LAST-1, THETA, COVARIATES, YB, YEAR_LAST_EXAM, true);
       grl = gr_survival(YEAR_FIRST, YEAR_LAST-1, THETA, COVARIATES, YB, YEAR_LAST_EXAM, ROTATED)/exp(logout);

       double tmp = hazard(OUTCOME, YEAR_LAST, THETA, COVARIATES, YB, true);
       logout += tmp;
       grl += gr_hazard(OUTCOME, YEAR_LAST, THETA, COVARIATES, YB, ROTATED)/exp(tmp);
     }
   }else if(OUTCOME==1){
     if(YEAR_LAST_EXAM>YEAR_LAST){
       logout = -10000; //return R_NegInf;
     } else{
       logout  = survival(YEAR_FIRST, YEAR_LAST-1, THETA, COVARIATES, YB, YEAR_LAST_EXAM, true);
       grl = gr_survival(YEAR_FIRST, YEAR_LAST-1, THETA, COVARIATES, YB, YEAR_LAST_EXAM, ROTATED)/exp(logout);

       double tmp = hazard(OUTCOME, YEAR_LAST, THETA, COVARIATES, YB, true);
       logout += tmp;
       grl += gr_hazard(OUTCOME, YEAR_LAST, THETA, COVARIATES, YB, ROTATED)/exp(tmp);
     }


   }

   return(grl);
 }
}
#endif
