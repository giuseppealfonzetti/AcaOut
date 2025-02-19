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
    const double ABILITY,
    const double SPEED,
    const int YB,
    const bool LOGFLAG = false
  ){

    Eigen::VectorXd latent(2); latent << ABILITY, SPEED;
    Eigen::VectorXd theta_cr  = THETA.tail(2*(YB+2)+1);
    double out;


    if((OUTCOME == 2 | OUTCOME == 3) && YEAR > YB) Rcpp::stop("`YEAR` larger than `YB`");

    if(OUTCOME == 2 | OUTCOME == 3){
      const double int_d = theta_cr(YEAR);
      const double int_t = theta_cr(YEAR+YB+2);
      const Eigen::VectorXd beta_d = theta_cr.segment(YB+1, 2);
      const Eigen::VectorXd beta_t = theta_cr.segment(2*YB+3, 2);
      const double eta_d = int_d + beta_d.dot(latent);
      const double eta_t = int_t + beta_t.dot(latent);

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
    const double ABILITY,
    const double SPEED,
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
      logout += log1mexp(-R::logspace_add(hazard(2, year, THETA, ABILITY, SPEED, YB, true),
                                          hazard(3, year, THETA, ABILITY, SPEED, YB, true)));
    }

  }else if(YEAR_LAST_EXAM <= YEAR_LAST){

    // Regime where graduation is possible from year `YEAR_LAST_EXAM`

    // Remain enrolled until YEAR_LAST_EXAM
    for(unsigned int year = YEAR_FIRST; year < YEAR_LAST_EXAM; year++){
      logout += log1mexp(-R::logspace_add(hazard(2, year, THETA, ABILITY, SPEED, YB, true),
                                          hazard(3, year, THETA, ABILITY, SPEED, YB, true)));
    }

    // Remain enrolled from YEAR_LAST_EXAM to YEAR_LAST
    for(unsigned int year = YEAR_LAST_EXAM; year <= YEAR_LAST; year++){
      logout += log1mexp(-hazard(1, year - YEAR_LAST_EXAM + 1, THETA, ABILITY, SPEED, YB, true));
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
      const int OUTCOME,
      const int YEAR_FIRST,
      const int YEAR_LAST,
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const double ABILITY,
      const double SPEED,
      const int YB,
      const int YEAR_LAST_EXAM = 100,
      const bool LOGFLAG = false
  ){
    double logout;

    if(OUTCOME==0){
      // The student is still enrolled
      logout = survival(YEAR_FIRST, YEAR_LAST, THETA, ABILITY, SPEED, YB, YEAR_LAST_EXAM, true);
    }else if(OUTCOME==2|OUTCOME==3){
      if((YEAR_LAST_EXAM<=YEAR_LAST) & !LOGFLAG) return 0;
      if((YEAR_LAST_EXAM<=YEAR_LAST) & LOGFLAG) return -10000;//return R_NegInf;

      //The student remain enrolled until YEAR_LAST-1 and the experience OUTCOME
      logout  = survival(YEAR_FIRST, YEAR_LAST-1, THETA, ABILITY, SPEED, YB, YEAR_LAST_EXAM, true);
      logout += hazard(OUTCOME, YEAR_LAST, THETA, ABILITY, SPEED, YB, true);
    }else if(OUTCOME==1){
      if((YEAR_LAST_EXAM>YEAR_LAST) & !LOGFLAG) return 0;
      if((YEAR_LAST_EXAM>YEAR_LAST) & LOGFLAG) return -10000;//return R_NegInf;

      //The student remain enrolled until YEAR_LAST-1 and the experience OUTCOME
      logout  = survival(YEAR_FIRST, YEAR_LAST-1, THETA, ABILITY, SPEED, YB, YEAR_LAST_EXAM, true);
      logout += hazard(1, YEAR_LAST-YEAR_LAST_EXAM+1, THETA, ABILITY, SPEED, YB, true);
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
     const double ABILITY,
     const double SPEED,
     const int YB,
     const bool LATPARFLAG
 ){

    const int n_cov   = COVARIATES.size();
    const int dim_cr  = 2*(YB+2)+1;
    const int dim_lat = 2+2*n_cov;
    const int dim_irt = THETA.size()-dim_cr-dim_lat;

    Eigen::VectorXd latent(2); latent << ABILITY, SPEED;
    Eigen::VectorXd theta_cr  = THETA.tail(dim_cr);
    Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA.size());
    double lpr;
    if((OUTCOME == 2 | OUTCOME == 3) && YEAR > YB) Rcpp::stop("`YEAR` larger than `YB`");

    if(OUTCOME == 2 | OUTCOME == 3){
       const double int_d = theta_cr(YEAR);
       const double int_t = theta_cr(YEAR+YB+2);
       const Eigen::VectorXd beta_d = theta_cr.segment(YB+1, 2);
       const Eigen::VectorXd beta_t = theta_cr.segment(2*YB+3, 2);
       const double eta_d = int_d + beta_d.dot(latent);
       const double eta_t = int_t + beta_t.dot(latent);

       double logConstNorm = log1pexp(R::logspace_add(eta_d, eta_t));
       double lpr_d = eta_d-logConstNorm;
       double lpr_t = eta_t-logConstNorm;
       double nprod = -exp(lpr_d+lpr_t);

     if(OUTCOME == 2){
       const double tmp = exp(lpr_d) - pow(exp(lpr_d), 2);
       gr(dim_irt+dim_lat+YEAR) = tmp;
       gr.segment(dim_irt+dim_lat+YB+1, 2) = tmp*latent;
       gr.segment(dim_irt+dim_lat+2*YB+3, 2) = nprod*latent;
       gr(dim_irt+dim_lat+YEAR+YB+2) = nprod;

       if(LATPARFLAG){
          const double d1 = ABILITY-(THETA.segment(dim_irt+2, n_cov).transpose()*COVARIATES);
          const double d2 = (SPEED -THETA(dim_irt)*d1-(THETA.segment(dim_irt+2+n_cov, n_cov).transpose()*COVARIATES))/THETA(dim_irt+1);

          // derivative wrt l1 (L(1,0) where L is the lower Cholesky decomp of lat covariance matrix)
          gr(dim_irt) = tmp*beta_d(1)*d1+nprod*beta_t(1)*d1;
          // derivative wrt l1 (L(1,1) where L is the lower Cholesky decomp of lat covariance matrix)
          gr(dim_irt+1) = tmp*beta_d(1)*d2+nprod*beta_t(1)*d2;

          // ability regression coefficients
          gr.segment(dim_irt+2, n_cov) = tmp*beta_d(0)*COVARIATES+nprod*beta_t(0)*COVARIATES;
          // speed regression coefficients
          gr.segment(dim_irt+2+n_cov, n_cov) = tmp*beta_d(1)*COVARIATES+nprod*beta_t(1)*COVARIATES;

       }


     }else if(OUTCOME == 3){
       const double tmp = exp(lpr_t) - pow(exp(lpr_t), 2);
       gr(dim_irt+dim_lat+YEAR+YB+2) = tmp;
       gr.segment(dim_irt+dim_lat+2*YB+3, 2) = tmp*latent;
       gr.segment(dim_irt+dim_lat+YB+1, 2) = nprod*latent;
       gr(dim_irt+dim_lat+YEAR) = nprod;

       if(LATPARFLAG){
          const double d1 = ABILITY-(THETA.segment(dim_irt+2, n_cov).transpose()*COVARIATES);
          const double d2 = (SPEED -THETA(dim_irt)*d1-(THETA.segment(dim_irt+2+n_cov, n_cov).transpose()*COVARIATES))/THETA(dim_irt+1);

          // derivative wrt l1 (L(1,0) where L is the lower Cholesky decomp of lat covariance matrix)
          gr(dim_irt) = tmp*beta_t(1)*d1+nprod*beta_d(1)*d1;
          // derivative wrt l1 (L(1,1) where L is the lower Cholesky decomp of lat covariance matrix)
          gr(dim_irt+1) = tmp*beta_t(1)*d2+nprod*beta_d(1)*d2;

          // ability regression coefficients
          gr.segment(dim_irt+2, n_cov) = tmp*beta_t(0)*COVARIATES+nprod*beta_d(0)*COVARIATES;
          // speed regression coefficients
          gr.segment(dim_irt+2+n_cov, n_cov) = tmp*beta_t(1)*COVARIATES+nprod*beta_d(1)*COVARIATES;

       }
     }


   }

    if(OUTCOME == 1){
     const double tmp = exp(theta_cr(0) - log1pexp(theta_cr(0)));
     gr(dim_irt+dim_lat) = tmp-pow(tmp, 2); // Probability of graduation when all exams are done
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
     const int YEAR_FIRST,
     const int YEAR_LAST,
     const Eigen::Ref<const Eigen::VectorXd> THETA,
     const Eigen::Ref<const Eigen::VectorXd> COVARIATES,
     const double ABILITY,
     const double SPEED,
     const int YB,
     const int YEAR_LAST_EXAM,
     const bool LATPARFLAG
 ){

    const int n_cov   = COVARIATES.size();
    const int dim_cr  = 2*(YB+2)+1;
    const int dim_lat = 2+2*n_cov;
    const int dim_irt = THETA.size()-dim_cr-dim_lat;

    Eigen::VectorXd latent(2); latent << ABILITY, SPEED;
    Eigen::VectorXd theta_cr  = THETA.tail(dim_cr);
    Eigen::VectorXd grl = Eigen::VectorXd::Zero(THETA.size());
    double logout = 0;


    if(YEAR_LAST_EXAM > YEAR_LAST){
     // Regime where graduation is not possible
     if(YEAR_LAST > YB) Rcpp::stop("Regime 1 mismatch: `YEAR_LAST` > `NYB`");

     // Remain enrolled until YEAR_LAST (conditioned on not having all exams)
     for(unsigned int year = YEAR_FIRST; year<=YEAR_LAST; year++){
       double logouty = log1mexp(-R::logspace_add(hazard(2, year, THETA, ABILITY, SPEED, YB, true), hazard(3, year, THETA, ABILITY, SPEED, YB, true)));
       logout += logouty;
       Eigen::VectorXd gry = - gr_hazard(2, year, THETA, COVARIATES, ABILITY, SPEED, YB, LATPARFLAG) - gr_hazard(3, year, THETA, COVARIATES, ABILITY, SPEED, YB, LATPARFLAG);
       grl += gry/exp(logouty);
     }
    }else if(YEAR_LAST_EXAM <= YEAR_LAST){

     // Regime where graduation is possible from year `YEAR_LAST_EXAM`

     // Remain enrolled until YEAR_LAST_EXAM
     for(unsigned int year = YEAR_FIRST; year < YEAR_LAST_EXAM; year++){
       double logouty = log1mexp(-R::logspace_add(hazard(2, year, THETA, ABILITY, SPEED, YB, true), hazard(3, year, THETA, ABILITY, SPEED, YB, true)));
       logout += logouty;
       Eigen::VectorXd gry = - gr_hazard(2, year, THETA, COVARIATES, ABILITY, SPEED, YB, LATPARFLAG) - gr_hazard(3, year, THETA, COVARIATES, ABILITY, SPEED, YB, LATPARFLAG);
       grl += gry/exp(logouty);
     }
     // Remain enrolled from YEAR_LAST_EXAM to YEAR_LAST
     for(unsigned int year = YEAR_LAST_EXAM; year <= YEAR_LAST; year++){
       double logouty = log1mexp(-hazard(1, year - YEAR_LAST_EXAM + 1, THETA, ABILITY, SPEED, YB, true));
       logout += logouty;
       grl(dim_irt+dim_lat)+=exp(logouty)-1;
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
     const double ABILITY,
     const double SPEED,
     const unsigned int YB,
     const unsigned int YEAR_LAST_EXAM,
     const bool LATPARFLAG
 ){
   Eigen::VectorXd grl = Eigen::VectorXd::Zero(THETA.size());
   double logout;

   if(OUTCOME==0){

     logout = survival(YEAR_FIRST, YEAR_LAST, THETA, ABILITY, SPEED, YB, YEAR_LAST_EXAM, true);
     grl = gr_survival(YEAR_FIRST, YEAR_LAST, THETA, COVARIATES, ABILITY, SPEED, YB, YEAR_LAST_EXAM, LATPARFLAG)/exp(logout);

   }else if(OUTCOME==2 | OUTCOME==3){
     if(YEAR_LAST_EXAM<=YEAR_LAST){
       logout = -10000; //return R_NegInf;
     } else {
       logout  = survival(YEAR_FIRST, YEAR_LAST-1, THETA, ABILITY, SPEED, YB, YEAR_LAST_EXAM, true);
       grl = gr_survival(YEAR_FIRST, YEAR_LAST-1, THETA, COVARIATES, ABILITY, SPEED, YB, YEAR_LAST_EXAM, LATPARFLAG)/exp(logout);

       double tmp = hazard(OUTCOME, YEAR_LAST, THETA, ABILITY, SPEED, YB, true);
       logout += tmp;
       grl += gr_hazard(OUTCOME, YEAR_LAST, THETA, COVARIATES, ABILITY, SPEED, YB, LATPARFLAG)/exp(tmp);
     }
   }else if(OUTCOME==1){
     if(YEAR_LAST_EXAM>YEAR_LAST){
       logout = -10000; //return R_NegInf;
     } else{
       logout  = survival(YEAR_FIRST, YEAR_LAST-1, THETA, ABILITY, SPEED, YB, YEAR_LAST_EXAM, true);
       grl = gr_survival(YEAR_FIRST, YEAR_LAST-1, THETA, COVARIATES, ABILITY, SPEED, YB, YEAR_LAST_EXAM, LATPARFLAG)/exp(logout);

       double tmp = hazard(OUTCOME, YEAR_LAST, THETA, ABILITY, SPEED, YB, true);
       logout += tmp;
       grl += gr_hazard(OUTCOME, YEAR_LAST, THETA, COVARIATES, ABILITY, SPEED, YB, LATPARFLAG)/exp(tmp);
     }


   }

   return(grl);
 }
}

namespace cr{
   // COMPETING RISK CONDITIONAL MODEL
   class CCR
   {
   private:
     Eigen::VectorXd _theta;
     int _outcome;
     Eigen::VectorXd _covariates;
     int _yb;
     int _year_first;
     int _year_last;
     int _year_last_exam;
     bool _latparflag;

   public:
      CCR(Eigen::VectorXd THETA,
            const int OUTCOME,
            Eigen::VectorXd COVARIATES,
            const int YB,
            const int YEAR_FIRST,
            const int YEAR_LAST,
            const int YEAR_LAST_EXAM,
            const bool LATPARFLAG):
     _theta(THETA),
     _outcome(OUTCOME),
     _covariates(COVARIATES),
     _yb(YB),
     _year_first(YEAR_FIRST),
     _year_last(YEAR_LAST),
     _year_last_exam(YEAR_LAST_EXAM),
     _latparflag(LATPARFLAG){}

     // conditional log-likelihood
     double ll(const double ABILITY, const double SPEED);

     //gradient of conditional log-likelihood
     Eigen::VectorXd grll(const double ABILITY, const double SPEED);

   };
   double CCR::ll(const double ABILITY, const double SPEED){

     double out = cr::outcomeLik(
       _outcome,
       _year_first,
       _year_last,
       _theta,
       ABILITY,
       SPEED,
       _yb,
       _year_last_exam,
       true);

     return out;
   }
   Eigen::VectorXd CCR::grll(const double ABILITY, const double SPEED){

     Eigen::VectorXd out=cr::grad::grl_outcomeLik(
       _outcome,
       _year_first,
       _year_last,
       _theta,
       _covariates,
       ABILITY,
       SPEED,
       _yb,
       _year_last_exam,
       _latparflag
     );

     return(out);
   }

   Rcpp::List ccr_sample(
         const Eigen::Ref<const Eigen::VectorXd> THETA,
         const Eigen::Ref<const Eigen::VectorXd> OUTCOME,
         const Eigen::Ref<const Eigen::MatrixXd> COVARIATES,
         const Eigen::Ref<const Eigen::VectorXd> YEAR_FIRST,
         const Eigen::Ref<const Eigen::VectorXd> YEAR_LAST,
         const Eigen::Ref<const Eigen::VectorXd> YEAR_LAST_EXAM,
         const Eigen::Ref<const Eigen::MatrixXd> LATMAT,
         const int YB,
         const bool GRFLAG = true
   ){
      const int n = OUTCOME.size();
      double ll = 0;
      Eigen::VectorXd grll = Eigen::VectorXd::Zero(THETA.size());

      for(int i = 0; i < n; i++){

        // Initialize conditional CR model
        CCR ccr_mod(THETA,
                    OUTCOME(i),
                    COVARIATES.row(i),
                    YB,
                    YEAR_FIRST(i),
                    YEAR_LAST(i),
                    YEAR_LAST_EXAM(i),
                    false);


        ll += ccr_mod.ll(LATMAT(i,0), LATMAT(i,1));
        if(GRFLAG){
          grll += ccr_mod.grll(LATMAT(i,0), LATMAT(i,1));
        }

      }

      Rcpp::List output =
         Rcpp::List::create(
            Rcpp::Named("gr") = grll,
            Rcpp::Named("ll") = ll
         );

      return output;
   }

}
#endif
