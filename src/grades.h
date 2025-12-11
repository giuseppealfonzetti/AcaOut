#ifndef grades_H
#define grades_H
#include "extractParams.h"

namespace exams{
   //' Evaluate the probability of grades greater or equal than the reference one
   //'
   //' @param GRADE Grade used as reference. Integers from 1 to N_GRADES.
   //' @param EXAM Exam of interest. Integers from 0 to N_EXAMS -1.
   //' @param THETA_IRT Portion of the parameter vector related to the IRT model
   //' @param N_GRADES Number of grades modelled.
   //' @param N_EXAMS Number of exams.
   //' @param ABILITY Ability value.
   //' @param LOGFLAG Set TRUE to return log value.
   //'
   //' @returns It returns the probability of obtaining grades higher than `GRADE` on exam `EXAM`.
   //'
   double pGreaterGrades(
       const unsigned int GRADE,
       const unsigned int EXAM,
       const Eigen::Ref<const Eigen::VectorXd> THETA_IRT,
       const unsigned int N_GRADES,
       const unsigned int N_EXAMS,
       const double ABILITY,
       const bool LOGFLAG = false
   ){
     if(EXAM > N_EXAMS) Rcpp::stop("`EXAM` larger than `N_EXAMS`");
     if(GRADE > N_GRADES) Rcpp::stop("`GRADE` larger than `N_GRADES`");

     double out;
     switch (GRADE)
     {
     case 0:
       out = 0;
       break;
     default:
       const double intercept = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 2, EXAM)(GRADE-1);
     const double coeff = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 1, EXAM)(0);
     out = -log1pexp(intercept - coeff*ABILITY);
     break;
     }

     if(!LOGFLAG) out = exp(out);

     return(out);
   }

    //' Evaluate the probability of getting a specific grade
    //'
    //' @param GRADE Grade used as reference
    //' @param EXAM Exam of interest. Integers from 0 to N_EXAMS -1.
    //' @param THETA_IRT Portion of the parameter vector related to the IRT model
    //' @param N_GRADES Number of grades modelled.
    //' @param N_EXAMS Number of exams.
    //' @param ABILITY Ability value.
    //' @param LOGFLAG Set TRUE to return log value.
    //'
    //' @returns It returns the probability of obtaining the grade `GRADE` on exam `EXAM`.
    double pGrade(
        const unsigned int GRADE,
        const unsigned int EXAM,
        const Eigen::Ref<const Eigen::VectorXd> THETA_IRT,
        const unsigned int N_GRADES,
        const unsigned int N_EXAMS,
        const double ABILITY,
        const bool LOGFLAG = false
    ){
      double out;

      if(LOGFLAG){
        if(GRADE==N_GRADES){
          out = pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true);
        }else if(GRADE==0){
          out = log1mexp(-pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true));
        }else if((GRADE<N_GRADES) & (GRADE >0)){
          out = R::logspace_sub(pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true), pGreaterGrades(GRADE+1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY, true));
        }

        // avoid log(0)
        out = std::max(-10000.0, out);
      }else{
        if(GRADE==N_GRADES){
          out = pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
        }else if(GRADE==0){
          out = 1 - pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
        }else if((GRADE<N_GRADES) & (GRADE >0)){
          out = pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY)-pGreaterGrades(GRADE+1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
        }
      }


      return(out);
    }
}

namespace exams::grad{

    //' Evaluate the gradient of the probability of grades greater or equal than the reference one
    //'
    //' @param GRADE Grade used as reference. Integers from 1 to N_GRADES.
    //' @param EXAM Exam of interest. Integers from 0 to N_EXAMS -1.
    //' @param THETA_IRT Portion of the parameter vector related to the IRT model
    //' @param THETA_LAT  Portion of the parameter vector related to the latent model
    //' @param COVARIATES Vector of external covariates
    //' @param N_GRADES Number of grades modelled.
    //' @param N_EXAMS Number of exams.
    //' @param ABILITY Ability value.
    //'
    //' @returns It returns the gradient of the probability of obtaining grades higher than `GRADE` on exam `EXAM`.
    //'
    Eigen::VectorXd gr_pGreaterGrades2(
         const unsigned int GRADE,
         const unsigned int EXAM,
         const Eigen::Ref<const Eigen::VectorXd> THETA,
         const Eigen::Ref<const Eigen::VectorXd> COVARIATES,
         const unsigned int N_GRADES,
         const unsigned int N_EXAMS,
         const double ABILITY,
         const bool LATPARFLAG
    ){
      if(EXAM > N_EXAMS) Rcpp::stop("`EXAM` larger than `N_EXAMS`");
      if(GRADE > N_GRADES) Rcpp::stop("`GRADE` larger than `N_GRADES`");
      const int dim_irt = N_EXAMS*(N_GRADES+3);
      const int n_cov   = COVARIATES.size();
      const int dim_lat = 2+2*n_cov;


      Eigen::VectorXd gr = Eigen::VectorXd::Zero(dim_irt+dim_lat);

      const double intercept = extract_params_irt(THETA, N_GRADES, N_EXAMS, 2, EXAM)(GRADE-1);
      const double coeff = extract_params_irt(THETA, N_GRADES, N_EXAMS, 1, EXAM)(0);

      double pr = exp(-log1pexp(intercept - coeff*ABILITY));
      std::vector<unsigned int> idx_i = extract_params_idx_irt(THETA, N_GRADES, N_EXAMS, 2, EXAM);
      std::vector<unsigned int> idx_c = extract_params_idx_irt(THETA, N_GRADES, N_EXAMS, 1, EXAM);

      double val_i = pr - pow(pr,2);
      gr(idx_c[0]) = val_i*ABILITY;

      Eigen::VectorXd tmpi = THETA.segment(idx_i[0], idx_i[1]).array().exp();
      tmpi(0) = 1;
      tmpi.tail(N_GRADES-GRADE) = Eigen::VectorXd::Zero(N_GRADES-GRADE);

      gr.segment(idx_i[0], idx_i[1]) = -tmpi*val_i;

      if(LATPARFLAG && n_cov>0) gr.segment(dim_irt+2, n_cov) = val_i*coeff*COVARIATES;

      return(gr);
     }



    //' Evaluate the gradient of the probability of getting a specific grade
    //'
    //' @param GRADE Grade used as reference
    //' @param EXAM Exam of interest. Integers from 0 to N_EXAMS -1.
    //' @param THETA_IRT Portion of the parameter vector related to the IRT model
    //' @param N_GRADES Number of grades modelled.
    //' @param N_EXAMS Number of exams.
    //' @param ABILITY Ability value.
    //'
    //' @returns It returns the gradient of the probability of obtaining the grade `GRADE` on exam `EXAM`.
    Eigen::VectorXd gr_pGrade2(
        const unsigned int GRADE,
        const unsigned int EXAM,
        const Eigen::Ref<const Eigen::VectorXd> THETA,
        const Eigen::Ref<const Eigen::VectorXd> COVARIATES,
        const unsigned int N_GRADES,
        const unsigned int N_EXAMS,
        const double ABILITY,
        const bool LATPARFLAG
    ){
      double out;
      Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA.size());

      if(GRADE==N_GRADES){
           gr = gr_pGreaterGrades2(GRADE, EXAM, THETA, COVARIATES, N_GRADES, N_EXAMS, ABILITY, LATPARFLAG);
      }else if(GRADE==0){
           gr = - gr_pGreaterGrades2(1, EXAM, THETA, COVARIATES, N_GRADES, N_EXAMS, ABILITY, LATPARFLAG);
      }else if((GRADE<N_GRADES) & (GRADE >0)){
           gr = gr_pGreaterGrades2(GRADE, EXAM, THETA, COVARIATES, N_GRADES, N_EXAMS, ABILITY, LATPARFLAG) - gr_pGreaterGrades2(GRADE+1, EXAM, THETA, COVARIATES, N_GRADES, N_EXAMS, ABILITY, LATPARFLAG);
      }

      return(gr);
    }

    //' Evaluate the gradient of the probability of grades greater or equal than the reference one
    //'
    //' @param GRADE Grade used as reference. Integers from 1 to N_GRADES.
    //' @param EXAM Exam of interest. Integers from 0 to N_EXAMS -1.
    //' @param THETA_IRT Portion of the parameter vector related to the IRT model
    //' @param N_GRADES Number of grades modelled.
    //' @param N_EXAMS Number of exams.
    //' @param ABILITY Ability value.
    //'
    //' @returns It returns the gradient of the probability of obtaining grades higher than `GRADE` on exam `EXAM`.
    //'
    Eigen::VectorXd gr_pGreaterGrades(
         const unsigned int GRADE,
         const unsigned int EXAM,
         const Eigen::Ref<const Eigen::VectorXd> THETA_IRT,
         const unsigned int N_GRADES,
         const unsigned int N_EXAMS,
         const double ABILITY
    ){
      if(EXAM > N_EXAMS) Rcpp::stop("`EXAM` larger than `N_EXAMS`");
      if(GRADE > N_GRADES) Rcpp::stop("`GRADE` larger than `N_GRADES`");
      Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA_IRT.size());

      const double intercept = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 2, EXAM)(GRADE-1);
      const double coeff = extract_params_irt(THETA_IRT, N_GRADES, N_EXAMS, 1, EXAM)(0);

      double pr = exp(-log1pexp(intercept - coeff*ABILITY));
      std::vector<unsigned int> idx_i = extract_params_idx_irt(THETA_IRT, N_GRADES, N_EXAMS, 2, EXAM);
      std::vector<unsigned int> idx_c = extract_params_idx_irt(THETA_IRT, N_GRADES, N_EXAMS, 1, EXAM);

      double val_i = pr - pow(pr,2);
      gr(idx_c[0]) = val_i*ABILITY;

      Eigen::VectorXd tmpi = THETA_IRT.segment(idx_i[0], idx_i[1]).array().exp();
      tmpi(0) = 1;
      tmpi.tail(N_GRADES-GRADE) = Eigen::VectorXd::Zero(N_GRADES-GRADE);

      gr.segment(idx_i[0], idx_i[1]) = -tmpi*val_i;

      return(gr);
     }



    //' Evaluate the gradient of the probability of getting a specific grade
    //'
    //' @param GRADE Grade used as reference
    //' @param EXAM Exam of interest. Integers from 0 to N_EXAMS -1.
    //' @param THETA_IRT Portion of the parameter vector related to the IRT model
    //' @param N_GRADES Number of grades modelled.
    //' @param N_EXAMS Number of exams.
    //' @param ABILITY Ability value.
    //'
    //' @returns It returns the gradient of the probability of obtaining the grade `GRADE` on exam `EXAM`.
    Eigen::VectorXd gr_pGrade(
        const unsigned int GRADE,
        const unsigned int EXAM,
        const Eigen::Ref<const Eigen::VectorXd> THETA_IRT,
        const unsigned int N_GRADES,
        const unsigned int N_EXAMS,
        const double ABILITY
    ){
      double out;
      Eigen::VectorXd gr = Eigen::VectorXd::Zero(THETA_IRT.size());

      if(GRADE==N_GRADES){
           gr = gr_pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
      }else if(GRADE==0){
           gr = - gr_pGreaterGrades(1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
      }else if((GRADE<N_GRADES) & (GRADE >0)){
           gr = gr_pGreaterGrades(GRADE, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY) - gr_pGreaterGrades(GRADE+1, EXAM, THETA_IRT, N_GRADES, N_EXAMS, ABILITY);
      }

      return(gr);
    }
}







#endif
