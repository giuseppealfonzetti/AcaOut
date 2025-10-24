#ifndef extractParams_H
#define extractParams_H
#include "thresholds.h"
#include <cmath>

//' Extract indices of parameters related to the IRT model
//'
//' `THETA_IRT` is assumed to be ordered as follows
//' c(exam1_int, exam1_slope, exam1_time_loc, exam1, time_sca, exam2_int, ...)
//'
//' @param THETA_IRT parameter vector related to irt model
//' @param N_GRADES number of grades modelled.
//' @param N_EXAMS number of exams.
//' @param OPTION Select parameters of interest. `1` for exam-grades slopes,
//' `2` for exam-grade intercepts, `3` for time-location parameter,
//' `4` for time-scale parameter.
//' @param EXAM exam of interest. Possible values in `0:(N_EXAMS-1)`.
//'
//' @returns It return a vector with two values representing the starting
//' index and the length of the segment of THETA_IRT corresponding to
//' the parameter of interest.
//'
//' @export
// [[Rcpp::export]]
std::vector<unsigned int> extract_params_idx_irt(
     Eigen::VectorXd THETA_IRT,
     const unsigned int N_GRADES,
     const unsigned int N_EXAMS,
     const unsigned int OPTION,
     const unsigned int EXAM
 ){
  std::vector<unsigned int> out(2);

   switch(OPTION){
   case 1: // slope
     out[0] = EXAM*(N_GRADES+3) + N_GRADES ;
     out[1] = 1;
     break;
   case 2: // intercept
     out[0] = EXAM*(N_GRADES+3);
     out[1] = N_GRADES;
     break;
   case 3: // time location
     out[0] = EXAM*(N_GRADES+3) + N_GRADES + 1;
     out[1] = 1;
     break;
   case 4: // time scale
     out[0] = EXAM*(N_GRADES+3) + N_GRADES + 2;
     out[1] = 1;
     break;
   }


   return(out);
 }


//' Extract parameters related to the IRT model
//'
//' See [extract_params_idx_irt]
//' @param THETA_IRT parameter vector related to irt model
//' @param N_GRADES number of grades modelled.
//' @param N_EXAMS number of exams.
//' @param OPTION Select parameters of interest. `1` for exam-grades slopes,
//' `2` for exam-grade intercepts, `3` for time-location parameter,
//' `4` for time-scale parameter.
//' @param EXAM exam of interest. Possible values in `1:N_EXAMS`.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd extract_params_irt(
    Eigen::VectorXd THETA_IRT,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const unsigned int OPTION,
    const unsigned int EXAM
){
  std::vector<unsigned int> idx = extract_params_idx_irt(THETA_IRT, N_GRADES, N_EXAMS, OPTION, EXAM);
  Eigen::VectorXd out = THETA_IRT.segment(idx[0], idx[1]);

  // Reparameterisations
  if(OPTION == 2){
    out = reparThr(out, false);
  }else if(OPTION == 4){
    out = out.array().exp();
  }
  return(out);
}


namespace params
{
 Eigen::MatrixXd get_L(
     const Eigen::Ref<const Eigen::VectorXd> THETA,
     const int N_EXAMS,
     const int N_GRADES
 ){
   const unsigned int dim_irt = N_EXAMS*(N_GRADES+3);
   Eigen::MatrixXd L(2,2);
   const double l21 = THETA(dim_irt);
   const double l22 = std::exp(THETA(dim_irt+1));
   L << 1, l21, 0, l22;
   return(L);
 }
}
#endif
