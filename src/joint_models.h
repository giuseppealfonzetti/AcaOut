#ifndef jointMods_H
#define jointMods_H
#include "grades.h"
#include "times.h"
#include "exams.h"
#include "extractParams.h"
#include "latent.h"
#include "cr.h"
#include "conditional_models.h"

//' IRT complete observation
//'
//' Evaluate the IRT complete model log likelihood of a single observation
//'
//' @param THETA Suitable parameter vector as provided by [parList2Vec]
//' @param EXAMS_GRADES Vector of grades.
//' @param EXAMS_DAYS Vector of times.
//' @param EXAMS_OBSFLAG Vector of booleans.`TRUE` elements represent observed exams. `FALSE` elements the unobserved ones.
//' @param EXAMS_SET Vector filled with booleans.`TRUE` elements represent exams in the study plan. `FALSE` elements non-relevant ones.
//' @param MAX_DAY Last day of observation
//' @param N_GRADES Number of grades modelled.
//' @param N_EXAMS Number of exams modelled
//' @param ABILITY Ability value.
//' @param SPEED Speed value.
//'
//' @return It returns the value of the integrand function,
//' given the parameters and the data of a single observation.
//'
//' @export
// [[Rcpp::export]]
double GRTCM_complete_obs(
    Eigen::VectorXd THETA,
    Eigen::VectorXd EXAMS_GRADES,
    Eigen::VectorXd EXAMS_DAYS,
    Eigen::VectorXd EXAMS_OBSFLAG,
    Eigen::VectorXd EXAMS_SET,
    Eigen::VectorXd COVARIATES,
    const unsigned int MAX_DAY,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY,
    const double SPEED
){

  const unsigned int dim_irt = N_EXAMS*(N_GRADES+3);

  // Initialize conditional IRT model
  GRTC_MOD irt_mod(THETA, EXAMS_GRADES, EXAMS_DAYS,
                  EXAMS_SET, EXAMS_OBSFLAG, COVARIATES, MAX_DAY, N_GRADES, N_EXAMS, false);

  double cll = irt_mod.cll(ABILITY, SPEED);

  return cll;
}


//' @export
// [[Rcpp::export]]
double CRGRTCM_complete_obs(
    Eigen::VectorXd THETA,
    Eigen::VectorXd EXAMS_GRADES,
    Eigen::VectorXd EXAMS_DAYS,
    Eigen::VectorXd EXAMS_OBSFLAG,
    Eigen::VectorXd EXAMS_SET,
    const unsigned int OUTCOME,
    const unsigned int YEAR_FIRST,
    const unsigned int YEAR_LAST,
    const unsigned int YEAR_LAST_EXAM,
    Eigen::VectorXd COVARIATES,
    const unsigned int MAX_DAY,
    const unsigned int YB,
    const unsigned int N_GRADES,
    const unsigned int N_EXAMS,
    const double ABILITY,
    const double SPEED
){

  const unsigned int dim_irt = N_EXAMS*(N_GRADES+3);

  // Initialize conditional IRT model
  GRTC_MOD irt_mod(THETA, EXAMS_GRADES, EXAMS_DAYS,
                  EXAMS_SET, EXAMS_OBSFLAG, COVARIATES, MAX_DAY, N_GRADES, N_EXAMS, false);

  // Initialize conditional CR model
  CR_MOD cr_mod(THETA, OUTCOME, COVARIATES, YB, YEAR_FIRST, YEAR_LAST, YEAR_LAST_EXAM, false);

  double cll = irt_mod.cll(ABILITY, SPEED) + cr_mod.ll(ABILITY, SPEED);

  return cll;
}
#endif
