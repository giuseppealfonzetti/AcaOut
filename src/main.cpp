#include <Rcpp.h>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>
#include <cmath>

#include "extractParams.h"
#include "latent.h"
#include "cr.h"
#include "testthat_wrappers.h"
// #include "conditional_models.h"
#include "joint_models.h"
#include "GRTCM.h"
#include "gq.h"


// [[Rcpp::export]]
Rcpp::List CCR(
    Eigen::Map<Eigen::VectorXd> THETA,
    Eigen::Map<Eigen::VectorXd> OUTCOME,
    Eigen::Map<Eigen::MatrixXd> COVARIATES,
    Eigen::Map<Eigen::VectorXd> YEAR_FIRST,
    Eigen::Map<Eigen::VectorXd> YEAR_LAST,
    Eigen::Map<Eigen::VectorXd> YEAR_LAST_EXAM,
    Eigen::Map<Eigen::MatrixXd> LATMAT,
    const int YB,
    const bool GRFLAG
){

  Rcpp::List output = cr::ccr_sample(
    THETA,
    OUTCOME,
    COVARIATES,
    YEAR_FIRST,
    YEAR_LAST,
    YEAR_LAST_EXAM,
    LATMAT,
    YB,
    GRFLAG
  );

  return output;
}



// [[Rcpp::export]]
Rcpp::List cpp_EM(
    Eigen::VectorXd THETA_START,
    Eigen::MatrixXd EXAMS_GRADES,
    Eigen::MatrixXd EXAMS_DAYS,
    Eigen::MatrixXd EXAMS_SET,
    Eigen::MatrixXd EXAMS_OBSFLAG,
    Eigen::VectorXd MAX_DAY,
    Eigen::VectorXd OUTCOME,
    Eigen::MatrixXd EXT_COVARIATES,
    Eigen::VectorXd YEAR_FIRST,
    Eigen::VectorXd YEAR_LAST,
    Eigen::VectorXd YEAR_LAST_EXAM,
    Eigen::MatrixXd GRID,
    Eigen::VectorXd WEIGHTS,
    const int YB,
    const int N_GRADES,
    const int N_EXAMS,
    const int M_MAX_ITER,
    const int MAX_ITER,
    const double TOL,
    const std::string MOD,
    const bool VERBOSE
){
  double enjll = 0;
  Eigen::VectorXd theta = THETA_START;

  const int n = EXAMS_GRADES.rows();
  const int nq = GRID.rows();
  const int dim_irt = N_EXAMS*(N_GRADES+3);
  const int dim_cr = 2*(YB+EXT_COVARIATES.cols()+2)+1;
  std::vector<double> path_enjll; path_enjll.push_back(std::numeric_limits<double>::infinity());
  std::vector<Eigen::VectorXd> path_theta; path_theta.push_back(theta);

  int last_iter = MAX_ITER;
  int convergence = 0;
  for(int iter=0; iter < MAX_ITER; iter++){
    Rcpp::checkUserInterrupt();
    // Eigen::MatrixXd L{{1,0},{theta(dim_irt), theta(dim_irt+1)}};
    // Eigen::MatrixXd grid = GRID * L.transpose();

    if(VERBOSE)Rcpp::Rcout << "Iter " << iter << ":\n";
    if(VERBOSE)Rcpp::Rcout << "- E-STEP...";
    Eigen::MatrixXd Ew = EM::Estep(theta,
                                   EXAMS_GRADES,
                                   EXAMS_DAYS,
                                   EXAMS_SET,
                                   EXAMS_OBSFLAG,
                                   MAX_DAY,
                                   OUTCOME,
                                   EXT_COVARIATES,
                                   YEAR_FIRST,
                                   YEAR_LAST,
                                   YEAR_LAST_EXAM,
                                   GRID,
                                   WEIGHTS,
                                   YB,
                                   N_GRADES,
                                   N_EXAMS,
                                   MOD);

    if(VERBOSE)Rcpp::Rcout << " Weights computed.\n";
    EM::EAPLOGJ eclass(
        EXAMS_GRADES,
        EXAMS_DAYS,
        EXAMS_SET,
        EXAMS_OBSFLAG,
        MAX_DAY,
        OUTCOME,
        EXT_COVARIATES,
        YEAR_FIRST,
        YEAR_LAST,
        YEAR_LAST_EXAM,
        YB,
        N_GRADES,
        N_EXAMS,
        MOD
    );

    eclass.update_quadrature(GRID, Ew);

    MCOUNT = 0;
    int status = optim_lbfgs(eclass, theta, enjll, M_MAX_ITER);
    double tol_check = (path_enjll.back() - enjll) / path_enjll.back();

    MCOUNT = 0;
    if(VERBOSE){
      if(status==0) Rcpp::Rcout <<  " | Converged.";
      Rcpp::Rcout <<  "\n- obj=" << enjll << ", obj_pdiff:" << tol_check  <<"|\n";
    }

    path_theta.push_back(theta);
    path_enjll.push_back(enjll);

    if(tol_check<TOL){
      last_iter = iter;
      convergence = 1;
      if(VERBOSE)Rcpp::Rcout << "EM converged correctly\n";
      break;
    }
  }

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("path_theta") = path_theta,
      Rcpp::Named("path_enjll") = path_enjll,
      Rcpp::Named("last_iter") = last_iter,
      Rcpp::Named("convergence") = convergence,
      Rcpp::Named("par") = theta
    );

  return output;
}


// [[Rcpp::export]]
Rcpp::List cpp_GQ(
    Eigen::VectorXd THETA,
    Eigen::MatrixXd EXAMS_GRADES,
    Eigen::MatrixXd EXAMS_DAYS,
    Eigen::MatrixXd EXAMS_SET,
    Eigen::MatrixXd EXAMS_OBSFLAG,
    Eigen::MatrixXd COVARIATES,
    Eigen::VectorXd MAX_DAY,
    Eigen::VectorXd OUTCOME,
    Eigen::VectorXd YEAR_FIRST,
    Eigen::VectorXd YEAR_LAST,
    Eigen::VectorXd YEAR_LAST_EXAM,
    Eigen::MatrixXd GRID,
    Eigen::VectorXd WEIGHTS,
    const int N_GRADES,
    const int N_EXAMS,
    const int YB,
    const std::string MOD,
    const bool GRFLAG = true,
    const bool LATPARFLAG = true,
  const bool HFLAG=false
){

  Rcpp::List output = gq::evaluate_grid(
    THETA,
    EXAMS_GRADES,
    EXAMS_DAYS,
    EXAMS_SET,
    EXAMS_OBSFLAG,
    COVARIATES,
    MAX_DAY,
    OUTCOME,
    YEAR_FIRST,
    YEAR_LAST,
    YEAR_LAST_EXAM,
    GRID,
    WEIGHTS,
    N_GRADES,
    N_EXAMS,
    YB,
    MOD,
    GRFLAG,
    true,
    HFLAG
  );

  return output;
}



 // [[Rcpp::export]]
 Rcpp::List cpp_EAP(
     Eigen::VectorXd THETA,
     Eigen::MatrixXd EXAMS_GRADES,
     Eigen::MatrixXd EXAMS_DAYS,
     Eigen::MatrixXd EXAMS_SET,
     Eigen::MatrixXd EXAMS_OBSFLAG,
     Eigen::VectorXd MAX_DAY,
     Eigen::VectorXd OUTCOME,
     Eigen::MatrixXd EXT_COVARIATES,
     Eigen::VectorXd YEAR_FIRST,
     Eigen::VectorXd YEAR_LAST,
     Eigen::VectorXd YEAR_LAST_EXAM,
     Eigen::MatrixXd GRID,
     Eigen::VectorXd WEIGHTS,
     const int YB,
     const int N_GRADES,
     const int N_EXAMS,
     const std::string MOD,
     const bool VERBOSE
 ){
   double enjll = 0;

   const int n = EXAMS_GRADES.rows();
   const int nq = GRID.rows();
   const int dim_irt = N_EXAMS*(N_GRADES+3);
   const int dim_cr = 2*(YB+EXT_COVARIATES.cols()+2)+1;
   const int n_cov = EXT_COVARIATES.cols();


  Eigen::MatrixXd Ew = EM::Estep(THETA,
                                 EXAMS_GRADES,
                                 EXAMS_DAYS,
                                 EXAMS_SET,
                                 EXAMS_OBSFLAG,
                                 MAX_DAY,
                                 OUTCOME,
                                 EXT_COVARIATES,
                                 YEAR_FIRST,
                                 YEAR_LAST,
                                 YEAR_LAST_EXAM,
                                 GRID,
                                 WEIGHTS,
                                 YB,
                                 N_GRADES,
                                 N_EXAMS,
                                 MOD);


  Eigen::MatrixXd grid = GRID;
  const double l21 = THETA(dim_irt);
  const double l22 = std::exp(THETA(dim_irt+1));
  Eigen::MatrixXd L{{1,0},{l21, l22}};
  grid = grid * L.transpose();

  Eigen::MatrixXd eap=Eigen::MatrixXd::Zero(n,2);
  for(unsigned int i = 0; i < n; i++){
    Rcpp::checkUserInterrupt();

    double mu_i_ability = 0;
    double mu_i_speed = 0;
    mu_i_ability = EXT_COVARIATES.row(i)*THETA.segment(dim_irt+2, n_cov);
    mu_i_speed = EXT_COVARIATES.row(i)*THETA.segment(dim_irt+2+n_cov, n_cov);

    for(unsigned int point = 0; point < nq; point++){

      Eigen::VectorXd lat = grid.row(point);
      lat(0) +=mu_i_ability; lat(1)+=mu_i_speed;
      eap.row(i) += lat*Ew(i, point);

    }
  }

   Rcpp::List output =
     Rcpp::List::create(
       Rcpp::Named("Ew") = Ew,
       Rcpp::Named("EAP")=eap
     );

   return output;
 }
