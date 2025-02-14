#include <Rcpp.h>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_DONT_PARALLELIZE
#include <RcppEigen.h>

#include "extractParams.h"
#include "latent.h"
#include "cr.h"
#include "testthat_wrappers.h"
#include "conditional_models.h"
#include "joint_models.h"

//' GRTCM computations with Gaussian quadrature
//'
//' Compute log-likelihood an gradient of the Graded response time-censored model
//'
//' @param THETA Parameter vector
//' @param EXAMS_GRADES Matrix of exam grades
//' @param EXAMS_DAYS Matrix of exam days
//' @param EXAMS_SET Matrix of binary values. 1 to include an exam in the study-plan, 0 otherwise
//' @param EXAMS_OBSFLAG Matrix of binary values. 1 if the exam is observed, 0 otherwise
//' @param MAX_DAY Vector with maximum day per observation
//' @param GRID Grid of quadrature points
//' @param WEIGHTS Weights of quadrature points
//' @param N_GRADES Number of possible grades modelled
//' @param N_EXAMS Number of possible exams modelled
//' @param ROTGRID TRUE to rotate the quadrature grid using the latent covariance matrix
//' @param GRFLAG Set to true to return the gradient along the log-likelihood
//'
//' @export
// [[Rcpp::export]]
Rcpp::List GRTCM_GH(
    Eigen::VectorXd THETA,
    Eigen::MatrixXd EXAMS_GRADES,
    Eigen::MatrixXd EXAMS_DAYS,
    Eigen::MatrixXd EXAMS_SET,
    Eigen::MatrixXd EXAMS_OBSFLAG,
    Eigen::MatrixXd COVARIATES,
    Eigen::VectorXd MAX_DAY,
    Eigen::MatrixXd GRID,
    Eigen::VectorXd WEIGHTS,
    const int N_GRADES,
    const int N_EXAMS,
    const bool GRFLAG = true,
    const bool ROTGRID = true
){
  double ll = 0;

  const int n = EXAMS_GRADES.rows();
  const int nq = GRID.rows();
  const int n_cov = COVARIATES.cols();
  const int dim_irt = N_EXAMS*(N_GRADES+3);
  const int dim_lat = 2+2*n_cov;
  Eigen::VectorXd grll = Eigen::VectorXd::Zero(THETA.size());


  Eigen::MatrixXd grid=GRID;
  if(ROTGRID){
    Eigen::MatrixXd L{{1,0},{THETA(dim_irt), THETA(dim_irt+1)}};
    grid = grid * L.transpose();

  }

  Eigen::MatrixXd llMat(n,nq);
  for(int i = 0; i < n; i++){

    double mu_i_ability = 0;
    double mu_i_speed = 0;
    if(ROTGRID){
      mu_i_ability = COVARIATES.row(i)*THETA.segment(dim_irt+2, n_cov);
      mu_i_speed = COVARIATES.row(i)*THETA.segment(dim_irt+2+n_cov, n_cov);
    }



    // Initialize conditional IRT model
    GRTC_MOD irt_mod(THETA,
                     EXAMS_GRADES.row(i),
                     EXAMS_DAYS.row(i),
                     EXAMS_SET.row(i),
                     EXAMS_OBSFLAG.row(i),
                     COVARIATES.row(i),
                     MAX_DAY(i),
                     N_GRADES,
                     N_EXAMS,
                     ROTGRID);

    Eigen::VectorXd f(nq);
    Eigen::MatrixXd gr = Eigen::MatrixXd::Zero(THETA.size(), nq);

    for(int point = 0; point < nq; point++){
      f(point) = exp(irt_mod.ll(grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed));
      llMat(i,point)=(irt_mod.ll(grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed));

      if(GRFLAG){
        Eigen::VectorXd gr_point = irt_mod.grll(grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed);
        gr_point *= f(point);
        gr.col(point) = gr_point;
      }

    }
    double lli = std::max(-10000.0, log(f.dot(WEIGHTS)));

    ll += lli;
    if(GRFLAG){
      grll += gr*WEIGHTS/exp(lli);
    }
  }

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("grid") = grid,
      Rcpp::Named("llMat") = llMat,
      Rcpp::Named("gr") = grll,
      Rcpp::Named("ll") = ll
    );

  return output;

}

//' @export
// [[Rcpp::export]]
Rcpp::List CRGRTCM_GH(
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
    const bool GRFLAG = true,
    const bool ROTGRID = true
){
  double ll = 0;
  Eigen::VectorXd grll = Eigen::VectorXd::Zero(THETA.size());

  const int n = EXAMS_GRADES.rows();
  const int nq = GRID.rows();
  const int dim_irt = N_EXAMS*(N_GRADES+3);



  if(ROTGRID){
    Eigen::MatrixXd L{{1,0},{THETA(dim_irt), THETA(dim_irt+1)}};
    GRID = GRID * L.transpose();
  }

  Eigen::MatrixXd irt_grid(n,nq);
  Eigen::MatrixXd cr_grid(n,nq);
for(int i = 0; i < n; i++){

    // Initialize conditional IRT model
    GRTC_MOD irt_mod(THETA,
                     EXAMS_GRADES.row(i),
                     EXAMS_DAYS.row(i),
                     EXAMS_SET.row(i),
                     EXAMS_OBSFLAG.row(i),
                     EXT_COVARIATES.row(i),
                     MAX_DAY(i),
                     N_GRADES,
                     N_EXAMS,
                     ROTGRID);

    // Initialize conditional CR model
    CR_MOD cr_mod(THETA,
                  OUTCOME(i),
                  EXT_COVARIATES.row(i),
                  YB,
                  YEAR_FIRST(i),
                  YEAR_LAST(i),
                  YEAR_LAST_EXAM(i),
                  ROTGRID);

    Eigen::VectorXd f(nq);
    Eigen::MatrixXd gr = Eigen::MatrixXd::Zero(THETA.size(), nq);

    for(int point = 0; point < nq; point++){
      double irt_cll = irt_mod.ll(GRID(point, 0), GRID(point, 1));
      double cr_cll  = cr_mod.ll( GRID(point, 0), GRID(point, 1));

      f(point) = exp(irt_cll+cr_cll);
      irt_grid(i,point)=irt_cll;
      cr_grid(i,point)=cr_cll;

      if(GRFLAG){
        Eigen::VectorXd gr_point = irt_mod.grll(GRID(point, 0), GRID(point, 1))+cr_mod.grll(GRID(point, 0), GRID(point, 1));
        gr_point *= f(point);
        gr.col(point) = gr_point;
      }

    }
    double lli = std::max(-100.0, log(f.dot(WEIGHTS)));

    ll += lli;
    if(GRFLAG){
      grll += gr*WEIGHTS/exp(lli);
    }
  }

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("irt_grid") = irt_grid,
      Rcpp::Named("cr_grid") = cr_grid,
      Rcpp::Named("gr") = grll,
      Rcpp::Named("ll") = ll
    );

  return output;
}

//' @export
// [[Rcpp::export]]
Rcpp::List CCR(
    Eigen::VectorXd THETA,
    Eigen::VectorXd OUTCOME,
    Eigen::MatrixXd EXT_COVARIATES,
    Eigen::VectorXd YEAR_FIRST,
    Eigen::VectorXd YEAR_LAST,
    Eigen::VectorXd YEAR_LAST_EXAM,
    Eigen::MatrixXd LATMAT,
    const int YB,
    const bool GRFLAG = true
){
  const int n = OUTCOME.size();
  double ll = 0;
  Eigen::VectorXd grll = Eigen::VectorXd::Zero(THETA.size());

  for(int i = 0; i < n; i++){

    // Initialize conditional CR model
    CR_MOD cr_mod(THETA,
                  OUTCOME(i),
                  EXT_COVARIATES.row(i),
                  YB,
                  YEAR_FIRST(i),
                  YEAR_LAST(i),
                  YEAR_LAST_EXAM(i),
                  false);

    ll += cr_mod.ll(LATMAT(i,0), LATMAT(i,1));
    if(GRFLAG){
      grll += cr_mod.grll(LATMAT(i,0), LATMAT(i,1));
    }

  }

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("gr") = grll,
      Rcpp::Named("ll") = ll
    );

  return output;
}


//' @export
// [[Rcpp::export]]
Rcpp::List CRGRTCM_EM(
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
    const std::string MOD
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
    Eigen::MatrixXd L{{1,0},{theta(dim_irt), theta(dim_irt+1)}};
    Eigen::MatrixXd grid = GRID * L.transpose();

    Rcpp::Rcout << "Iter " << iter << " | E-STEP...";
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
                                   grid,
                                   WEIGHTS,
                                   YB,
                                   N_GRADES,
                                   N_EXAMS,
                                   MOD);

    Rcpp::Rcout << " done! | M_STEP...";
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

    eclass.update_quadrature(grid, Ew);

    int status = optim_lbfgs(eclass, theta, enjll, M_MAX_ITER);
    double tol_check = (path_enjll.back() - enjll) / path_enjll.back();

    Rcpp::Rcout << " status="<<status<<  " | obj=" << enjll << ", obj_pdiff:" << tol_check  <<"|\n";

    path_theta.push_back(theta);
    path_enjll.push_back(enjll);

    if(tol_check<TOL){
      last_iter = iter;
      convergence = 1;
      Rcpp::Rcout << "Converged\n";
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


//' @export
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
    Eigen::MatrixXd L{{1,0},{theta(dim_irt), theta(dim_irt+1)}};
    Eigen::MatrixXd grid = GRID * L.transpose();

    if(VERBOSE)Rcpp::Rcout << "Iter " << iter << " | E-STEP...";
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
                                   grid,
                                   WEIGHTS,
                                   YB,
                                   N_GRADES,
                                   N_EXAMS,
                                   MOD);

    if(VERBOSE)Rcpp::Rcout << " done! | M_STEP...";
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

    eclass.update_quadrature(grid, Ew);

    int status = optim_lbfgs(eclass, theta, enjll, M_MAX_ITER);
    double tol_check = (path_enjll.back() - enjll) / path_enjll.back();

    if(VERBOSE)Rcpp::Rcout << " status="<<status<<  " | obj=" << enjll << ", obj_pdiff:" << tol_check  <<"|\n";

    path_theta.push_back(theta);
    path_enjll.push_back(enjll);

    if(tol_check<TOL){
      last_iter = iter;
      convergence = 1;
      if(VERBOSE)Rcpp::Rcout << "Converged\n";
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



//' @export
// [[Rcpp::export]]
Rcpp::List CRGRTCM_GH2(
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
    const bool GRFLAG = true,
    const bool ROTGRID = true
){
  double ll = 0;
  Eigen::VectorXd grll = Eigen::VectorXd::Zero(THETA.size());

  const int n = EXAMS_GRADES.rows();
  const int nq = GRID.rows();
  const int dim_irt = N_EXAMS*(N_GRADES+3);



  if(ROTGRID){
    Eigen::MatrixXd L{{1,0},{THETA(dim_irt), THETA(dim_irt+1)}};
    GRID = GRID * L.transpose();
  }

  Eigen::MatrixXd irt_grid(n,nq);
  Eigen::MatrixXd cr_grid(n,nq);
for(int i = 0; i < n; i++){

    // Initialize conditional IRT model
    GRTC_MOD irt_mod(THETA,
                     EXAMS_GRADES.row(i),
                     EXAMS_DAYS.row(i),
                     EXAMS_SET.row(i),
                     EXAMS_OBSFLAG.row(i),
                     EXT_COVARIATES.row(i),
                     MAX_DAY(i),
                     N_GRADES,
                     N_EXAMS,
                     ROTGRID);

    // Initialize conditional CR model
    CR_MOD cr_mod(THETA,
                  OUTCOME(i),
                  EXT_COVARIATES.row(i),
                  YB,
                  YEAR_FIRST(i),
                  YEAR_LAST(i),
                  YEAR_LAST_EXAM(i),
                  ROTGRID);

    Eigen::VectorXd f(nq);
    Eigen::MatrixXd gr = Eigen::MatrixXd::Zero(THETA.size(), nq);

    for(int point = 0; point < nq; point++){
      double irt_cll = irt_mod.ll(GRID(point, 0), GRID(point, 1));
      double cr_cll  = cr_mod.ll( GRID(point, 0), GRID(point, 1));

      f(point) = exp(irt_cll+cr_cll);
      irt_grid(i,point)=irt_cll;
      cr_grid(i,point)=cr_cll;

      if(GRFLAG){
        Eigen::VectorXd gr_point = irt_mod.grll(GRID(point, 0), GRID(point, 1))+cr_mod.grll(GRID(point, 0), GRID(point, 1));
        gr_point *= f(point);
        gr.col(point) = gr_point;
      }

    }
    double lli = std::max(-10000.0, log(f.dot(WEIGHTS)));

    ll += lli;
    if(GRFLAG){
      grll += gr*WEIGHTS/exp(lli);
    }
  }

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("irt_grid") = irt_grid,
      Rcpp::Named("cr_grid") = cr_grid,
      Rcpp::Named("gr") = grll,
      Rcpp::Named("ll") = ll
    );

  return output;
}
