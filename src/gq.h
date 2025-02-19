#ifndef gq_H
#define gq_H

#include "grades.h"
#include "times.h"
#include "exams.h"
#include "extractParams.h"
#include "latent.h"
#include "GRTCM.h"
#include "cr.h"

namespace gq{
  // Evaluate marginal loglikelihood and its gradient on a quadrature grid
  Rcpp::List evaluate_grid(
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const Eigen::Ref<const Eigen::MatrixXd> EXAMS_GRADES,
      const Eigen::Ref<const Eigen::MatrixXd> EXAMS_DAYS,
      const Eigen::Ref<const Eigen::MatrixXd> EXAMS_SET,
      const Eigen::Ref<const Eigen::MatrixXd> EXAMS_OBSFLAG,
      const Eigen::Ref<const Eigen::MatrixXd> COVARIATES,
      const Eigen::Ref<const Eigen::VectorXd> MAX_DAY,
      const Eigen::Ref<const Eigen::VectorXd> OUTCOME,
      const Eigen::Ref<const Eigen::VectorXd> YEAR_FIRST,
      const Eigen::Ref<const Eigen::VectorXd> YEAR_LAST,
      const Eigen::Ref<const Eigen::VectorXd> YEAR_LAST_EXAM,
      const Eigen::Ref<const Eigen::MatrixXd> GRID,
      const Eigen::Ref<const Eigen::VectorXd> WEIGHTS,
      const int N_GRADES,
      const int N_EXAMS,
      const int YB,
      const std::string MOD,
      const bool GRFLAG = true,
      const bool LATPARFLAG = true
  ){
    double ll = 0;

    const int n = EXAMS_GRADES.rows();
    const int nq = GRID.rows();
    const int n_cov = COVARIATES.cols();
    const int dim_irt = N_EXAMS*(N_GRADES+3);
    const int dim_lat = 2+2*n_cov;
    Eigen::VectorXd grll = Eigen::VectorXd::Zero(THETA.size());


    Eigen::MatrixXd grid=GRID;
    if(LATPARFLAG){
      Eigen::MatrixXd L{{1,0},{THETA(dim_irt), THETA(dim_irt+1)}};
      grid = grid * L.transpose();
    }

    Eigen::MatrixXd llMat(n,nq);
    for(int i = 0; i < n; i++){

      double mu_i_ability = 0;
      double mu_i_speed = 0;
      if(LATPARFLAG){
        mu_i_ability = COVARIATES.row(i)*THETA.segment(dim_irt+2, n_cov);
        mu_i_speed = COVARIATES.row(i)*THETA.segment(dim_irt+2+n_cov, n_cov);
      }



      // Initialize conditional IRT model
      grtcm::GRTC grtc(THETA,
                       EXAMS_GRADES.row(i),
                       EXAMS_DAYS.row(i),
                       EXAMS_SET.row(i),
                       EXAMS_OBSFLAG.row(i),
                       COVARIATES.row(i),
                       MAX_DAY(i),
                       N_GRADES,
                       N_EXAMS,
                       LATPARFLAG);

      // Initialize conditional CR model
      cr::CCR ccrm(THETA,
                   OUTCOME(i),
                   COVARIATES.row(i),
                   YB,
                   YEAR_FIRST(i),
                   YEAR_LAST(i),
                   YEAR_LAST_EXAM(i),
                   LATPARFLAG);

      Eigen::VectorXd f(nq);
      Eigen::MatrixXd gr = Eigen::MatrixXd::Zero(THETA.size(), nq);

      for(int point = 0; point < nq; point++){

        if(MOD=="full"){

          double irt_cll = grtc.ll(GRID(point, 0)+mu_i_ability, GRID(point, 1)+mu_i_speed);
          double cr_cll  = ccrm.ll( GRID(point, 0)+mu_i_ability, GRID(point, 1)+mu_i_speed);

          f(point) = exp(irt_cll+cr_cll);
          llMat(i,point) = irt_cll+cr_cll;

        }else{
          f(point) = exp(grtc.ll(grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed));
          llMat(i,point)=(grtc.ll(grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed));
        }


        if(GRFLAG){
          if(MOD=="full"){
            Eigen::VectorXd gr_point = grtc.grll(GRID(point, 0)+mu_i_ability, GRID(point, 1)+mu_i_speed)+ccrm.grll(GRID(point, 0)+mu_i_ability, GRID(point, 1)+mu_i_speed);
            gr_point *= f(point);
            gr.col(point) = gr_point;
          }else{
            Eigen::VectorXd gr_point = grtc.grll(grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed);
            gr_point *= f(point);
            gr.col(point) = gr_point;
          }

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
}
#endif
