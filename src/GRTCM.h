#ifndef grtcm_H
#define grtcm_H

#include "grades.h"
#include "times.h"
#include "exams.h"
#include "extractParams.h"
#include "latent.h"

namespace grtcm{


  // GRADED RESPONSE TIME CENSORED CLASS
  class GRTC
  {
  private:
    Eigen::VectorXd _theta;
    Eigen::VectorXd _exams_grades;
    Eigen::VectorXd _exams_days;
    Eigen::VectorXd _exams_set;
    Eigen::VectorXd _exams_obsflag;
    Eigen::VectorXd _covariates;
    unsigned int _max_day;
    unsigned int _n_grades;
    unsigned int _n_exams;
    unsigned int _dim_irt;
    unsigned int _dim_lat;
    bool _latparflag;

  public:
    GRTC(Eigen::VectorXd THETA,
             Eigen::VectorXd EXAMS_GRADES,
             Eigen::VectorXd EXAMS_DAYS,
             Eigen::VectorXd EXAMS_SET,
             Eigen::VectorXd EXAMS_OBSFLAG,
             Eigen::VectorXd COVARIATES,
             const int MAX_DAY,
             const int N_GRADES,
             const int N_EXAMS,
             const bool LATPARFLAG):
    _theta(THETA),
    _exams_grades(EXAMS_GRADES),
    _exams_days(EXAMS_DAYS),
    _exams_set(EXAMS_SET),
    _exams_obsflag(EXAMS_OBSFLAG),
    _covariates(COVARIATES),
    _max_day(MAX_DAY),
    _n_grades(N_GRADES),
    _n_exams(N_EXAMS),
    _latparflag(LATPARFLAG)
    {
      _dim_irt = N_EXAMS*(N_GRADES+3);
      _dim_lat = 2+2*(COVARIATES.size());
    }

    // conditional log-likelihood
    double ll(const double ABILITY, const double SPEED);

    //gradient of conditional log-likelihood
    Eigen::VectorXd grll(const double ABILITY, const double SPEED);


    // complete log-likelihood
    double cll(const double ABILITY, const double SPEED);

    //gradient of complete log-likelihood
    Eigen::VectorXd grcll(const double ABILITY, const double SPEED);


  };

  double GRTC::ll(const double ABILITY, const double SPEED) {

    double out = 0;

    for(unsigned int exam = 0; exam < _n_exams; exam++){

      if(_exams_set[exam]){
        double ell = exams::examLik(exam,
                                    _exams_grades(exam),
                                    _exams_days(exam),
                                    _max_day,
                                    _exams_obsflag(exam),
                                    _theta,
                                    _n_grades,
                                    _n_exams,
                                    ABILITY, SPEED, 1);
        out+=ell;

      }
    }

    return out;
  }
  double GRTC::cll(const double ABILITY, const double SPEED) {

    latent::LAT lat(_theta.segment(_dim_irt, _dim_lat), _covariates);
    double out = lat.ll(ABILITY, SPEED);

    for(unsigned int exam = 0; exam < _n_exams; exam++){

      if(_exams_set[exam]){
        out += exams::examLik(exam,
                              _exams_grades(exam),
                              _exams_days(exam),
                              _max_day,
                              _exams_obsflag(exam),
                              _theta,
                              _n_grades,
                              _n_exams,
                              ABILITY, SPEED, 1);

      }
    }

    return out;
  }
  Eigen::VectorXd GRTC::grll(const double ABILITY, const double SPEED){

    Eigen::VectorXd gr = Eigen::VectorXd::Zero(_theta.size());
    Eigen::VectorXd gr_irt = Eigen::VectorXd::Zero(_dim_irt+_dim_lat);

    for(unsigned int exam = 0; exam < _n_exams; exam++){

      if(_exams_set[exam]){
        gr_irt += exams::grad::grl_examLik(exam,
                                           _exams_grades(exam),
                                           _exams_days(exam),
                                           _max_day,
                                           _exams_obsflag(exam),
                                           _theta,
                                           _covariates,
                                           _n_grades,
                                           _n_exams,
                                           ABILITY, SPEED,
                                           _latparflag);
      }
    }

    gr.segment(0, _dim_irt+_dim_lat) = gr_irt;
    return gr;
  }

  Eigen::VectorXd GRTC::grcll(const double ABILITY, const double SPEED){

    Eigen::VectorXd gr = Eigen::VectorXd::Zero(_theta.size());
    Eigen::VectorXd gr_irt = Eigen::VectorXd::Zero(_dim_irt+_dim_lat);

    latent::LAT lat(_theta.segment(_dim_irt, _dim_lat), _covariates);


    for(unsigned int exam = 0; exam < _n_exams; exam++){

      if(_exams_set[exam]){
        gr_irt += exams::grad::grl_examLik(exam,
                                           _exams_grades(exam),
                                           _exams_days(exam),
                                           _max_day,
                                           _exams_obsflag(exam),
                                           _theta,
                                           _covariates,
                                           _n_grades,
                                           _n_exams,
                                           ABILITY, SPEED,
                                           0);
      }
    }

    gr.segment(0, _dim_irt+_dim_lat) = gr_irt;
    gr.segment(_dim_irt, _dim_lat) += lat.grll(ABILITY, SPEED);
    return gr;
  }


  // namespace gq{
  //
  //   // Evaluate marginal loglikelihood and its gradient
  //   Rcpp::List evaluate_grid(
  //       const Eigen::Ref<const Eigen::VectorXd> THETA,
  //       const Eigen::Ref<const Eigen::MatrixXd> EXAMS_GRADES,
  //       const Eigen::Ref<const Eigen::MatrixXd> EXAMS_DAYS,
  //       const Eigen::Ref<const Eigen::MatrixXd> EXAMS_SET,
  //       const Eigen::Ref<const Eigen::MatrixXd> EXAMS_OBSFLAG,
  //       const Eigen::Ref<const Eigen::MatrixXd> COVARIATES,
  //       const Eigen::Ref<const Eigen::VectorXd> MAX_DAY,
  //       const Eigen::Ref<const Eigen::MatrixXd> GRID,
  //       const Eigen::Ref<const Eigen::VectorXd> WEIGHTS,
  //       const int N_GRADES,
  //       const int N_EXAMS,
  //       const bool GRFLAG = true,
  //       const bool ROTGRID = true
  //   ){
  //     double ll = 0;
  //
  //     const int n = EXAMS_GRADES.rows();
  //     const int nq = GRID.rows();
  //     const int n_cov = COVARIATES.cols();
  //     const int dim_irt = N_EXAMS*(N_GRADES+3);
  //     const int dim_lat = 2+2*n_cov;
  //     Eigen::VectorXd grll = Eigen::VectorXd::Zero(THETA.size());
  //
  //
  //     Eigen::MatrixXd grid=GRID;
  //     if(ROTGRID){
  //       Eigen::MatrixXd L{{1,0},{THETA(dim_irt), THETA(dim_irt+1)}};
  //       grid = grid * L.transpose();
  //     }
  //
  //     Eigen::MatrixXd llMat(n,nq);
  //     for(int i = 0; i < n; i++){
  //
  //       double mu_i_ability = 0;
  //       double mu_i_speed = 0;
  //       if(ROTGRID){
  //         mu_i_ability = COVARIATES.row(i)*THETA.segment(dim_irt+2, n_cov);
  //         mu_i_speed = COVARIATES.row(i)*THETA.segment(dim_irt+2+n_cov, n_cov);
  //       }
  //
  //
  //
  //       // Initialize conditional IRT model
  //       GRTC grtc(THETA,
  //                 EXAMS_GRADES.row(i),
  //                 EXAMS_DAYS.row(i),
  //                 EXAMS_SET.row(i),
  //                 EXAMS_OBSFLAG.row(i),
  //                 COVARIATES.row(i),
  //                 MAX_DAY(i),
  //                 N_GRADES,
  //                 N_EXAMS,
  //                 ROTGRID);
  //
  //       Eigen::VectorXd f(nq);
  //       Eigen::MatrixXd gr = Eigen::MatrixXd::Zero(THETA.size(), nq);
  //
  //       for(int point = 0; point < nq; point++){
  //         f(point) = exp(grtc.ll(grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed));
  //         llMat(i,point)=(grtc.ll(grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed));
  //
  //         if(GRFLAG){
  //           Eigen::VectorXd gr_point = grtc.grll(grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed);
  //           gr_point *= f(point);
  //           gr.col(point) = gr_point;
  //         }
  //
  //       }
  //       double lli = std::max(-10000.0, log(f.dot(WEIGHTS)));
  //
  //       ll += lli;
  //       if(GRFLAG){
  //         grll += gr*WEIGHTS/exp(lli);
  //       }
  //     }
  //
  //     Rcpp::List output =
  //       Rcpp::List::create(
  //         Rcpp::Named("grid") = grid,
  //         Rcpp::Named("llMat") = llMat,
  //         Rcpp::Named("gr") = grll,
  //         Rcpp::Named("ll") = ll
  //       );
  //
  //     return output;
  //
  //   }
  // }
  //
  //
  // namespace em{
  //
  // }


}



#endif
