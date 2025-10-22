#ifndef em_H
#define em_H
#include "grades.h"
#include "times.h"
#include "exams.h"
#include "extractParams.h"
#include "latent.h"
#include "cr.h"
#include "conditional_models.h"
#include "GRTCM.h"
#include <RcppNumerical.h>

// Define global variable used to count gradient evaluations in M-Step
int MCOUNT;

namespace EM
{
  Eigen::MatrixXd Estep(
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
      const unsigned int YB,
      const unsigned int N_GRADES,
      const unsigned int N_EXAMS,
      std::string MOD,
      const bool LATPARFLAG = true
  ){

    const unsigned int n = EXAMS_GRADES.rows();
    const int n_cov = EXT_COVARIATES.cols();
    const unsigned int nq = GRID.rows();
    const unsigned int dim_irt = N_EXAMS*(N_GRADES+3);

    Eigen::MatrixXd Emat(n,nq);

    Eigen::MatrixXd grid=GRID;
    if(LATPARFLAG){
      Eigen::MatrixXd L{{1,0},{THETA(dim_irt), THETA(dim_irt+1)}};
      grid = grid * L.transpose();
    }

    for(unsigned int i = 0; i < n; i++){

      double mu_i_ability = 0;
      double mu_i_speed = 0;
      if(LATPARFLAG){
        mu_i_ability = EXT_COVARIATES.row(i)*THETA.segment(dim_irt+2, n_cov);
        mu_i_speed = EXT_COVARIATES.row(i)*THETA.segment(dim_irt+2+n_cov, n_cov);
      }

      // Initialize conditional IRT model
      grtcm::GRTC grtc(THETA,
                       EXAMS_GRADES.row(i),
                       EXAMS_DAYS.row(i),
                       EXAMS_SET.row(i),
                       EXAMS_OBSFLAG.row(i),
                       EXT_COVARIATES.row(i),
                       MAX_DAY(i),
                       N_GRADES,
                       N_EXAMS,
                       false);

      // Initialize conditional CR model
      cr::CCR ccrm(THETA,
                   OUTCOME(i),
                   EXT_COVARIATES.row(i),
                   YB,
                   YEAR_FIRST(i),
                   YEAR_LAST(i),
                   YEAR_LAST_EXAM(i),
                   false);


      for(unsigned int point = 0; point < nq; point++){
        double logd = grtc.ll(grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed);

        if(MOD=="full"){
          logd  += ccrm.ll( grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed);
        }

        Emat(i, point) = exp(logd)*WEIGHTS(point);
      }

      Emat.row(i) /= Emat.row(i).sum();
    }
    return Emat;
  }






  class EAPLOGJ: public Numer::MFuncGrad
  {
  private:
    Eigen::MatrixXd _exams_grades;
    Eigen::MatrixXd _exams_days;
    Eigen::MatrixXd _exams_set;
    Eigen::MatrixXd _exams_obsflag;
    Eigen::VectorXd _max_day;
    Eigen::VectorXd _outcome;
    Eigen::MatrixXd _ext_covariates;
    Eigen::VectorXd _year_first;
    Eigen::VectorXd _year_last;
    Eigen::VectorXd _year_last_exam;


    unsigned int _yb;
    unsigned int _n_grades;
    unsigned int _n_exams;
    unsigned int _n;
    unsigned int _nq;

    std::string _mod;

  public:
    Eigen::MatrixXd _grid;
    Eigen::MatrixXd _eweights;
    EAPLOGJ(Eigen::MatrixXd EXAMS_GRADES,
          Eigen::MatrixXd EXAMS_DAYS,
          Eigen::MatrixXd EXAMS_SET,
          Eigen::MatrixXd EXAMS_OBSFLAG,
          Eigen::VectorXd MAX_DAY,
          Eigen::VectorXd OUTCOME,
          Eigen::MatrixXd EXT_COVARIATES,
          Eigen::VectorXd YEAR_FIRST,
          Eigen::VectorXd YEAR_LAST,
          Eigen::VectorXd YEAR_LAST_EXAM,
          const unsigned int YB,
          const unsigned int N_GRADES,
          const unsigned int N_EXAMS,
          const std::string MOD):
      _exams_grades(EXAMS_GRADES),
      _exams_days(EXAMS_DAYS),
      _exams_set(EXAMS_SET),
      _exams_obsflag(EXAMS_OBSFLAG),
      _max_day(MAX_DAY),
      _outcome(OUTCOME),
      _ext_covariates(EXT_COVARIATES),
      _year_first(YEAR_FIRST),
      _year_last(YEAR_LAST),
      _year_last_exam(YEAR_LAST_EXAM),
      _yb(YB),
      _n_grades(N_GRADES),
      _n_exams(N_EXAMS),
      _mod(MOD){
      _n = EXAMS_GRADES.rows();
    }

    void update_quadrature(Eigen::MatrixXd GRID,
                           Eigen::MatrixXd EWEIGHTS){
      _grid    = GRID;
      _eweights= EWEIGHTS;
      _nq      = GRID.rows();
    }

    double irt_cll(Eigen::VectorXd theta, double ABILITY, double SPEED){

      double nll=0;
      for(unsigned int i = 0; i < _n; i++){

        // Initialize conditional IRT model
        grtcm::GRTC grtcm(theta,
                         _exams_grades.row(i),
                         _exams_days.row(i),
                         _exams_set.row(i),
                         _exams_obsflag.row(i),
                         _ext_covariates.row(i),
                         _max_day(i),
                         _n_grades,
                         _n_exams,
                         false);


        nll-= grtcm.cll(ABILITY, SPEED);




      }
      return nll;

      }

    double f_grad(Numer::Constvec& theta, Numer::Refvec grad){

      double nll = 0;
      const int dim_irt = _n_exams*(_n_grades+3);
      const int n_cov = _ext_covariates.cols();

      Eigen::MatrixXd grid=_grid;
      Eigen::MatrixXd L{{1,0},{theta(dim_irt), theta(dim_irt+1)}};
      grid = grid * L.transpose();

      Eigen::VectorXd gr = Eigen::VectorXd::Zero(theta.size());
      for(unsigned int i = 0; i < _n; i++){
        Rcpp::checkUserInterrupt();

        double mu_i_ability = 0;
        double mu_i_speed = 0;
        mu_i_ability = _ext_covariates.row(i)*theta.segment(dim_irt+2, n_cov);
        mu_i_speed = _ext_covariates.row(i)*theta.segment(dim_irt+2+n_cov, n_cov);

        // Initialize conditional IRT model
        grtcm::GRTC grtcm(theta,
                         _exams_grades.row(i),
                         _exams_days.row(i),
                         _exams_set.row(i),
                         _exams_obsflag.row(i),
                         _ext_covariates.row(i),
                         _max_day(i),
                         _n_grades,
                         _n_exams,
                         false);

        // Initialize conditional CR model
        cr::CCR ccrm(theta,
                      _outcome(i),
                      _ext_covariates.row(i),
                      _yb,
                      _year_first(i),
                      _year_last(i),
                      _year_last_exam(i),
                      false);


        for(unsigned int point = 0; point < _nq; point++){
          double logd = grtcm.cll(grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed);
          Eigen::VectorXd grlogd = grtcm.grcll( grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed);

          if(_mod=="full"){
            logd  += ccrm.ll( grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed);
            grlogd += ccrm.grll( grid(point, 0)+mu_i_ability, grid(point, 1)+mu_i_speed);
          }

          nll -= logd*_eweights(i, point);
          gr -= grlogd*_eweights(i, point);

        }
      }

      MCOUNT++;
      Rcpp::Rcout << "\r- M-STEP... Grad evals:"<<MCOUNT <<"| nll:"<< nll/_n;

      grad = gr/_n;
      return nll/_n;


    }






  };
}

#endif




















