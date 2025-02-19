#ifndef latent_H
#define latent_H

const double pi = 3.1415926535897;
const double neglog2pi = -1.837877;

class LAT_DISTR{
  private:
    Eigen::VectorXd _theta;

  public:
    LAT_DISTR(Eigen::VectorXd THETA):
    _theta(THETA){}

  double ll(const double ABILITY, const double SPEED);

  Eigen::VectorXd grll(const double ABILITY, const double SPEED);
};

double LAT_DISTR::ll(const double ABILITY, const double SPEED){
  double ll;
  Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
  Eigen::MatrixXd L{{1,0},{_theta(0), _theta(1)}};
  Eigen::MatrixXd cov = L*L.transpose();
  double det = cov.determinant();
  Eigen::MatrixXd inv_cov = cov.inverse();
  ll = neglog2pi-0.5*log(det) - 0.5 * double(lat.transpose()*inv_cov*lat);
  return(ll);
}

Eigen::VectorXd LAT_DISTR::grll(const double ABILITY, const double SPEED){
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(2);
  Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
  Eigen::MatrixXd L{{1,0},{_theta(0), _theta(1)}};
  Eigen::MatrixXd iL = L.inverse();
  Eigen::MatrixXd Z1 = Eigen::MatrixXd::Zero(2,2); Z1(1,0) = 1;
  Eigen::MatrixXd Z2 = Eigen::MatrixXd::Zero(2,2); Z2(1,1) = 1;

  Eigen::MatrixXd M1 = -iL.transpose()*Z1.transpose()*iL.transpose()*iL;
  Eigen::MatrixXd M2 = -iL.transpose()*Z2.transpose()*iL.transpose()*iL;

  gr(0) = -(iL*Z1).trace() - .5*lat.transpose()*(M1 + M1.transpose())*lat;
  gr(1) = -(iL*Z2).trace() - .5*lat.transpose()*(M2 + M2.transpose())*lat;
  return(gr);
}

//' Latent distribution
//'
//' Provides a wrapper for the cpp class `LAT_DISTR`, which computes
//' log density and gradient of the joint distribution of ability and speed.
//'
//' @param ABILITY ability value.
//' @param SPEED speed value.
//' @param PARS Two elements, used as (`L[2,1]`, `L[2,2]`), where L is the lower triangular Cholesky of the latent covariance matrix.
//' @param GRFLAG `TRUE` to compute the gradient.
//'
//' @return It returns a list with:
//' \itemize{
//'   \item ll - The log-likelihood of the latent distribution.
//'   \item gr - The gradient of the log-likelihood.
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List internal_lat(
    const double ABILITY,
    const double SPEED,
    Eigen::VectorXd PARS,
    const bool GRFLAG = true
){

  LAT_DISTR lat(PARS);

  double ll = lat.ll(ABILITY, SPEED);
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(2);
  if(GRFLAG){
    gr = lat.grll(ABILITY, SPEED);
  }

  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("gr") = gr,
      Rcpp::Named("ll") = ll
    );

  return output;
 }


class LAT{
private:
  Eigen::VectorXd _theta;
  unsigned int _dim_irt;
  bool _rotated;

public:
  LAT(Eigen::VectorXd THETA,
      const unsigned int DIM_IRT,
      const bool ROTATED):
  _theta(THETA), _dim_irt(DIM_IRT), _rotated(ROTATED){}

  double ll(const double ABILITY, const double SPEED);

  Eigen::VectorXd grll(const double ABILITY, const double SPEED);
};

double LAT::ll(const double ABILITY, const double SPEED){
  double ll;
  Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
  Eigen::MatrixXd L{{1,0},{_theta(_dim_irt), _theta(_dim_irt+1)}};
  Eigen::MatrixXd cov = L*L.transpose();
  double det = cov.determinant();
  Eigen::MatrixXd inv_cov = cov.inverse();
  ll = neglog2pi-0.5*log(det) - 0.5 * double(lat.transpose()*inv_cov*lat);
  return(ll);
}

Eigen::VectorXd LAT::grll(const double ABILITY, const double SPEED){
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(_theta.size());
  Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
  Eigen::MatrixXd L{{1,0},{_theta(_dim_irt), _theta(_dim_irt+1)}};
  Eigen::MatrixXd iL = L.inverse();
  Eigen::MatrixXd iS = iL.transpose()*iL;

  Eigen::MatrixXd Z1 = Eigen::MatrixXd::Zero(2,2); Z1(1,0) = 1;
  Eigen::MatrixXd Z2 = Eigen::MatrixXd::Zero(2,2); Z2(1,1) = 1;

  // M1 + M1.transpose is dS^-1 wrt l1
  Eigen::MatrixXd M1 = -iL.transpose()*Z1.transpose()*iL.transpose()*iL;
  Eigen::MatrixXd M2 = -iL.transpose()*Z2.transpose()*iL.transpose()*iL;

  gr(_dim_irt)   = -(iL*Z1).trace() - .5*lat.transpose()*(M1 + M1.transpose())*lat;
  gr(_dim_irt+1) = -(iL*Z2).trace() - .5*lat.transpose()*(M2 + M2.transpose())*lat;

  if(_rotated){
    double speed_raw = (SPEED -_theta(_dim_irt)*ABILITY)/_theta(_dim_irt+1);

    Eigen::VectorXd dlat1(2); dlat1 << 0, ABILITY;
    Eigen::VectorXd dlat2(2); dlat2 << 0, speed_raw;

    gr(_dim_irt) -= dlat1.transpose()*iS*lat;
    gr(_dim_irt+1) -= dlat2.transpose()*iS*lat;


  }


  return(gr);
}



class LAT2{
private:
  Eigen::VectorXd _theta;
  Eigen::VectorXd _covariates;
  int _dim_irt;
  int _dim_lat;
  int _n_cov;
  bool _latparflag;
  Eigen::MatrixXd _L = Eigen::MatrixXd::Zero(2,2);
  Eigen::MatrixXd _B;
  Eigen::VectorXd _mean;

public:
  LAT2(const Eigen::Ref<const Eigen::VectorXd> THETA,
       const Eigen::Ref<const Eigen::VectorXd> COVARIATES,
       const int DIM_IRT,
       const bool LATPARFLAG):
  _theta(THETA), _covariates(COVARIATES), _dim_irt(DIM_IRT),  _latparflag(LATPARFLAG){
    _n_cov = _covariates.size();
    _dim_lat = 2+2*_n_cov;
    _L(0,0)= 1;
    _L(1,0)= _theta(_dim_irt);
    _L(1,1)= _theta(_dim_irt+1);
    _B = THETA.segment(_dim_irt+2, 2*_n_cov).reshaped(2,_n_cov);
    _mean = _B * COVARIATES;
  }

  double ll(const double ABILITY, const double SPEED);

  Eigen::VectorXd grll(const double ABILITY, const double SPEED);
};

double LAT2::ll(const double ABILITY, const double SPEED){
  double ll;
  Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
  Eigen::MatrixXd cov = _L*_L.transpose();
  double det = cov.determinant();
  Eigen::MatrixXd inv_cov = cov.inverse();
  ll = neglog2pi-0.5*log(det) - 0.5 * double((lat-_mean).transpose()*inv_cov*(lat-_mean));
  return(ll);
}

Eigen::VectorXd LAT2::grll(const double ABILITY, const double SPEED){
  Eigen::VectorXd gr = Eigen::VectorXd::Zero(_theta.size());
  Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
  Eigen::MatrixXd L{{1,0},{_theta(_dim_irt), _theta(_dim_irt+1)}};
  Eigen::MatrixXd iL = L.inverse();
  Eigen::MatrixXd iS = iL.transpose()*iL;

  Eigen::MatrixXd Z1 = Eigen::MatrixXd::Zero(2,2); Z1(1,0) = 1;
  Eigen::MatrixXd Z2 = Eigen::MatrixXd::Zero(2,2); Z2(1,1) = 1;

  // M1 + M1.transpose is dS^-1 wrt l1
  Eigen::MatrixXd M1 = -iL.transpose()*Z1.transpose()*iL.transpose()*iL;
  Eigen::MatrixXd M2 = -iL.transpose()*Z2.transpose()*iL.transpose()*iL;

  gr(_dim_irt)   = -(iL*Z1).trace() - .5*lat.transpose()*(M1 + M1.transpose())*lat;
  gr(_dim_irt+1) = -(iL*Z2).trace() - .5*lat.transpose()*(M2 + M2.transpose())*lat;

  if(_latparflag){
    double ability_raw = ABILITY-(_theta.segment(_dim_irt+2, _n_cov).transpose()*-_covariates);
    double speed_raw = (SPEED -_theta(_dim_irt)*ability_raw-(_theta.segment(_dim_irt+2+_n_cov, _n_cov).transpose()*_covariates))/_theta(_dim_irt+1);


    Eigen::VectorXd dlat1(2); dlat1 << 0, ability_raw;
    Eigen::VectorXd dlat2(2); dlat2 << 0, speed_raw;

    gr(_dim_irt) -= dlat1.transpose()*iS*lat;
    gr(_dim_irt+1) -= dlat2.transpose()*iS*lat;


  }


  return(gr);
}






namespace latent{


  class LAT{
  private:
    Eigen::VectorXd theta_lat;
    Eigen::VectorXd covariates;
    const int n_cov;

  public:
    LAT(const Eigen::Ref<const Eigen::VectorXd> THETA_LAT_,
        const Eigen::Ref<const Eigen::VectorXd> COVARIATES_):
    theta_lat(THETA_LAT_), covariates(COVARIATES_), n_cov(COVARIATES_.size()){}

    double ll(const double ABILITY, const double SPEED);

    Eigen::VectorXd grll(const double ABILITY, const double SPEED);
  };

  double LAT::ll(const double ABILITY, const double SPEED){
    double ll;
    Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
    Eigen::MatrixXd L{{1,0},{theta_lat(0), theta_lat(1)}};
    Eigen::MatrixXd cov = L*L.transpose();
    double det = cov.determinant();
    Eigen::MatrixXd inv_cov = cov.inverse();
    Eigen::MatrixXd B(2, n_cov);
    B.row(0) = theta_lat.segment(2, n_cov);
    B.row(1) = theta_lat.segment(2+n_cov, n_cov);
    Eigen::VectorXd mu = B*covariates;
    ll = neglog2pi-0.5*log(det) - 0.5 * double((lat-mu).transpose()*inv_cov*(lat-mu));
    return(ll);
  }

  Eigen::VectorXd LAT::grll(const double ABILITY, const double SPEED){
    Eigen::VectorXd gr = Eigen::VectorXd::Zero(2+2*n_cov);
    Eigen::VectorXd lat(2); lat << ABILITY, SPEED;
    Eigen::MatrixXd L{{1,0},{theta_lat(0), theta_lat(1)}};
    Eigen::MatrixXd B(2, n_cov);
    B.row(0) = theta_lat.segment(2, n_cov);
    B.row(1) = theta_lat.segment(2+n_cov, n_cov);
    Eigen::VectorXd mu = B*covariates;
    Eigen::VectorXd latc = lat-mu;

    Eigen::MatrixXd iL = L.inverse();
    Eigen::MatrixXd Z1 = Eigen::MatrixXd::Zero(2,2); Z1(1,0) = 1;
    Eigen::MatrixXd Z2 = Eigen::MatrixXd::Zero(2,2); Z2(1,1) = 1;

    Eigen::MatrixXd M1 = -iL.transpose()*Z1.transpose()*iL.transpose()*iL;
    Eigen::MatrixXd M2 = -iL.transpose()*Z2.transpose()*iL.transpose()*iL;

    gr(0) = -(iL*Z1).trace() - .5*latc.transpose()*(M1 + M1.transpose())*latc;
    gr(1) = -(iL*Z2).trace() - .5*latc.transpose()*(M2 + M2.transpose())*latc;
    gr.segment(2, n_cov) = ((1+pow(L(1,0),2)/pow(L(1,1),2))*latc(0) - (L(1,0)/pow(L(1,1),2))*latc(1) ) *covariates;
    gr.segment(2+n_cov, n_cov) = ((1/pow(L(1,1),2))*(latc(1))- (L(1,0)/pow(L(1,1),2))*latc(0)) *covariates;

    return(gr);
  }

}
#endif





























