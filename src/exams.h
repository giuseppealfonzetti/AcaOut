#ifndef exams_H
#define exams_H
#include "grades.h"
#include "times.h"
#include "extractParams.h"

namespace exams{
  //' Evaluate exam specific likelihood
  //'
  //' @param EXAM Exam of interest. Integers from 0 to N_EXAMS -1.
  //' @param GRADE Grade used as reference. Integers from 1 to N_GRADES.
  //' @param DAY Day of interest.
  //' @param MAX_DAY Last day of observation.
  //' @param OBSFLAG TRUE for observed, FALSE for not-observed.
  //' @param THETA_IRT Portion of the parameter vector related to the IRT model
  //' @param N_GRADES Number of grades modelled.
  //' @param N_EXAMS Number of exams.
  //' @param ABILITY ability value.
  //' @param SPEED speed value.
  //' @param LOGFLAG Set TRUE to return log value.
  //'
  //' @returns It returns the probability of observing or not a specific
  //' grade on a given exam before a given day conditioned on ability and speed.
  double examLik(
      const unsigned int EXAM,
      const unsigned int GRADE,
      const double DAY,
      const int MAX_DAY,
      const bool OBSFLAG,
      const Eigen::Ref<const Eigen::VectorXd> THETA,
      const unsigned int N_GRADES,
      const unsigned int N_EXAMS,
      const double ABILITY,
      const double SPEED,
      const bool LOGFLAG = false
  ){
    double out, logpExam, logpTime;


    if(OBSFLAG){
      logpExam = exams::pGrade(GRADE, EXAM, THETA, N_GRADES, N_EXAMS, ABILITY, true);
      logpTime = exams::pTimeExam(EXAM, DAY, THETA, N_GRADES,  N_EXAMS, SPEED, false, true);
      out = logpExam+logpTime;
    }else{
      logpExam = exams::pGreaterGrades(1, EXAM, THETA, N_GRADES, N_EXAMS, ABILITY, true);
      logpTime = exams::pTimeExam(EXAM, MAX_DAY, THETA, N_GRADES,  N_EXAMS, SPEED, true, true);
      // Rcpp::Rcout<<"\nexamLik |:"<< ABILITY <<", "<< SPEED<<"\n";
      // Rcpp::Rcout<<"MAX_DAY |:"<< MAX_DAY <<"\n";
      // Rcpp::Rcout << "logpExam "<< logpExam<< " | logpTime "<< logpTime<<"\n";

      out = log1mexp(-logpExam-logpTime);
      out = std::max(-10000.0, out);  // avoid log(0)
    }

    if(LOGFLAG){
      return(out);
    }else{
      return(exp(out));
    }

  }
}

namespace exams::grad{
 //' Evaluate the gradient of exam-specific log-likelihood
 //'
 //' @param EXAM Exam of interest. Integers from 0 to N_EXAMS -1.
 //' @param GRADE Grade used as reference. Integers from 1 to N_GRADES.
 //' @param DAY Day of interest.
 //' @param MAX_DAY Last day of observation.
 //' @param OBSFLAG TRUE for observed, FALSE for not-observed.
 //' @param THETA_IRT Portion of the parameter vector related to the IRT model
 //' @param N_GRADES Number of grades modelled.
 //' @param N_EXAMS Number of exams.
 //' @param ABILITY ability value.
 //' @param SPEED speed value.
 //' @param ROTATED Have the latent scores been rotated using their covariance matrix?
 //'
 //' @returns It returns the probability of observing or not a specific
 //' grade on a given exam before a given day conditioned on ability and speed.
 //'
 Eigen::VectorXd grl_examLik(
     const unsigned int EXAM,
     const unsigned int GRADE,
     const double DAY,
     const double MAX_DAY,
     const bool OBSFLAG,
     const Eigen::Ref<const Eigen::VectorXd> THETA,
     const Eigen::Ref<const Eigen::VectorXd> COVARIATES,
     const unsigned int N_GRADES,
     const unsigned int N_EXAMS,
     const double ABILITY,
     const double SPEED,
     const bool LATPARFLAG
 ){
   double logp, logpGrade, logpTime;
   const int dim_irt = N_EXAMS*(N_GRADES+3);
   const int dim_lat = 2+2*COVARIATES.size();

   Eigen::VectorXd gr = Eigen::VectorXd::Zero(dim_irt+dim_lat);


   if(OBSFLAG){
     logpGrade = exams::pGrade(GRADE, EXAM, THETA, N_GRADES, N_EXAMS, ABILITY, true);
     logpTime = exams::pTimeExam(EXAM, DAY, THETA, N_GRADES,  N_EXAMS, SPEED, false, true);
     Eigen::VectorXd grlTime = exams::grad::gr_pTimeExam2(EXAM, DAY, THETA, COVARIATES, N_GRADES, N_EXAMS, SPEED, ABILITY, false, LATPARFLAG, true);
     // Rcpp::Rcout << "Speed:"<< SPEED<<", Ability:"<<ABILITY<<", l1:"<<THETA(dim_irt)<<", l2:"<< THETA(dim_irt+1);
     // Rcpp::Rcout << ", grl1:"<< grlTime(dim_irt)<< ", grl2:"<< grlTime(dim_irt+1)<< "\n";
     Eigen::VectorXd grGrade = exams::grad::gr_pGrade2(GRADE,EXAM, THETA, COVARIATES, N_GRADES, N_EXAMS, ABILITY, LATPARFLAG);
     // Rcpp::Rcout << "Grades:"<< grGrade.segment(dim_irt,2)<< "\n";

     gr = grlTime + grGrade /exp(logpGrade);
     logp = logpGrade+logpTime;
   }else{
     logpGrade = exams::pGreaterGrades(1, EXAM, THETA, N_GRADES, N_EXAMS, ABILITY, true);
     logpTime = exams::pTimeExam(EXAM, MAX_DAY, THETA, N_GRADES,  N_EXAMS, SPEED, true, true);
     Eigen::VectorXd grGrade = exams::grad::gr_pGreaterGrades2(1,EXAM, THETA, COVARIATES, N_GRADES, N_EXAMS, ABILITY, LATPARFLAG);
     Eigen::VectorXd grTime = exams::grad::gr_pTimeExam2(EXAM, MAX_DAY, THETA, COVARIATES, N_GRADES, N_EXAMS, SPEED, ABILITY, true, LATPARFLAG, false);
     gr = -(grTime*exp(logpGrade) + exp(logpTime)*grGrade);
     logp = log1mexp(-logpGrade-logpTime);
     logp = std::max(-1000.0, logp);  // avoid log(0)
     gr/= exp(logp);
   }

   // gr/=std::max(1e-20, exp(logp));
   // logp/= 10;

   // Eigen::VectorXi is_selected = (gr.array() != 0).cast<int>();
   // gr/= exp(logp);


   return(gr);

 }
}


#endif
