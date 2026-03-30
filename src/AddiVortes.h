#ifndef SOFT_ADDI_VORTES_H
#define SOFT_ADDI_VORTES_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#include <RcppEigen.h>
#pragma GCC diagnostic pop

#include <vector>
#include <string>
#include "functions.h"

#ifndef M_2PI
#define M_2PI 6.283185307179586476925286766559
#endif

// Mirrors softBART's Hypers struct
struct Hypers {
    double sigma;
    double sigma_mu;
    double omega;
    double lambda_poisson;
    double p_shape;
    double p_rate;
    double p_sd;
    double nu;
    double lambda_invgamma;
    int num_covariates;

    void UpdateSigma(const Eigen::VectorXd& res);
};

// Mirrors softBART's Node/Tree class
class Tessellation {
public:
    std::vector<double> tess;
    std::vector<int> dim;
    int n_centres;
    int n_dims;
    double p_val;
    Eigen::VectorXd mu;

    Eigen::MatrixXd cached_sqdist;
    Eigen::MatrixXd cached_weights;

    Tessellation();
    ~Tessellation();

    void GetW(const Eigen::MatrixXd& X, double epsilon = 1e-10);
    void RecomputeW(double epsilon = 1e-10);
};

// Core Mathematical Functions
void GetSuffStats(Tessellation* t, const Eigen::VectorXd& y, const Eigen::MatrixXd& X, 
                  const Hypers& hypers, Eigen::VectorXd& mu_hat_out, Eigen::MatrixXd& Omega_inv_out);

double LogLT(Tessellation* t, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X, const Hypers& hypers);

void UpdateMu(Tessellation* t, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X, const Hypers& hypers);

void UpdateP(Tessellation* t, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X, const Hypers& hypers);

double TransRatio(int old_c, int new_c, int old_d, int new_d, const std::string& mod, const Hypers& hypers);

void perturb_decision_rule(Tessellation* tree, const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, 
                           const Hypers& hypers, double sd_val, const double* mu_val);

void TessBackfit(std::vector<Tessellation*>& forest, Eigen::VectorXd& Y_hat, Hypers& hypers, 
                 const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, double sd_val, const double* mu_val);

#endif