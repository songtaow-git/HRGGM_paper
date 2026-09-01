#ifndef _MATRIX_COMPUTE_H_
#define _MATRIX_COMPUTE_H_
// VERSION: 2026-02-12 01

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <stdexcept>
#include <algorithm>
#include <cstddef>
#include "Class_define.h"

double Logdet_sparse(const Eigen::SparseMatrix<double>& X);
double Logdet_sparse(const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& chol);
double Logdet_dense(const Eigen::MatrixXd& X);

double Trace_dense_sparse(const Eigen::MatrixXd& S, const Eigen::SparseMatrix<double>& X);
double Trace_dense_dense(const Eigen::MatrixXd& S, const Eigen::MatrixXd& X);

double Likelihood_f(const Eigen::MatrixXd& S, const Eigen::SparseMatrix<double>& X);
double Likelihood_f(const Eigen::MatrixXd& S, const Eigen::SparseMatrix<double>& X
, const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& chol);

double Offdiag_l1_norm_ll(
    const Eigen::SparseMatrix<double>& X, const int Nl);

double Block_l2_norm_lg(const Eigen::SparseMatrix<double>& X,
    const int Nl, const int Nh, const int Ng);

double Block_l1_norm_lh(const Eigen::SparseMatrix<double>& X, const int Nl, const int Nh);

double Block_l2_norm_lh(const Eigen::SparseMatrix<double>& X, const int Nl, const int Nh);

RowMat Inverse_spd_sparse(const Eigen::SparseMatrix<double>& X,
	const int blockCols = 1024);

RowMat Inverse_spd_sparse_from_chol(
    const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& chol,
    const int n, const int blockCols = 1024);

std::vector<std::pair<int, int>> Find_update_entries_ll(const Eigen::SparseMatrix<double>& X_t,
    const Eigen::MatrixXd& Nabla_g, const double Lambda, const int Nl);

std::vector<std::pair<int, int>> Find_update_entries_lh(
    const Eigen::SparseMatrix<double>& X_t, const Eigen::MatrixXd& Nabla_g,
    const double Lambda, const double Gamma, const int Nl, const int Nh);

void Find_update_entries_lh_lg_diag(
    const std::vector<std::pair<int, int>>& pairs_base,
    const int Nl, const int Nh, const int Ng,
    std::vector<std::pair<int, int>>& pairs_lg,
    std::vector<std::pair<int, int>>& pairs_diag,
    std::vector<std::pair<int, int>>& pairs_all);

void Find_update_entries_ll_diag(
    const std::vector<std::pair<int, int>>& pairs_base,
    const int Nl, std::vector<std::pair<int, int>>& pairs_diag,
    std::vector<std::pair<int, int>>& pairs_all);

static int Find_value_index_in_col(const Eigen::SparseMatrix<double>& A,
    const int row, const int col);

void D_index_pattern(Eigen::SparseMatrix<double>& D_t, const int n,
    const std::vector<std::pair<int, int>>& new_pairs, FixedDindex& out);

void Update_mu_ll_offdiag(
    const Eigen::SparseMatrix<double>& X_t,
    const int i, const int j, const std::size_t k,
    const double Lambda,
    FixedDindex& Didx,
    const RowMat& W_t,
    const RowMat& W_T,
    RowMat& U_T,
    const Eigen::MatrixXd& S);

void Update_mu_lh(
    const Eigen::SparseMatrix<double>& X_t,
    const int i, const int j, const std::size_t k,
    const double Lambda, const double Alpha, const double Gamma,
    FixedDindex& Didx,
    const RowMat& W_t,
    const RowMat& W_T,
    RowMat& U_T,
    const Eigen::MatrixXd& S);

void Update_mu_lg(
    const Eigen::SparseMatrix<double>& X_t,
    const int i, const int j, const std::size_t k,
    const double Alpha, FixedDindex& Didx,
    const RowMat& W_t,
    const RowMat& W_T,
    RowMat& U_T,
    const Eigen::MatrixXd& S);

void Update_mu_diag(
    const int i, const std::size_t k,
    FixedDindex& Didx,
    const RowMat& W_t,
    const RowMat& W_T,
    RowMat& U_T,
    const Eigen::MatrixXd& S);

double sgn(const double x);

std::size_t tri_index_upper(const int i, const int j, const int NL);

KKT_residual KKTResidual_lhg_log(const Eigen::SparseMatrix<double>& Theta,
    const Eigen::MatrixXd& S, const RowMat& W,
    const int Nl, const int Nh, const int Ng,
    const double lambda, const double alpha, const double gamma);

bool KKTResidual_lhg_convergence(const Eigen::SparseMatrix<double>& Theta,
    const Eigen::MatrixXd& S, const RowMat& W,
    const int Nl, const int Nh, const int Ng, const double Scale,
    const double lambda, const double alpha, const double gamma);

bool KKTResidual_l_convergence(const Eigen::SparseMatrix<double>& Theta,
    const Eigen::MatrixXd& S, const RowMat& W, const int Nl,
    const double Scale, const double lambda);

Eigen::MatrixXd Nabla_LL_only(const Eigen::MatrixXd& S, const RowMat& W_t, const int Nl);

Eigen::MatrixXd Nabla_lhg_LH(const Eigen::MatrixXd& S, const RowMat& W_t,
    const Eigen::SparseMatrix<double>& X_t, const int Nl, const int Nh,
    const double Alpha, const double Gamma);

void add_to_UT_col_from_Wrow(RowMat& U_T, const int col, const RowMat& W,
    const int rowW, const double mu);

Eigen::SparseMatrix<double> Expand_theta(const Eigen::SparseMatrix<double>& A, const int Nf);

#endif // _MATRIX_COMPUTE_H_
