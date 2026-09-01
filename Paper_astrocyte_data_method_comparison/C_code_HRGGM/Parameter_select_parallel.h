#ifndef _PARAMETER_SELECT_PARALLEL_H_
#define _PARAMETER_SELECT_PARALLEL_H_

#include <Eigen/Dense>
#include <vector>
#include <stdexcept>
#include <omp.h>
#include <chrono>
#include "Matrix_io.h"
#include "Class_define.h"
#include "LHG_update.h"
#include "Newton_method_offdiag.h"

using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

Input_format Copy_input_lambda(const Input_format& Base_input, const double Lam);

Input_format Copy_input_alpha(const Input_format& Base_input, const double Alp);

void Compute_fold_lhg_cov_parallel(const KFoldPlan& plan, std::vector<Eigen::MatrixXd>& LHG_train,
    std::vector<Eigen::MatrixXd>& LHG_test, std::vector<Eigen::MatrixXd>& Cov_train,
    std::vector<Eigen::MatrixXd>& Cov_test, const std::vector<int>& valid_idx_h,
    const std::vector<int>& valid_idx_g, const int Nh, const int Ng, const int Blockcols = 4096);

Lambda_metric Compute_fold_lambda_parallel(const KFoldPlan& plan, const std::vector<double>& lambda_list,
    const Input_format& Input_train_base, const Input_format& Input_test_base,
    const std::vector<Eigen::MatrixXd>& Cov_train, const std::vector<Eigen::MatrixXd>& Cov_test,
    const std::vector<Eigen::MatrixXd>& LHG_train, const std::vector<Eigen::MatrixXd>& LHG_test,
    const std::vector<int>& valid_idx_h, const std::vector<int>& valid_idx_g,
    std::vector<Eigen::SparseMatrix<double>>& Theta_pre, const int N_thread = 1);

Alpha_metric Compute_fold_alpha_parallel(const KFoldPlan& plan, const std::vector<double>& alpha_list,
    const Input_format& Input_train_base, const std::vector<Eigen::MatrixXd>& LHG_train,
    const std::vector<Eigen::MatrixXd>& Cov_train, const std::vector<Eigen::MatrixXd>& Cov_test,
    const std::vector<int>& valid_idx_h, const std::vector<int>& valid_idx_g,
    std::vector<Eigen::SparseMatrix<double>>& Theta_pre, const int N_thread = 1);

Lambda_metric Compute_orig_lambda_parallel(const KFoldPlan& plan, const std::vector<double>& lambda_list,
    const Input_format& Input_train_base, const std::vector<Eigen::MatrixXd>& Cov_train,
    const Input_format& Input_test_base, const std::vector<Eigen::MatrixXd>& Cov_test,
    std::vector<Eigen::SparseMatrix<double>>& Theta_pre, const int N_thread = 1);

#endif // _PARAMETER_SELECT_PARALLEL_H_