#ifndef _LHG_UPDATE_H_
#define _LHG_UPDATE_H_

#include <Eigen/Dense>
#include <vector>
#include <functional>
#include <random>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>
#include "Class_define.h"
#include "Matrix_compute.h"
#include "Newton_method_offdiag.h"

using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using BlockReader = std::function<void(int col0, int bs, Eigen::Ref<RowMat> out)>;

RowMat Vstack(const RowMat& A, const RowMat& B);

Eigen::MatrixXd Cov_scdata(const RowMat& Data);

RowMat Topfactor_eigcov_pca(const BlockReader& read_block_raw,
    const int Rsel, const int C, const int nGlobals, const int BlockCols = 4096);

void Zscore_DataRM(RowMat& DataRM);

RowMat Topfactor_zscore_out(RowMat& Data, const int nHidden, const int nGlobal,
    const std::vector<int>& Idx_h, const std::vector<int>& Idx_g, const int BlockCols = 4096);

void LHG_zscore_out(RowMat& Data, Input_format& Input_l, const std::vector<int>& Idx_h,
    const std::vector<int>& Idx_g, const int BlockCols = 4096);

RowMat Update_hg_factor(const RowMat& Data, Input_format& Input_l);

bool SparseMatricesClose(const Eigen::SparseMatrix<double>& A,
    const Eigen::SparseMatrix<double>& B, const double D_theta);

RowMat Solution_theta_and_hg_factor(RowMat& Data, Input_format& Input_l,
    const std::vector<int>& Idx_h, const std::vector<int>& Idx_g, const int Max_loop = 200,
    const double Stop_threshold = 1e-4, const int Blockcol_inv = 1024, const int BlockCol_data = 4096);

double Likelihood_marginal_theta(Input_format& Input_l, const Eigen::MatrixXd& Cov_test);

std::pair<double, double> Likelihood_marginal_theta_and_global_magnitude(
    Input_format& Input_l, const Eigen::MatrixXd& Cov_test);

double Entropy_theta_lg(Input_format& Input_l);

double Global_dominance(const RowMat& Data_local, const RowMat& Factor_hg,
    Input_format& Input_l, const int BlockCol_data = 4096);

std::pair<double, double> Score_local_nd(const Eigen::SparseMatrix<double>& Theta_train,
    const Eigen::SparseMatrix<double>& Theta_test, const int Nl);

Eigen::MatrixXd Expandcov_lhg_reader(const Eigen::MatrixXd& S_LL, const BlockReader& Datareader,
    const int Nl, const int C, const RowMat& H, const int Blockcols = 4096);

KFoldPlan MakeKfoldplan_reorderdata(const RowMat& Data, const int Fold_k, unsigned int seed = 1);

SumStat Compute_sumsq_onerange(const RowMat& Data_perm, const int c0, const int c1);

SumStat Mergesumstat(const SumStat& a, const SumStat& b);

ZScoreParams Build_zscoreparams(const SumStat& st);

void Precompute_foldzscores(KFoldPlan& plan);

BlockReader Make_trainreader(const RowMat& Data_perm, const FoldRanges& fr, const ZScoreParams& zp);

BlockReader Make_testreader(const RowMat& Data_perm, const FoldRanges& fr, const ZScoreParams& zp);

Eigen::MatrixXd Cov_scdata_reader(const BlockReader& Reader, const int Nl,
    const int C, const int Blockcols = 4096);

BlockReader Make_trainreader_selectedrows(const RowMat& Data_perm,
    const FoldRanges& fr, const ZScoreParams& zp, const std::vector<int>& selected_rows);

BlockReader Make_testreader_selectedrows(const RowMat& Data_perm,
    const FoldRanges& fr, const ZScoreParams& zp, const std::vector<int>& selected_rows);

Eigen::MatrixXd HG_out_fold_train(const RowMat& Data_perm,
    const BlockReader& fullTrainReader, const FoldRanges& fr, const ZScoreParams& train_z,
    const std::vector<int>& Idx_h, const std::vector<int>& Idx_g,
    const Eigen::MatrixXd& Cov_train, const int Nh, const int Ng, const int Blockcols = 4096);

Eigen::MatrixXd HG_out_fold_test(const RowMat& Data_perm,
    const BlockReader& fullTrainReader, const FoldRanges& fr, const ZScoreParams& test_z,
    const std::vector<int>& Idx_h, const std::vector<int>& Idx_g,
    const Eigen::MatrixXd& Cov_test, const int Nh, const int Ng, const int Blockcols = 4096);

void LHG_zscore_out_fold_train(const RowMat& Data_perm, const FoldRanges& fr,
    const ZScoreParams& train_z, Input_format& Input_l, const Eigen::MatrixXd& S_LL,
    const std::vector<int>& Idx_h, const std::vector<int>& Idx_g, int Blockcols = 4096);

void LHG_zscore_out_fold_test(const RowMat& Data_perm, const FoldRanges& fr,
    const ZScoreParams& test_z, Input_format& Input_l, const Eigen::MatrixXd& S_LL,
    const std::vector<int>& Idx_h, const std::vector<int>& Idx_g, const int Blockcols = 4096);

RowMat Update_hg_factor_reader(const BlockReader& FullReader, const int C,
    Input_format& Input_l, const Eigen::MatrixXd& S_LL, const int Blockcols = 4096);

RowMat Solution_theta_and_hg_factor_fold(const BlockReader& FullReader, const int C,
    Input_format& Input_l, const Eigen::MatrixXd& S_LL, const int Max_loop = 200,
    const double Stop_threshold = 1e-4, const int Blockcol_inv = 1024, const int Blockcol_data = 4096);

#endif // _LHG_UPDATE_H_