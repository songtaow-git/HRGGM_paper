#ifndef _CLASS_DEFINE_H_
#define _CLASS_DEFINE_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <cstddef> 
#include <utility>

using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ColMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

class Input_format {
private:
	int Nlocal;
	int Nhidden;
	int Nglobal;
	double Lambda;
	double Alpha;
	double Gamma;
	Eigen::MatrixXd Observed_cov;
	Eigen::SparseMatrix<double> Int_theta;
	int Max_iteration;
	double Convergence_scale;
	bool If_save;

public:
	Input_format(const Eigen::MatrixXd& O_cov,
		const Eigen::SparseMatrix<double>& I_theta,
		const int N_l, const int N_h = 0, const int N_g = 0,
		const double Lam = 0.1, const double Alp = 0, const double Gam = 0, 
		const int M_iteration = 2000, 
		const double C_scale = 5e-5, 
		const bool I_save = false);
	Input_format(const Input_format& I_format);
	void Set_Observed_cov(const Eigen::MatrixXd& O_cov);
	void Set_Int_theta(const Eigen::SparseMatrix<double>& I_theta);
	void Set_Lambda(double Lam);
	void Set_Alpha(double Alp);
	void Print_input_format(const int P_idx);
	int Get_Nlocal();
	int Get_Nhidden();
	int Get_Nglobal();
	const Eigen::SparseMatrix<double>& Get_theta() const;
	const Eigen::MatrixXd& Get_observed_cov() const;
    
	friend void Newton_method_offdiag_l(Input_format& Input_l, const int Blockcol_inv);
	friend void Newton_method_offdiag_lhg(Input_format& Input_l, const int Blockcol_inv);
};

class MarkAndPush {
private:
	int B_;
	std::vector<uint8_t>& Mark_;
	std::vector<std::pair<int, int>>& Updates_;
public:
	MarkAndPush(const int B, std::vector<uint8_t>& Mark,
		std::vector<std::pair<int, int>>& Updates);
	void operator()(const int i, const int j);
};

class KKT_residual {
private:
	double KKT_r_diag;
	double KKT_r_ll_offdiag_nz;
	double KKT_r_ll_offdiag_z;
	double KKT_r_lh_nz;
	double KKT_r_lh_z;
	double KKT_r_lg;

public:
	double KKT_error;
	KKT_residual(const double Diag, const double LL_offdiag_nz, 
		const double LL_offdiag_z, const double LH_nz,
		const double LH_z, const double LG);
	KKT_residual(const KKT_residual& KKT_r);
	void KKT_print();
};

struct FixedDindex {
	Eigen::SparseMatrix<double>* D;             
	std::vector<std::pair<int, int>> pairs; 
	std::vector<int> idx_ij;                    
	std::vector<int> idx_ji;                   
	int n;

	FixedDindex() : D(0), n(0) {}
};

struct FoldSplit {
	std::vector<int> Train_idx; // 0-based indices
	std::vector<int> Test_idx;  // 0-based indices
};

struct FoldRanges
{
	int test_begin = 0;   // inclusive, in Data_perm
	int test_end = 0;   // exclusive, in Data_perm
	int n_used = 0;
	int fold_size = 0;
	int n_train() const { return n_used - fold_size; }
	int n_test()  const { return fold_size; }
};

struct ZScoreParams
{
	Eigen::VectorXd mu;        // size = nRows
	Eigen::VectorXd invsigma;  // size = nRows
};

struct KFoldPlan
{
	RowMat Data_perm;                    // Data(:, used_idx)
	std::vector<int> used_idx;           // original sample indices kept
	std::vector<int> dropped_idx;        // original sample indices dropped
	std::vector<FoldRanges> folds;       // one per fold
	std::vector<ZScoreParams> train_z;   // one per fold
	std::vector<ZScoreParams> test_z;    // one per fold
};

struct SumStat
{
	Eigen::VectorXd sum;
	Eigen::VectorXd sqsum;
	int count = 0;
};

struct Lambda_metric
{
	RowMat Share_mat;   // Fold_k x N_Alpha
	RowMat Difference_mat;
	RowMat Score_ratio;
	RowMat Score_norm;
	Lambda_metric(const int N_f, const int N_l);
	void Compute_score();
};

struct Alpha_metric
{
	RowMat Likelihood_mat;   // Fold_k x N_Alpha
	RowMat Entropy_mat;
	RowMat Score;
	Alpha_metric(const int N_f, const int N_a);
	void Compute_score();
};

#endif // _CLASS_DEFINE_H_

