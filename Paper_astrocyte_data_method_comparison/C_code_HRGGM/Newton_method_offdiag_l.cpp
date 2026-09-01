#include "Newton_method_offdiag.h"

void Newton_method_offdiag_l(Input_format& Input_l, const int Blockcol_inv) 
{
	int Size = Input_l.Nlocal;
	double Likelihood_0 = Likelihood_f(Input_l.Observed_cov, Input_l.Int_theta);
	double Norm_l1_ll_0 = Input_l.Lambda * Offdiag_l1_norm_ll(Input_l.Int_theta, Input_l.Nlocal);
	double Objective_0 = Likelihood_0 + Norm_l1_ll_0;
	std::vector<double> Likelihood;
	std::vector<double> Norm_l1_ll;
	std::vector<double> Objective;
	std::vector<double> D_objective;
	Likelihood.push_back(Likelihood_0);
	Norm_l1_ll.push_back(Norm_l1_ll_0);
	Objective.push_back(Objective_0);
	D_objective.push_back(0);
	Eigen::SparseMatrix<double> Theta_t = Input_l.Int_theta;
	Eigen::SparseMatrix<double> D_t(Size, Size);
	RowMat U_T = RowMat::Zero(Size, Size);
	RowMat W_t = Inverse_spd_sparse(Theta_t, Blockcol_inv);
	RowMat W_T = W_t.transpose().eval();
	Eigen::MatrixXd Nabla_ll = Eigen::MatrixXd(Input_l.Nlocal, Input_l.Nlocal);
	std::vector<std::pair<int, int>> Update_vec_ll;
	std::vector<std::pair<int, int>> Update_vec_diag;
	std::vector<std::pair<int, int>> Update_vec_ll_diag;
	bool If_converge(false);

	for (int i = 0; i < Input_l.Max_iteration; ++i) {
		// Update local-local(offdiag)
		Nabla_ll = Nabla_LL_only(Input_l.Observed_cov, W_t, Input_l.Nlocal);
		Update_vec_ll = Find_update_entries_ll(Theta_t, Nabla_ll, Input_l.Lambda, Input_l.Nlocal);
		Find_update_entries_ll_diag(Update_vec_ll, Input_l.Nlocal, Update_vec_diag, Update_vec_ll_diag);
		FixedDindex Didx_ll_diag;
		D_index_pattern(D_t, Size, Update_vec_ll_diag, Didx_ll_diag);
		for (size_t k = 0; k < Update_vec_ll.size(); ++k) {
			int C_x = Update_vec_ll[k].first;
			int C_y = Update_vec_ll[k].second;
			Update_mu_ll_offdiag(Theta_t, C_x, C_y, k, Input_l.Lambda, Didx_ll_diag,
				W_t, W_T, U_T, Input_l.Observed_cov);
		}
		// Update diag
		for (size_t k = 0; k < Update_vec_diag.size(); ++k) {
			int C_x = Update_vec_diag[k].first;
			Update_mu_diag(C_x, (k + Update_vec_ll.size()),Didx_ll_diag, 
				W_t, W_T, U_T, Input_l.Observed_cov);
		}
		int h = 0;
		double Likelihood_h = 0;
		double Norm_l1_ll_h = 0;
		double Objective_h = 0;
		double Delta_f = 0;
		Eigen::SparseMatrix<double> Theta_old = Theta_t;
		double Beta_h = 1.0;
		Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol_l;
		chol_l.analyzePattern(Theta_old + D_t);
		while (true) {
			Theta_t = Theta_old + Beta_h * D_t;
			chol_l.factorize(Theta_t);
			if (chol_l.info() == Eigen::Success) {
				Likelihood_h = Likelihood_f(Input_l.Observed_cov, Theta_t, chol_l);
				Norm_l1_ll_h = Input_l.Lambda * Offdiag_l1_norm_ll(Theta_t, Input_l.Nlocal);
				Delta_f = (Likelihood_h - Likelihood[i]) +(Norm_l1_ll_h - Norm_l1_ll[i]);
				if (Delta_f <= 0) {
					Objective_h = Likelihood_h + Norm_l1_ll_h;
					W_t = Inverse_spd_sparse_from_chol(chol_l, Size, Blockcol_inv);
					W_T = W_t.transpose().eval();
					U_T.setZero();
					Update_vec_ll.clear();
					Update_vec_diag.clear();
					Update_vec_ll_diag.clear();
					Likelihood.push_back(Likelihood_h);
					Norm_l1_ll.push_back(Norm_l1_ll_h);
					Objective.push_back(Objective_h);
					D_objective.push_back(Delta_f);
					
					const int sd = Delta_f;
					break;
				}
			}
			Beta_h *= 0.5;
			if (Beta_h < 1e-6) {
				break;
			}
		}
		If_converge = KKTResidual_l_convergence(Theta_t, Input_l.Observed_cov,
			W_t, Input_l.Nlocal, Input_l.Convergence_scale, Input_l.Lambda);
		if (If_converge) {
			Input_l.Set_Int_theta(Theta_t);
			break;
		}
	}
	Input_l.Set_Int_theta(Theta_t);
}