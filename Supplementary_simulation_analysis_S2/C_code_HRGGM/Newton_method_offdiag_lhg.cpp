#include "Newton_method_offdiag.h"

void Newton_method_offdiag_lhg(Input_format& Input_l, const int Blockcol_inv)
{
	int Size = Input_l.Nlocal + Input_l.Nhidden + Input_l.Nglobal;
	double Likelihood_0 = Likelihood_f(Input_l.Observed_cov, Input_l.Int_theta);
	double Norm_l1_ll_0 = Input_l.Lambda * Offdiag_l1_norm_ll(Input_l.Int_theta,
		Input_l.Nlocal);
	double Norm_l1_lh_0 = Input_l.Gamma * Input_l.Lambda * Block_l1_norm_lh(Input_l.Int_theta,
		Input_l.Nlocal, Input_l.Nhidden);
	double Norm_l2_lh_0 = (1 - Input_l.Gamma) * Input_l.Alpha * Block_l2_norm_lh(Input_l.Int_theta,
		Input_l.Nlocal, Input_l.Nhidden);
	double Norm_l2_lg_0 = Input_l.Alpha * Block_l2_norm_lg(Input_l.Int_theta,
		Input_l.Nlocal, Input_l.Nhidden, Input_l.Nglobal);
	double Objective_0 = Likelihood_0 + Norm_l1_ll_0 + Norm_l1_lh_0 + Norm_l2_lh_0 + Norm_l2_lg_0;

	std::vector<double> Likelihood;
	std::vector<double> Norm_l1_ll;
	std::vector<double> Norm_l1_lh;
	std::vector<double> Norm_l2_lh;
	std::vector<double> Norm_l2_lg;
	std::vector<double> Objective;
	std::vector<double> D_objective;

	Likelihood.push_back(Likelihood_0);
	Norm_l1_ll.push_back(Norm_l1_ll_0);
	Norm_l1_lh.push_back(Norm_l1_lh_0);
	Norm_l2_lh.push_back(Norm_l2_lh_0);
	Norm_l2_lg.push_back(Norm_l2_lg_0);
	Objective.push_back(Objective_0);
	D_objective.push_back(0);

	Eigen::SparseMatrix<double> Theta_t = Input_l.Int_theta;
	Eigen::SparseMatrix<double> D_t(Size, Size);
	RowMat U_T = RowMat::Zero(Size, Size);
	RowMat W_t = Inverse_spd_sparse(Theta_t, Blockcol_inv);
	RowMat W_T = W_t.transpose();

	Eigen::MatrixXd Nabla_ll = Eigen::MatrixXd(Input_l.Nlocal, Input_l.Nlocal);
	Eigen::MatrixXd Nabla_lh = Eigen::MatrixXd(Input_l.Nlocal, Input_l.Nhidden);

	std::vector<std::pair<int, int>> Update_vec_ll;
	std::vector<std::pair<int, int>> Update_vec_lh;
	std::vector<std::pair<int, int>> Update_vec_lg;
	std::vector<std::pair<int, int>> Update_vec_diag;
	std::vector<std::pair<int, int>> Update_vec_lh_lg_diag;

	bool If_converge(false);

	for (int i = 0; i < Input_l.Max_iteration; ++i) {

		// ============================================================
		// Minimal fix:
		// use the last successfully saved objective-history index.
		// ============================================================
		std::size_t last = Likelihood.size() - 1;

		// Update local-local(offdiag)
		Nabla_ll = Nabla_LL_only(Input_l.Observed_cov, W_t, Input_l.Nlocal);

		Update_vec_ll = Find_update_entries_ll(
			Theta_t,
			Nabla_ll,
			Input_l.Lambda,
			Input_l.Nlocal
		);

		FixedDindex Didx_ll;
		D_index_pattern(D_t, Size, Update_vec_ll, Didx_ll);

		for (size_t k = 0; k < Update_vec_ll.size(); ++k) {
			int C_x = Update_vec_ll[k].first;
			int C_y = Update_vec_ll[k].second;

			Update_mu_ll_offdiag(
				Theta_t,
				C_x,
				C_y,
				k,
				Input_l.Lambda,
				Didx_ll,
				W_t,
				W_T,
				U_T,
				Input_l.Observed_cov
			);
		}

		// Update X_t/W_t
		int h = 0;
		double Likelihood_h1 = 0;
		double Norm_l1_ll_h1 = 0;
		double Norm_l1_lh_h1 = 0;
		double Norm_l2_lh_h1 = 0;
		double Norm_l2_lg_h1 = 0;
		double Delta_f1 = 0;

		Eigen::SparseMatrix<double> Theta_old = Theta_t;
		Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol_ll;
		chol_ll.analyzePattern(Theta_old + D_t);

		double Beta_h1 = 1.0;

		while (true) {
			Theta_t = Theta_old + Beta_h1 * D_t;
			chol_ll.factorize(Theta_t);

			if (chol_ll.info() == Eigen::Success) {
				Likelihood_h1 = Likelihood_f(Input_l.Observed_cov, Theta_t, chol_ll);

				Norm_l1_ll_h1 = Input_l.Lambda * Offdiag_l1_norm_ll(
					Theta_t,
					Input_l.Nlocal
				);

				Norm_l1_lh_h1 = Input_l.Gamma * Input_l.Lambda * Block_l1_norm_lh(
					Theta_t,
					Input_l.Nlocal,
					Input_l.Nhidden
				);

				Norm_l2_lh_h1 = (1 - Input_l.Gamma) * Input_l.Alpha * Block_l2_norm_lh(
					Theta_t,
					Input_l.Nlocal,
					Input_l.Nhidden
				);

				Norm_l2_lg_h1 = Input_l.Alpha * Block_l2_norm_lg(
					Theta_t,
					Input_l.Nlocal,
					Input_l.Nhidden,
					Input_l.Nglobal
				);

				Delta_f1 = (Likelihood_h1 - Likelihood[last]) + (Norm_l1_ll_h1 - Norm_l1_ll[last]) +
					(Norm_l1_lh_h1 - Norm_l1_lh[last]) + (Norm_l2_lh_h1 - Norm_l2_lh[last]) +
					(Norm_l2_lg_h1 - Norm_l2_lg[last]);

				if (Delta_f1 <= 0) {
					W_t = Inverse_spd_sparse_from_chol(chol_ll, Size, Blockcol_inv);
					W_T = W_t.transpose();
					U_T.setZero();
					break;
				}
			}

			Beta_h1 *= 0.5;

			if (Beta_h1 < 1e-6) {
				break;
			}
		}

		// Update local-hidden/hidden-local
		Nabla_lh = Nabla_lhg_LH(
			Input_l.Observed_cov,
			W_t,
			Theta_t,
			Input_l.Nlocal,
			Input_l.Nhidden,
			Input_l.Alpha,
			Input_l.Gamma
		);

		Update_vec_lh = Find_update_entries_lh(
			Theta_t,
			Nabla_lh,
			Input_l.Lambda,
			Input_l.Gamma,
			Input_l.Nlocal,
			Input_l.Nhidden
		);

		Find_update_entries_lh_lg_diag(
			Update_vec_lh,
			Input_l.Nlocal,
			Input_l.Nhidden,
			Input_l.Nglobal,
			Update_vec_lg,
			Update_vec_diag,
			Update_vec_lh_lg_diag
		);

		FixedDindex Didx_lh_lg_diag;
		D_index_pattern(D_t, Size, Update_vec_lh_lg_diag, Didx_lh_lg_diag);

		for (std::size_t k = 0; k < Update_vec_lh.size(); ++k) {
			const int C_x = Update_vec_lh[k].first;
			const int C_y = Update_vec_lh[k].second;

			Update_mu_lh(
				Theta_t,
				C_x,
				C_y,
				k,
				Input_l.Lambda,
				Input_l.Alpha,
				Input_l.Gamma,
				Didx_lh_lg_diag,
				W_t,
				W_T,
				U_T,
				Input_l.Observed_cov
			);
		}

		// Update local-global/global-local
		for (size_t k = 0; k < Update_vec_lg.size(); ++k) {
			int C_x = Update_vec_lg[k].first;
			int C_y = Update_vec_lg[k].second;

			Update_mu_lg(
				Theta_t,
				C_x,
				C_y,
				(k + Update_vec_lh.size()),
				Input_l.Alpha,
				Didx_lh_lg_diag,
				W_t,
				W_T,
				U_T,
				Input_l.Observed_cov
			);
		}

		// Update diag
		for (size_t k = 0; k < Update_vec_diag.size(); ++k) {
			int C_x = Update_vec_diag[k].first;

			Update_mu_diag(
				C_x,
				(k + Update_vec_lh.size() + Update_vec_lg.size()),
				Didx_lh_lg_diag,
				W_t,
				W_T,
				U_T,
				Input_l.Observed_cov
			);
		}

		// Update X_t/W_t
		h = 0;
		double Likelihood_h2 = 0;
		double Norm_l1_ll_h2 = 0;
		double Norm_l1_lh_h2 = 0;
		double Norm_l2_lh_h2 = 0;
		double Norm_l2_lg_h2 = 0;
		double Delta_f2 = 0;
		double Objective_h = 0;
		double Delta_f = 0;

		Theta_old = Theta_t;
		double Beta_h2 = 1.0;

		Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol_lhg;
		chol_lhg.analyzePattern(Theta_old + D_t);

		while (true) {
			Theta_t = Theta_old + Beta_h2 * D_t;
			chol_lhg.factorize(Theta_t);

			if (chol_lhg.info() == Eigen::Success) {
				Likelihood_h2 = Likelihood_f(Input_l.Observed_cov, Theta_t, chol_lhg);

				Norm_l1_ll_h2 = Input_l.Lambda * Offdiag_l1_norm_ll(
					Theta_t,
					Input_l.Nlocal
				);

				Norm_l1_lh_h2 = Input_l.Gamma * Input_l.Lambda * Block_l1_norm_lh(
					Theta_t,
					Input_l.Nlocal,
					Input_l.Nhidden
				);

				Norm_l2_lh_h2 = (1 - Input_l.Gamma) * Input_l.Alpha * Block_l2_norm_lh(
					Theta_t,
					Input_l.Nlocal,
					Input_l.Nhidden
				);

				Norm_l2_lg_h2 = Input_l.Alpha * Block_l2_norm_lg(
					Theta_t,
					Input_l.Nlocal,
					Input_l.Nhidden,
					Input_l.Nglobal
				);

				Delta_f2 = (Likelihood_h2 - Likelihood_h1) + (Norm_l1_ll_h2 - Norm_l1_ll_h1) +
					(Norm_l1_lh_h2 - Norm_l1_lh_h1) + (Norm_l2_lh_h2 - Norm_l2_lh_h1) +
					(Norm_l2_lg_h2 - Norm_l2_lg_h1);

				if (Delta_f2 <= 0) {
					Delta_f = (Likelihood_h2 - Likelihood[last]) + (Norm_l1_ll_h2 - Norm_l1_ll[last]) +
						(Norm_l1_lh_h2 - Norm_l1_lh[last]) + (Norm_l2_lh_h2 - Norm_l2_lh[last]) +
						(Norm_l2_lg_h2 - Norm_l2_lg[last]);

					Objective_h = Likelihood_h2 + Norm_l1_ll_h2 + Norm_l1_lh_h2 +
						Norm_l2_lh_h2 + Norm_l2_lg_h2;

					W_t = Inverse_spd_sparse_from_chol(chol_lhg, Size, Blockcol_inv);
					W_T = W_t.transpose();
					U_T.setZero();

					Update_vec_ll.clear();
					Update_vec_lh.clear();
					Update_vec_lg.clear();
					Update_vec_diag.clear();
					Update_vec_lh_lg_diag.clear();

					Likelihood.push_back(Likelihood_h2);
					Norm_l1_ll.push_back(Norm_l1_ll_h2);
					Norm_l1_lh.push_back(Norm_l1_lh_h2);
					Norm_l2_lh.push_back(Norm_l2_lh_h2);
					Norm_l2_lg.push_back(Norm_l2_lg_h2);
					Objective.push_back(Objective_h);
					D_objective.push_back(Delta_f);

					break;
				}
			}

			Beta_h2 *= 0.5;

			if (Beta_h2 < 1e-6) {
				break;
			}
		}

		If_converge = KKTResidual_lhg_convergence(
			Theta_t,
			Input_l.Observed_cov,
			W_t,
			Input_l.Nlocal,
			Input_l.Nhidden,
			Input_l.Nglobal,
			Input_l.Convergence_scale,
			Input_l.Lambda,
			Input_l.Alpha,
			Input_l.Gamma
		);

		if (If_converge) {
			Input_l.Set_Int_theta(Theta_t);
			break;
		}
	}

	Input_l.Set_Int_theta(Theta_t);
}