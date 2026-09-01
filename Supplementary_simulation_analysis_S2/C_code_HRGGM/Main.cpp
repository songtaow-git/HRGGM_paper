#include <iostream>
#include <string>
#include "Matrix_io.h"
#include "Class_define.h"
#include "Newton_method_offdiag.h"
#include "LHG_update.h"
#include "Parameter_select_parallel.h"
#include <filesystem>

int main(int argc, char* argv[]) {
	// Adjustable parameters.
	const int Max_iter = 10000;              // Maximum iterations for each optimization call.
	const double Converge_thre = 1e-3;       // Convergence tolerance for precision-matrix optimization.
	const int N_thread = 8;                   // Number of parallel worker threads used in parameter selection.
	const int Blockcols = 4096;               // Number of sample columns processed per covariance-computation block.
	const int Fold_k = 3;                     // Number of cross-validation folds.
	const int Nhidden = 1;                    // Number of explicitly modeled hidden features.
	const int Nglobal = 1;                    // Number of explicitly modeled global features.
	const double Gamma = 0.9;                 // Penalty weight used for local-hidden regularization.
	const double Alpha_retention_tau = 0.20; // Retention threshold used in the alpha-selection score.

	Eigen::setNbThreads(0);

	// Locate the nearest parent directory containing the Data folder.
	std::filesystem::path exe_path = std::filesystem::absolute(argv[0]);
	std::filesystem::path exe_dir = exe_path.parent_path();
	std::filesystem::path base_dir = exe_dir;
	while (!std::filesystem::exists(base_dir / "Data")) {
		auto parent = base_dir.parent_path();
		if (parent == base_dir) {
			throw std::runtime_error(
				"Could not locate 'Data' directory above: " + exe_dir.string());
		}
		base_dir = parent;
	}

	// Input files are read from Data. Results are written to the sibling Result folder.
	std::filesystem::path data_dir = base_dir / "Data";
	std::filesystem::path out_dir = base_dir / "Result";
	std::filesystem::path out_dir_1 = out_dir / "lambda_ratio";
	std::filesystem::path out_dir_2 = out_dir / "Alpha_metric";

	std::filesystem::create_directories(out_dir);
	std::filesystem::create_directories(out_dir_1);
	std::filesystem::create_directories(out_dir_2);

	std::string Path_data = (data_dir / "Data_whole.bin").string();
	std::string Path_lambda = (data_dir / "Lambda_list.txt").string();
	std::string Path_alpha = (data_dir / "Alpha_list.txt").string();
	std::string Path_idx_g = (data_dir / "Idx_g_list.txt").string();
	std::string Path_idx_h = (data_dir / "Idx_h_list.txt").string();

	std::string filename_1;
	std::string filename_2;

	// Nlocal is obtained from the number of rows in Data_whole.bin.
	RowMat Data = Load_bin_f64_data(Path_data);
	const int Nlocal = Data.rows();
	const int Nfeature = Nlocal + Nhidden + Nglobal;

	std::vector<double> Lambda_list = ReadParameterFromTxt(Path_lambda);
	std::vector<Eigen::SparseMatrix<double>> Theta_preserve_l(Lambda_list.size());
	Eigen::SparseMatrix<double> I_theta(Nlocal, Nlocal);
	I_theta.setIdentity();
	for (int i = 0; i < static_cast<int>(Lambda_list.size()); ++i) {
		Theta_preserve_l[i] = I_theta;
	}

	std::vector<double> Alpha_list = ReadParameterFromTxt(Path_alpha);
	std::vector<Eigen::SparseMatrix<double>> Theta_preserve_a(Alpha_list.size());

	double Lambda_fix;
	double Alpha_fix;

	std::vector<int> valid_idx_h = ReadidxFromTxt(Path_idx_h);
	std::vector<int> valid_idx_g = ReadidxFromTxt(Path_idx_g);

	auto t1 = std::chrono::high_resolution_clock::now();

	// Construct cross-validation folds and compute fold-specific covariance inputs.
	std::cout << System_time() << '\t' << "Compute fold value" << std::endl;
	KFoldPlan plan = MakeKfoldplan_reorderdata(Data, Fold_k);
	Precompute_foldzscores(plan);

	std::vector<Eigen::MatrixXd> LHG_train(Fold_k);
	std::vector<Eigen::MatrixXd> LHG_test(Fold_k);
	std::vector<Eigen::MatrixXd> Cov_train(Fold_k);
	std::vector<Eigen::MatrixXd> Cov_test(Fold_k);

	Compute_fold_lhg_cov_parallel(
		plan, LHG_train, LHG_test, Cov_train, Cov_test,
		valid_idx_h, valid_idx_g, Nhidden, Nglobal, Blockcols
	);

	auto t2 = std::chrono::high_resolution_clock::now();
	double elapsed = std::chrono::duration<double>(t2 - t1).count();
	std::cout << System_time() << '\t'
		<< std::fixed << std::setprecision(3)
		<< "Compute fold value (Elapsed time, all tasks): "
		<< elapsed << " s\n\n";

	t1 = std::chrono::high_resolution_clock::now();

	// Select the initial lambda value using local-only cross-validation.
	std::cout << System_time() << '\t' << "Search initial lambda" << std::endl;

	Eigen::MatrixXd I_cov(Nlocal, Nlocal);
	Input_format Input_train_orig(I_cov, I_theta, Nlocal, Nhidden, Nglobal,
		1, 1, Gamma, Max_iter, Converge_thre);
	Input_format Input_test_orig(Input_train_orig);

	Lambda_metric Score_orig = Compute_orig_lambda_parallel(
		plan, Lambda_list, Input_train_orig, Cov_train,
		Input_test_orig, Cov_test, Theta_preserve_l, N_thread
	);
	Score_orig.Compute_score();

	Eigen::Index m_idx_0;
	double max_ratio_0 = Score_orig.Score_norm.row(0).maxCoeff(&m_idx_0);
	int max_idx_0 = static_cast<int>(m_idx_0);
	Lambda_fix = Lambda_list[max_idx_0];

	t2 = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration<double>(t2 - t1).count();
	std::cout << System_time() << '\t'
		<< std::fixed << std::setprecision(3)
		<< "Search initial lambda (Elapsed time, all tasks): "
		<< elapsed << " s\n\n";

	// Initial lambda outputs:
	// Result/lambda_ratio/Score_norm_0.txt stores normalized lambda-selection scores.
	// Result/lambda_ratio/Score_ratio_0.txt stores the corresponding score ratios.
	filename_1 = "Score_norm_" + std::to_string(0) + ".txt";
	Write_metric(Score_orig.Score_norm, out_dir_1, filename_1);

	filename_1 = "Score_ratio_" + std::to_string(0) + ".txt";
	Write_metric(Score_orig.Score_ratio, out_dir_1, filename_1);

	// Alternate alpha and lambda selection until the selected value is unchanged.
	int loop = 1;
	double Lambda_old = 0;
	double Alpha_old = 0;
	I_cov.resize(Nfeature, Nfeature);
	I_cov.setIdentity();
	I_theta = Expand_theta(Theta_preserve_l[max_idx_0], Nfeature);
	for (int j = 0; j < static_cast<int>(Alpha_list.size()); ++j) {
		Theta_preserve_a[j] = I_theta;
	}

	while (1) {
		std::cout << System_time() << '\t'
			<< "Search optimal alpha: "
			<< "loop " << loop << std::endl;

		t1 = std::chrono::high_resolution_clock::now();

		Input_format Input_train_a(I_cov, I_theta, Nlocal, Nhidden, Nglobal,
			Lambda_fix, 1, Gamma, Max_iter, Converge_thre);

		Alpha_metric Alp_metric = Compute_fold_alpha_parallel(
			plan, Alpha_list, Input_train_a, LHG_train, Cov_train, Cov_test,
			valid_idx_h, valid_idx_g, Theta_preserve_a, N_thread
		);
		Alp_metric.Compute_score(Alpha_retention_tau);

		Eigen::Index min_idx;
		double min_score = Alp_metric.Score.row(0).minCoeff(&min_idx);
		int Min_idx = static_cast<int>(min_idx);
		Alpha_fix = Alpha_list[Min_idx];

		t2 = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration<double>(t2 - t1).count();
		std::cout << System_time() << '\t'
			<< std::fixed << std::setprecision(3)
			<< "Search optimal alpha (Elapsed time, all tasks): "
			<< elapsed << " s\n\n";

		// Likelihood_mat_<loop>.txt stores Alp_metric.Likelihood_mat.
		filename_2 = "Likelihood_mat_" + std::to_string(loop) + ".txt";
		Write_metric(Alp_metric.Likelihood_mat, out_dir_2, filename_2);

		// Entropy_mat_<loop>.txt stores Alp_metric.Entropy_mat.
		filename_2 = "Entropy_mat_" + std::to_string(loop) + ".txt";
		Write_metric(Alp_metric.Entropy_mat, out_dir_2, filename_2);

		// Global_magnitude_mat_<loop>.txt stores Alp_metric.Global_magnitude_mat.
		filename_2 = "Global_magnitude_mat_" + std::to_string(loop) + ".txt";
		Write_metric(Alp_metric.Global_magnitude_mat, out_dir_2, filename_2);

		// Retention_vec_<loop>.txt stores Alp_metric.Retention_vec.
		filename_2 = "Retention_vec_" + std::to_string(loop) + ".txt";
		Write_metric(Alp_metric.Retention_vec, out_dir_2, filename_2);

		// Adjusted_entropy_vec_<loop>.txt stores Alp_metric.Adjusted_entropy_vec.
		filename_2 = "Adjusted_entropy_vec_" + std::to_string(loop) + ".txt";
		Write_metric(Alp_metric.Adjusted_entropy_vec, out_dir_2, filename_2);

		// Score_vec_<loop>.txt stores the alpha-selection score vector.
		filename_2 = "Score_vec_" + std::to_string(loop) + ".txt";
		Write_metric(Alp_metric.Score, out_dir_2, filename_2);

		if (Alpha_fix == Alpha_old) {
			I_theta = Theta_preserve_a[Min_idx];
			break;
		}

		if (loop == 1) {
			for (int i = 0; i < static_cast<int>(Lambda_list.size()); ++i) {
				Theta_preserve_l[i] = Theta_preserve_a[Min_idx];
			}
		}

		std::cout << System_time() << '\t'
			<< "Search optimal lambda: "
			<< "loop " << loop << std::endl;

		t1 = std::chrono::high_resolution_clock::now();

		Input_format Input_train_l(I_cov, I_theta, Nlocal, Nhidden, Nglobal,
			1, Alpha_fix, Gamma, Max_iter, Converge_thre);
		Input_format Input_test_l(Input_train_l);

		Lambda_metric Lam_metric = Compute_fold_lambda_parallel(
			plan, Lambda_list, Input_train_l, Input_test_l,
			Cov_train, Cov_test, LHG_train, LHG_test,
			valid_idx_h, valid_idx_g, Theta_preserve_l, N_thread
		);
		Lam_metric.Compute_score();

		Eigen::Index max_idx;
		double max_ratio = Lam_metric.Score_norm.row(0).maxCoeff(&max_idx);
		int Max_idx = static_cast<int>(max_idx);
		Lambda_fix = Lambda_list[Max_idx];

		t2 = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration<double>(t2 - t1).count();
		std::cout << System_time() << '\t'
			<< std::fixed << std::setprecision(3)
			<< "Search optimal lambda (Elapsed time, all tasks): "
			<< elapsed << " s\n\n";

		// Score_norm<loop>.txt stores Lam_metric.Score_norm.
		filename_1 = "Score_norm" + std::to_string(loop) + ".txt";
		Write_metric(Lam_metric.Score_norm, out_dir_1, filename_1);

		// Score_ratio<loop>.txt stores Lam_metric.Score_ratio.
		filename_1 = "Score_ratio" + std::to_string(loop) + ".txt";
		Write_metric(Lam_metric.Score_ratio, out_dir_1, filename_1);

		for (int j = 0; j < static_cast<int>(Alpha_list.size()); ++j) {
			Theta_preserve_a[j] = Theta_preserve_l[Max_idx];
		}

		if (Lambda_fix == Lambda_old) {
			I_theta = Theta_preserve_l[Max_idx];
			break;
		}

		Lambda_old = Lambda_fix;
		Alpha_old = Alpha_fix;
		++loop;
		if (loop > 6) {
			std::cout << System_time() << '\t'
				<< "Warning: Maximum number of loops reached. "
				<< "The optimal parameters may not have converged." << std::endl;
			break;
		}
	}

	std::cout << System_time() << '\t' << "Optimal parameters:" << std::endl;
	std::cout << System_time() << '\t' << "Lambda:" << '\t' << Lambda_fix << std::endl;
	std::cout << System_time() << '\t' << "Alpha:" << '\t' << Alpha_fix << std::endl;
	std::cout << System_time() << '\t' << "Gamma:" << '\t' << Gamma << std::endl;

	// Fit the model to the complete dataset using the selected parameters.
	Input_format Input_whole(I_cov, I_theta, Nlocal, Nhidden, Nglobal,
		Lambda_fix, Alpha_fix, Gamma, Max_iter, Converge_thre);

	Solution_theta_and_hg_factor(Data, Input_whole, valid_idx_h, valid_idx_g);

	// Result/Theta_final.txt stores the final estimated sparse precision matrix.
	std::filesystem::path out_file = out_dir / "Theta_final.txt";
	WriteSparseForMatlab(Input_whole.Get_theta(), out_file.string());

	return 0;
}
