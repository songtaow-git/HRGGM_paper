#include <iostream>
#include <string>
#include "Matrix_io.h"
#include "Class_define.h"
#include "Newton_method_offdiag.h"
#include "LHG_update.h"
#include "Parameter_select_parallel.h"
#include <filesystem>

int main(int argc, char* argv[]) {
	// 0) Original data input
	const int Max_iter = 1000;
	const double Converge_thre = 1e-3;
	const int N_thread = 96;
	const int Blockcols = 4096;

	// Assume executable is in: Project/build/C_GRN
	// Then project root is:     Project/
	std::filesystem::path exe_path = std::filesystem::absolute(argv[0]);
	std::filesystem::path exe_dir = exe_path.parent_path();
	std::filesystem::path base_dir = exe_dir.parent_path();

	std::filesystem::path data_dir = base_dir / "Data";
	std::filesystem::path out_dir = base_dir / "Results";
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

	RowMat Data = Load_bin_f64_data(Path_data);
	const int Fold_k = 3;
	const int Nlocal = Data.rows();
	const int Nhidden = 1;
	const int Nglobal = 1;
	const int Nfeature = Nlocal + Nhidden + Nglobal;
	const double Gamma = 0.8;

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

	// 1) Build fold plan and reorder data
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

	// 3) Find initial lambda
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

	filename_1 = "Score_norm_" + std::to_string(0) + ".txt";
	Write_metric(Score_orig.Score_norm, out_dir_1, filename_1);

	filename_1 = "Score_ratio_" + std::to_string(0) + ".txt";
	Write_metric(Score_orig.Score_ratio, out_dir_1, filename_1);

	// 4) Loop find optimal lambda and alpha
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
		Alp_metric.Compute_score();

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

		filename_2 = "Likelihood_mat_" + std::to_string(loop) + ".txt";
		Write_metric(Alp_metric.Likelihood_mat, out_dir_2, filename_2);

		filename_2 = "Entropy_mat_" + std::to_string(loop) + ".txt";
		Write_metric(Alp_metric.Entropy_mat, out_dir_2, filename_2);

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

		filename_1 = "Score_norm" + std::to_string(loop) + ".txt";
		Write_metric(Lam_metric.Score_norm, out_dir_1, filename_1);

		filename_1 = "Score_ratio" + std::to_string(loop) + ".txt";
		Write_metric(Lam_metric.Score_ratio, out_dir_1, filename_1);

		for (int j = 0; j < static_cast<int>(Alpha_list.size()); ++j) {
			Theta_preserve_a[j] = Theta_preserve_l[Max_idx];
		}

		if (Lambda_fix == Lambda_old) {
			I_theta = Theta_preserve_l[Max_idx];
			break;
		}

		// update optimal parameters
		Lambda_old = Lambda_fix;
		Alpha_old = Alpha_fix;
		++loop;
	}

	std::cout << System_time() << '\t' << "Optimal parameters:" << std::endl;
	std::cout << System_time() << '\t' << "Lambda:" << '\t' << Lambda_fix << std::endl;
	std::cout << System_time() << '\t' << "Alpha:" << '\t' << Alpha_fix << std::endl;
	std::cout << System_time() << '\t' << "Gamma:" << '\t' << Gamma << std::endl;

	Input_format Input_whole(I_cov, I_theta, Nlocal, Nhidden, Nglobal,
		Lambda_fix, Alpha_fix, Gamma, Max_iter, Converge_thre);

	Solution_theta_and_hg_factor(Data, Input_whole, valid_idx_h, valid_idx_g);

	std::filesystem::path out_file = out_dir / "Theta_final.txt";
	WriteSparseForMatlab(Input_whole.Get_theta(), out_file.string());

	return 0;
}
