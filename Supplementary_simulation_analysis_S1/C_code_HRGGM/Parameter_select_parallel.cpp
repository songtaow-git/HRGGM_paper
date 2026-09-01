#include "Parameter_select_parallel.h"
#include <filesystem>

Input_format Copy_input_lambda(const Input_format& Base_input, const double Lam)
{
    Input_format out = Base_input;
    out.Set_Lambda(Lam);
    return out;
}

Input_format Copy_input_alpha(const Input_format& Base_input, const double Alp)
{
    Input_format out = Base_input;
    out.Set_Alpha(Alp);
    return out;
}

void Compute_fold_lhg_cov_parallel(const KFoldPlan& plan, std::vector<Eigen::MatrixXd>& LHG_train,
    std::vector<Eigen::MatrixXd>& LHG_test, std::vector<Eigen::MatrixXd>& Cov_train,
    std::vector<Eigen::MatrixXd>& Cov_test, const std::vector<int>& valid_idx_h, 
    const std::vector<int>& valid_idx_g, const int Nh, const int Ng, const int Blockcols)
{
    const int Fold_k = static_cast<int>(plan.folds.size());
    if (Fold_k <= 0) {
        throw std::invalid_argument("Fold_k must be > 0.");
    }
    if (static_cast<int>(plan.train_z.size()) != Fold_k ||
        static_cast<int>(plan.test_z.size()) != Fold_k) {
        throw std::invalid_argument("plan.train_z / plan.test_z size mismatch with plan.folds.");
    }
    for (int k = 0; k < Fold_k; ++k) {
        std::cout << System_time() << '\t' << "Fold " << k + 1 << " start calculation."<< "\n";
        auto t1 = std::chrono::high_resolution_clock::now();
        const FoldRanges& fr = plan.folds[k];
        const ZScoreParams& train_z = plan.train_z[k];
        const ZScoreParams& test_z = plan.test_z[k];
        // build readers
        BlockReader Train_reader = Make_trainreader(plan.Data_perm, fr, train_z);
        BlockReader Test_reader = Make_testreader(plan.Data_perm, fr, test_z);
        const int C_train = fr.n_train();
        const int C_test = fr.n_test();
        // compute Cov_obs
        Eigen::MatrixXd S_obs_train = Cov_scdata_reader(Train_reader,
            static_cast<int>(plan.Data_perm.rows()), fr.n_train());
        Eigen::MatrixXd S_obs_test = Cov_scdata_reader(Test_reader,
            static_cast<int>(plan.Data_perm.rows()), fr.n_test());
        Cov_train[k] = S_obs_train;
        Cov_test[k] = S_obs_test;
        // hidden and global and initial out
        LHG_train[k] = HG_out_fold_train(plan.Data_perm, Train_reader, fr, train_z,
            valid_idx_h, valid_idx_g, S_obs_train, Nh, Ng, Blockcols);
        LHG_test[k] = HG_out_fold_test(plan.Data_perm, Test_reader, fr, test_z,
            valid_idx_h, valid_idx_g, S_obs_test, Nh, Ng, Blockcols);
        auto t2 = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(t2 - t1).count();
        std::cout << System_time() << '\t' << std::fixed << std::setprecision(3) <<
            "Fold " << k + 1 << " : " << elapsed << " s\n";
    }
}

Lambda_metric Compute_fold_lambda_parallel(const KFoldPlan& plan, const std::vector<double>& lambda_list,
    const Input_format& Input_train_base, const Input_format& Input_test_base,
    const std::vector<Eigen::MatrixXd>& Cov_train, const std::vector<Eigen::MatrixXd>& Cov_test,
    const std::vector<Eigen::MatrixXd>& LHG_train, const std::vector<Eigen::MatrixXd>& LHG_test,
    const std::vector<int>& valid_idx_h, const std::vector<int>& valid_idx_g, 
    std::vector<Eigen::SparseMatrix<double>>& Theta_pre, const int N_thread)
{
    const int Fold_k = static_cast<int>(plan.folds.size());
    const int N_Lambda = static_cast<int>(lambda_list.size());
    if (Fold_k <= 0) {
        throw std::invalid_argument("Fold_k must be > 0.");
    }
    if (N_Lambda <= 0) {
        throw std::invalid_argument("N_Lambda must be > 0.");
    }
    if (N_thread <= 0) {
        throw std::invalid_argument("N_thread must be > 0.");
    }
    if (static_cast<int>(plan.train_z.size()) != Fold_k ||
        static_cast<int>(plan.test_z.size()) != Fold_k) {
        throw std::invalid_argument("plan.train_z / plan.test_z size mismatch with plan.folds.");
    }
    Lambda_metric Score(Fold_k, N_Lambda);
    omp_set_num_threads(N_thread);
    // ------------------------------------------------------------
    // Stage 1:
    // Compute fold 0 for all lambdas in parallel.
    // Only fold 0 is allowed to update Theta_pre[lam_id].
    // ------------------------------------------------------------
#pragma omp parallel for schedule(dynamic)
    for (int lam_id = 0; lam_id < N_Lambda; ++lam_id) {
        int tid = omp_get_thread_num();
        const int k = 0;
        const double lambda_val = lambda_list[lam_id];
#pragma omp critical
        {
            std::cout << System_time() << '\t'
                << "Stage 1 | lambda: " << lam_id + 1 << '\t'
                << "fold: " << k + 1
                << " on thread " << tid + 1 << "\n";
        }
        auto t1 = std::chrono::high_resolution_clock::now();

        const FoldRanges& fr = plan.folds[k];
        const ZScoreParams& train_z = plan.train_z[k];
        const ZScoreParams& test_z = plan.test_z[k];
        Input_format Input_train = Copy_input_lambda(Input_train_base, lambda_val);
        Input_train.Set_Int_theta(Theta_pre[lam_id]);
        BlockReader Train_reader = Make_trainreader(plan.Data_perm, fr, train_z);
        const int C_train = fr.n_train();
        Input_train.Set_Observed_cov(LHG_train[k]);
        Solution_theta_and_hg_factor_fold(Train_reader, C_train, Input_train, Cov_train[k]);
        Theta_pre[lam_id] = Input_train.Get_theta();

        Input_format Input_test = Copy_input_lambda(Input_test_base, lambda_val);
        Input_test.Set_Int_theta(Theta_pre[lam_id]);
        BlockReader Test_reader = Make_testreader(plan.Data_perm, fr, test_z);
        const int C_test = fr.n_test();
        Input_test.Set_Observed_cov(LHG_test[k]);
        Solution_theta_and_hg_factor_fold(Test_reader, C_test, Input_test, Cov_test[k]);
        // metrics
        const std::pair<double, double> Pair_l = Score_local_nd(Input_train.Get_theta(),
            Input_test.Get_theta(), Input_train.Get_Nlocal());
        Score.Share_mat(k, lam_id) = Pair_l.first;
        Score.Difference_mat(k, lam_id) = Pair_l.second;
        auto t2 = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(t2 - t1).count();
#pragma omp critical
        {
            std::cout << System_time() << '\t'
                << std::fixed << std::setprecision(3)
                << "Stage 1 | lambda: " << lam_id + 1 << '\t'
                << "fold " << k + 1
                << " : " << elapsed << " s\n";
        }
    }
    // ------------------------------------------------------------
    // Stage 2:
    // After all fold-0 tasks are completed, compute folds 1..K-1.
    // ------------------------------------------------------------
    if (Fold_k > 1) {
        const int total_jobs = N_Lambda * (Fold_k - 1);
#pragma omp parallel for schedule(dynamic)
        for (int job = 0; job < total_jobs; ++job) {
            const int lam_id = job / (Fold_k - 1);
            const int k = 1 + (job % (Fold_k - 1));
            const int tid = omp_get_thread_num();
            const double lambda_val = lambda_list[lam_id];
#pragma omp critical
            {
                std::cout << System_time() << '\t'
                    << "Stage 2 | lambda: " << lam_id + 1 << '\t'
                    << "fold: " << k + 1
                    << " on thread " << tid + 1 << "\n";
            }
            auto t1 = std::chrono::high_resolution_clock::now();

            const FoldRanges& fr = plan.folds[k];
            const ZScoreParams& train_z = plan.train_z[k];
            const ZScoreParams& test_z = plan.test_z[k];
            Input_format Input_train = Copy_input_lambda(Input_train_base, lambda_val);
            Input_train.Set_Int_theta(Theta_pre[lam_id]);
            BlockReader Train_reader = Make_trainreader(plan.Data_perm, fr, train_z);
            const int C_train = fr.n_train();
            Input_train.Set_Observed_cov(LHG_train[k]);
            Solution_theta_and_hg_factor_fold(Train_reader, C_train, Input_train, Cov_train[k]);

            Input_format Input_test = Copy_input_lambda(Input_test_base, lambda_val);
            Input_test.Set_Int_theta(Theta_pre[lam_id]);
            BlockReader Test_reader = Make_testreader(plan.Data_perm, fr, test_z);
            const int C_test = fr.n_test();
            Input_test.Set_Observed_cov(LHG_test[k]);
            Solution_theta_and_hg_factor_fold(Test_reader, C_test, Input_test, Cov_test[k]);
            // metrics
            const std::pair<double, double> Pair_l = Score_local_nd(Input_train.Get_theta(),
                Input_test.Get_theta(), Input_train.Get_Nlocal());
            Score.Share_mat(k, lam_id) = Pair_l.first;
            Score.Difference_mat(k, lam_id) = Pair_l.second;
            auto t2 = std::chrono::high_resolution_clock::now();
            double elapsed = std::chrono::duration<double>(t2 - t1).count();

#pragma omp critical
            {
                std::cout << System_time() << '\t'
                    << std::fixed << std::setprecision(3)
                    << "Stage 2 | lambda: " << lam_id + 1 << '\t'
                    << "fold " << k + 1
                    << " : " << elapsed << " s\n";
            }
        }
    }
    return Score;
}

Alpha_metric Compute_fold_alpha_parallel(const KFoldPlan& plan, const std::vector<double>& alpha_list,
    const Input_format& Input_train_base, const std::vector<Eigen::MatrixXd>& LHG_train,
    const std::vector<Eigen::MatrixXd>& Cov_train, const std::vector<Eigen::MatrixXd>& Cov_test,
    const std::vector<int>& valid_idx_h, const std::vector<int>& valid_idx_g, 
    std::vector<Eigen::SparseMatrix<double>>& Theta_pre, const int N_thread)
{
    const int Fold_k = static_cast<int>(plan.folds.size());
    const int N_Alpha = static_cast<int>(alpha_list.size());
    if (Fold_k <= 0) {
        throw std::invalid_argument("Fold_k must be > 0.");
    }
    if (N_Alpha <= 0) {
        throw std::invalid_argument("N_Alpha must be > 0.");
    }
    if (N_thread <= 0) {
        throw std::invalid_argument("N_thread must be > 0.");
    }
    if (static_cast<int>(plan.train_z.size()) != Fold_k ||
        static_cast<int>(plan.test_z.size()) != Fold_k) {
        throw std::invalid_argument("plan.train_z / plan.test_z size mismatch with plan.folds.");
    }
    Alpha_metric Metric(Fold_k, N_Alpha);

    const int total_jobs = Fold_k * N_Alpha;
    std::vector<omp_lock_t> locks(N_Alpha);
    for (int i = 0; i < N_Alpha; ++i) {
        omp_init_lock(&locks[i]);
    }
    omp_set_num_threads(N_thread);
#pragma omp parallel for schedule(dynamic)
    for (int job = 0; job < total_jobs; ++job) {
        int tid = omp_get_thread_num();
        const int k = job / N_Alpha;
        const int alp_id = job % N_Alpha;
        const double alpha_val = alpha_list[alp_id];
#pragma omp critical
        {
            std::cout << System_time() << '\t'
                << "Procedure | alpha: " << alp_id + 1 << '\t'
                << "fold: " << k + 1
                << " on thread " << tid + 1 << "\n";
        }
        auto t1 = std::chrono::high_resolution_clock::now();

        const FoldRanges& fr = plan.folds[k];
        const ZScoreParams& train_z = plan.train_z[k];
        const ZScoreParams& test_z = plan.test_z[k];
        // Copy the input
        Input_format Input_train = Copy_input_alpha(Input_train_base, alpha_val);
        omp_set_lock(&locks[alp_id]);
        Input_train.Set_Int_theta(Theta_pre[alp_id]);
        omp_unset_lock(&locks[alp_id]);
        // Build readers
        BlockReader Train_reader = Make_trainreader(plan.Data_perm, fr, train_z);
        const int C_train = fr.n_train();
        // z-score and initial out
        Input_train.Set_Observed_cov(LHG_train[k]);
        Eigen::MatrixXd S_test_obs = Cov_test[k];
        // solve
        Solution_theta_and_hg_factor_fold(Train_reader, C_train, Input_train, Cov_train[k]);
        omp_set_lock(&locks[alp_id]);
        Theta_pre[alp_id] = Input_train.Get_theta();
        omp_unset_lock(&locks[alp_id]);
        // metrics: likelihood, entropy, and magnitude used in E_adj(alpha)
        const std::pair<double, double> Alpha_metric_pair =
            Likelihood_marginal_theta_and_global_magnitude(Input_train, S_test_obs);
        const double Likelihood_m = Alpha_metric_pair.first;
        const double Global_magnitude = Alpha_metric_pair.second;
        const double Entropy_lg = Entropy_theta_lg(Input_train);
        Metric.Likelihood_mat(k, alp_id) = Likelihood_m;
        Metric.Entropy_mat(k, alp_id) = Entropy_lg;
        Metric.Global_magnitude_mat(k, alp_id) = Global_magnitude;
        auto t2 = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(t2 - t1).count();
#pragma omp critical
        {
            std::cout << System_time() << '\t'
                << std::fixed << std::setprecision(3)
                << "Procedure | alpha: " << alp_id + 1 << '\t'
                << "fold " << k + 1
                << " : " << elapsed << " s\n";
        }
    }
    return Metric;
}

Lambda_metric Compute_orig_lambda_parallel(
    const KFoldPlan& plan,
    const std::vector<double>& lambda_list,
    const Input_format& Input_train_base,
    const std::vector<Eigen::MatrixXd>& Cov_train,
    const Input_format& Input_test_base,
    const std::vector<Eigen::MatrixXd>& Cov_test,
    std::vector<Eigen::SparseMatrix<double>>& Theta_pre,
    const int N_thread)
{
    const int Fold_k = static_cast<int>(plan.folds.size());
    const int N_Lambda = static_cast<int>(lambda_list.size());

    if (Fold_k <= 0) {
        throw std::invalid_argument("Fold_k must be > 0.");
    }
    if (N_Lambda <= 0) {
        throw std::invalid_argument("N_Lambda must be > 0.");
    }
    if (N_thread <= 0) {
        throw std::invalid_argument("N_thread must be > 0.");
    }
    if (static_cast<int>(plan.train_z.size()) != Fold_k ||
        static_cast<int>(plan.test_z.size()) != Fold_k) {
        throw std::invalid_argument("plan.train_z / plan.test_z size mismatch with plan.folds.");
    }
    if (static_cast<int>(Theta_pre.size()) != N_Lambda) {
        throw std::invalid_argument("Theta_pre size mismatch with lambda_list.");
    }

    Lambda_metric Score(Fold_k, N_Lambda);
    omp_set_num_threads(N_thread);

    // ------------------------------------------------------------
    // Stage 1:
    // Compute fold 0 for all lambdas in parallel.
    // Only fold 0 is allowed to update Theta_pre[lam_id].
    // ------------------------------------------------------------
#pragma omp parallel for schedule(dynamic)
    for (int lam_id = 0; lam_id < N_Lambda; ++lam_id) {
        const int k = 0;
        const int tid = omp_get_thread_num();
        const double lambda_val = lambda_list[lam_id];
#pragma omp critical
        {
            std::cout << System_time() << '\t'
                << "Stage 1 | lambda: " << lam_id + 1 <<'\t'
                << "fold: " << k + 1
                << " on thread " << tid + 1 << "\n";
        }

        auto t1 = std::chrono::high_resolution_clock::now();

        Input_format Input_train = Copy_input_lambda(Input_train_base, lambda_val);
        Input_format Input_test = Copy_input_lambda(Input_test_base, lambda_val);
        Input_train.Set_Observed_cov(Cov_train[k]);
        Input_test.Set_Observed_cov(Cov_test[k]);
        // Use the current shared warm start for this lambda.
        Input_train.Set_Int_theta(Theta_pre[lam_id]);
        Newton_method_offdiag_l(Input_train);
        Theta_pre[lam_id] = Input_train.Get_theta();
        Input_test.Set_Int_theta(Theta_pre[lam_id]);
        Newton_method_offdiag_l(Input_test);
        const std::pair<double, double> Pair_l = Score_local_nd(Input_train.Get_theta(),
            Input_test.Get_theta(), Input_train.Get_Nlocal());
        Score.Share_mat(k, lam_id) = Pair_l.first;
        Score.Difference_mat(k, lam_id) = Pair_l.second;
        auto t2 = std::chrono::high_resolution_clock::now();
        const double elapsed = std::chrono::duration<double>(t2 - t1).count();
#pragma omp critical
        {
            std::cout << System_time() << '\t'
                << std::fixed << std::setprecision(3)
                << "Stage 1 | lambda: " << lam_id  + 1 << '\t'
                << "fold " << k + 1
                << " : " << elapsed << " s\n";
        }
    }
    // ------------------------------------------------------------
    // Stage 2:
    // After all fold-0 tasks are completed, compute folds 1..K-1.
    // ------------------------------------------------------------
    if (Fold_k > 1) {
        const int total_jobs = N_Lambda * (Fold_k - 1);

#pragma omp parallel for schedule(dynamic)
        for (int job = 0; job < total_jobs; ++job) {
            const int lam_id = job / (Fold_k - 1);
            const int k = 1 + (job % (Fold_k - 1));
            const int tid = omp_get_thread_num();
            const double lambda_val = lambda_list[lam_id];
#pragma omp critical
            {
                std::cout << System_time() << '\t'
                    << "Stage 2 | lambda: " << lam_id + 1 << '\t'
                    << "fold: " << k + 1
                    << " on thread " << tid + 1 << "\n";
            }
            auto t1 = std::chrono::high_resolution_clock::now();
            Input_format Input_train = Copy_input_lambda(Input_train_base, lambda_val);
            Input_format Input_test = Copy_input_lambda(Input_test_base, lambda_val);
            Input_train.Set_Observed_cov(Cov_train[k]);
            Input_test.Set_Observed_cov(Cov_test[k]);
            // Use the theta produced by fold 0 as a read-only warm start.
            Input_train.Set_Int_theta(Theta_pre[lam_id]);
            Input_test.Set_Int_theta(Theta_pre[lam_id]);
            Newton_method_offdiag_l(Input_train);
            Newton_method_offdiag_l(Input_test);
            const std::pair<double, double> Pair_l = Score_local_nd(Input_train.Get_theta(), 
                Input_test.Get_theta(), Input_train.Get_Nlocal());
            Score.Share_mat(k, lam_id) = Pair_l.first;
            Score.Difference_mat(k, lam_id) = Pair_l.second;
            auto t2 = std::chrono::high_resolution_clock::now();
            const double elapsed = std::chrono::duration<double>(t2 - t1).count();
#pragma omp critical
            {
                std::cout << System_time() << '\t'
                    << std::fixed << std::setprecision(3)
                    << "Stage 2 | lambda: " << lam_id + 1 << '\t'
                    << "fold " << k + 1
                    << " : " << elapsed << " s\n";
            }
        }
    }
    return Score;
}
