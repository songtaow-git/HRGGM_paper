#include "Class_define.h"

Input_format::Input_format(const Eigen::MatrixXd& O_cov,
    const Eigen::SparseMatrix<double>& I_theta,
    const int N_l, const int N_h, const int N_g,
    const double Lam, const double Alp, const double Gam,    
    const int M_iteration, const double C_scale, const bool I_save): 
    Nlocal(N_l),
    Nglobal(N_g),
    Nhidden(N_h),
    Lambda(Lam),
    Alpha(Alp),
    Gamma(Gam),
    Observed_cov(O_cov),
    Int_theta(I_theta),
    Max_iteration(M_iteration),
    Convergence_scale(C_scale),
    If_save(I_save)
{
}

Input_format::Input_format(const Input_format& I_format)
    : Nlocal(I_format.Nlocal),
    Nglobal(I_format.Nglobal),
    Nhidden(I_format.Nhidden),
    Lambda(I_format.Lambda),
    Alpha(I_format.Alpha),
    Gamma(I_format.Gamma),
    Observed_cov(I_format.Observed_cov),
    Int_theta(I_format.Int_theta),
    Max_iteration(I_format.Max_iteration),
    Convergence_scale(I_format.Convergence_scale),
    If_save(I_format.If_save)
{
}

void Input_format::Set_Observed_cov(const Eigen::MatrixXd& O_cov)
{
	Observed_cov = O_cov;
}

void Input_format::Set_Int_theta(const Eigen::SparseMatrix<double>& I_theta)
{
	Int_theta = I_theta;
}

void Input_format::Set_Lambda(double Lam) 
{
    Lambda = Lam;
}

void Input_format::Set_Alpha(double Alp)
{
    Alpha = Alp;
}

void Input_format::Print_input_format(const int P_idx)
{
    switch (P_idx) {
    case 0:
        std::cout << "Observed_cov = " << Observed_cov << std::endl;
        break;

    case 1:
        std::cout << "Nlocal = " << Nlocal << std::endl;
        std::cout << "Nglobal = " << Nglobal << std::endl;
        std::cout << "Nhidden = " << Nhidden << std::endl;
        break;

    case 2:
        std::cout << "Lambda = " << Lambda << std::endl;
        std::cout << "Alpha = " << Alpha << std::endl;
        std::cout << "Gamma = " << Gamma << std::endl;
        break;

    case 3:
        std::cout << "Int_theta = " << Int_theta << std::endl;
        break;

    case 4:
        std::cout << "Max_iteration = " << Max_iteration << std::endl;
        break;

    case 5:
        std::cout << "Convergence_scale = " << Convergence_scale << std::endl;
        break;

    case 6:
        std::cout << "If_save = " << If_save << std::endl;
        break;

    default:
        std::cout << "Wrong print idx!" << std::endl;
        break;
    }
}

int Input_format::Get_Nlocal()
{
    return Nlocal;
}

int Input_format::Get_Nhidden()
{
    return Nhidden;
}

int Input_format::Get_Nglobal()
{
    return Nglobal;
}

const Eigen::SparseMatrix<double>& Input_format::Get_theta() const
{
    return Int_theta;
}
const Eigen::MatrixXd& Input_format::Get_observed_cov() const
{
    return Observed_cov;
}

MarkAndPush::MarkAndPush(const int B, std::vector<uint8_t>& Mark,
    std::vector<std::pair<int, int>>& Updates)
    : B_(B), Mark_(Mark), Updates_(Updates) 
{
}

void MarkAndPush::operator()(const int i, const int j)
{
    const size_t idx =
        static_cast<size_t>(i) + static_cast<size_t>(j) * static_cast<size_t>(B_); // col-major
    if (!Mark_[idx]) {
        Mark_[idx] = 1;
        Updates_.emplace_back(i, j); 
    }
}

KKT_residual::KKT_residual(const double Diag, const double LL_offdiag_nz,
    const double LL_offdiag_z, const double LH_nz,
    const double LH_z, const double LG):
    KKT_r_diag(Diag),
    KKT_r_ll_offdiag_nz(LL_offdiag_nz),
    KKT_r_ll_offdiag_z(LL_offdiag_z),
    KKT_r_lh_nz(LH_nz),
    KKT_r_lh_z(LH_z),
    KKT_r_lg(LG)
{ 
    KKT_error = std::max({ KKT_r_diag, KKT_r_ll_offdiag_nz, KKT_r_ll_offdiag_z,
        KKT_r_lh_nz, KKT_r_lh_z, KKT_r_lg });
}

KKT_residual::KKT_residual(const KKT_residual& KKT_r):
    KKT_r_diag(KKT_r.KKT_r_diag),
    KKT_r_ll_offdiag_nz(KKT_r.KKT_r_ll_offdiag_nz),
    KKT_r_ll_offdiag_z(KKT_r.KKT_r_ll_offdiag_z),
    KKT_r_lh_nz(KKT_r.KKT_r_lh_nz),
    KKT_r_lh_z(KKT_r.KKT_r_lh_z),
    KKT_r_lg(KKT_r.KKT_r_lg)
{
    KKT_error = std::max({ KKT_r_diag, KKT_r_ll_offdiag_nz, KKT_r_ll_offdiag_z,
        KKT_r_lh_nz, KKT_r_lh_z, KKT_r_lg });
}

void KKT_residual::KKT_print()
{
    std::cout << "KKT:" << std::endl;
    std::cout << "Diag = " << '\t' << KKT_r_diag << std::endl;
    std::cout << "LL_offdiag_nz= " << '\t' << KKT_r_ll_offdiag_nz << std::endl;
    std::cout << "LL_offdiag_z =" << '\t' << KKT_r_ll_offdiag_z << std::endl;
    std::cout << "LH_nz = " << '\t' << KKT_r_lh_nz << std::endl;
    std::cout << "LH_z = " << '\t' << KKT_r_lh_z << std::endl;
    std::cout << "LG = " << '\t' << KKT_r_lg << std::endl;
}

Lambda_metric::Lambda_metric(const int N_f, const int N_l)
{
    Share_mat = RowMat::Zero(N_f, N_l);
    Difference_mat = RowMat::Zero(N_f, N_l);
}

void Lambda_metric::Compute_score()
{
    const double Zero_r = 1e-6;
    RowMat Share_vec = Share_mat.colwise().mean();
    RowMat Difference_vec = Difference_mat.colwise().mean();
    double Share_max = Share_vec.maxCoeff();
    double Difference_max = Difference_vec.maxCoeff();
    RowMat Share_norm = Share_vec / Share_max;
    RowMat Difference_norm = Difference_vec / Difference_max;
    Score_norm = Share_norm - Difference_norm;
    Score_ratio = (Share_norm.array() * Share_norm.array()) / (Difference_norm.array() + Zero_r);
}

namespace
{
    void Zscore_row_safe(RowMat& X, const double Epsilon)
    {
        const int C = static_cast<int>(X.cols());
        if (C <= 1) {
            X.setZero();
            return;
        }
        double* p = X.data();
        double sum = 0.0;
        for (int j = 0; j < C; ++j) {
            sum += p[j];
        }
        const double mu = sum / static_cast<double>(C);
        double ss = 0.0;
        for (int j = 0; j < C; ++j) {
            const double d = p[j] - mu;
            ss += d * d;
        }
        const double sigma = std::sqrt(ss / static_cast<double>(C - 1));
        if (sigma > Epsilon && std::isfinite(sigma)) {
            const double inv = 1.0 / sigma;
            for (int j = 0; j < C; ++j) {
                p[j] = (p[j] - mu) * inv;
            }
        }
        else {
            X.setZero();
        }
    }
}

Alpha_metric::Alpha_metric(const int N_f, const int N_a)
{
    Likelihood_mat = RowMat::Zero(N_f, N_a);
    Entropy_mat = RowMat::Zero(N_f, N_a);
    Global_magnitude_mat = RowMat::Zero(N_f, N_a);
    Retention_vec = RowMat::Zero(1, N_a);
    Adjusted_entropy_vec = RowMat::Zero(1, N_a);
    Score = RowMat::Zero(1, N_a);
}

void Alpha_metric::Compute_score(const double Tau, const double Epsilon)
{
    if (!(Tau > 0.0 && Tau <= 1.0)) {
        throw std::invalid_argument("Alpha_metric::Compute_score: Tau must be in (0, 1].");
    }
    if (!(Epsilon > 0.0)) {
        throw std::invalid_argument("Alpha_metric::Compute_score: Epsilon must be positive.");
    }

    RowMat Likelihood_vec = Likelihood_mat.colwise().mean();
    const RowMat Entropy_vec = Entropy_mat.colwise().mean();
    const RowMat Magnitude_vec = Global_magnitude_mat.colwise().mean();
    const int C = static_cast<int>(Likelihood_vec.cols());

    const double M_ref = Magnitude_vec.maxCoeff();
    Retention_vec.setZero(1, C);
    if (std::isfinite(M_ref) && M_ref > Epsilon) {
        const double threshold = Tau * M_ref + Epsilon;
        for (int j = 0; j < C; ++j) {
            const double raw = Magnitude_vec(0, j) / threshold;
            Retention_vec(0, j) = std::min(1.0, std::max(0.0, raw));
        }
    }

    Adjusted_entropy_vec = Entropy_vec.array() * Retention_vec.array();

    RowMat Z_likelihood = Likelihood_vec;
    RowMat Z_adjusted_entropy = Adjusted_entropy_vec;
    Zscore_row_safe(Z_likelihood, Epsilon);
    Zscore_row_safe(Z_adjusted_entropy, Epsilon);

    // Smaller values indicate a better alpha.
    Score = 0.5 * Z_likelihood - 0.5 * Z_adjusted_entropy;
}
