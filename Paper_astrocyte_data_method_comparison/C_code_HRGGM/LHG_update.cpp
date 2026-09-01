#include "LHG_update.h"

RowMat Vstack(const RowMat& A, const RowMat& B)
{
    if (A.cols() != B.cols())
        throw std::runtime_error("Column mismatch.");
    RowMat C(A.rows() + B.rows(), A.cols());
    C.topRows(A.rows()) = A;
    C.bottomRows(B.rows()) = B;
    return C;
}

Eigen::MatrixXd Cov_scdata(const RowMat& Data)
{
    const int R = static_cast<int>(Data.rows());
    const int C = static_cast<int>(Data.cols());
    if (C <= 1)
        throw std::runtime_error("Need at least 2 observations (cols) for cov.");
    // 1) mu = row means
    Eigen::VectorXd mu = Data.rowwise().mean();   // (R)
    // 2) G = Data * Data^T   (ColMajor for speed)
    Eigen::MatrixXd G = Data * Data.transpose();  // (R x R)
    // 3) subtract C * mu*mu^T
    G.noalias() -= double(C) * (mu * mu.transpose());
    // 4) divide by (C-1)
    G /= double(C - 1);
    // 5) convert to RowMajor at return
    return G;
}

RowMat Topfactor_eigcov_pca(const BlockReader& read_block_raw,
    const int Rsel, const int C, const int nGlobals, const int BlockCols)
{
    if (Rsel <= 0) throw std::runtime_error("Rsel must be > 0.");
    if (C <= 1)    throw std::runtime_error("C must be > 1.");
    if (nGlobals <= 0) throw std::runtime_error("nGlobals must be > 0.");
    const int k = std::min(nGlobals, Rsel);
    RowMat Bfull(Rsel, BlockCols);
    // ---------- 1) mean ----------
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(Rsel);
    for (int col0 = 0; col0 < C; col0 += BlockCols) {
        const int bs = std::min(BlockCols, C - col0);
        auto B = Bfull.leftCols(bs);               // view: (Rsel x bs)
        read_block_raw(col0, bs, B);
        sum.noalias() += B.rowwise().sum();
    }
    const Eigen::VectorXd mu = sum / double(C);
    // ---------- 2) covariance (ColMajor) ----------
    Eigen::MatrixXd Cov = Eigen::MatrixXd::Zero(Rsel, Rsel);
    for (int col0 = 0; col0 < C; col0 += BlockCols) {
        const int bs = std::min(BlockCols, C - col0);
        auto B = Bfull.leftCols(bs);
        read_block_raw(col0, bs, B);
        // center in-place
        B.colwise() -= mu;
        Cov.noalias() += B * B.transpose();
    }
    Cov /= double(C - 1);
    // ---------- 3) eigen-decompose ----------
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Cov);
    if (es.info() != Eigen::Success) throw std::runtime_error("Eigen decomposition failed.");
    Eigen::MatrixXd V_k(Rsel, k);
    for (int j = 0; j < k; ++j)
        V_k.col(j) = es.eigenvectors().col(Rsel - 1 - j);
    // ---------- 4) TF (RowMajor) ----------
    RowMat TF(k, C);
    for (int col0 = 0; col0 < C; col0 += BlockCols) {
        const int bs = std::min(BlockCols, C - col0);
        auto B = Bfull.leftCols(bs);
        read_block_raw(col0, bs, B);
        B.colwise() -= mu;
        // TF_block = V_k^T * B  => (k x bs)
        TF.middleCols(col0, bs).noalias() = V_k.transpose() * B;
    }
    // ---------- 5) deterministic sign fix ----------
    for (int i = 0; i < k; ++i) {
        Eigen::Index idx;
        TF.row(i).cwiseAbs().maxCoeff(&idx);
        if (TF(i, idx) < 0) TF.row(i) *= -1.0;
    }
    return TF;
}

void Zscore_DataRM(RowMat& DataRM)
{
    const int k = static_cast<int>(DataRM.rows());
    const int C = static_cast<int>(DataRM.cols());
    if (C <= 1) return;
    for (int i = 0; i < k; ++i) {
        double* p = DataRM.row(i).data(); // RowMajor
        // pass 1: mean
        double sum = 0.0;
        for (int j = 0; j < C; ++j) sum += p[j];
        const double mu = sum / double(C);
        // pass 2: sample variance
        double ss = 0.0;
        for (int j = 0; j < C; ++j) {
            const double d = p[j] - mu;
            ss += d * d;
        }
        const double sigma = std::sqrt(ss / double(C - 1));
        if (sigma > 0.0 && std::isfinite(sigma)) {
            const double inv = 1.0 / sigma;
            for (int j = 0; j < C; ++j) {
                p[j] = (p[j] - mu) * inv;
            }
        }
        else {
            constexpr double nanv = std::numeric_limits<double>::quiet_NaN();
            for (int j = 0; j < C; ++j) p[j] = nanv;
        }
    }
}

RowMat Topfactor_zscore_out(RowMat& Data, const int nHidden, const int nGlobal,
    const std::vector<int>& Idx_h, const std::vector<int>& Idx_g, const int BlockCols)
{
    Zscore_DataRM(Data);
    int C = (int)Data.cols();
    int Rsel_h = (int)Idx_h.size();
    BlockReader Reader_h = [&](int col0, int bs, Eigen::Ref<RowMat> out) {
        // out: (Rsel x bs) row-major
        for (int i = 0; i < out.rows(); ++i) {
            int g = Idx_h[i];
            out.row(i) = Data.row(g).segment(col0, bs);
        }
        };
    RowMat Hidden_factor = Topfactor_eigcov_pca(Reader_h, Rsel_h, C, nHidden, BlockCols);
    int Rsel_g = (int)Idx_g.size();
    BlockReader Reader_g = [&](int col0, int bs, Eigen::Ref<RowMat> out) {
        // out: (Rsel x bs) row-major
        for (int i = 0; i < out.rows(); ++i) {
            int g = Idx_g[i];
            out.row(i) = Data.row(g).segment(col0, bs);
        }
        };
    RowMat Global_factor = Topfactor_eigcov_pca(Reader_g, Rsel_g, C, nGlobal, BlockCols);
    RowMat HG_factor = Vstack(Hidden_factor, Global_factor);
    Zscore_DataRM(HG_factor);
    return HG_factor;
}

void LHG_zscore_out(RowMat& Data, Input_format& Input_l, const std::vector<int>& Idx_h,
    const std::vector<int>& Idx_g, const int BlockCols)
{
    Zscore_DataRM(Data);
    int C = (int)Data.cols();
    int Rsel_h = (int)Idx_h.size();
    BlockReader Reader_h = [&](int col0, int bs, Eigen::Ref<RowMat> out) {
        // out: (Rsel x bs) row-major
        for (int i = 0; i < out.rows(); ++i) {
            int g = Idx_h[i];
            out.row(i) = Data.row(g).segment(col0, bs);
        }
        };
    RowMat Hidden_factor = Topfactor_eigcov_pca(Reader_h, Rsel_h, C, Input_l.Get_Nhidden(), BlockCols);
    int Rsel_g = (int)Idx_g.size();
    BlockReader Reader_g = [&](int col0, int bs, Eigen::Ref<RowMat> out) {
        // out: (Rsel x bs) row-major
        for (int i = 0; i < out.rows(); ++i) {
            int g = Idx_g[i];
            out.row(i) = Data.row(g).segment(col0, bs);
        }
        };
    RowMat Global_factor = Topfactor_eigcov_pca(Reader_g, Rsel_g, C, Input_l.Get_Nglobal(), BlockCols);
    RowMat HG_factor = Vstack(Hidden_factor, Global_factor);
    Zscore_DataRM(HG_factor);
    RowMat LHG_Data = Vstack(Data, HG_factor);
    Input_l.Set_Observed_cov(Cov_scdata(LHG_Data));
}

RowMat Update_hg_factor(const RowMat& Data, Input_format& Input_l)
{
    const int Nl = Input_l.Get_Nlocal();
    const int Nhg = Input_l.Get_Nhidden() + Input_l.Get_Nglobal();
    const int N = Nl + Nhg;
    const int nSamples = static_cast<int>(Data.cols());
    RowMat E(Nhg, nSamples);
    E.setZero();
    std::vector<double> invDiag(Nhg, 0.0);
    for (int c = Nl; c < N; ++c) {
        double diag = 0.0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(Input_l.Get_theta(), c); it; ++it) {
            if (it.row() == c) { diag = it.value(); break; }
        }
        if (diag == 0.0) throw std::runtime_error("Theta(c,c) is zero.");
        invDiag[c - Nl] = 1.0 / diag;
    }
    Eigen::RowVectorXd acc(nSamples);
    for (int c = Nl; c < N; ++c) {
        acc.setZero();
        for (Eigen::SparseMatrix<double>::InnerIterator it(Input_l.Get_theta(), c); it; ++it) {
            const int r = it.row();
            if (r < Nl) {
                // Theta(c,r) = Theta(r,c) by symmetry (Upper stored)
                const double val = it.value();
                acc.noalias() += val * Data.row(r);
            }
        }
        E.row(c - Nl) = -(invDiag[c - Nl]) * acc;
    }
    RowMat E_zscore = E;
    Zscore_DataRM(E_zscore);
    RowMat Data_com = Vstack(Data, E_zscore);
    Input_l.Set_Observed_cov(Cov_scdata(Data_com));
    return E_zscore;
}

bool SparseMatricesClose(const Eigen::SparseMatrix<double>& A,
    const Eigen::SparseMatrix<double>& B, const double D_theta)
{
    if (A.rows() != B.rows() || A.cols() != B.cols())
        throw std::runtime_error("Matrix size mismatch.");
    if (A.rows() != A.cols())
        throw std::runtime_error("A must be square.");
    if (D_theta < 0)
        throw std::runtime_error("D_theta must be non-negative.");
    for (int outer = 0; outer < A.outerSize(); ++outer)
    {
        Eigen::SparseMatrix<double>::InnerIterator itA(A, outer);
        Eigen::SparseMatrix<double>::InnerIterator itB(B, outer);
        while (itA || itB)
        {
            int row, col;

            if (itA && (!itB || itA.index() < itB.index())) {
                row = itA.index();
                col = outer; // ColMajor: outer
                if (row <= col) {
                    if (std::abs(itA.value()) >= D_theta) return false;
                }
                ++itA;
            }
            else if (itB && (!itA || itB.index() < itA.index())) {
                row = itB.index();
                col = outer;
                if (row <= col) {
                    if (std::abs(itB.value()) >= D_theta) return false;
                }
                ++itB;
            }
            else {
                row = itA.index(); // == itB.index()
                col = outer;
                if (row <= col) {
                    if (std::abs(itA.value() - itB.value()) >= D_theta) return false;
                }
                ++itA;
                ++itB;
            }
        }
    }
    return true;
}

RowMat Solution_theta_and_hg_factor(RowMat& Data, Input_format& Input_l, 
    const std::vector<int>& Idx_h, const std::vector<int>& Idx_g, const int Max_loop, 
    const double Stop_threshold, const int Blockcol_inv, const int BlockCol_data)
{
    LHG_zscore_out(Data, Input_l, Idx_h, Idx_g, BlockCol_data);
    RowMat  Updated_hg;
    Eigen::SparseMatrix<double> Old_theta = Input_l.Get_theta();
    for (int i = 0; i < Max_loop; ++i) {
        Newton_method_offdiag_lhg(Input_l, Blockcol_inv);
        Updated_hg = Update_hg_factor(Data, Input_l);
        bool If_break = SparseMatricesClose(Old_theta, Input_l.Get_theta(), Stop_threshold);
        if (If_break) {
            std::cout << "Loop is completed:" << '\t' << i + 1 << '\t' << "times;" << std::endl;
            break;
        }
        Old_theta = Input_l.Get_theta();
    }
    return Updated_hg;
}

#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <stdexcept>

double Likelihood_marginal_theta(Input_format& Input_l, const Eigen::MatrixXd& Cov_test)
{
    const int Nl = Input_l.Get_Nlocal();
    const int N = Nl + Input_l.Get_Nhidden() + Input_l.Get_Nglobal();
    const int Nr = N - Nl;
    Eigen::SparseMatrix<double> Theta_full = Input_l.Get_theta();
    if (Theta_full.rows() != N || Theta_full.cols() != N)
        throw std::runtime_error("Theta_full size mismatch.");
    const int R0 = Nl;
    const int R1 = N;
    Eigen::MatrixXd LR = Eigen::MatrixXd::Zero(Nl, Nr);
    for (int jj = 0; jj < Nr; ++jj)
    {
        const int col = R0 + jj;
        for (Eigen::SparseMatrix<double>::InnerIterator it(Theta_full, col); it; ++it)
        {
            const int row = it.row();
            if (row < Nl)
            {
                LR(row, jj) = it.value();
            }
        }
    }
    Eigen::Matrix<double, 5, 5> RRbuf;
    RRbuf.setZero();
    for (int jj = 0; jj < Nr; ++jj)
    {
        const int col = R0 + jj;
        for (Eigen::SparseMatrix<double>::InnerIterator it(Theta_full, col); it; ++it)
        {
            const int row = it.row();
            if (row >= R0 && row < R1)
            {
                const int ii = row - R0;
                // Only consume upper triangle in RR (ii <= jj) and mirror
                if (ii <= jj)
                {
                    const double v = it.value();
                    RRbuf(ii, jj) = v;
                    RRbuf(jj, ii) = v;
                }
            }
        }
    }
    Eigen::MatrixXd RR = RRbuf.topLeftCorner(Nr, Nr);
    Eigen::MatrixXd LL = Eigen::MatrixXd::Zero(Nl, Nl);
    for (int col = 0; col < Nl; ++col)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Theta_full, col); it; ++it)
        {
            const int row = it.row();
            if (row < Nl)
            {
                // Only consume upper triangle in LL (row <= col) and mirror
                if (row <= col)
                {
                    const double v = it.value();
                    LL(row, col) = v;
                    LL(col, row) = v;
                }
            }
        }
    }
    Eigen::LDLT<Eigen::MatrixXd> ldlt(RR);
    if (ldlt.info() != Eigen::Success)
        throw std::runtime_error("LDLT factorization of RR failed.");
    /* Solve X = RR^{-1} * LRᵀ  (Nr × Nl) */
    Eigen::MatrixXd X = ldlt.solve(LR.transpose());
    /* Compute Schur complement: Theta_marg_LL = LL - LR * X */
    Eigen::MatrixXd Theta_marg_LL = LL;
    Theta_marg_LL.noalias() -= LR * X;
    double Logdet = Logdet_dense(Theta_marg_LL);
    double Trace = Trace_dense_dense(Cov_test, Theta_marg_LL);
    return  Trace - Logdet;
}

double Entropy_theta_lg(Input_format& Input_l)
{
    Eigen::SparseMatrix<double> Theta_full = Input_l.Get_theta();
    // LG columns: [c0, c1)
    const int Nl = Input_l.Get_Nlocal();
    const int Nh = Input_l.Get_Nhidden();
    const int Ng = Input_l.Get_Nglobal();
    const int c0 = Nl + Nh;
    const int c1 = c0 + Ng;
    // n = numel(x) where x = Theta_LG(:) and Theta_LG is Nl×Ng
    const double n = static_cast<double>(Nl) * static_cast<double>(Ng);
    if (n <= 1.0) {
        // log(0) or log(1) makes normalization undefined; return 0 as a safe convention
        return 0.0;
    }
    // Accumulate:
    // s = sum(v)
    // sum_v_logv = sum(v * log(v))
    // where v = abs(Theta(i,j)) for (i in L rows, j in G cols)
    double s = 0.0;
    double sum_v_logv = 0.0;
    for (int col = c0; col < c1; ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Theta_full, col); it; ++it) {
            const int row = it.row();
            if (row >= 0 && row < Nl) {
                const double v = std::abs(it.value());
                if (v > 0.0) { // skip exact zeros (shouldn't appear in sparse, but safe)
                    s += v;
                    sum_v_logv += v * std::log(v);
                }
            }
        }
    }
    // if s == 0, E = 0
    if (s == 0.0)
        return 0.0;
    // Entropy:
    // H = -sum(p * log(p)), p = v/s over v>0
    // Compute stably without storing p:
    // H = log(s) - (1/s) * sum(v * log(v))
    const double H = std::log(s) - (sum_v_logv / s);
    // Normalized entropy: E = H / log(n)
    const double denom = std::log(n);
    if (denom == 0.0) return 0.0;
    return H / denom;
}

double Global_dominance(const RowMat& Data_local, const RowMat& Factor_hg,
    Input_format& Input_l, const int BlockCol_data)
{
    const int Nl = static_cast<int>(Data_local.rows());
    const int Nh = Input_l.Get_Nhidden();
    const int Ng = Input_l.Get_Nglobal();
    const int nSamples = static_cast<int>(Data_local.cols());
    constexpr double eps = std::numeric_limits<double>::epsilon();
    // Accumulate row-wise squared norms for Data_local
    Eigen::VectorXd xnorm_sq = Eigen::VectorXd::Zero(Nl);
    // Accumulate squared norms for each G row (Ng is tiny)
    Eigen::VectorXd gnorm_sq = Eigen::VectorXd::Zero(Ng);
    // Accumulate dot products num(:,g) = Data_local * G(g,:)^T  (Nl×Ng)
    Eigen::MatrixXd num = Eigen::MatrixXd::Zero(Nl, Ng);
    // Temporary vector to hold a G row block as a column vector
    Eigen::VectorXd gvec;
    // Blocked pass over columns to improve cache locality and reduce peak temporaries
    for (int j0 = 0; j0 < nSamples; j0 += BlockCol_data){
        const int B = std::min(BlockCol_data, nSamples - j0);
        // View of Data_local block: Nl × B (no copy)
        const auto Xblk = Data_local.block(0, j0, Nl, B);
        // Update xnorm_sq: sum of squares per row
        xnorm_sq.noalias() += Xblk.array().square().rowwise().sum().matrix();
        // Update gnorm_sq and num(:,g)
        for (int g = 0; g < Ng; ++g)
        {
            const auto Grow = Factor_hg.block(Nh + g, j0, 1, B); // 1 × B view
            // Update gnorm_sq
            gnorm_sq(g) += Grow.array().square().sum();
            // Compute and accumulate num(:,g) += Xblk * Grow^T
            gvec = Grow.transpose(); // B×1
            num.col(g).noalias() += Xblk * gvec; // GEMV: (Nl×B) * (B×1) -> (Nl×1)
        }
    }
    // Final norms
    Eigen::VectorXd xnorm = xnorm_sq.array().sqrt();
    Eigen::VectorXd gnorm = gnorm_sq.array().sqrt();
    xnorm.array() += eps;
    gnorm.array() += eps;
    // Compute dominance score D = max_{i,g} | num(i,g) / (xnorm(i) * gnorm(g)) |
    double D = 0.0;
    for (int g = 0; g < Ng; ++g)
    {
        const double inv_g = 1.0 / gnorm(g);
        for (int i = 0; i < Nl; ++i)
        {
            const double val = std::abs(num(i, g) * inv_g / xnorm(i));
            if (val > D) D = val;
        }
    }
    return D;
}

std::pair<double, double> Score_local_nd(const Eigen::SparseMatrix<double>& Theta_train,
    const Eigen::SparseMatrix<double>& Theta_test, const int Nl)
{
    if (Theta_train.rows() != Theta_train.cols() || Theta_test.rows() != Theta_test.cols())
        throw std::runtime_error("Inputs must be square.");
    if (Theta_train.rows() < Nl || Theta_test.rows() < Nl)
        throw std::runtime_error("Nl exceeds matrix size.");
    double num_ut = 0.0, den_ut = 0.0;
    for (int col = 0; col < Nl; ++col)
    {
        Eigen::SparseMatrix<double>::InnerIterator itA(Theta_train, col);
        Eigen::SparseMatrix<double>::InnerIterator itB(Theta_test, col);
        while (itA && itA.row() < col && itB && itB.row() < col)
        {
            const int ra = itA.row();
            const int rb = itB.row();
            if (ra == rb)
            {
                const double a = itA.value();
                const double b = itB.value();
                const double d = a - b;
                den_ut += d * d;
                if (a * b > 0.0)
                {
                    const double s = std::min(std::abs(a), std::abs(b));
                    num_ut += s * s;
                }
                ++itA; ++itB;
            }
            else if (ra < rb)
            {
                const double a = itA.value();
                den_ut += a * a;
                ++itA;
            }
            else
            {
                const double b = itB.value();
                den_ut += b * b;
                ++itB;
            }
        }
        while (itA && itA.row() < col)
        {
            const double a = itA.value();
            den_ut += a * a;
            ++itA;
        }
        while (itB && itB.row() < col)
        {
            const double b = itB.value();
            den_ut += b * b;
            ++itB;
        }
    }
    return std::pair((2.0 * num_ut), (2.0 * den_ut));
}

Eigen::MatrixXd Expandcov_lhg_reader(const Eigen::MatrixXd& S_LL, const BlockReader& Datareader,
    const int Nl, const int C, const RowMat& H, const int Blockcols)
{
    const int Nhg = static_cast<int>(H.rows());
    if (S_LL.rows() != Nl || S_LL.cols() != Nl)
        throw std::runtime_error("Expandcov_lhg_reader: S_LL size mismatch.");
    if (H.cols() != C)
        throw std::runtime_error("Expandcov_lhg_reader: H must have same number of cols as train data.");
    if (C <= 1)
        throw std::runtime_error("Expandcov_lhg_reader: need at least 2 train samples.");
    if (Blockcols <= 0)
        throw std::runtime_error("Expandcov_lhg_reader: Blockcols must be > 0.");
    // muX and muH are kept for exact consistency with Cov_scdata(...)
    Eigen::VectorXd sumX = Eigen::VectorXd::Zero(Nl);
    Eigen::VectorXd muH = H.rowwise().mean();
    Eigen::MatrixXd GLH = Eigen::MatrixXd::Zero(Nl, Nhg);
    Eigen::MatrixXd GHH = Eigen::MatrixXd::Zero(Nhg, Nhg);
    RowMat Xblk(Nl, std::min(Blockcols, C));
    for (int j0 = 0; j0 < C; j0 += Blockcols)
    {
        const int B = std::min(Blockcols, C - j0);
        // Read z-scored train block: Xblk(:, 0:B-1)
        Datareader(j0, B, Xblk.leftCols(B));
        const auto Xview = Xblk.leftCols(B);
        const auto Hblk = H.block(0, j0, Nhg, B);
        // Accumulate row sums of X for exact Cov_scdata-style centering
        sumX.noalias() += Xview.rowwise().sum();
        // Accumulate cross Gram and HG Gram
        GLH.noalias() += Xview * Hblk.transpose(); // (Nl×B)*(B×Nhg)
        GHH.noalias() += Hblk * Hblk.transpose();  // (Nhg×B)*(B×Nhg)
    }
    const Eigen::VectorXd muX = sumX / double(C);
    const double denom = double(C - 1);
    Eigen::MatrixXd S_LH =
        (GLH - double(C) * (muX * muH.transpose())) / denom;
    Eigen::MatrixXd S_HH =
        (GHH - double(C) * (muH * muH.transpose())) / denom;
    Eigen::MatrixXd S_LHG(Nl + Nhg, Nl + Nhg);
    S_LHG.topLeftCorner(Nl, Nl) = S_LL;
    S_LHG.topRightCorner(Nl, Nhg) = S_LH;
    S_LHG.bottomLeftCorner(Nhg, Nl) = S_LH.transpose();
    S_LHG.bottomRightCorner(Nhg, Nhg) = S_HH;
    return S_LHG;
}

KFoldPlan MakeKfoldplan_reorderdata(const RowMat& Data,
    const int Fold_k,unsigned int seed)
{
    const int nRows = static_cast<int>(Data.rows());
    const int nSamples = static_cast<int>(Data.cols());
    if (Fold_k < 2)
        throw std::runtime_error("Fold_k must be >= 2.");
    if (Fold_k > nSamples)
        throw std::runtime_error("Fold_k cannot exceed number of samples.");
    KFoldPlan plan;
    // 1) randperm
    std::vector<int> perm(nSamples);
    std::iota(perm.begin(), perm.end(), 0);
    std::mt19937 rng(seed);
    std::shuffle(perm.begin(), perm.end(), rng);
    // 2) equal fold size, drop remainder
    const int fold_size = nSamples / Fold_k;
    const int n_used = fold_size * Fold_k;
    plan.used_idx.assign(perm.begin(), perm.begin() + n_used);
    plan.dropped_idx.assign(perm.begin() + n_used, perm.end());
    // 3) reorder Data -> Data_perm = Data(:, used_idx)
    plan.Data_perm.resize(nRows, n_used);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < nRows; ++i) {
        const double* src = Data.row(i).data();
        double* dst = plan.Data_perm.row(i).data();
        for (int j = 0; j < n_used; ++j) {
            dst[j] = src[plan.used_idx[j]];
        }
    }
    // 4) fold ranges
    plan.folds.resize(Fold_k);
    for (int k = 0; k < Fold_k; ++k) {
        FoldRanges fr;
        fr.test_begin = k * fold_size;
        fr.test_end = (k + 1) * fold_size;
        fr.n_used = n_used;
        fr.fold_size = fold_size;
        plan.folds[k] = fr;
    }
    return plan;
}

SumStat Compute_sumsq_onerange(const RowMat& Data_perm, const int c0, const int c1)
{
    const int R = static_cast<int>(Data_perm.rows());
    if (c0 < 0 || c1 < c0 || c1 > Data_perm.cols())
        throw std::runtime_error("ComputeSumSq_OneRange: invalid range.");
    SumStat st;
    st.sum = Eigen::VectorXd::Zero(R);
    st.sqsum = Eigen::VectorXd::Zero(R);
    st.count = c1 - c0;
    for (int i = 0; i < R; ++i) {
        const double* p = Data_perm.row(i).data();
        double s = 0.0, ss = 0.0;
        for (int j = c0; j < c1; ++j) {
            const double x = p[j];
            s += x;
            ss += x * x;
        }
        st.sum(i) = s;
        st.sqsum(i) = ss;
    }
    return st;
}

SumStat Mergesumstat(const SumStat& a, const SumStat& b)
{
    if (a.sum.size() != b.sum.size())
        throw std::runtime_error("MergeSumStat: size mismatch.");
    SumStat out;
    out.sum = a.sum + b.sum;
    out.sqsum = a.sqsum + b.sqsum;
    out.count = a.count + b.count;
    return out;
}

ZScoreParams Build_zscoreparams(const SumStat& st)
{
    const int R = static_cast<int>(st.sum.size());
    if (st.count <= 1)
        throw std::runtime_error("Build_zscoreparams: count must be > 1.");
    ZScoreParams zp;
    zp.mu.resize(R);
    zp.invsigma.resize(R);
    const double C = static_cast<double>(st.count);
    for (int i = 0; i < R; ++i) {
        const double mu = st.sum(i) / C;
        const double centered_ss = st.sqsum(i) - C * mu * mu;
        const double var = centered_ss / double(st.count - 1);
        zp.mu(i) = mu;
        if (var > 0.0 && std::isfinite(var)) {
            const double sigma = std::sqrt(var);
            zp.invsigma(i) = 1.0 / sigma;
        }
        else {
            zp.invsigma(i) = 0.0;
        }
    }
    return zp;
}

void Precompute_foldzscores(KFoldPlan& plan)
{
    const int K = static_cast<int>(plan.folds.size());
    plan.train_z.resize(K);
    plan.test_z.resize(K);
    for (int k = 0; k < K; ++k) {
        const FoldRanges& fr = plan.folds[k];
        // test = [test_begin, test_end)
        SumStat st_test = Compute_sumsq_onerange(
            plan.Data_perm, fr.test_begin, fr.test_end);
        // train = [0, test_begin) U [test_end, n_used)
        SumStat st_left = Compute_sumsq_onerange(
            plan.Data_perm, 0, fr.test_begin);
        SumStat st_right = Compute_sumsq_onerange(
            plan.Data_perm, fr.test_end, fr.n_used);
        SumStat st_train = Mergesumstat(st_left, st_right);
        plan.test_z[k] = Build_zscoreparams(st_test);
        plan.train_z[k] = Build_zscoreparams(st_train);
    }
}

BlockReader Make_trainreader(const RowMat& Data_perm, const FoldRanges& fr, 
    const ZScoreParams& zp)
{
    const int R = static_cast<int>(Data_perm.rows());

    return [&Data_perm, &fr, &zp, R](int col0, int bs, Eigen::Ref<RowMat> out)
        {
            if (out.rows() != R || out.cols() != bs)
                throw std::runtime_error("Make_trainreader: out size mismatch.");
            for (int i = 0; i < R; ++i)
            {
                const double* src = Data_perm.row(i).data();
                double* dst = out.row(i).data();
                const double mu = zp.mu(i);
                const double inv = zp.invsigma(i);
                if (inv <= 0.0 || !std::isfinite(inv))
                {
                    constexpr double nanv = std::numeric_limits<double>::quiet_NaN();
                    for (int j = 0; j < bs; ++j) dst[j] = nanv;
                    continue;
                }
                for (int j = 0; j < bs; ++j)
                {
                    const int tcol = col0 + j;
                    const int real_col =
                        (tcol < fr.test_begin) ? tcol : (tcol + fr.fold_size);
                    dst[j] = (src[real_col] - mu) * inv;
                }
            }
        };
}

BlockReader Make_testreader(const RowMat& Data_perm, const FoldRanges& fr, 
    const ZScoreParams& zp)
{
    const int R = static_cast<int>(Data_perm.rows());
    return [&Data_perm, &fr, &zp, R](int col0, int bs, Eigen::Ref<RowMat> out)
        {
            if (out.rows() != R || out.cols() != bs)
                throw std::runtime_error("MakeTestReader: out size mismatch.");
            for (int i = 0; i < R; ++i) {
                const double* src = Data_perm.row(i).data();
                double* dst = out.row(i).data();
                const double mu = zp.mu(i);
                const double inv = zp.invsigma(i);
                if (inv <= 0.0 || !std::isfinite(inv)) {
                    constexpr double nanv = std::numeric_limits<double>::quiet_NaN();
                    for (int j = 0; j < bs; ++j) dst[j] = nanv;
                    continue;
                }
                for (int j = 0; j < bs; ++j) {
                    const int real_col = fr.test_begin + col0 + j;
                    dst[j] = (src[real_col] - mu) * inv;
                }
            }
        };
}

Eigen::MatrixXd Cov_scdata_reader(const BlockReader& Reader, const int Nl, 
    const int C, const int Blockcols)
{
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(Nl);
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(Nl, Nl);
    RowMat Xblk(Nl, std::min(Blockcols, C));
    for (int j0 = 0; j0 < C; j0 += Blockcols)
    {
        const int B = std::min(Blockcols, C - j0);
        Reader(j0, B, Xblk.leftCols(B));
        const auto Xview = Xblk.leftCols(B);
        // Accumulate row sums
        sum.noalias() += Xview.rowwise().sum();
        // Accumulate Gram matrix
        G.noalias() += Xview * Xview.transpose();
    }
    const Eigen::VectorXd mu = sum / double(C);
    G.noalias() -= double(C) * (mu * mu.transpose());
    G /= double(C - 1);
    return G;
}

BlockReader Make_trainreader_selectedrows(const RowMat& Data_perm,
    const FoldRanges& fr, const ZScoreParams& zp, const std::vector<int>& selected_rows)
{
    return [&Data_perm, &fr, &zp, &selected_rows](int col0, int bs, Eigen::Ref<RowMat> out)
        {
            const int Rsel = static_cast<int>(selected_rows.size());
            if (out.rows() != Rsel || out.cols() != bs)
                throw std::runtime_error("Make_trainreader_selectedrows: out size mismatch.");
            for (int i = 0; i < Rsel; ++i) {
                const int r = selected_rows[i];
                const double* src = Data_perm.row(r).data();
                double* dst = out.row(i).data();
                const double mu = zp.mu(r);
                const double inv = zp.invsigma(r);
                if (inv <= 0.0 || !std::isfinite(inv)) {
                    constexpr double nanv = std::numeric_limits<double>::quiet_NaN();
                    for (int j = 0; j < bs; ++j) dst[j] = nanv;
                    continue;
                }
                for (int j = 0; j < bs; ++j) {
                    const int tcol = col0 + j;
                    const int real_col =
                        (tcol < fr.test_begin) ? tcol : (tcol + fr.fold_size);
                    dst[j] = (src[real_col] - mu) * inv;
                }
            }
        };
}

BlockReader Make_testreader_selectedrows(const RowMat& Data_perm,
    const FoldRanges& fr, const ZScoreParams& zp, const std::vector<int>& selected_rows)
{
    return [&Data_perm, &fr, &zp, &selected_rows](int col0, int bs, Eigen::Ref<RowMat> out)
        {
            const int Rsel = static_cast<int>(selected_rows.size());

            if (out.rows() != Rsel || out.cols() != bs)
                throw std::runtime_error("Make_testreader_selectedrows: out size mismatch.");
            for (int i = 0; i < Rsel; ++i) {
                const int r = selected_rows[i];
                const double* src = Data_perm.row(r).data();
                double* dst = out.row(i).data();
                const double mu = zp.mu(r);
                const double inv = zp.invsigma(r);
                if (inv <= 0.0 || !std::isfinite(inv)) {
                    constexpr double nanv = std::numeric_limits<double>::quiet_NaN();
                    for (int j = 0; j < bs; ++j) dst[j] = nanv;
                    continue;
                }
                for (int j = 0; j < bs; ++j) {
                    const int real_col = fr.test_begin + col0 + j;
                    dst[j] = (src[real_col] - mu) * inv;
                }
            }
        };
}

Eigen::MatrixXd HG_out_fold_train(const RowMat& Data_perm, 
    const BlockReader& fullTrainReader,const FoldRanges& fr, const ZScoreParams& train_z, 
    const std::vector<int>& Idx_h, const std::vector<int>& Idx_g, 
    const Eigen::MatrixXd& Cov_train, const int Nh, const int Ng, const int Blockcols)
{
    const int C_train = fr.n_train();
    // Reader for hidden rows from z-scored train data
    int Rsel_h = static_cast<int>(Idx_h.size());
    BlockReader Reader_h = [&](int col0, int bs, Eigen::Ref<RowMat> out) {
        RowMat full_block(Data_perm.rows(), bs);
        fullTrainReader(col0, bs, full_block);
        for (int i = 0; i < Rsel_h; ++i) {
            out.row(i) = full_block.row(Idx_h[i]);
        }
        };
    RowMat Hidden_factor = Topfactor_eigcov_pca(
        Reader_h, Rsel_h, C_train, Nh, Blockcols);
    int Rsel_g = static_cast<int>(Idx_g.size());
    BlockReader Reader_g = [&](int col0, int bs, Eigen::Ref<RowMat> out) {
        RowMat full_block(Data_perm.rows(), bs);
        fullTrainReader(col0, bs, full_block);
        for (int i = 0; i < Rsel_g; ++i) {
            out.row(i) = full_block.row(Idx_g[i]);
        }
        };
    RowMat Global_factor = Topfactor_eigcov_pca(
        Reader_g, Rsel_g, C_train, Ng, Blockcols);
    RowMat HG_factor = Vstack(Hidden_factor, Global_factor);
    Zscore_DataRM(HG_factor);
    const int Nl = static_cast<int>(Data_perm.rows());
    Eigen::MatrixXd S_LHG = Expandcov_lhg_reader(
        Cov_train, fullTrainReader, Nl, C_train, HG_factor, Blockcols);
    return S_LHG;
}

Eigen::MatrixXd HG_out_fold_test(const RowMat& Data_perm,
    const BlockReader& fullTrainReader, const FoldRanges& fr, const ZScoreParams& test_z,
    const std::vector<int>& Idx_h, const std::vector<int>& Idx_g,
    const Eigen::MatrixXd& Cov_test, const int Nh, const int Ng, const int Blockcols)
{
    const int C_test = fr.n_test();
    // Selected-row readers for PCA
    const int Rsel_h = static_cast<int>(Idx_h.size());
    BlockReader Reader_h = Make_testreader_selectedrows(
        Data_perm, fr, test_z, Idx_h);
    RowMat Hidden_factor = Topfactor_eigcov_pca(
        Reader_h, Rsel_h, C_test, Nh, Blockcols);
    const int Rsel_g = static_cast<int>(Idx_g.size());
    BlockReader Reader_g = Make_testreader_selectedrows(
        Data_perm, fr, test_z, Idx_g);
    RowMat Global_factor = Topfactor_eigcov_pca(
        Reader_g, Rsel_g, C_test, Ng, Blockcols);
    RowMat HG_factor = Vstack(Hidden_factor, Global_factor);
    Zscore_DataRM(HG_factor);
    const int Nl = static_cast<int>(Data_perm.rows());
    Eigen::MatrixXd S_LHG = Expandcov_lhg_reader(
        Cov_test, fullTrainReader, Nl, C_test, HG_factor, Blockcols);
    return S_LHG;
}

void LHG_zscore_out_fold_train(const RowMat& Data_perm, const FoldRanges& fr,
    const ZScoreParams& train_z, Input_format& Input_l, const Eigen::MatrixXd& S_LL,
    const std::vector<int>& Idx_h, const std::vector<int>& Idx_g, const int Blockcols)
{
    const int C_train = fr.n_train();
    // Reader for hidden rows from z-scored train data
    int Rsel_h = static_cast<int>(Idx_h.size());
    BlockReader fullTrainReader = Make_trainreader(Data_perm, fr, train_z);
    BlockReader Reader_h = [&](int col0, int bs, Eigen::Ref<RowMat> out) {
        RowMat full_block(Data_perm.rows(), bs);
        fullTrainReader(col0, bs, full_block);
        for (int i = 0; i < Rsel_h; ++i) {
            out.row(i) = full_block.row(Idx_h[i]);
        }
        };
    RowMat Hidden_factor = Topfactor_eigcov_pca(
        Reader_h, Rsel_h, C_train, Input_l.Get_Nhidden(), Blockcols);
    int Rsel_g = static_cast<int>(Idx_g.size());
    BlockReader Reader_g = [&](int col0, int bs, Eigen::Ref<RowMat> out) {
        RowMat full_block(Data_perm.rows(), bs);
        fullTrainReader(col0, bs, full_block);
        for (int i = 0; i < Rsel_g; ++i) {
            out.row(i) = full_block.row(Idx_g[i]);
        }
        };
    RowMat Global_factor = Topfactor_eigcov_pca(
        Reader_g, Rsel_g, C_train, Input_l.Get_Nglobal(), Blockcols);
    RowMat HG_factor = Vstack(Hidden_factor, Global_factor);
    Zscore_DataRM(HG_factor);
    const int Nl = static_cast<int>(Data_perm.rows());
    Eigen::MatrixXd S_LHG = Expandcov_lhg_reader(S_LL, fullTrainReader, Nl, 
        C_train, HG_factor, Blockcols);
    Input_l.Set_Observed_cov(S_LHG);
}

void LHG_zscore_out_fold_test(const RowMat& Data_perm, const FoldRanges& fr,
    const ZScoreParams& test_z, Input_format& Input_l, const Eigen::MatrixXd& S_LL,
    const std::vector<int>& Idx_h, const std::vector<int>& Idx_g, const int Blockcols)
{
    const int C_test = fr.n_test();
    // Full local-test reader for covariance expansion
    BlockReader fullTestReader = Make_testreader(Data_perm, fr, test_z);
    // Selected-row readers for PCA
    const int Rsel_h = static_cast<int>(Idx_h.size());
    BlockReader Reader_h = Make_testreader_selectedrows(
        Data_perm, fr, test_z, Idx_h);
    RowMat Hidden_factor = Topfactor_eigcov_pca(
        Reader_h, Rsel_h, C_test, Input_l.Get_Nhidden(), Blockcols);
    const int Rsel_g = static_cast<int>(Idx_g.size());
    BlockReader Reader_g = Make_testreader_selectedrows(
        Data_perm, fr, test_z, Idx_g);
    RowMat Global_factor = Topfactor_eigcov_pca(
        Reader_g, Rsel_g, C_test, Input_l.Get_Nglobal(), Blockcols);
    RowMat HG_factor = Vstack(Hidden_factor, Global_factor);
    Zscore_DataRM(HG_factor);
    const int Nl = static_cast<int>(Data_perm.rows());
    Eigen::MatrixXd S_LHG = Expandcov_lhg_reader(
        S_LL, fullTestReader, Nl, C_test, HG_factor, Blockcols);
    Input_l.Set_Observed_cov(S_LHG);
}

RowMat Update_hg_factor_reader(const BlockReader& FullReader, const int C, 
    Input_format& Input_l, const Eigen::MatrixXd& S_LL, const int Blockcols)
{
    const int Nhg = Input_l.Get_Nhidden() + Input_l.Get_Nglobal();
    const int Nl = Input_l.Get_Nlocal();
    const int N = Nl + Nhg;
    RowMat E(Nhg, C);
    E.setZero();
    // Precompute inverse diagonal of Theta(hg,hg)
    std::vector<double> invDiag(Nhg, 0.0);
    const auto& Theta = Input_l.Get_theta();
    for (int c = Nl; c < N; ++c) {
        double diag = 0.0;
        for (Eigen::SparseMatrix<double>::InnerIterator it(Theta, c); it; ++it) {
            if (it.row() == c) {
                diag = it.value();
                break;
            }
        }
        if (diag == 0.0)
            throw std::runtime_error("Update_hg_factor_reader: Theta(c,c) is zero.");
        invDiag[c - Nl] = 1.0 / diag;
    }
    // Temporary local-data block
    RowMat Xblk(Nl, std::min(Blockcols, C));
    // Fill E block by block
    for (int j0 = 0; j0 < C; j0 += Blockcols) {
        const int B = std::min(Blockcols, C - j0);
        FullReader(j0, B, Xblk.leftCols(B));
        const auto Xview = Xblk.leftCols(B);
        for (int c = Nl; c < N; ++c) {
            auto Eseg = E.block(c - Nl, j0, 1, B);
            Eseg.setZero();
            const double scale_diag = -invDiag[c - Nl];
            for (Eigen::SparseMatrix<double>::InnerIterator it(Theta, c); it; ++it) {
                const int r = it.row();
                if (r < Nl) {
                    const double val = it.value();
                    Eseg.noalias() += (scale_diag * val) * Xview.row(r);
                }
            }
        }
    }
    RowMat E_zscore = E;
    Zscore_DataRM(E_zscore);
    Eigen::MatrixXd S_LHG = Expandcov_lhg_reader(S_LL, FullReader, Nl, C, E_zscore, Blockcols);
    Input_l.Set_Observed_cov(S_LHG);
    return E_zscore;
}

RowMat Solution_theta_and_hg_factor_fold(const BlockReader& FullReader,const int C, 
    Input_format& Input_l, const Eigen::MatrixXd& S_LL, const int Max_loop,
    const double Stop_threshold, const int Blockcol_inv, const int Blockcol_data)
{
    RowMat  Updated_hg;
    Eigen::SparseMatrix<double> Old_theta = Input_l.Get_theta();
    for (int i = 0; i < Max_loop; ++i) {
        Newton_method_offdiag_lhg(Input_l, Blockcol_inv);
        Updated_hg = Update_hg_factor_reader(FullReader, C, Input_l, S_LL, Blockcol_data);
        bool If_break = SparseMatricesClose(Old_theta, Input_l.Get_theta(), Stop_threshold);
        if (If_break) {
            //std::cout << "Loop is completed:" << '\t' << i + 1 << '\t' << "times;" << std::endl;
            break;
        }
        Old_theta = Input_l.Get_theta();
    }
    return Updated_hg;
}
