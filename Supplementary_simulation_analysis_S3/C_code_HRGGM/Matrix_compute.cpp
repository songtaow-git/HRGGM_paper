#include "Matrix_compute.h"

double Logdet_sparse(const Eigen::SparseMatrix<double>& X)
{
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol;
    chol.compute(X);
    if (chol.info() != Eigen::Success)
        throw std::runtime_error("Cholesky failed");
    // log(det(X)) = sum(log(diag(L)))
    const Eigen::SparseMatrix<double>& L = chol.matrixL();
    double logdet = 2 * L.diagonal().array().log().sum();
    return logdet;
}

double Logdet_sparse(const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& chol)
{
    if (chol.info() != Eigen::Success)
        throw std::runtime_error("Cholesky failed");
    // log(det(X)) = sum(log(diag(L)))
    const Eigen::SparseMatrix<double>& L = chol.matrixL();
    double logdet = 2 * L.diagonal().array().log().sum();
    return logdet;
}

double Logdet_dense(const Eigen::MatrixXd& X)
{
    if (X.rows() != X.cols())
        throw std::runtime_error("Logdet_dense_spd: X must be square.");
    // LLT assumes symmetric positive definite.
    Eigen::LLT<Eigen::MatrixXd> llt;
    llt.compute(X);
    if (llt.info() != Eigen::Success)
        throw std::runtime_error("Logdet_dense_spd: LLT failed (matrix not SPD or numerical issue).");
    const Eigen::MatrixXd& L = llt.matrixL();
    // log(det(X)) = 2 * sum(log(diag(L)))
    const double logdet = 2.0 * L.diagonal().array().log().sum();
    return logdet;
}

double Trace_dense_sparse(const Eigen::MatrixXd& S, const Eigen::SparseMatrix<double>& X)
{
    // trace(S*X) = sum_{r,c} S(c,r) * X(r,c)
    double tr = 0.0;
    for (int k = 0; k < X.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(X, k); it; ++it) {
            tr += S(it.col(), it.row()) * it.value();
        }
    }
    return tr;
}

double Trace_dense_dense(const Eigen::MatrixXd& S, const Eigen::MatrixXd& X)
{
    if (S.rows() != X.rows() || S.cols() != X.cols())
        throw std::runtime_error("Trace_dense_dense: size mismatch.");
    // trace(S*X) = sum(S^T .* X)
    return (S.transpose().cwiseProduct(X)).sum();
}

double Likelihood_f(const Eigen::MatrixXd& S,
    const Eigen::SparseMatrix<double>& X)
{
    double Likelihood_val = Trace_dense_sparse(S, X) - Logdet_sparse(X);
    return Likelihood_val;
}

double Likelihood_f(const Eigen::MatrixXd& S, const Eigen::SparseMatrix<double>& X
    , const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& chol)
{
    double Likelihood_val = Trace_dense_sparse(S, X) - Logdet_sparse(chol);
    return Likelihood_val;
}

double Offdiag_l1_norm_ll(const Eigen::SparseMatrix<double>& X, const int Nl)
{
    double sum = 0.0;
    // Eigen::SparseMatrix 
    // outerSize() == N_col
    for (int j = 0; j < X.outerSize(); ++j) {
        if (j >= Nl) break;
        for (Eigen::SparseMatrix<double>::InnerIterator it(X, j); it; ++it) {
            const int i = it.row();
            // off-diag
            if (i < Nl && i != j) {
                sum += std::abs(it.value());
            }
        }
    }
    return sum;
}

double Block_l1_norm_lh(const Eigen::SparseMatrix<double>& X, const int Nl, const int Nh)
{
    double sum_abs = 0.0;
    const int col_start = Nl;
    const int col_end = Nl + Nh; // exclusive
    // Eigen::SparseMatrix is column-major by default → loop over columns
    for (int j = col_start; j < col_end; ++j) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(X, j); it; ++it) {
            const int i = it.row();
            if (i < Nl) {
                sum_abs += std::abs(it.value());
            }
        }
    }
    return 2.0 * sum_abs;
}

double Block_l2_norm_lh(const Eigen::SparseMatrix<double>& X, const int Nl, const int Nh)
{
    double sum_sq = 0.0;
    const int col_start = Nl;
    const int col_end = Nl + Nh;   // exclusive
    // SparseMatrix is column-major by default → iterate columns
    for (int j = col_start; j < col_end; ++j) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(X, j); it; ++it) {
            const int i = it.row();
            if (i < Nl) {
                const double v = it.value();
                sum_sq += v * v;
            }
        }
    }
    return 2.0 * sum_sq;
}

double Block_l2_norm_lg(const Eigen::SparseMatrix<double>& X, 
    const int Nl, const int Nh, const int Ng)
{
    double sum_sq = 0.0;
    int col_start = Nl + Nh;
    int col_end = Nl + Nh + Ng;  // exclusive
    for (int k = col_start; k < col_end; ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(X, k); it; ++it) {
            int i = it.row();
            if (i < Nl) {
                double v = it.value();
                sum_sq += v * v;
            }
        }
    }
    return 2.0 * sum_sq;
}

RowMat Inverse_spd_sparse(const Eigen::SparseMatrix<double>& X,
    const int blockCols)
{
    const int n = static_cast<int>(X.rows());
    if (X.cols() != n) {
        throw std::runtime_error("inverse_spd_sparse_to_dense: X must be square.");
    }
    if (blockCols <= 0) {
        throw std::runtime_error("inverse_spd_sparse_to_dense: blockCols must be > 0.");
    }
    // SPD sparse Cholesky
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> chol;
    chol.analyzePattern(X);
    chol.factorize(X);
    if (chol.info() != Eigen::Success) {
        throw std::runtime_error("inverse_spd_sparse_to_dense: Cholesky failed (X not SPD or numerical issue).");
    }
    // Output dense inverse
    RowMat Xinv(n, n);
    // Blocked RHS and solution
    Eigen::MatrixXd RHS(n, blockCols);
    Eigen::MatrixXd Y(n, blockCols);
    for (int j0 = 0; j0 < n; j0 += blockCols) {
        const int k = std::min(blockCols, n - j0);

        RHS.setZero(n, k);
        RHS.block(j0, 0, k, k).setIdentity();   // columns j0..j0+k-1 of I

        Y.noalias() = chol.solve(RHS);
        if (chol.info() != Eigen::Success) {
            throw std::runtime_error("inverse_spd_sparse_to_dense: solve failed.");
        }

        Xinv.block(0, j0, n, k) = Y.leftCols(k);
    }
    return Xinv;
}

RowMat Inverse_spd_sparse_from_chol(
    const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& chol,
    const int n, const int blockCols)
{
    if (n <= 0) {
        throw std::runtime_error("Inverse_spd_from_chol: n must be > 0.");
    }
    if (blockCols <= 0) {
        throw std::runtime_error("Inverse_spd_from_chol: blockCols must be > 0.");
    }
    if (chol.info() != Eigen::Success) {
        throw std::runtime_error("Inverse_spd_from_chol: invalid factorization (chol.info != Success).");
    }
    RowMat Xinv(n, n);
    // Allocate once; conservativeResize inside loop to avoid reallocation of capacity
    Eigen::MatrixXd RHS(n, blockCols);
    Eigen::MatrixXd Y(n, blockCols);
    for (int j0 = 0; j0 < n; j0 += blockCols) {
        const int k = std::min(blockCols, n - j0);
        // Resize views without reallocating too often
        RHS.conservativeResize(n, k);
        Y.conservativeResize(n, k);
        RHS.setZero();
        RHS.block(j0, 0, k, k).setIdentity(); // columns j0..j0+k-1 of I
        // Solve: X * Y = RHS  =>  Y = X^{-1} * RHS
        Y.noalias() = chol.solve(RHS);
        // Copy into output
        Xinv.block(0, j0, n, k) = Y;
    }
    return Xinv;
}

std::vector<std::pair<int, int>> Find_update_entries_ll(
    const Eigen::SparseMatrix<double>& X_t, const Eigen::MatrixXd& Nabla_g,
    const double Lambda, const int Nl)
{
    const int B = Nl;
    if (B <= 1) {
        return std::vector<std::pair<int, int>>();
    }
    std::vector<std::pair<int, int>> Updates;
    const size_t cap_tri = (static_cast<size_t>(B) * static_cast<size_t>(B - 1)) / 2u;
    const size_t cap_nnz = static_cast<size_t>(X_t.nonZeros());
    Updates.reserve(std::min(cap_tri, cap_nnz) + 64u);
    // Mark for the full BxB (simple & fast); MarkAndPush uses col-major indexing
    std::vector<uint8_t> Mark(static_cast<size_t>(B) * static_cast<size_t>(B), 0);
    MarkAndPush Mark_and_push(B, Mark, Updates);
    // (1) X_t != 0 inside local block, but keep only upper triangle (row < col)
    for (int col = 0; col < X_t.outerSize() && col < B; ++col) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(X_t, col); it; ++it) {
            const int row = it.row();
            if (row < col && row < B) {
                Mark_and_push(row, col);
            }
        }
    }
    // (2) abs(Nabla_g) > Lambda inside local block, upper triangle only (i < j)
    for (int j = 1; j < B; ++j) {
        for (int i = 0; i < j; ++i) {
            if (std::abs(Nabla_g(i, j)) > Lambda) {
                Mark_and_push(i, j);
            }
        }
    }
    return Updates; 
}

std::vector<std::pair<int, int>> Find_update_entries_lh(
    const Eigen::SparseMatrix<double>& X_t, const Eigen::MatrixXd& Nabla_lh,
    const double Lambda, const double Gamma, const int Nl, const int Nh)
{
    const double thr = Gamma * Lambda;
    std::vector<std::pair<int, int>> Updates;
    const size_t cap_block = static_cast<size_t>(Nl) * static_cast<size_t>(Nh);
    const size_t cap_nnz = static_cast<size_t>(X_t.nonZeros());
    Updates.reserve(std::min(cap_block, cap_nnz) + 64u);
    // Mark on (Nl x Nh) block using col-major indexing, stride B_=Nl
    std::vector<uint8_t> Mark(static_cast<size_t>(Nl) * static_cast<size_t>(Nh), 0);
    MarkAndPush Mark_and_push(Nl, Mark, Updates);
    int colH0 = Nl;
    int colH1 = Nh + Nl; // exclusive
    // (1) X_t != 0 in (L,H): scan only H columns (Nh columns)
    for (int col = colH0; col < X_t.outerSize() && col < colH1; ++col) {
        const int j_off = col - Nl; // [0, Nh)
        for (Eigen::SparseMatrix<double>::InnerIterator it(X_t, col); it; ++it) {
            const int row = it.row();
            if (row >= 0 && row < Nl) {
                // store as (rowL, colHoffset) first
                Mark_and_push(row, j_off);
            }
        }
    }
    colH0 = 0;
    colH1 = Nh; // exclusive
    // (2) abs(Nabla_g) > Gamma*Lambda in (L,H):
    for (int col = colH0; col < colH1; ++col) {
        for (int row = 0; row < Nl; ++row) {
            if (std::abs(Nabla_lh(row, col)) > thr) {
                Mark_and_push(row, col);
            }
        }
    }
    const size_t m = Updates.size();
    for (size_t k = 0; k < m; ++k) {
        Updates[k].second += Nl; // colH = Nl + colHoffset
    }

    return Updates;
}

void Find_update_entries_lh_lg_diag(
    const std::vector<std::pair<int, int>>& pairs_base,
    const int Nl, const int Nh, const int Ng,
    std::vector<std::pair<int, int>>& pairs_lg,
    std::vector<std::pair<int, int>>& pairs_diag,
    std::vector<std::pair<int, int>>& pairs_all
)
{
    if (Nl < 0 || Nh < 0 || Ng < 0) {
        throw std::runtime_error("Build_lg_diag_all_pairs: Nl/Nh/Ng must be >= 0.");
    }
    const int N = Nl + Nh + Ng;
    if (N <= 0) {
        pairs_lg.clear();
        pairs_diag.clear();
        pairs_all.clear();
        return;
    }
    const int colG0 = Nl + Nh;
    const int colG1 = Nl + Nh + Ng; // exclusive
    const std::size_t lg_sz =
        static_cast<std::size_t>(Nl) * static_cast<std::size_t>(Ng);
    const std::size_t diag_sz =
        static_cast<std::size_t>(N);
    // ---- build lg
    pairs_lg.clear();
    pairs_lg.reserve(lg_sz);
    for (int col = colG0; col < colG1; ++col) {
        for (int row = 0; row < Nl; ++row) {
            pairs_lg.emplace_back(row, col);
        }
    }
    // ---- build diag
    pairs_diag.clear();
    pairs_diag.reserve(diag_sz);
    for (int d = 0; d < N; ++d) {
        pairs_diag.emplace_back(d, d);
    }
    // ---- build all = base + lg + diag
    pairs_all.clear();
    pairs_all.reserve(pairs_base.size() + pairs_lg.size() + pairs_diag.size());
    pairs_all.insert(pairs_all.end(), pairs_base.begin(), pairs_base.end());
    pairs_all.insert(pairs_all.end(), pairs_lg.begin(), pairs_lg.end());
    pairs_all.insert(pairs_all.end(), pairs_diag.begin(), pairs_diag.end());
}

void Find_update_entries_ll_diag(
    const std::vector<std::pair<int, int>>& pairs_base,
    const int Nl, std::vector<std::pair<int, int>>& pairs_diag,
    std::vector<std::pair<int, int>>& pairs_all)
{
    
    const int N = Nl;
    if (N <= 0) {
        pairs_diag.clear();
        pairs_all.clear();
        return;
    }
    const std::size_t diag_sz = static_cast<std::size_t>(N);
    // ---- build diag
    pairs_diag.clear();
    pairs_diag.reserve(diag_sz);
    for (int d = 0; d < N; ++d) {
        pairs_diag.emplace_back(d, d);
    }
    // ---- build all = base + diag
    pairs_all.clear();
    pairs_all.reserve(pairs_base.size() + pairs_diag.size());
    pairs_all.insert(pairs_all.end(), pairs_base.begin(), pairs_base.end());
    pairs_all.insert(pairs_all.end(), pairs_diag.begin(), pairs_diag.end());
}

static int Find_value_index_in_col(const Eigen::SparseMatrix<double>& A,
    const int row, const int col) 
{
    const int* outer = A.outerIndexPtr();
    const int* inner = A.innerIndexPtr();
    const int begin = outer[col];
    const int end = outer[col + 1];
    const int* it = std::lower_bound(inner + begin, inner + end, row);
    if (it == inner + end || *it != row) {
        return -1;
    }
    return static_cast<int>(it - inner);
}

void D_index_pattern(Eigen::SparseMatrix<double>& D_t, const int n,
    const std::vector<std::pair<int, int>>& new_pairs, FixedDindex& out)
{
    if (n <= 0) {
        throw std::runtime_error("RebuildPatternAndIndexDT: n must be > 0.");
    }

    if (D_t.rows() != n || D_t.cols() != n) {
        D_t.resize(n, n);
    }
    std::vector<Eigen::Triplet<double> > trips;
    trips.reserve(new_pairs.size() * 2u);
    for (std::size_t k = 0; k < new_pairs.size(); ++k) {
        const int i = new_pairs[k].first;
        const int j = new_pairs[k].second;
        if (i < 0 || j < 0 || i >= n || j >= n) {
            throw std::runtime_error("RebuildPatternAndIndexDT: pair index out of range.");
        }

        trips.push_back(Eigen::Triplet<double>(i, j, 0.0));
        if (i != j) {
            trips.push_back(Eigen::Triplet<double>(j, i, 0.0));
        }
    }
    D_t.setFromTriplets(trips.begin(), trips.end());
    D_t.makeCompressed();
    // Fill out mapping
    out.D = &D_t;
    out.n = n;
    out.pairs = new_pairs;
    out.idx_ij.resize(new_pairs.size());
    out.idx_ji.resize(new_pairs.size());
    // Precompute indices into value array
    for (std::size_t k = 0; k < new_pairs.size(); ++k) {
        const int i = new_pairs[k].first;
        const int j = new_pairs[k].second;
        int idx1 = Find_value_index_in_col(D_t, i, j);
        if (idx1 < 0) {
            throw std::runtime_error("RebuildPatternAndIndexDT: (i,j) not found in D_t pattern.");
        }
        out.idx_ij[k] = idx1;
        int idx2 = Find_value_index_in_col(D_t, j, i);
        if (idx2 < 0) {
            throw std::runtime_error("RebuildPatternAndIndexDT: (j,i) not found in D_t pattern.");
        }
        out.idx_ji[k] = idx2;
    }
}

void Update_mu_ll_offdiag(
    const Eigen::SparseMatrix<double>& X_t,
    const int i, const int j, const std::size_t k,
    const double Lambda,
    FixedDindex& Didx,      
    const RowMat& W_t,
    const RowMat& W_T,
    RowMat& U_T,
    const Eigen::MatrixXd& S)
{
    Eigen::SparseMatrix<double>& D_t = *Didx.D;
    double* dvals = D_t.valuePtr();
    const int idx_ij = Didx.idx_ij[k];
    const int idx_ji = Didx.idx_ji[k];
    // a
    const double Wij = W_t(i, j);
    double a = Wij * Wij;
    if (i != j) {
        a += W_t(i, i) * W_t(j, j);
    }
    // b = S(i,j) - W(i,j) + (W(:,i))' * U(:,j)
    //const double b = S(i, j) - Wij + W_t.col(i).dot(U_t.col(j));
    const double b = S(i, j) - Wij + W_T.row(i).dot(U_T.row(j));
    double mu_ij = 0.0;
    if (i != j) {
        // c = X_t(i,j) + D_t(i,j)  (D_t(i,j) read in O(1))
        const double Dij = dvals[idx_ij];
        const double c = X_t.coeff(i, j) + Dij;
        mu_ij = -c;
        const double ba = b / a;
        const double la = Lambda / a;
        const double t = c - ba;
        if (std::abs(t) > std::abs(la)) {
            if (t > 0.0) {
                mu_ij = -(b + Lambda) / a;
            }
            else {
                mu_ij = -(b - Lambda) / a;
            }
        }
    }
    else {
        // i == j
        mu_ij = -b / a;
    }
    // D_t(i,j) += mu_ij; D_t(j,i) += mu_ij (symmetric) -- O(1) updates
    if (i != j) {
        dvals[idx_ij] += mu_ij;
        dvals[idx_ji] += mu_ij;
    }
    else {
        dvals[idx_ij] += mu_ij;
    }
    // U_t(i,:) += mu_ij * W_t(j,:)
    //U_t.row(i).noalias() += mu_ij * W_t.row(j);
    add_to_UT_col_from_Wrow(U_T, i, W_t, j, mu_ij);
    // if i != j: U_t(j,:) += mu_ij * W_t(i,:)
    if (i != j) {
        //U_t.row(j).noalias() += mu_ij * W_t.row(i);
        add_to_UT_col_from_Wrow(U_T, j, W_t, i, mu_ij);
    }
}


void Update_mu_lh(
    const Eigen::SparseMatrix<double>& X_t,
    const int i, const int j, const std::size_t k,
    const double Lambda, const double Alpha, const double Gamma,
    FixedDindex& Didx,
    const RowMat& W_t,
    const RowMat& W_T,
    RowMat& U_T,
    const Eigen::MatrixXd& S)
{
    Eigen::SparseMatrix<double>& D_t = *Didx.D;
    double* dvals = D_t.valuePtr();
    const int idx_ij = Didx.idx_ij[k];
    const int idx_ji = Didx.idx_ji[k];
    const double Dij = dvals[idx_ij];
    const double Wij = W_t(i, j);
    const double a = Wij * Wij + W_t(i, i) * W_t(j, j) + 2.0 * Alpha * (1.0 - Gamma);
    const double c = X_t.coeff(i, j) + Dij;
    //const double dotWi_Uj = W_t.col(i).dot(U_t.col(j));
    const double dotWi_Uj = W_T.row(i).dot(U_T.row(j));
    const double b = S(i, j) - Wij + dotWi_Uj + 2.0 * Alpha * (1.0 - Gamma) * c;
    double mu_ij = -c;
    const double ba = b / a;
    const double thr = (Gamma * Lambda) / a;
    if (std::abs(c - ba) > std::abs(thr)) {
        if ((c - ba) > 0.0) {
            mu_ij = -(b + Gamma * Lambda) / a;
        }
        else {
            mu_ij = -(b - Gamma * Lambda) / a;
        }
    }
    dvals[idx_ij] += mu_ij;
    dvals[idx_ji] += mu_ij;
    //U_t.row(i).noalias() += mu_ij * W_t.row(j);
    //U_t.row(j).noalias() += mu_ij * W_t.row(i);
    add_to_UT_col_from_Wrow(U_T, i, W_t, j, mu_ij);
    add_to_UT_col_from_Wrow(U_T, j, W_t, i, mu_ij);
}

void Update_mu_lg(
    const Eigen::SparseMatrix<double>& X_t,
    const int i, const int j, const std::size_t k,
    const double Alpha, FixedDindex& Didx,
    const RowMat& W_t,
    const RowMat& W_T,
    RowMat& U_T,
    const Eigen::MatrixXd& S)
{
    Eigen::SparseMatrix<double>& D_t = *Didx.D;
    double* dvals = D_t.valuePtr();
    const int idx_ij = Didx.idx_ij[k];
    const int idx_ji = Didx.idx_ji[k];
    // a = W(i,j)^2 + W(i,i)*W(j,j) + 2*Alpha
    const double Wij = W_t(i, j);
    const double a = Wij * Wij + W_t(i, i) * W_t(j, j) + 2.0 * Alpha;
    // c = X(i,j) + D(i,j)
    const double c = X_t.coeff(i, j) + dvals[idx_ij];
    // b = S(i,j) - W(i,j) + W(:,i)' * U(:,j) + 2*Alpha*c
    const double b = S(i, j)
        - Wij
        //+ W_t.col(i).dot(U_t.col(j))
        + W_T.row(i).dot(U_T.row(j))
        + 2.0 * Alpha * c;
    // mu_ij = -b / a
    const double mu_ij = -b / a;
    // D_t(i,j) += mu_ij;
    // D_t(j,i) += mu_ij;
    dvals[idx_ij] += mu_ij;
    dvals[idx_ji] += mu_ij;
    // U_t(i,:) += mu_ij * W_t(j,:);
    //U_t.row(i).noalias() += mu_ij * W_t.row(j);
    // U_t(j,:) += mu_ij * W_t(i,:);
    //U_t.row(j).noalias() += mu_ij * W_t.row(i);
    add_to_UT_col_from_Wrow(U_T, i, W_t, j, mu_ij);
    add_to_UT_col_from_Wrow(U_T, j, W_t, i, mu_ij);
}

void Update_mu_diag(
    const int i, const std::size_t k,
    FixedDindex& Didx,
    const RowMat& W_t,
    const RowMat& W_T,
    RowMat& U_T,
    const Eigen::MatrixXd& S)
{
    Eigen::SparseMatrix<double>& D_t = *Didx.D;
    double* dvals = D_t.valuePtr();
    // For diagonal update, idx_ij[k] should correspond to (i,i)
    const int idx_ii = Didx.idx_ij[k];
    // a = W(i,i)^2
    const double Wii = W_t(i, i);
    const double a = Wii * Wii;
    // b = S(i,i) - W(i,i) + W(:,i)' * U(:,i)
    //const double b = S(i, i) - Wii + W_t.col(i).dot(U_t.col(i));
    const double b = S(i, i) - Wii + W_T.row(i).dot(U_T.row(i));
    // mu_ii = -b / a
    const double mu_ii = -b / a;
    // D_t(i,i) += mu_ii  (O(1) sparse update)
    dvals[idx_ii] += mu_ii;
    // U_t(i,:) += mu_ii * W_t(i,:)
   // U_t.row(i).noalias() += mu_ii * W_t.row(i);
    add_to_UT_col_from_Wrow(U_T, i, W_t, i, mu_ii);
}

double sgn(const double x) 
{
    // caller guarantees x != 0
    return (x > 0.0) ? 1.0 : -1.0;
}

std::size_t tri_index_upper(const int i, const int j, const int NL) 
{
    const std::size_t base =
        static_cast<std::size_t>(i) *
        static_cast<std::size_t>(2 * NL - i - 1) / 2;
    return base + static_cast<std::size_t>(j - i - 1);
}

KKT_residual KKTResidual_lhg_log(const Eigen::SparseMatrix<double>& Theta,
    const Eigen::MatrixXd& S, const RowMat& W,
    const int Nl, const int Nh, const int Ng,
    const double lambda, const double alpha, const double gamma)
{
    if (gamma < 0.0 || gamma > 1.0) {
        throw std::invalid_argument("gamma must be in [0,1].");
    }
    const int Size = Nl + Nh + Ng;
    if (S.rows() != Size || S.cols() != Size) {
        throw std::invalid_argument("S dimension mismatch.");
    }
    if (W.rows() != Size || W.cols() != Size) {
        throw std::invalid_argument("W dimension mismatch.");
    }
    if (Theta.rows() != Size || Theta.cols() != Size) {
        throw std::invalid_argument("Theta dimension mismatch.");
    }
    if (Nl <= 0 || Nh < 0 || Ng < 0 || Nl + Nh + Ng <= 0) {
        throw std::invalid_argument("Invalid block sizes.");
    }
    if (Nh + Ng > 5) {
        // not required, but you mentioned this regime; remove if undesired
    }
    const int L0 = 0;
    const int H0 = Nl;
    const int G0 = Nl + Nh;
    const double tau = gamma * lambda;                 // L1 on LH
    const double c2 = 2.0 * (1.0 - gamma) * alpha;    // L2 grad multiplier on LH
    const double cLG = 2.0 * alpha;                    // L2 grad multiplier on LG

    // -------------------- 1) Diagonal residual (all diagonal entries are free) --------------------
    double r_diag = 0.0;
    for (int i = 0; i < Size; ++i) {
        const double grad = S(i, i) - W(i, i);
        r_diag = std::max(r_diag, std::abs(grad));
    }
    // -------------------- 2) LL off-diagonal (upper triangle): exact active/inactive KKT --------------------
    const std::size_t tri_size = static_cast<std::size_t>(Nl) * static_cast<std::size_t>(Nl - 1) / 2;
    std::vector<uint8_t> ll_active(tri_size, uint8_t{ 0 });
    std::vector<double>  ll_value(tri_size, 0.0); // store Theta_ij for active entries only
    ll_value.clear(); ll_value.shrink_to_fit();
    double r_LL_nz = 0.0;
    // Iterate Theta nonzeros once: mark LL active (upper-triangle unique) and compute active residual.
    for (int k = 0; k < Theta.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Theta, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            const double val = it.value();
            if (val == 0.0) continue;
            if (i == j) continue;
            // LL block only
            if (i < Nl && j < Nl) {
                int a = i, b = j;
                if (a > b) std::swap(a, b); // ensure a<b
                const std::size_t idx = tri_index_upper(a, b, Nl);
                if (ll_active[idx]) continue; // already processed this pair
                ll_active[idx] = 1;
                const double grad = S(a, b) - W(a, b);
                // active KKT: |grad + lambda*sign(theta)| (theta is val with correct sign)
                r_LL_nz = std::max(r_LL_nz, std::abs(grad + lambda * sgn(val)));
            }
        }
    }
    // Inactive LL: scan all i<j in LL and consider those NOT active.
    double r_LL_z = 0.0;
    for (int i = 0; i < Nl; ++i) {
        const std::size_t base = static_cast<std::size_t>(i) *
            static_cast<std::size_t>(2 * Nl - i - 1) / 2;
        for (int j = i + 1; j < Nl; ++j) {
            const std::size_t idx = base + static_cast<std::size_t>(j - i - 1);
            if (ll_active[idx]) continue;
            const double grad = S(i, j) - W(i, j);
            const double viol = std::abs(grad) - lambda;
            if (viol > r_LL_z) r_LL_z = viol;
        }
    }
    if (r_LL_z < 0.0) r_LL_z = 0.0;
    // -------------------- 3) LH block (L x H): elastic-net, H small --------------------
    double r_LH_nz = 0.0, r_LH_z = 0.0;
    if (Nh > 0) {
        for (int i = 0; i < Nl; ++i) {
            for (int j = 0; j < Nh; ++j) {
                const int col = H0 + j; // global index in Theta
                const double grad = S(i, col) - W(i, col);

                const double th = Theta.coeff(i, col);
                if (th != 0.0) {
                    const double res = std::abs(grad + c2 * th + tau * sgn(th));
                    if (res > r_LH_nz) r_LH_nz = res;
                }
                else {
                    const double viol = std::abs(grad) - tau;
                    if (viol > r_LH_z) r_LH_z = viol;
                }
            }
        }
        if (r_LH_z < 0.0) r_LH_z = 0.0;
    }
    // -------------------- 4) LG block (L x G): L2 only, G small --------------------
    double r_LG = 0.0;
    if (Ng > 0) {
        for (int i = 0; i < Nl; ++i) {
            for (int j = 0; j < Ng; ++j) {
                const int col = G0 + j;
                const double grad = S(i, col) - W(i, col);
                const double th = Theta.coeff(i, col); 
                const double res = std::abs(grad + cLG * th);
                if (res > r_LG) r_LG = res;
            }
        }
    }
    KKT_residual KKT_out(r_diag, r_LL_nz, r_LL_z, r_LH_nz, r_LH_z, r_LG);
    return KKT_out;
}

bool KKTResidual_lhg_convergence(const Eigen::SparseMatrix<double>& Theta,
    const Eigen::MatrixXd& S, const RowMat& W,
    const int Nl, const int Nh, const int Ng, const double Scale,
    const double lambda, const double alpha, const double gamma)
{
    if (gamma < 0.0 || gamma > 1.0) {
        throw std::invalid_argument("gamma must be in [0,1].");
    }
    const int Size = Nl + Nh + Ng;
    if (S.rows() != Size || S.cols() != Size) {
        throw std::invalid_argument("S dimension mismatch.");
    }
    if (W.rows() != Size || W.cols() != Size) {
        throw std::invalid_argument("W dimension mismatch.");
    }
    if (Theta.rows() != Size || Theta.cols() != Size) {
        throw std::invalid_argument("Theta dimension mismatch.");
    }
    if (Nl <= 0 || Nh < 0 || Ng < 0 || Nl + Nh + Ng <= 0) {
        throw std::invalid_argument("Invalid block sizes.");
    }
    if (Nh + Ng > 5) {
        // not required, but you mentioned this regime; remove if undesired
    }
    const int L0 = 0;
    const int H0 = Nl;
    const int G0 = Nl + Nh;
    const double tau = gamma * lambda;                 // L1 on LH
    const double c2 = 2.0 * (1.0 - gamma) * alpha;    // L2 grad multiplier on LH
    const double cLG = 2.0 * alpha;                    // L2 grad multiplier on LG

    // -------------------- 1) Diagonal residual (all diagonal entries are free) --------------------
    for (int i = 0; i < Size; ++i) {
        const double grad = S(i, i) - W(i, i);
        if (std::abs(grad) > Scale) return false;
    }
    // -------------------- 2) LL off-diagonal (upper triangle): exact active/inactive KKT --------------------
    const std::size_t tri_size = static_cast<std::size_t>(Nl) * static_cast<std::size_t>(Nl - 1) / 2;
    std::vector<uint8_t> ll_active(tri_size, uint8_t{ 0 });
    std::vector<double>  ll_value(tri_size, 0.0); // store Theta_ij for active entries only
    ll_value.clear(); ll_value.shrink_to_fit();
    // Iterate Theta nonzeros once: mark LL active (upper-triangle unique) and compute active residual.
    for (int k = 0; k < Theta.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Theta, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            const double val = it.value();
            if (val == 0.0) continue;
            if (i == j) continue;
            // LL block only
            if (i < Nl && j < Nl) {
                int a = i, b = j;
                if (a > b) std::swap(a, b); // ensure a<b
                const std::size_t idx = tri_index_upper(a, b, Nl);
                if (ll_active[idx]) continue; // already processed this pair
                ll_active[idx] = 1;
                const double grad = S(a, b) - W(a, b);
                // active KKT: |grad + lambda*sign(theta)| (theta is val with correct sign)
                const double res = std::abs(grad + lambda * sgn(val));
                if (res > Scale) return false;
            }
        }
    }
    // Inactive LL: scan all i<j in LL and consider those NOT active.
    for (int i = 0; i < Nl; ++i) {
        const std::size_t base = static_cast<std::size_t>(i) *
            static_cast<std::size_t>(2 * Nl - i - 1) / 2;
        for (int j = i + 1; j < Nl; ++j) {
            const std::size_t idx = base + static_cast<std::size_t>(j - i - 1);
            if (ll_active[idx]) continue;
            const double grad = S(i, j) - W(i, j);
            const double viol = std::abs(grad) - lambda;
            if (viol > Scale) return false;
        }
    }
    // -------------------- 3) LH block (L x H): elastic-net, H small --------------------
    if (Nh > 0) {
        for (int i = 0; i < Nl; ++i) {
            for (int j = 0; j < Nh; ++j) {
                const int col = H0 + j; // global index in Theta
                const double grad = S(i, col) - W(i, col);

                const double th = Theta.coeff(i, col);
                if (th != 0.0) {
                    const double res = std::abs(grad + c2 * th + tau * sgn(th));
                    if (res > Scale) return false;
                }
                else {
                    const double viol = std::abs(grad) - tau;
                    if (viol > Scale) return false;
                }
            }
        }
    }
    // -------------------- 4) LG block (L x G): L2 only, G small --------------------
    if (Ng > 0) {
        for (int i = 0; i < Nl; ++i) {
            for (int j = 0; j < Ng; ++j) {
                const int col = G0 + j;
                const double grad = S(i, col) - W(i, col);
                const double th = Theta.coeff(i, col);
                const double res = std::abs(grad + cLG * th);
                if (res > Scale) return false;
            }
        }
    }
    return true;
}

bool KKTResidual_l_convergence(const Eigen::SparseMatrix<double>& Theta,
    const Eigen::MatrixXd& S, const RowMat& W, const int Nl, 
    const double Scale, const double lambda)
{
    const int Size = Nl;
    if (S.rows() != Size || S.cols() != Size) {
        throw std::invalid_argument("S dimension mismatch.");
    }
    if (W.rows() != Size || W.cols() != Size) {
        throw std::invalid_argument("W dimension mismatch.");
    }
    if (Theta.rows() != Size || Theta.cols() != Size) {
        throw std::invalid_argument("Theta dimension mismatch.");
    }
    const int L0 = 0;
    const int H0 = Nl;

    // -------------------- 1) Diagonal residual (all diagonal entries are free) --------------------
    for (int i = 0; i < Size; ++i) {
        const double grad = S(i, i) - W(i, i);
        if (std::abs(grad) > Scale) return false;
    }
    // -------------------- 2) LL off-diagonal (upper triangle): exact active/inactive KKT --------------------
    const std::size_t tri_size = static_cast<std::size_t>(Nl) * static_cast<std::size_t>(Nl - 1) / 2;
    std::vector<uint8_t> ll_active(tri_size, uint8_t{ 0 });
    std::vector<double>  ll_value(tri_size, 0.0); // store Theta_ij for active entries only
    ll_value.clear(); ll_value.shrink_to_fit();
    // Iterate Theta nonzeros once: mark LL active (upper-triangle unique) and compute active residual.
    for (int k = 0; k < Theta.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Theta, k); it; ++it) {
            int i = it.row();
            int j = it.col();
            const double val = it.value();
            if (val == 0.0) continue;
            if (i == j) continue;
            // LL block only
            if (i < Nl && j < Nl) {
                int a = i, b = j;
                if (a > b) std::swap(a, b); // ensure a<b
                const std::size_t idx = tri_index_upper(a, b, Nl);
                if (ll_active[idx]) continue; // already processed this pair
                ll_active[idx] = 1;
                const double grad = S(a, b) - W(a, b);
                // active KKT: |grad + lambda*sign(theta)| (theta is val with correct sign)
                const double res = std::abs(grad + lambda * sgn(val));
                if (res > Scale) return false;
            }
        }
    }
    // Inactive LL: scan all i<j in LL and consider those NOT active.
    for (int i = 0; i < Nl; ++i) {
        const std::size_t base = static_cast<std::size_t>(i) *
            static_cast<std::size_t>(2 * Nl - i - 1) / 2;
        for (int j = i + 1; j < Nl; ++j) {
            const std::size_t idx = base + static_cast<std::size_t>(j - i - 1);
            if (ll_active[idx]) continue;
            const double grad = S(i, j) - W(i, j);
            const double viol = std::abs(grad) - lambda;
            if (viol > Scale) return false;
        }
    }
    return true;
}

Eigen::MatrixXd Nabla_LL_only(const Eigen::MatrixXd& S, const RowMat& W_t, const int Nl)
{
    return S.topLeftCorner(Nl, Nl) - W_t.topLeftCorner(Nl, Nl);
}

Eigen::MatrixXd Nabla_lhg_LH(const Eigen::MatrixXd& S, const RowMat& W_t,
    const Eigen::SparseMatrix<double>& X_t,const int Nl, const int Nh,
    const double Alpha, const double Gamma)
{
    const int H0 = Nl;
    const int H1 = Nl + Nh;

    // base: Nl × Nh
    Eigen::MatrixXd LH = S.block(0, H0, Nl, Nh) - W_t.block(0, H0, Nl, Nh);
    const double scaleLH = 2.0 * Alpha * (1.0 - Gamma);

    // 只扫 H 列：j ∈ [H0, H1)
    for (int j = H0; j < H1; ++j) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(X_t, j); it; ++it) {
            const int i = it.row();
            if (i >= 0 && i < Nl) {
                LH(i, j - H0) += scaleLH * it.value();
            }
        }
    }
    return LH;
}

void add_to_UT_col_from_Wrow(RowMat& U_T, const int col, const RowMat& W,
    const int rowW, const double mu)
{
    const int n = U_T.rows();
    const int ld = U_T.cols();
    double* ut = U_T.data() + col;
    const double* w = W.data() + rowW * W.cols();
    for (int r = 0; r < n; ++r) ut[r * ld] += mu * w[r];
}

Eigen::SparseMatrix<double> Expand_theta(const Eigen::SparseMatrix<double>& A, const int Nf)
{
    const int Nl = A.rows();
    const int Nextra = Nf - Nl;
    const int N = Nl + Nextra;
    Eigen::SparseMatrix<double> B(N, N);
    // total nnz = original nnz + new diagonal ones
    B.reserve(A.nonZeros() + Nextra);
    // copy A into top-left block
    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            B.insert(it.row(), it.col()) = it.value();
        }
    }
    // add identity on the new diagonal block
    for (int i = 0; i < Nextra; ++i) {
        B.insert(Nl + i, Nl + i) = 1.0;
    }
    B.makeCompressed();
    return B;
}