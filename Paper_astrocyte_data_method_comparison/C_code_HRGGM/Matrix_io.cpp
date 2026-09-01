#include "Matrix_io.h"

std::string System_time()
{
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);

    std::tm local_tm{};
#ifdef _WIN32
    localtime_s(&local_tm, &now_c);
#else
    localtime_r(&now_c, &local_tm);
#endif

    std::ostringstream oss;
    oss << std::put_time(&local_tm, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

void Write_metric(const RowMat& A, const std::filesystem::path Out_dir, 
    const std::string Filename)
{
    std::filesystem::path out_file = Out_dir / Filename;
    std::ofstream fout(out_file);
    if (!fout) {
        throw std::runtime_error("Cannot open file: " + out_file.string());
    }
    fout << std::scientific << std::setprecision(17);

    for (int r = 0; r < A.rows(); ++r) {
        for (int c = 0; c < A.cols(); ++c) {
            fout << A(r, c);
            if (c + 1 < A.cols()) fout << ' ';
        }
        fout << '\n';
    }
}


RowMat Load_bin_f64_data(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("open failed: " + path);

    uint64_t r = 0, c = 0;
    in.read(reinterpret_cast<char*>(&r), sizeof(r));
    in.read(reinterpret_cast<char*>(&c), sizeof(c));
    if (!in) throw std::runtime_error("header read failed");

    RowMat A((Eigen::Index)r, (Eigen::Index)c);
    in.read(reinterpret_cast<char*>(A.data()), sizeof(double) * r * c);
    if (!in) throw std::runtime_error("data read failed");

    return A;
}

void WriteSparseForMatlab(const Eigen::SparseMatrix<double>& A, const std::string& filename)
{
    std::ofstream fout(filename);
    if (!fout) throw std::runtime_error("open failed: " + filename);

    fout << std::scientific << std::setprecision(17);

    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            fout << A.coeff(i, j);
            if (j + 1 < A.cols()) fout << ' ';
        }
        fout << '\n';
    }
}

void WriteDenseForMatlab(const Eigen::MatrixXd& A, const std::string& filename)
{
    std::ofstream fout(filename);
    if (!fout) throw std::runtime_error("open failed: " + filename);

    fout << std::scientific << std::setprecision(16);

    const int rows = A.rows();
    const int cols = A.cols();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            fout << A(i, j);
            if (j + 1 < cols) fout << ' ';
        }
        fout << '\n';
    }
}

void WriteVectorToTxt(const std::vector<int>& vec, const std::string& filename)
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    for (int x : vec) {
        out << x << '\n';
    }
}

std::vector<double> ReadParameterFromTxt(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string content = buffer.str();
    std::vector<double> values;
    std::stringstream ss(content);
    std::string token;
    while (std::getline(ss, token, ',')) {
        std::stringstream token_stream(token);
        double val;
        if (token_stream >> val) {
            values.push_back(val);
        }
    }
    return values;
}

std::vector<int> ReadidxFromTxt(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    std::stringstream buffer;
    buffer << in.rdbuf();
    std::string content = buffer.str();
    std::vector<int> values;
    std::stringstream ss(content);
    std::string token;
    while (std::getline(ss, token, ',')) {
        std::stringstream token_stream(token);
        int val;
        if (token_stream >> val) {
            values.push_back(val);
        }
    }
    return values;
}

void WritePairVectorToFile(const std::vector<std::pair<int, int>>& vec,
    const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    for (const auto& p : vec) {
        out << p.first << ' ' << p.second << '\n';
    }
}