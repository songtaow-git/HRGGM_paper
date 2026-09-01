#ifndef _MATRIX_IO_H_
#define _MATRIX_IO_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <filesystem>

using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

std::string System_time();

void Write_metric(const RowMat& A, const std::filesystem::path Out_dir, 
	const std::string Filename);

RowMat Load_bin_f64_data(const std::string& path);

void WriteSparseForMatlab(const Eigen::SparseMatrix<double>& A, const std::string& filename);

void WriteDenseForMatlab(const Eigen::MatrixXd& A, const std::string& filename);

void WriteVectorToTxt(const std::vector<int>& vec, const std::string& filename);

std::vector<double> ReadParameterFromTxt(const std::string& filename);

std::vector<int> ReadidxFromTxt(const std::string& filename);

void WritePairVectorToFile(const std::vector<std::pair<int, int>>& vec,
	const std::string& filename);

#endif // _MATRIX_IO_H_