#ifndef _NEWTON_METHOD_OFFDIAG_H_
#define _NEWTON_METHOD_OFFDIAG_H_

#include <string>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>	
#include "Matrix_compute.h"
#include "Class_define.h"

void Newton_method_offdiag_l(Input_format& Input_l, const int Blockcol_inv = 1024);

void Newton_method_offdiag_lhg(Input_format& Input_l, const int Blockcol_inv = 1024);

#endif // _NEWTON_METHOD_OFFDIAG_H_