#ifndef PECSB2_H
#define PECSB2_H
#include <complex>
#include <cmath>
void pecsb2(double f_MHz, double res_pe, int X, int Y, double Z,double dz,std::complex<double>* pYZ,std::complex<double>* pXZ,std::complex<double>* PXY);
//void pecsb2(const char* filename, const char * dataset_r, const char* dataset_i);
#endif