// #ifndef PEBLOCK_H
// #define PEBLOCK_H
// #include"DyArray.h"
// #include<complex>
// #ifdef STD_COMPLEX
// template<typename T>
// using Complex = std::complex<T>;
// #else
// #include "MyComplex.h"
// template<typename T>
// using Complex = MyComplex<T>;
// #endif
// //DyArray2 �����ͨ�����������ͷ��ڴ�
// void peblock(DyArray2<Complex<double> > U, double step1, double dx, double dy, double dz, double Z2, int Xmax, int Ymax, DyArray2<double> filx, DyArray2<double> fily, double K0, double front, double rwy, DyArray2<Complex<double> > &U2,float*** N, DyArray2<Complex<double> > &ddm);
// #endif

#ifndef PEBLOCK_H
#define PEBLOCK_H
#include"DyArray.h"
#include<complex>
#ifdef STD_COMPLEX
template<typename T>
using Complex = std::complex<T>;
#else
#include "MyComplex.h"
template<typename T>
using Complex = MyComplex<T>;
#endif
//DyArray2 �����ͨ�����������ͷ��ڴ�
void peblock(DyArray2<Complex<double> > U, double dx, double dy, double dz, int Xmax, int Ymax, DyArray2<double> filx, DyArray2<double> fily, double K0, double front, DyArray2<Complex<double> > &U2,DyArray2<Complex<double> > &yz,DyArray2<Complex<double> > &xz,DyArray2<Complex<double> > &xy);
#endif
