#ifndef COMPLEX_H
#define COMPLEX_H
#include<iostream>
#include<fstream>
#include<complex>
template<typename T>
class MyComplex{
	T real, imag;
public:
	MyComplex(const T&r, const T&i) :real(r), imag(i){}
	MyComplex(const T&r) :real(r), imag(0){}
	MyComplex() :real(0), imag(0){}
	std::complex<T> toStd(){ return std::complex<T>(real, imag); }
	friend MyComplex<T> operator+(const MyComplex<T>&c1, const MyComplex<T>&c2){
		return MyComplex<T>(c1.real + c2.real, c1.imag + c2.imag);
	}
	friend MyComplex<T> operator-(const MyComplex<T>&c1, const MyComplex<T>&c2){
		return MyComplex<T>(c1.real - c2.real, c1.imag - c2.imag);
	}
	friend MyComplex<T> operator*(const MyComplex<T>&c1, const MyComplex<T>&c2){
		return MyComplex<T>(c1.real*c2.real - c1.imag*c2.imag, c1.real*c2.imag + c1.imag*c2.real);
	}
	//template<typename T>
	friend MyComplex<T> specialMul(const MyComplex<T>&c1, const MyComplex<T>&c2){
		return MyComplex<T>(c1.real*c2.real - c1.imag*c2.imag, c1.real*c2.imag + c1.imag*c2.real);
	}
	friend MyComplex<T> operator/(const MyComplex<T>&c1, const MyComplex<T>&c2){
		T deno = c2.real*c2.real + c2.imag*c2.imag;
		return MyComplex<T>((c1.real*c2.real + c1.imag*c2.imag) / deno, (c1.imag*c2.real - c1.real*c2.imag) / deno);
	}
	friend std::ostream& operator<<(std::ostream&out, const MyComplex<T>&c){
		out << '(' << c.real << ',' << c.imag << ')';
		return out;
	}
	friend std::ofstream& operator<<(std::ofstream&fout, const MyComplex<T>&c){
		fout << '(' << c.real << ',' << c.imag << ')';
		return fout;
	}
	friend bool operator==(const MyComplex&c1, const MyComplex&c2){
		return c1.real == c2.real&&c1.imag == c2.imag;
	}
	MyComplex operator-()const{
		return MyComplex(-real, -imag);
	}
};
#endif
