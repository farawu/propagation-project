#ifndef COMPLEX_H
#define COMPLEX_H
#include<iostream>
#include<fstream>
template<typename T>
class Complex{
	T real, imag;
public:
	Complex(const T&r, const T&i) :real(r), imag(i){}
	Complex(const T&r) :real(r), imag(0){}
	Complex() :real(0), imag(0){}
	friend Complex<T> operator+(const Complex<T>&c1, const Complex<T>&c2){
		return Complex<T>(c1.real + c2.real, c1.imag + c2.imag);
	}
	friend Complex<T> operator-(const Complex<T>&c1, const Complex<T>&c2){
		return Complex<T>(c1.real - c2.real, c1.imag - c2.imag);
	}
	friend Complex<T> operator*(const Complex<T>&c1, const Complex<T>&c2){
		return Complex<T>(c1.real*c2.real - c1.imag*c2.imag, c1.real*c2.imag + c1.imag*c2.real);
	}
	template<typename T>
	friend Complex<T> specialMul(const Complex<T>&c1, const Complex<T>&c2){
		return Complex<T>(c1.real*c2.real - c1.imag*c2.imag, c1.real*c2.imag + c1.imag*c2.real);
	}
	friend
	friend Complex<T> operator/(const Complex<T>&c1, const Complex<T>&c2){
		T deno = c2.real*c2.real + c2.imag*c2.imag;
		return Complex<T>((c1.real*c2.real + c1.imag*c2.imag) / deno, (c1.imag*c2.real - c1.real*c2.imag) / deno);
	}
	friend std::ostream& operator<<(std::ostream&out, const Complex<T>&c){
		out << '(' << c.real << ',' << c.imag << ')';
		return out;
	}
	friend std::ofstream& operator<<(std::ofstream&fout, const Complex<T>&c){
		fout << '(' << c.real << ',' << c.imag << ')';
		return fout;
	}
	friend bool operator==(const Complex&c1, const Complex&c2){
		return c1.real == c2.real&&c1.imag == c2.imag;
	}
	Complex operator-()const{
		return Complex(-real, -imag);
	}
};
#endif