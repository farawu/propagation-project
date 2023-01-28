#include <complex>
#include "DyArray.h"

#include "peblock.h"
#include <cmath>
#include "Triangle.h"
#include <iostream>
#include <iomanip>
//#include "arraytest.h"



void pecsb2( double f_MHz,double res_pe,int X,int Y,double Z,double dz,std::complex<double>* pYZ,std::complex<double>* pXZ,std::complex<double>* PXY)

{
	const double pi = 3.141592653589793;

	double cLight = 3e8;
	double epsilon = 30;

	double hh=Y;

	double dx = 1 / res_pe;
	double dy = 1 / res_pe;
	double front1=Z;
	double frq = f_MHz*1e6;
	double K0 = 2 * pi*frq / cLight;
	double lamda = cLight / frq;

	int X_max = 1 * X;
	int Xmax = floor(X_max / dx);
	int Ymax2;



		double h=hh;
		int Y_max = h;

		int Ymax = floor(Y_max / dy);
		DyArray2<Complex<double> > Ux(Xmax, Ymax, 0);
		DyArray2<Complex<double> > Ux1(Xmax, Ymax, 0);
		DyArray2<Complex<double> > yzmap(Ymax,(int)(front1/dz),0);
		DyArray2<Complex<double> > xzmap(Xmax,(int)(front1/dz),0);
		DyArray2<Complex<double> > xymap(Xmax,Ymax,0);
	
		

				int hcen,hcen2,hcen3,xcen;
	double beta1,beta2,theta;
	DyArray2<Complex<double>> uver(Xmax,Ymax,0);
	DyArray2<double> uvercos(Xmax,Ymax,0);
	DyArray2<double> uversin(Xmax,Ymax,0);
	DyArray2<double> uvercos2(Xmax,Ymax,0);
	DyArray2<double> uversin2(Xmax,Ymax,0);
	DyArray2<Complex<double>> uhor(Xmax,Ymax,0);
	DyArray2<double> uhorcos(Xmax,Ymax,0);
	DyArray2<double> uhorsin(Xmax,Ymax,0);

	theta=0;
	xcen=0.5*Xmax;
	hcen2=30;
	beta1=15*pi/180;
	beta2=15*pi/180;
///////////////////////////////////////////////

	for(int i=1;i<Xmax;i++){
		for(int j=1;j<Ymax;j++){
			uversin[i][j]=sqrt(2*log(2))/(K0*sin(beta1/2));
			uver[i][j]=exp(-(((i-xcen)*dx)/uversin[i][j])*(((i-xcen)*dx)/uversin[i][j]));
			uhorsin[i][j]=sqrt(2*log(2))/(K0*sin(beta2/2));
			uhor[i][j]=exp(-(((j-hcen2)*dy)/uhorsin[i][j])*(((j-hcen2)*dy)/uhorsin[i][j]));
			Ux[i][j]=uver[i][j]*uhor[i][j];
		}
	}




//////////////////////////////////////////////////////

		double fx = floor(Xmax / 16);
		double fy = floor(Ymax / 16);
		DyArray2<double > filt(1, Ymax, 1);
		DyArray2<double > filttemp(1, Xmax, 1);

		for (int j = 1; j <= Xmax; ++j) {
			if (j <= Xmax / 16)
				filttemp.getData()[j - 1] = 0.5 + 0.5*cos(pi*(Xmax / 16 - j + 1) / fx);
			else if (j > 0.9375*Xmax)
				filttemp.getData()[j - 1] = 0.5 + 0.5*cos(pi*(j - 0.9375*Xmax) / fx);

		}
		for (int j = 1; j <= Ymax; ++j) {
			if (j > 0.9375*Ymax)
				filt.getData()[j - 1] = 0.5 + 0.5*cos(pi*(j - 0.9375*Ymax) / fy);
		}
		//DyArray2<double> fil(Xmax, Ymax,0);
		DyArray2<double> filx(Xmax, Ymax, 0);
		DyArray2<double> fily(Xmax, Ymax, 0);
		for (int i = 1; i < Xmax; ++i) {
			fily.copyRow(i - 1, filt.getData(), Ymax);
		}
		for (int i = 1; i <= Ymax; ++i) {
			filx.copyCol(i - 1, filttemp.getData(), Xmax);
		}


		DyArray2<double> fil(filx);
		fil.HadmardProduct(fily);
		DyArray2<Complex<double> > xytemp(Xmax,Ymax,0);
		DyArray2<Complex<double> > yztemp(Ymax,(int)(front1/dz),0);
        DyArray2<Complex<double> > xztemp(Xmax,(int)(front1/dz),0);
		
		peblock(Ux, dx, dy, dz, Xmax, Ymax, filx, fily,  K0, front1, Ux1, yztemp,xztemp,xytemp);
		
	

		yzmap.copyArray(0, yzmap.getDims(0), 0, front1 / dz, yztemp);
		xzmap.copyArray(0, xzmap.getDims(0), 0, front1 / dz, xztemp);
		xymap.copyArray(0, xymap.getDims(0), 0,Ymax,xytemp);



	for (int i=0;i<yzmap.getDims(0)*yzmap.getDims(1);++i)
		pYZ[i]=yzmap.getData()[i].toStd();
	for (int i=0;i<xzmap.getDims(0)*xzmap.getDims(1);++i)
		pXZ[i]=xzmap.getData()[i].toStd();	
	for (int i=0;i<xymap.getDims(0)*xymap.getDims(1);++i)
		PXY[i]=xymap.getData()[i].toStd();
	// successfully terminated
}
