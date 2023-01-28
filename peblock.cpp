#include "DyArray.h" 
#include "Solve_Tridiagonal.h"
#include "omp.h"
#include "peblock.h"
#include <complex>
#include <cmath>
#include <iostream>
#include <iomanip>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctime>



void peblock(DyArray2<Complex<double> > U, double dx, double dy, double dz, int Xmax, int Ymax, DyArray2<double> filx, DyArray2<double> fily, double K0, double front, DyArray2<Complex<double> > &U2, DyArray2<Complex<double> > &yz,DyArray2<Complex<double> > &xz,DyArray2<Complex<double> > &xy){
	const double pi = 3.141592653589793;
	int Zmax = front / dz;
    double cLight = 3e8;
	double f_MHz;
    double lamda = 2*pi/K0;
    double Order; 
	Order=0.0833333333333;//1/12
    int Pmax=3; 
	int pol=0;
	double epsilon_g = 30;
    double sigma_g = 0.01;       
	int x=0;
	int n = 1;
	int y=0;
	int m = 2;
	int l = 2;
	printf("x=%d\n",x);    
    int hcen2=30;    

     /////////////////////////////////////////////////////////////////
	 double para=60;
	 complex<double> epsilon_e;
	 epsilon_e=epsilon_g+complex<double>(0,1)*para*lamda*sigma_g;
	 
	 complex<double> sqrepsilon_e;
	 sqrepsilon_e=sqrt(epsilon_e);
	 double a;
     a=sqrepsilon_e.real();
	 double b;
     b=sqrepsilon_e.imag();
     Complex<double> gd;
	 gd=a+Complex<double>(0,1)*b;
	 Complex<double> g=1/gd;
     printf("a=%\n",a);  
	
	 DyArray2<Complex<double> > beta(1, 1, 0);
	

     if (pol==0){
	     beta[0][0]=Complex<double>(0,1)*K0/g;
		  cout<<beta[0][0]<<endl;
	 }
     else   {   
		 Complex<double> epsilon_e;                                                                
         beta[0][0]== gd/epsilon_e;
	     }
     


	DyArray2<Complex<double> > aa (1, 3, 0);
    DyArray2<Complex<double> > bb (1, 3, 0);


	
    DyArray2<Complex<double> > Rx(Xmax, Ymax, 0);
	DyArray2<Complex<double> > Ry(Xmax, Ymax, 0);
	DyArray2<Complex<double> > main(Xmax, Ymax, 0);
	DyArray2<Complex<double> > Utemp(Xmax, Ymax, 0);
	DyArray2<Complex<double> > su(Xmax, Ymax-1, 0);
	DyArray2<Complex<double> > su2(Xmax-1, Ymax, 0);
    DyArray2<Complex<double> > rhs(Xmax, Ymax, 0);
	


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
for (int n = 1; n <= Zmax; ++n) {
   for (int p = 1; p <= Pmax; ++p) {	
	    #pragma omp parallel for private(x,y)
		for (x = 0;x < Xmax;++x) {
			for (y = 0;y < Ymax;++y) {
                 Rx[x][y]=1/(dx*K0*dx*K0);
				 Ry[x][y]=1/(dy*K0*dy*K0);
                 main[x][y]=1-2*(Order+bb[0][p-1]*Ry[x][y]);


			}
		}
		
		#pragma omp parallel for private(x,y)
		for (x=0;x<Xmax;++x){
			for (y=0;y<Ymax-1;++y){
				su[x][y]=Order+bb[0][p-1]*Ry[x][y];
			}
		}
		#pragma omp parallel for private(x,y)
		for (x=0;x<Xmax-1;++x){
			for (y=0;y<Ymax;++y){
				su2[x][y]=Order+bb[0][p-1]*Ry[x][y];
			}
		}


        #pragma omp parallel for private(x,y)
		for (x = 2; x <= Xmax - 1; ++x) {
			//#pragma omp simd
			for (y = 2; y <= Ymax - 1; ++y) {
				rhs[x - 1][y - 1] = (Order+aa[0][p-1]*Ry[x-1][y-1])*U[x - 1][y - 2] + (1-2*(Order+aa[0][p-1]*Ry[x-1][y-1]))*U[x - 1][y - 1] +(Order+aa[0][p-1]*Ry[x-1][y-1])*U[x - 1][y];
			}
            rhs[x - 1][0] = 2*(Order+aa[0][p-1]*Ry[x-1][0])*U[x - 1][1] + 2*(1-(Order+aa[0][p-1]*Ry[x-1][0])*(1-beta[0][0]*dy))*U[x - 1][0];//IBC
		}
               

		DyArray2<Complex<double> > x_1(main.getDims(0), 1, 0);
		DyArray2<Complex<double> > rhs_col_1(rhs.getDims(0), 1, 0);
		DyArray2<Complex<double> > main_col_1(main.getDims(0), 1, 0);
		DyArray2<Complex<double> > su_col(su.getDims(0), 1, 0);
		rhs.subCol(rhs_col_1, 0);
		main.subCol(main_col_1, 0);
		su.subCol(su_col, 0);
		//cout << "Before Solve:" <<"rhs:"<< rhs[1300][0] <<"U:"<<U[1300][0]<< endl;
		Solve_Tridiagonal(main_col_1.getData(), su_col.getData(), su_col.getData(), rhs_col_1.getData(), x_1.getData(), Xmax);
		//cout<<"After Solve:" << "main:"<<main.getData()[100]<<"   sup:"<<sup.getData()[100]<<"   sub:"<<sub.getData()[100]<<"   rhs_col_1:"<<rhs_col_1.getData()[1298]<<"   x:"<<x.getData()[10] << endl;

		Utemp.copyCol(0, x_1.getData(), x_1.getDims(0)*x_1.getDims(1));

        #pragma omp parallel for private(l)
		for (l = 2; l <= Ymax - 1; ++l) {
			DyArray2<Complex<double> > rhs_col(rhs.getDims(0), 1, 0);
			DyArray2<Complex<double> > main_col(main.getDims(0), 1, 0);
			DyArray2<Complex<double> > su_col(su.getDims(0), 1, 0);
			rhs.subCol(rhs_col, l - 1);
			main.subCol(main_col, l - 1);
			su.subCol(su_col, l - 1);
			DyArray2<Complex<double> > xx(main.getDims(0), 1, 0);
			Solve_Tridiagonal(main_col.getData(), su_col.getData(), su_col.getData(), rhs_col.getData(), xx.getData(), Xmax);
			Utemp.copyCol(l - 1, xx.getData(), xx.getDims(0)*xx.getDims(1));
		}


		#pragma omp parallel for private(x)
		for (x = 0;x < Xmax;++x) {
                 main[x][0]=2*(1-(Order+bb[0][p-1]*Ry[x-1][0])*(1-beta[0][0]*dy));//LBC
		}



        #pragma omp parallel for private(x,y)
		for (x = 2; x <= Xmax - 1; ++x) {
			//#pragma omp simd
			for (y = 1; y <= Ymax ; ++y) {
				rhs[x - 1][y - 1] = (Order+aa[0][p-1]*Ry[x-1][y-1])*Utemp[x - 2][y - 1] + (1-2*(Order+aa[0][p-1]*Ry[x-1][y-1]))*Utemp[x - 1][y - 1] + (Order+aa[0][p-1]*Ry[x-1][y-1])*Utemp[x][y - 1];
			}
		}
        #pragma omp parallel for private(m)
		for (m = 2; m <= Xmax - 1; ++m){
			Solve_Tridiagonal(main.getData()+(m-1)*main.getDims(1), su2.getData()+(m-1)*su2.getDims(1), su2.getData()+(m-1)*su2.getDims(1), rhs.getData() + (m - 1)*rhs.getDims(1), U.getData() + (m - 1)*U.getDims(1), Ymax);
		}

	
//////////////////////////////////////////////////////////////////////////////////


		// 	for (x = 1; x <= Xmax - 1; ++x) {//pec
		// 	//#pragma omp simd
		// 	for (y = 1; y <= floor(10/dy); ++y) {
		// 		U[x - 1][y-1] =0;
		// 	}
		// }
		
    }
//////////////////////////////save data/////////////////////////////////////////////////

		U.HadmardProduct(fily);
		U.HadmardProduct(filx);

		double Rwy;
		Rwy=round(0.5*Xmax);
		for (int i=0;i<yz.getDims(0);++i){
			yz[i][n-1]=U[Rwy-1][i];
		}

        int Rwx;
		Rwx=round(hcen2/dy);
		for (int i=0;i<Xmax;++i){
			xz[i][n-1]=U[i][Rwx-1];
		}




		if (n==Zmax){
			xy=U;
		}
		printf("obstacle:%d\n", n);
	}
	U2 = U;
   
}









