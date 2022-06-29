#include<stdio.h>
#include <iostream>
#include<stdlib.h>
#include<string.h>
#include<ctime>
#include "pecsb2.h"
#include<hdf5.h>
using namespace std;

/////////////////////////////data save in sbopaperout.h5///////////////
#define SAVE "sbopaperout.h5"

#define YZMAPR "yzmapr"
#define YZMAPI "yzmapi"

#define XZMAPR "xzmapr"
#define XZMAPI "xzmapi"

#define XYMAPR "xymapr"
#define XYMAPI "xymapi"


#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif
#ifdef WIN32
#endif
#ifdef WIN32

/////////////////////////////time///////////////////////////
int gettimeofday(struct timeval *tp, void *tzp)
{
	time_t clock;
	struct tm tm;
	SYSTEMTIME wtm;
	GetLocalTime(&wtm);
	tm.tm_year = wtm.wYear - 1900;
	tm.tm_mon = wtm.wMonth - 1;
	tm.tm_mday = wtm.wDay;
	tm.tm_hour = wtm.wHour;
	tm.tm_min = wtm.wMinute;
	tm.tm_sec = wtm.wSecond;
	tm.tm_isdst = -1;
	clock = mktime(&tm);
	tp->tv_sec = clock;
	tp->tv_usec = wtm.wMilliseconds * 1000;
	return (0);
}
#endif
long getCurrentTime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}


//////////////////////////////////////////////////////////////
int main(){
	clock_t sum_t = getCurrentTime();

	double f_MHz = 300;
	double res_pe = 1; 
	int X = 1000; 
	int Y = 1000; 
	double Z =6000; 
	double dz = 5;   
    int hcen2=30;
	int dimyz=Y*res_pe*Z/dz;
	int dimxz=Y*res_pe*Z/dz;
	int dimxy=X*Y*res_pe*res_pe;
	
	complex<double>* pYZ=new complex<double>[dimyz], *pXY=new complex<double>[dimxy];
    complex<double>* pXZ=new complex<double>[dimxz];
	double* pyzr=new double[dimyz];
	double* pyzi=new double[dimyz];
	double* pxzr=new double[dimxz];
	double* pxzi=new double[dimxz];
	double* pxyr=new double[dimxy];
	double* pxyi=new double[dimxy];

	

	hid_t save_id = H5Fcreate(SAVE, H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT);
	herr_t status;
	unsigned long  long dimsyz=Y*res_pe*Z/dz;
	unsigned long  long dimsxz=X*res_pe*Z/dz;
	unsigned long long dimsxy=X*Y*res_pe*res_pe;

	hsize_t dim5[1]={dimsyz};
	hsize_t dim5x[1]={dimsxz};
	hsize_t dim6[1]={dimsxy};
	
	hid_t yzdataspace=H5Screate_simple(1,dim5,NULL);
	hid_t xzdataspace=H5Screate_simple(1,dim5x,NULL);
	hid_t xydataspace=H5Screate_simple(1,dim6,NULL);
	
	hid_t yzdatasetr=H5Dcreate(save_id,YZMAPR,H5T_NATIVE_DOUBLE,yzdataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t yzdataseti=H5Dcreate(save_id,YZMAPI,H5T_NATIVE_DOUBLE,yzdataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	hid_t xzdatasetr=H5Dcreate(save_id,XZMAPR,H5T_NATIVE_DOUBLE,xzdataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t xzdataseti=H5Dcreate(save_id,XZMAPI,H5T_NATIVE_DOUBLE,xzdataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	hid_t xydatasetr=H5Dcreate(save_id,XYMAPR,H5T_NATIVE_DOUBLE,xydataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t xydataseti=H5Dcreate(save_id,XYMAPI,H5T_NATIVE_DOUBLE,xydataspace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	pecsb2(f_MHz,res_pe,X,Y,Z,dz,pYZ,pXZ,pXY);
	
	for (int i=0;i<dimyz;++i){
		pyzr[i]=pYZ[i].real();
		pyzi[i]=pYZ[i].imag();
	}

	for (int i=0;i<dimxz;++i){
		pxzr[i]=pXZ[i].real();
		pxzi[i]=pXZ[i].imag();
	}

	for (int i=0;i<dimxy;++i){
		pxyr[i]=pXY[i].real();
		pxyi[i]=pXY[i].imag();
	}

	
	
	status=H5Dwrite(yzdatasetr,H5T_NATIVE_DOUBLE,H5S_ALL, H5S_ALL, H5P_DEFAULT, pyzr);
	status=H5Dwrite(yzdataseti,H5T_NATIVE_DOUBLE,H5S_ALL, H5S_ALL, H5P_DEFAULT, pyzi);

	status=H5Dwrite(xzdatasetr,H5T_NATIVE_DOUBLE,H5S_ALL, H5S_ALL, H5P_DEFAULT, pxzr);
	status=H5Dwrite(xzdataseti,H5T_NATIVE_DOUBLE,H5S_ALL, H5S_ALL, H5P_DEFAULT, pxzi);

	status=H5Dwrite(xydatasetr,H5T_NATIVE_DOUBLE,H5S_ALL, H5S_ALL, H5P_DEFAULT, pxyr);
	status=H5Dwrite(xydataseti,H5T_NATIVE_DOUBLE,H5S_ALL, H5S_ALL, H5P_DEFAULT, pxyi);



	status=H5Sclose(yzdataspace);
	status=H5Sclose(xzdataspace);
	status=H5Sclose(xydataspace);

	status=H5Dclose(yzdatasetr);
	status=H5Dclose(yzdataseti);
	status=H5Dclose(xzdatasetr);
	status=H5Dclose(xzdataseti);
	status=H5Dclose(xydatasetr);
	status=H5Dclose(xydataseti);

	status=H5Fclose(save_id);

	cout << "time cost(ms):" << getCurrentTime() - sum_t << endl;
	//system("pause");
}
