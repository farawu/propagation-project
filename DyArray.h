#ifndef DYARRAY_H
#define DYARRAY_H
#include <cstring>
#include "MyComplex.h"
#include "hdf5.h"
#include <complex>
//
template<typename T>
class DyArray2{
	int dims[2];
	T *data;
public:
	T* getData(){
		return data;
	}
	const T* getData()const{
		return data;
	}
	int getDims(int i)const{
		return dims[i];
	}
	DyArray2(const DyArray2& a){
		dims[0] = a.dims[0];
		dims[1] = a.dims[1];
		data = new T[dims[0] * dims[1]];
		memcpy(data, a.data, sizeof(T)*dims[0] * dims[1]);
	}
	DyArray2(int x, int y, T init){
		dims[0] = x;
		dims[1] = y;
		data = new T[dims[0] * dims[1]];
		if (init==T(0))
		    memset(data, 0, sizeof(T)*dims[0] * dims[1]);
		else{
			for (int i = 0; i < dims[0] * dims[1]; ++i)
				data[i] = init;
		}
	}
	/*
	DyArray2(T start, T step, T end){
		dims[0] = 1;
		dims[1] = (end - start) / step + 1;
		data = new T[dims[1]];
		for (int i = 0; i < dims[1]; ++i)
			data[i] = start + i*step;
	}
	*/
	/*DyArray2(T* pData, int x,int y){
		dims[0] = x;
		dims[1] = y;
		data = pData;
	}*/
	T* operator[](int index){
		if (!data){
			throw "operator[]:data pointer = 0!";
		}
		return data + index*dims[1];
	}
	void subArray(DyArray2& dest,int x1, int x2, int y1, int y2){
		if (!(x1 >= 0 && x2 <= dims[0] && y1 >= 0 && y2 <= dims[1]))
			throw "destArray: dimensions doesn't equal!";
		for (int i = 0; i < dest.dims[0]; ++i)
			memcpy(dest.data + i*dest.dims[1], data + (x1 + i)*dims[1]+y1, sizeof(T)*dest.dims[1]);
	}
	void subCol(DyArray2&dest, int col){
		for (int i = 0; i < dims[0]; ++i)
			dest.data[i] = data[i*dims[1] + col];
	}
	void copyArray(int x1, int x2, int y1, int y2, const DyArray2<T>& src){
		if (!(x2 - x1 == src.getDims(0) && y2 - y1 == src.getDims(1)))
			throw "copyArray:dimensions doesn't equal!";
		for (int i = 0; i < src.dims[0]; ++i){
			memcpy(data + (x1 + i)*dims[1] + y1, src.data + i*src.dims[1], sizeof(T)*src.dims[1]);
		}
	}
	void copyRow(int row, const T* src, int len){
		if(!(row<dims[0]&&dims[1] == len))
			throw "copyRow:dimensions doesn't equal!";
		memcpy(data + row*dims[1], src, sizeof(T)*dims[1]);
	}
	void copyCol(int col, const T* src, int len){
		if (!(col < dims[1] && dims[0] == len))
			throw "copyCol:dimensions doesn't equal!";
		for (int i = 0; i < dims[0]; ++i)
			data[i*dims[1] + col] = src[i];
	}
	/*
	friend DyArray2 HadmardProduct(const DyArray2& a1, const DyArray2& a2){
		if (!(a1.dims[0] == a2.dims[0] && a1.dims[1] == a2.dims[1]))
			throw "HadmardProduct:dimensions doesn't equal!";
		DyArray2 mul(a1.dims[0], a1.dims[1], 0);
		for (int i = 0; i < a1.dims[0] * a1.dims[1]; ++i){
			mul.data[i] = a1.data[i] * a2.data[i];
		}
		return mul;
	}
	*/
    template<typename S>
	void HadmardProduct(const DyArray2<S>& a){
		if (!(dims[0] == a.getDims(0) && dims[1] == a.getDims(1)))
			throw "HadmardProduct(one argument):dimensions doesn't equal!";
		int i;
		#pragma omp parallel for private(i)
		for (i = 0; i < a.getDims(0) * a.getDims(1); ++i)
			data[i] = data[i]*a.getData()[i];
	}
	
	/*
	void free(){
		delete[]data;
		dims[0] = dims[1] = dims[2] = 0;
		data = 0;
	}
	*/

	~DyArray2(){
		if (data)
			delete[]data;
		dims[0] = dims[1] = 0;
		data = 0;
	}
	DyArray2& operator=(const DyArray2& a){
		if (data&&dims[0]*dims[1]!=a.dims[0]*a.dims[1]){
			delete[]data;
			data = 0;
		}
		if (!data)
			data = new T[a.dims[0] * a.dims[1]];
			
		//assert(data&&dims[0] == a.dims[0] && dims[1] == a.dims[1]);
		dims[0] = a.dims[0];
		dims[1] = a.dims[1];
		memcpy(data, a.data, sizeof(T)*dims[0] * dims[1]);
		return *this;
	}
	/*
	DyArray2 clone(){
		DyArray2 c(dims[0], dims[1], 0);
		memcpy(c.data, data, sizeof(T)*dims[0] * dims[1]);
		return c;
	}
	*/

};
/*
template<typename T>
DyArray2<T> read2(const char *file_name, const char *dataset_name){
	try
	{
		H5::H5File file(file_name, H5F_ACC_RDONLY);
		H5::DataSet dataset = file.openDataSet(dataset_name);
		H5::DataSpace dataspace = dataset.getSpace();
		int rank = dataspace.getSimpleExtentNdims();
		if (rank != 2){
			std::cout << "The rank of " << dataset_name << " is not 2!" << std::endl;
			exit(-1);
		}
		hsize_t dims_out[2];
		int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
		hsize_t  offset[2] = { 0, 0 };   // hyperslab offset in the file
		hsize_t  count[2] = { dims_out[0], dims_out[1] };    // size of the hyperslab in the file
		dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
		H5::DataSpace memspace(rank, dims_out);
		memspace.selectHyperslab(H5S_SELECT_SET, count, offset);
		DyArray2<T> data_out((int)dims_out[0], (int)dims_out[1],0);
		dataset.read(data_out.getData(), typeid(T) == typeid(int) ? H5::PredType::NATIVE_INT:(typeid(T) == typeid(float) ? H5::PredType::NATIVE_FLOAT : H5::PredType::NATIVE_DOUBLE), memspace, dataspace);
		return data_out;
	}
	catch (H5::FileIException error)
	{
		error.printErrorStack();
		exit(-1);
	}
	// catch failure caused by the DataSet operations
	catch (H5::DataSetIException error)
	{
		error.printErrorStack();
		exit(-1);
	}
	// catch failure caused by the H5::DataSpace operations
	catch (H5::DataSpaceIException error)
	{
		error.printErrorStack();
		exit(-1);
	}
	// catch failure caused by the H5::DataSpace operations
	catch (H5::DataTypeIException error)
	{
		error.printErrorStack();
		exit(-1);
	}
	return DyArray2<T>(0,0,0);
}
*/

template<typename T>
DyArray2<T> read2(const char *file_name, const char *dataset_name){

	// ��HDF5�ļ�
	hid_t file=H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);

	hid_t dataset = H5Dopen(file, dataset_name, H5P_DEFAULT);
	hid_t dataspace = H5Dget_space(dataset);
	int rank = H5Sget_simple_extent_ndims(dataspace);
	hsize_t dims_out[2];
	int status = H5Sget_simple_extent_dims(dataspace, dims_out,NULL);
	DyArray2<double> data_out(dims_out[0], dims_out[1], 0);
	hid_t memspace = H5Screate_simple(rank, dims_out, NULL);
	
	hsize_t offset[2] = { 0, 0 };
	hsize_t stride[2] = { 1, 1 };
	status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset,
		NULL, dims_out, NULL);


	status = H5Dread(dataset, typeid(T) == typeid(int) ? H5T_NATIVE_INT : (typeid(T) == typeid(float) ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE), memspace, dataspace, H5P_DEFAULT, data_out.getData());

	// 
	status = H5Dclose(dataset);
	status = H5Sclose(dataspace);

	// 
	status = H5Fclose(file);
	return data_out;
}




/*
template<typename T>
DyArray2<MyComplex<T> > toComplex(DyArray2<T> a1, DyArray2<T> a2){
	DyArray2<MyComplex<T> > c(a1.getDims(0), a2.getDims(1),0);
	for (int i = 0; i < a1.getDims(0); ++i){
		for (int j = 0; j < a1.getDims(1); ++j){
			c[i][j] = MyComplex<T>(a1[i][j], a2[i][j]);
		}
	}
	return c;
}
*/
template<typename T>
DyArray2<std::complex<T> > tocomplex(DyArray2<T> a1, DyArray2<T> a2){
	DyArray2<std::complex<T> > c(a1.getDims(0), a2.getDims(1), 0);
	for (int i = 0; i < a1.getDims(0); ++i){
		for (int j = 0; j < a1.getDims(1); ++j){
			c[i][j] = std::complex<T>(a1[i][j], a2[i][j]);
		}
	}
	return c;
}

template<typename T>
DyArray2<MyComplex<T> > toComplex(DyArray2<T> a1, DyArray2<T> a2){
	DyArray2<MyComplex<T> > c(a1.getDims(0), a2.getDims(1), 0);
	for (int i = 0; i < a1.getDims(0); ++i){
		for (int j = 0; j < a1.getDims(1); ++j){
			c[i][j] = MyComplex<T>(a1[i][j], a2[i][j]);
		}
	}
	return c;
}
template<typename T>
DyArray2<MyComplex<T> > toComplex_nega(DyArray2<T> a1, DyArray2<T> a2){
	DyArray2<MyComplex<T> > c(a1.getDims(0), a2.getDims(1), 0);
	for (int i = 0; i < a1.getDims(0); ++i){
		for (int j = 0; j < a1.getDims(1); ++j){
			c[i][j] = MyComplex<T>(a1[i][j], -a2[i][j]);
		}
	}
	return c;
}
/*
template<typename T>
class DyArray3{
int dims[3];
T *data;
public:
T* getData(){
return data;
}
int getDims(int i){
return dims[i];
}
DyArray3(int x,int y,int z){
dims[0] = x;
dims[1] = y;
dims[2] = z;
data = new T[dims[0] * dims[1] * dims[2]];
memset(data,0, sizeof(T)*dims[0]*dims[1]*dims[2]);
}
DyArray3(T *pData, int x, int y, int z){
dims[0] = x;
dims[1] = y;
dims[2] = z;
data = pData;
}
DyArray2<T> operator[](int index){
return DyArray2<T>(data + index*dims[1] * dims[2], dims[1], dims[2]);
}
void free(){
delete[]data;
data = dims[0] = dims[1] = dims[2] = 0;
}
};
*/
/*
template<typename T>
T** malloc2(int x, int y){
	T** dest = (T**)malloc(sizeof(T*)*x);
	for (int i = 0; i < x; ++i){
		dest[i] = (T*)malloc(sizeof(T)*y);
		memset(dest[i], 0, sizeof(T)*y);
	}
	return dest;
}
template<typename T>
T*** malloc3(int x, int y, int z){
	T***dest = (T***)malloc(sizeof(T**)*x);
	for (int i = 0; i < x; ++i){
		dest[i] = (T**)(malloc(sizeof(T*)*y));
		for (int j = 0; j < y; ++j){
			dest[i][j] = (T*)(malloc(sizeof(T)*z));
			memset(dest[i][j], 0, sizeof(T)*z);
		}
	}
	return dest;
}
template<typename T>
void free2(T** target, int x){
	for (int i = 0; i < x; ++i)
		free(target[x]);
}
template<typename T>
void free3(T***target, int x, int y){
	for (int i = 0; i < x; ++i){
		for (int j = 0; j < y; ++j)
			free(target[i][j]);
		free(target[i]);
	}
	free(target);
}
*/
/*
template<typename T>
DyArray3<T> read3(const char *file_name, const char *dataset_name){
	try{
		H5File file(file_name, H5F_ACC_RDONLY);
		DataSet dataset = file.openDataSet(dataset_name);
		H5::DataSpace dataspace = dataset.getSpace();
		int rank = dataspace.getSimpleExtentNdims();
		if (rank != 3){
			std::cout << "The rank of " << dataset_name << " is not 3!" << std::endl;
			exit(-1);
		}
		hsize_t dims_out[3];
		int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
		hsize_t   offset[3] = { 0, 0, 0 };   // hyperslab offset in the file
		hsize_t   count[3] = { dims_out[0], dims_out[1], dims_out[2] };    // size of the hyperslab in the file
		dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
		H5::DataSpace memspace(rank, dims_out);
		memspace.selectHyperslab(H5S_SELECT_SET, count, offset);
		DyArray3<int> data_out(dims_out[0], dims_out[1], dims_out[3]);
		dataset.read(data_out.getData(), H5::PredType::NATIVE_INT, memspace, dataspace);
		return data_out;
	}
	catch (FileIException error)
	{
		error.printErrorStack();
		exit(-1);
	}
	// catch failure caused by the DataSet operations
	catch (DataSetIException error)
	{
		error.printErrorStack();
		exit(-1);
	}
	// catch failure caused by the H5::DataSpace operations
	catch (DataSpaceIException error)
	{
		error.printErrorStack();
		exit(-1);
	}
	// catch failure caused by the H5::DataSpace operations
	catch (DataTypeIException error)
	{
		error.printErrorStack();
		exit(-1);
	}
	return DyArray2<T>(0, 0);
}

*/
#endif
