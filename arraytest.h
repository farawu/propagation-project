#include"hdf5.h"
template<typename T>
void read2disperse(T** mem_data,hsize_t dims[2], const char *file_name, const char *dataset_name){
	hid_t file = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);
	hid_t dataset = H5Dopen(file, dataset_name, H5P_DEFAULT);
	hid_t dataspace = H5Dget_space(dataset);
	int rank = H5Sget_simple_extent_ndims(dataspace);
	hsize_t dims_out[2];
	int status = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

	for (int i = 0; i < dims[0]; ++i){
		hsize_t offset[2] = { i, 0 };
		hsize_t count[2] = { 1, dims[1] };
		hid_t tmp_memspace = H5Screate_simple(2, count, NULL);
		status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset,
			NULL, dims_out, NULL);
		status = H5Dread(dataset, typeid(T) == typeid(int) ? H5T_NATIVE_INT : (typeid(T) == typeid(float) ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE), tmp_memspace, dataspace, H5P_DEFAULT, mem_data[i]);
	}
	// 关闭dataset相关对象
	status = H5Dclose(dataset);
	status = H5Sclose(dataspace);

	// 关闭文件对象
	status = H5Fclose(file);
}
template<typename T>
void write2disperse(T** mem_data, hsize_t dims[2], const char *file_name, const char *dataset_name){
	hid_t file_id = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);
	hid_t dataspace = H5Screate_simple(2, dims, NULL);
	int status;
	// 创建数据集中的数据本身
	hid_t dataset_id = H5Dcreate(file_id, dataset_name, typeid(T) == typeid(int) ? H5T_NATIVE_INT : (typeid(T) == typeid(float) ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE), dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	for (int i = 0; i < dims[0]; ++i){
		hsize_t offset[2] = { i, 0 };
		hsize_t count[2] = { 1, dims[1] };
		hid_t tmp_memspace = H5Screate_simple(2, count, NULL);
		status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset,
			NULL, count, NULL);
		status = H5Dwrite(dataset_id, typeid(T) == typeid(int) ? H5T_NATIVE_INT : (typeid(T) == typeid(float) ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE), tmp_memspace, dataspace, H5P_DEFAULT, mem_data[i]);
	}
	//status = H5Dwrite(dataset_id, typeid(T) == typeid(int) ? H5T_NATIVE_INT : (typeid(T) == typeid(float) ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE), H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix[0]);

	// 关闭dataset相关对象
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace);

	// 关闭文件对象
	status = H5Fclose(file_id);
}