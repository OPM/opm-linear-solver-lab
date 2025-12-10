//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#pragma once
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include "hip/hip_runtime.h"
#include <rocsparse/rocsparse.h>
#include <vector>
#include <rocblas/rocblas.h>
#include <hip/hip_version.h>
#include <bits/stdc++.h> 
#include <fstream>
#include <iterator>
#include <string>
#include <bits/stdc++.h> 
#include <sstream>


#define HIP_CHECK(STAT)                                  \
    do {                                                 \
        const hipError_t stat = (STAT);                  \
        if(stat != hipSuccess)                           \
        {                                                \
            std::ostringstream oss;                      \
            oss << "rocsparseSolverBackend::hip ";       \
            oss << "error: " << hipGetErrorString(stat); \
            std::cout<<oss.str()<<std::endl  ;           \
        }                                                \
    } while(0)

#define ROCSPARSE_CHECK(STAT)                            \
    do {                                                 \
        const rocsparse_status stat = (STAT);            \
        if(stat != rocsparse_status_success)             \
        {                                                \
            std::ostringstream oss;                      \
            oss << "rocsparseSolverBackend::rocsparse "; \
            oss << "error: " << stat ;                   \
            std::cout<<oss.str()<<std::endl   ;          \
        }                                                \
    } while(0)

#define ROCBLAS_CHECK(STAT)                              \
    do {                                                 \
        const rocblas_status stat = (STAT);              \
        if(stat != rocblas_status_success)               \
        {                                                \
            std::ostringstream oss;                      \
            oss << "rocsparseSolverBackend::rocblas ";   \
            oss << "error: " << stat;                    \
            std::cout<<oss.str()<<std::endl   ;          \
        }                                                \
    } while(0)





struct solver_data{


    int mb; //rows
    int block_size;
    int nnzb; //nnz blocks
    int nnz; // total nnz = nnzb * block_size^2
    int verbosity=1;
    int Nb; //rows
    int nnzbs_prec;
    int N; // size of b
   // int num_subdomains;
   // int rows_per_subdomain;
    
    //bsr on host;
   	std::vector<double> hbsr_val;
  	std::vector<int>    hbsr_row_ptr;
  	std::vector<int>    hbsr_col_ind;
    std::vector<double> hbsr_bval;
    
    //bsr preconditioner on host 
    std::vector<double> hbsr_pval;
  	std::vector<int>    hbsr_prow_ptr;
  	std::vector<int>    hbsr_pcol_ind;
    
    /*
    //graph_vectors
    std::vector <int> subdomain_offsets;
    std::vector <int> llevels_per_subdomain;
    std::vector <int> level_start_subdomain;
    std::vector <int> hmap2;
    std::vector <int> llevelsu_per_subdomain;
    std::vector <int> levelu_start_subdomain;
    std::vector <int> hmap3;
      
    //graph data
    std::vector<int> hdl2graph_rows, hdl2graph_cols, hdl2graph_levptr, hdl2graph_data_offsets;
    std::vector<int> hdu2graph_rows, hdu2graph_cols, hdu2graph_levptr, hdu2graph_data_offsets;
 
    
    //lower and upper vals
    double* bsrl_val;
    double* bsru_val;
    double* dval;
    
    //device pointers 
    int *dlevel_start_subdomain, *dhdl2graph_rows, *dhdl2graph_cols, *dhdl2graph_levptr, *dhdl2graph_data_offsets, *dllevels_per_subdomain, *dsubdomain_offsets;
    int *dlevelu_start_subdomain, *dhdu2graph_rows, *dhdu2graph_cols, *dhdu2graph_levptr, *dhdu2graph_data_offsets, *dllevelsu_per_subdomain;
    */
    rocsparse_direction dir = rocsparse_direction_row;
    rocsparse_operation operation = rocsparse_operation_none;
    rocsparse_handle handle;
    rocblas_handle blas_handle;
    rocsparse_mat_descr descr_A, descr_M, descr_L,descr_U;
    rocsparse_mat_info ilu_info;
#if HIP_VERSION >= 50400000
    rocsparse_mat_info spmv_info;
#endif
    hipStream_t stream;

    rocsparse_int *d_Arows, *d_Mrows;
    rocsparse_int *d_Acols, *d_Mcols;
    double *d_Avals, *d_Mvals;
    double *d_x, *d_b, *d_r, *d_rw, *d_p;     // vectors, used during linear solve
    double *d_pw, *d_s, *d_t, *d_v;
    void *d_buffer; // buffer space, used by rocsparse ilu0 analysis
    int  ver;
    char rev[64]; 
    
};



struct solver_data2{


    int mb; //rows
    int block_size;
    int nnzb; //nnz blocks
    int nnz; // total nnz = nnzb * block_size^2
    int verbosity=1;
    int Nb; //rows
    int nnzbs_prec;
    int N; // size of b
    int num_subdomains;
    int rows_per_subdomain;
    
    //bsr on host;
   	std::vector<double> hbsr_val;
  	std::vector<int>    hbsr_row_ptr;
  	std::vector<int>    hbsr_col_ind;
    std::vector<double> hbsr_bval;
    
    //bsr preconditioner on host 
    std::vector<double> hbsr_plval;
  	std::vector<int>    hbsr_plrow_ptr;
  	std::vector<int>    hbsr_plcol_ind;
   
    std::vector<double> hbsr_puval;
  	std::vector<int>    hbsr_purow_ptr;
  	std::vector<int>    hbsr_pucol_ind;
   
    //Pc diag
    std::vector <double> hdval;
    
    //graph_vectors
    std::vector <int> subdomain_offsets;
    std::vector <int> llevels_per_subdomain;
    std::vector <int> level_start_subdomain;
    std::vector <int> hmap2;
    std::vector <int> llevelsu_per_subdomain;
    std::vector <int> levelu_start_subdomain;
    std::vector <int> hmap3;
      
    //graph data
    std::vector<int> hdl2graph_rows, hdl2graph_cols, hdl2graph_levptr, hdl2graph_data_offsets;
    std::vector<int> hdu2graph_rows, hdu2graph_cols, hdu2graph_levptr, hdu2graph_data_offsets;
 
    
    //lower and upper vals
    double* bsrl_val;
    double* bsru_val;
    double* dval;
    
    //device pointers 
    int *dlevel_start_subdomain, *dhdl2graph_rows, *dhdl2graph_cols, *dhdl2graph_levptr, *dhdl2graph_data_offsets, *dllevels_per_subdomain, *dsubdomain_offsets;
    int *dlevelu_start_subdomain, *dhdu2graph_rows, *dhdu2graph_cols, *dhdu2graph_levptr, *dhdu2graph_data_offsets, *dllevelsu_per_subdomain;
    
    rocsparse_direction dir = rocsparse_direction_row;
    rocsparse_operation operation = rocsparse_operation_none;
    rocsparse_handle handle;
    rocblas_handle blas_handle;
    rocsparse_mat_descr descr_A, descr_M, descr_L,descr_U;
    rocsparse_mat_info ilu_info;
#if HIP_VERSION >= 50400000
    rocsparse_mat_info spmv_info;
#endif
    hipStream_t stream;

    rocsparse_int *d_Arows, *d_Mlrows, *d_Murows;
    rocsparse_int *d_Acols, *d_Mlcols, *d_Mucols;
    double *d_Avals, *d_Mlvals, *d_Muvals;
    double *d_x, *d_b, *d_r, *d_rw, *d_p;     // vectors, used during linear solve
    double *d_pw, *d_s, *d_t, *d_v;
    void *d_buffer; // buffer space, used by rocsparse ilu0 analysis
    int  ver;
    char rev[64]; 
    
};

struct solver_data_vc{


    int mb; //rows
    int block_size;
    int nnzb; //nnz blocks
    int nnz; // total nnz = nnzb * block_size^2
    int verbosity=1;
    int Nb; //rows
    int nnzbs_prec;
    int N; // size of b
    int num_subdomains;
    int rows_per_subdomain;
    
    //bsr on host;
   	std::vector<double> hbsr_val;
  	std::vector<int>    hbsr_row_ptr;
  	std::vector<int>    hbsr_col_ind;
    std::vector<double> hbsr_bval;
    
    //bsr preconditioner on host 
    std::vector<double> hbsr_plval;
  	std::vector<int>    hbsr_plrow_ptr;
  	std::vector<int>    hbsr_plcol_ind;
   
    std::vector<double> hbsr_puval;
  	std::vector<int>    hbsr_purow_ptr;
  	std::vector<int>    hbsr_pucol_ind;
   
    //Pc diag
    std::vector <double> hdval;
    
    //graph_vectors
    std::vector <int> subdomain_offsets;
    std::vector <int> llevels_per_subdomain;
    std::vector <int> level_start_subdomain;
    std::vector <int> hmap2;
    std::vector <int> llevelsu_per_subdomain;
    std::vector <int> levelu_start_subdomain;
    std::vector <int> hmap3;
      
    //graph data
    std::vector<int> hdl2graph_rows,  hdl2graph_levptr;
    std::vector<int> hdu2graph_rows,  hdu2graph_levptr;
 
    
    //lower and upper vals
    double* bsrl_val;
    double* bsru_val;
    double* dval;
    
    //device pointers 
    int *dlevel_start_subdomain, *dhdl2graph_rows,  *dhdl2graph_levptr, *dllevels_per_subdomain, *dsubdomain_offsets;
    int *dlevelu_start_subdomain, *dhdu2graph_rows,  *dhdu2graph_levptr,  *dllevelsu_per_subdomain;
    
    rocsparse_direction dir = rocsparse_direction_row;
    rocsparse_operation operation = rocsparse_operation_none;
    rocsparse_handle handle;
    rocblas_handle blas_handle;
    rocsparse_mat_descr descr_A, descr_M, descr_L,descr_U;
    rocsparse_mat_info ilu_info;
#if HIP_VERSION >= 50400000
    rocsparse_mat_info spmv_info;
#endif
    hipStream_t stream;

    rocsparse_int *d_Arows, *d_Mlrows, *d_Murows;
    rocsparse_int *d_Acols, *d_Mlcols, *d_Mucols;
    double *d_Avals, *d_Mlvals, *d_Muvals;
    double *d_x, *d_b, *d_r, *d_rw, *d_p;     // vectors, used during linear solve
    double *d_pw, *d_s, *d_t, *d_v;
    void *d_buffer; // buffer space, used by rocsparse ilu0 analysis
    int  ver;
    char rev[64]; 
    
};
