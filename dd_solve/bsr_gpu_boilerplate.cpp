//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "baseline_standalone.hpp"
#include "gen_graph.hpp"


void rocsparse_boilerplate_vc(solver_data_vc *bsr) 
{

    std::size_t d_bufferSize_M, d_bufferSize_L, d_bufferSize_U, d_bufferSize;
    int Nb = bsr->mb;
    int nnzbs_prec = bsr->nnzb;
    //Timer t;
    
    ROCSPARSE_CHECK(rocsparse_set_pointer_mode(bsr->handle, rocsparse_pointer_mode_host));
    
    ROCSPARSE_CHECK(rocsparse_create_mat_info(&bsr->ilu_info));
#if HIP_VERSION >= 50400000
    ROCSPARSE_CHECK(rocsparse_create_mat_info(&bsr->spmv_info));
#endif

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&bsr->descr_A));
    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&bsr->descr_M));
    
    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&bsr->descr_L));
    ROCSPARSE_CHECK(rocsparse_set_mat_fill_mode(bsr->descr_L, rocsparse_fill_mode_lower));
    ROCSPARSE_CHECK(rocsparse_set_mat_diag_type(bsr->descr_L, rocsparse_diag_type_unit));

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&bsr->descr_U));
    ROCSPARSE_CHECK(rocsparse_set_mat_fill_mode(bsr->descr_U, rocsparse_fill_mode_upper));
    ROCSPARSE_CHECK(rocsparse_set_mat_diag_type(bsr->descr_U, rocsparse_diag_type_non_unit));
    
    //ROCSPARSE_CHECK(rocsparse_dbsrilu0_buffer_size(bsr->handle, bsr->dir, Nb, nnzbs_prec,
    //                             bsr->descr_M, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, &d_bufferSize_M));
    
    //ROCSPARSE_CHECK(rocsparse_dbsrsv_buffer_size(bsr->handle, bsr->dir, bsr->operation, Nb, nnzbs_prec,
    //                           bsr->descr_L, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, &d_bufferSize_L));
    //ROCSPARSE_CHECK(rocsparse_dbsrsv_buffer_size(bsr->handle, bsr->dir, bsr->operation, Nb, nnzbs_prec,
    //                           bsr->descr_U, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, &d_bufferSize_U));

    //d_bufferSize = std::max(d_bufferSize_M, std::max(d_bufferSize_L, d_bufferSize_U));

    //HIP_CHECK(hipMalloc((void**)&bsr->d_buffer, d_bufferSize));

    // analysis of ilu LU decomposition
    //ROCSPARSE_CHECK(rocsparse_dbsrilu0_analysis(bsr->handle, bsr->dir, \
    //                           Nb, nnzbs_prec, bsr->descr_M, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, \
    //                           bsr->block_size, bsr->ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, bsr->d_buffer));

    int zero_position = 0;
    //rocsparse_status status = rocsparse_bsrilu0_zero_pivot(bsr->handle, bsr->ilu_info, &zero_position);
    //if (rocsparse_status_success != status) {
    //    printf("L has structural and/or numerical zero at L(%d,%d)\n", zero_position, zero_position);
    //   // return false;
    //    return;
    //}

    // analysis of ilu apply
    //ROCSPARSE_CHECK(rocsparse_dbsrsv_analysis(bsr->handle, bsr->dir, bsr->operation, \
    //                         Nb, nnzbs_prec, bsr->descr_L, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, \
    //                         bsr->block_size, bsr->ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, bsr->d_buffer));
    //ROCSPARSE_CHECK(rocsparse_dbsrsv_analysis(bsr->handle, bsr->dir, bsr->operation, \
    //                         Nb, nnzbs_prec, bsr->descr_U, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, \
    //                         bsr->block_size, bsr->ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, bsr->d_buffer));

#if HIP_VERSION >= 50400000
    std::cout<<"Analysis for rocpsarse dbsrmv"<<std::endl;
    ROCSPARSE_CHECK(rocsparse_dbsrmv_ex_analysis(bsr->handle, bsr->dir, bsr->operation, \
        Nb, Nb, bsr->nnzb,                                                              \
        bsr->descr_A, bsr->d_Avals, bsr->d_Arows, bsr->d_Acols,                         \
        bsr->block_size, bsr->spmv_info));
#endif

    /*
    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(stream));
        std::ostringstream out;
        out << "rocsparseSolver::analyze_matrix(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
    analysis_done = true;
    */
    //return true;

}
void rocsparse_boilerplate_g(solver_data2 *bsr) 
{

    std::size_t d_bufferSize_M, d_bufferSize_L, d_bufferSize_U, d_bufferSize;
    int Nb = bsr->mb;
    int nnzbs_prec = bsr->nnzb;
    //Timer t;
    
    ROCSPARSE_CHECK(rocsparse_set_pointer_mode(bsr->handle, rocsparse_pointer_mode_host));
    
    ROCSPARSE_CHECK(rocsparse_create_mat_info(&bsr->ilu_info));
#if HIP_VERSION >= 50400000
    ROCSPARSE_CHECK(rocsparse_create_mat_info(&bsr->spmv_info));
#endif

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&bsr->descr_A));
    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&bsr->descr_M));
    
    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&bsr->descr_L));
    ROCSPARSE_CHECK(rocsparse_set_mat_fill_mode(bsr->descr_L, rocsparse_fill_mode_lower));
    ROCSPARSE_CHECK(rocsparse_set_mat_diag_type(bsr->descr_L, rocsparse_diag_type_unit));

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&bsr->descr_U));
    ROCSPARSE_CHECK(rocsparse_set_mat_fill_mode(bsr->descr_U, rocsparse_fill_mode_upper));
    ROCSPARSE_CHECK(rocsparse_set_mat_diag_type(bsr->descr_U, rocsparse_diag_type_non_unit));
    
    //ROCSPARSE_CHECK(rocsparse_dbsrilu0_buffer_size(bsr->handle, bsr->dir, Nb, nnzbs_prec,
    //                             bsr->descr_M, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, &d_bufferSize_M));
    
    //ROCSPARSE_CHECK(rocsparse_dbsrsv_buffer_size(bsr->handle, bsr->dir, bsr->operation, Nb, nnzbs_prec,
    //                           bsr->descr_L, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, &d_bufferSize_L));
    //ROCSPARSE_CHECK(rocsparse_dbsrsv_buffer_size(bsr->handle, bsr->dir, bsr->operation, Nb, nnzbs_prec,
    //                           bsr->descr_U, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, &d_bufferSize_U));

    //d_bufferSize = std::max(d_bufferSize_M, std::max(d_bufferSize_L, d_bufferSize_U));

    //HIP_CHECK(hipMalloc((void**)&bsr->d_buffer, d_bufferSize));

    // analysis of ilu LU decomposition
    //ROCSPARSE_CHECK(rocsparse_dbsrilu0_analysis(bsr->handle, bsr->dir, \
    //                           Nb, nnzbs_prec, bsr->descr_M, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, \
    //                           bsr->block_size, bsr->ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, bsr->d_buffer));

    int zero_position = 0;
    //rocsparse_status status = rocsparse_bsrilu0_zero_pivot(bsr->handle, bsr->ilu_info, &zero_position);
    //if (rocsparse_status_success != status) {
    //    printf("L has structural and/or numerical zero at L(%d,%d)\n", zero_position, zero_position);
    //   // return false;
    //    return;
    //}

    // analysis of ilu apply
    //ROCSPARSE_CHECK(rocsparse_dbsrsv_analysis(bsr->handle, bsr->dir, bsr->operation, \
    //                         Nb, nnzbs_prec, bsr->descr_L, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, \
    //                         bsr->block_size, bsr->ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, bsr->d_buffer));
    //ROCSPARSE_CHECK(rocsparse_dbsrsv_analysis(bsr->handle, bsr->dir, bsr->operation, \
    //                         Nb, nnzbs_prec, bsr->descr_U, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, \
    //                         bsr->block_size, bsr->ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, bsr->d_buffer));

#if HIP_VERSION >= 50400000
    std::cout<<"Analysis for rocpsarse dbsrmv"<<std::endl;
    ROCSPARSE_CHECK(rocsparse_dbsrmv_ex_analysis(bsr->handle, bsr->dir, bsr->operation, \
        Nb, Nb, bsr->nnzb,                                                              \
        bsr->descr_A, bsr->d_Avals, bsr->d_Arows, bsr->d_Acols,                         \
        bsr->block_size, bsr->spmv_info));
#endif

    /*
    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(stream));
        std::ostringstream out;
        out << "rocsparseSolver::analyze_matrix(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
    analysis_done = true;
    */
    //return true;

}

void rocsparse_boilerplate(solver_data *bsr) 
{

    std::size_t d_bufferSize_M, d_bufferSize_L, d_bufferSize_U, d_bufferSize;
    int Nb = bsr->mb;
    int nnzbs_prec = bsr->nnzb;
    //Timer t;
    
    ROCSPARSE_CHECK(rocsparse_set_pointer_mode(bsr->handle, rocsparse_pointer_mode_host));
    
    ROCSPARSE_CHECK(rocsparse_create_mat_info(&bsr->ilu_info));
#if HIP_VERSION >= 50400000
    ROCSPARSE_CHECK(rocsparse_create_mat_info(&bsr->spmv_info));
#endif

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&bsr->descr_A));
    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&bsr->descr_M));
    
    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&bsr->descr_L));
    ROCSPARSE_CHECK(rocsparse_set_mat_fill_mode(bsr->descr_L, rocsparse_fill_mode_lower));
    ROCSPARSE_CHECK(rocsparse_set_mat_diag_type(bsr->descr_L, rocsparse_diag_type_unit));

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&bsr->descr_U));
    ROCSPARSE_CHECK(rocsparse_set_mat_fill_mode(bsr->descr_U, rocsparse_fill_mode_upper));
    ROCSPARSE_CHECK(rocsparse_set_mat_diag_type(bsr->descr_U, rocsparse_diag_type_non_unit));
    
    ROCSPARSE_CHECK(rocsparse_dbsrilu0_buffer_size(bsr->handle, bsr->dir, Nb, nnzbs_prec,
                                 bsr->descr_M, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, &d_bufferSize_M));
    
    ROCSPARSE_CHECK(rocsparse_dbsrsv_buffer_size(bsr->handle, bsr->dir, bsr->operation, Nb, nnzbs_prec,
                               bsr->descr_L, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, &d_bufferSize_L));
    ROCSPARSE_CHECK(rocsparse_dbsrsv_buffer_size(bsr->handle, bsr->dir, bsr->operation, Nb, nnzbs_prec,
                               bsr->descr_U, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, &d_bufferSize_U));

    d_bufferSize = std::max(d_bufferSize_M, std::max(d_bufferSize_L, d_bufferSize_U));

    HIP_CHECK(hipMalloc((void**)&bsr->d_buffer, d_bufferSize));

    // analysis of ilu LU decomposition
    ROCSPARSE_CHECK(rocsparse_dbsrilu0_analysis(bsr->handle, bsr->dir, \
                               Nb, nnzbs_prec, bsr->descr_M, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, \
                               bsr->block_size, bsr->ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, bsr->d_buffer));

    int zero_position = 0;
    rocsparse_status status = rocsparse_bsrilu0_zero_pivot(bsr->handle, bsr->ilu_info, &zero_position);
    if (rocsparse_status_success != status) {
        printf("L has structural and/or numerical zero at L(%d,%d)\n", zero_position, zero_position);
       // return false;
        return;
    }

    // analysis of ilu apply
    ROCSPARSE_CHECK(rocsparse_dbsrsv_analysis(bsr->handle, bsr->dir, bsr->operation, \
                             Nb, nnzbs_prec, bsr->descr_L, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, \
                             bsr->block_size, bsr->ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, bsr->d_buffer));
    ROCSPARSE_CHECK(rocsparse_dbsrsv_analysis(bsr->handle, bsr->dir, bsr->operation, \
                             Nb, nnzbs_prec, bsr->descr_U, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, \
                             bsr->block_size, bsr->ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, bsr->d_buffer));

#if HIP_VERSION >= 50400000
    std::cout<<"Analysis for rocpsarse dbsrmv"<<std::endl;
    ROCSPARSE_CHECK(rocsparse_dbsrmv_ex_analysis(bsr->handle, bsr->dir, bsr->operation, \
        Nb, Nb, bsr->nnzb,                                                              \
        bsr->descr_A, bsr->d_Avals, bsr->d_Arows, bsr->d_Acols,                         \
        bsr->block_size, bsr->spmv_info));
#endif

    /*
    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(stream));
        std::ostringstream out;
        out << "rocsparseSolver::analyze_matrix(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
    analysis_done = true;
    */
    //return true;

}




void data_xfers(solver_data *bsr) 
{

    int N = bsr->hbsr_bval.size(); //size of b
    int Nb = bsr->mb;
    int nnzb = bsr->nnzb;
    int nnz = bsr->nnz;
    
  
    hipStream_t stream;
    HIP_CHECK(hipStreamCreate (&stream));  
    HIP_CHECK(hipMemcpyAsync(bsr->d_Arows, bsr->hbsr_row_ptr.data(), sizeof(rocsparse_int) * (Nb + 1), hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_Acols, bsr->hbsr_col_ind.data(), sizeof(rocsparse_int) * nnzb, hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_Avals, bsr->hbsr_val.data(), sizeof(double) * nnz, hipMemcpyHostToDevice, stream));
    //compute ILU0 on duplicate of A so we retain the actual A in d_Avals
    //HIP_CHECK(hipMemcpyAsync(bsr->d_Mvals, bsr->d_Avals, sizeof(double) * nnz, hipMemcpyDeviceToDevice, stream));
    
    
    HIP_CHECK(hipMemcpyAsync(bsr->d_Mrows, bsr->hbsr_prow_ptr.data(), sizeof(rocsparse_int) *  bsr->hbsr_prow_ptr.size(), hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_Mcols, bsr->hbsr_pcol_ind.data(), sizeof(rocsparse_int) *  bsr->hbsr_pcol_ind.size(), hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_Mvals, bsr->hbsr_pval.data(), sizeof(double) *  bsr->hbsr_pval.size(), hipMemcpyHostToDevice, stream));
    
    HIP_CHECK(hipMemsetAsync(bsr->d_x, 0, sizeof(double) * N, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_b, bsr->hbsr_bval.data(), sizeof(double) * N, hipMemcpyHostToDevice, stream));
   
    
}

void data_xfers_g(solver_data2 *bsr) 
{

    int N = bsr->hbsr_bval.size(); //size of b
    int Nb = bsr->mb;
    int nnzb = bsr->nnzb;
    int nnz = bsr->nnz;
    
  
    hipStream_t stream;
    HIP_CHECK(hipStreamCreate (&stream));  
    HIP_CHECK(hipMemcpyAsync(bsr->d_Arows, bsr->hbsr_row_ptr.data(), sizeof(rocsparse_int) * (Nb + 1), hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_Acols, bsr->hbsr_col_ind.data(), sizeof(rocsparse_int) * nnzb, hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_Avals, bsr->hbsr_val.data(), sizeof(double) * nnz, hipMemcpyHostToDevice, stream));
    //compute ILU0 on duplicate of A so we retain the actual A in d_Avals
    //HIP_CHECK(hipMemcpyAsync(bsr->d_Mvals, bsr->d_Avals, sizeof(double) * nnz, hipMemcpyDeviceToDevice, stream));
    
    
    HIP_CHECK(hipMemcpyAsync(bsr->d_Mlrows, bsr->hbsr_plrow_ptr.data(), sizeof(rocsparse_int) *  bsr->hbsr_plrow_ptr.size(), hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_Mlcols, bsr->hbsr_plcol_ind.data(), sizeof(rocsparse_int) *  bsr->hbsr_plcol_ind.size(), hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_Mlvals, bsr->hbsr_plval.data(), sizeof(double) *  bsr->hbsr_plval.size(), hipMemcpyHostToDevice, stream));
    
    
    HIP_CHECK(hipMemcpyAsync(bsr->d_Murows, bsr->hbsr_purow_ptr.data(), sizeof(rocsparse_int) *  bsr->hbsr_purow_ptr.size(), hipMemcpyHostToDevice, stream));
    //std::cout<<" size of s"<<  bsr->hbsr_pucol_ind.size()<<"\n\n\n";
    HIP_CHECK(hipMemcpyAsync(bsr->d_Mucols, bsr->hbsr_pucol_ind.data(), sizeof(rocsparse_int) *  bsr->hbsr_pucol_ind.size(), hipMemcpyHostToDevice, stream));
    
    HIP_CHECK(hipMemcpyAsync(bsr->d_Muvals, bsr->hbsr_puval.data(), sizeof(double) *  bsr->hbsr_puval.size(), hipMemcpyHostToDevice, stream));
    
    
    
    HIP_CHECK(hipMemsetAsync(bsr->d_x, 0, sizeof(double) * N, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_b, bsr->hbsr_bval.data(), sizeof(double) * N, hipMemcpyHostToDevice, stream));
    
    
    hipMemcpy(bsr->dlevel_start_subdomain, bsr->level_start_subdomain.data(), sizeof(int) * bsr->level_start_subdomain.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dhdl2graph_rows, bsr->hdl2graph_rows.data(), sizeof(int) * bsr->hdl2graph_rows.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dhdl2graph_cols,bsr->hdl2graph_cols.data(), sizeof(int) * bsr->hdl2graph_cols.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dhdl2graph_levptr, bsr->hdl2graph_levptr.data(), sizeof(int) * bsr->hdl2graph_levptr.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dhdl2graph_data_offsets, bsr->hdl2graph_data_offsets.data(), sizeof(int) * bsr->hdl2graph_data_offsets.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dllevels_per_subdomain, bsr->llevels_per_subdomain.data(), sizeof(int) * bsr->llevels_per_subdomain.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dsubdomain_offsets, bsr->subdomain_offsets.data(), sizeof(int) * bsr->subdomain_offsets.size(), hipMemcpyHostToDevice);
    
    
    hipMemcpy(bsr->dlevelu_start_subdomain, bsr->levelu_start_subdomain.data(), sizeof(int) * bsr->levelu_start_subdomain.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dhdu2graph_rows, bsr->hdu2graph_rows.data(), sizeof(int) * bsr->hdu2graph_rows.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dhdu2graph_cols,bsr->hdu2graph_cols.data(), sizeof(int) * bsr->hdu2graph_cols.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dhdu2graph_levptr, bsr->hdu2graph_levptr.data(), sizeof(int) * bsr->hdu2graph_levptr.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dhdu2graph_data_offsets, bsr->hdu2graph_data_offsets.data(), sizeof(int) * bsr->hdu2graph_data_offsets.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dllevelsu_per_subdomain, bsr->llevelsu_per_subdomain.data(), sizeof(int) * bsr->llevelsu_per_subdomain.size(), hipMemcpyHostToDevice);
   
    hipMemcpy(bsr->dval, bsr->hdval.data(), sizeof(double) * bsr->hdval.size(), hipMemcpyHostToDevice);
    
    
}


void data_xfers_vc(solver_data_vc *bsr) 
{

    int N = bsr->hbsr_bval.size(); //size of b
    int Nb = bsr->mb;
    int nnzb = bsr->nnzb;
    int nnz = bsr->nnz;
    
  
    hipStream_t stream;
    HIP_CHECK(hipStreamCreate (&stream));  
    HIP_CHECK(hipMemcpyAsync(bsr->d_Arows, bsr->hbsr_row_ptr.data(), sizeof(rocsparse_int) * (Nb + 1), hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_Acols, bsr->hbsr_col_ind.data(), sizeof(rocsparse_int) * nnzb, hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_Avals, bsr->hbsr_val.data(), sizeof(double) * nnz, hipMemcpyHostToDevice, stream));
    //compute ILU0 on duplicate of A so we retain the actual A in d_Avals
    //HIP_CHECK(hipMemcpyAsync(bsr->d_Mvals, bsr->d_Avals, sizeof(double) * nnz, hipMemcpyDeviceToDevice, stream));
    
    
    HIP_CHECK(hipMemcpyAsync(bsr->d_Mlrows, bsr->hbsr_plrow_ptr.data(), sizeof(rocsparse_int) *  bsr->hbsr_plrow_ptr.size(), hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_Mlcols, bsr->hbsr_plcol_ind.data(), sizeof(rocsparse_int) *  bsr->hbsr_plcol_ind.size(), hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_Mlvals, bsr->hbsr_plval.data(), sizeof(double) *  bsr->hbsr_plval.size(), hipMemcpyHostToDevice, stream));
    
    
    HIP_CHECK(hipMemcpyAsync(bsr->d_Murows, bsr->hbsr_purow_ptr.data(), sizeof(rocsparse_int) *  bsr->hbsr_purow_ptr.size(), hipMemcpyHostToDevice, stream));
    //std::cout<<" size of s"<<  bsr->hbsr_pucol_ind.size()<<"\n\n\n";
    HIP_CHECK(hipMemcpyAsync(bsr->d_Mucols, bsr->hbsr_pucol_ind.data(), sizeof(rocsparse_int) *  bsr->hbsr_pucol_ind.size(), hipMemcpyHostToDevice, stream));
    
    HIP_CHECK(hipMemcpyAsync(bsr->d_Muvals, bsr->hbsr_puval.data(), sizeof(double) *  bsr->hbsr_puval.size(), hipMemcpyHostToDevice, stream));
    
    
    
    HIP_CHECK(hipMemsetAsync(bsr->d_x, 0, sizeof(double) * N, stream));
    HIP_CHECK(hipMemcpyAsync(bsr->d_b, bsr->hbsr_bval.data(), sizeof(double) * N, hipMemcpyHostToDevice, stream));
    
    
    hipMemcpy(bsr->dlevel_start_subdomain, bsr->level_start_subdomain.data(), sizeof(int) * bsr->level_start_subdomain.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dhdl2graph_rows, bsr->hdl2graph_rows.data(), sizeof(int) * bsr->hdl2graph_rows.size(), hipMemcpyHostToDevice);
   // hipMemcpy(bsr->dhdl2graph_cols,bsr->hdl2graph_cols.data(), sizeof(int) * bsr->hdl2graph_cols.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dhdl2graph_levptr, bsr->hdl2graph_levptr.data(), sizeof(int) * bsr->hdl2graph_levptr.size(), hipMemcpyHostToDevice);
   // hipMemcpy(bsr->dhdl2graph_data_offsets, bsr->hdl2graph_data_offsets.data(), sizeof(int) * bsr->hdl2graph_data_offsets.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dllevels_per_subdomain, bsr->llevels_per_subdomain.data(), sizeof(int) * bsr->llevels_per_subdomain.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dsubdomain_offsets, bsr->subdomain_offsets.data(), sizeof(int) * bsr->subdomain_offsets.size(), hipMemcpyHostToDevice);
    
    
    hipMemcpy(bsr->dlevelu_start_subdomain, bsr->levelu_start_subdomain.data(), sizeof(int) * bsr->levelu_start_subdomain.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dhdu2graph_rows, bsr->hdu2graph_rows.data(), sizeof(int) * bsr->hdu2graph_rows.size(), hipMemcpyHostToDevice);
    //hipMemcpy(bsr->dhdu2graph_cols,bsr->hdu2graph_cols.data(), sizeof(int) * bsr->hdu2graph_cols.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dhdu2graph_levptr, bsr->hdu2graph_levptr.data(), sizeof(int) * bsr->hdu2graph_levptr.size(), hipMemcpyHostToDevice);
   // hipMemcpy(bsr->dhdu2graph_data_offsets, bsr->hdu2graph_data_offsets.data(), sizeof(int) * bsr->hdu2graph_data_offsets.size(), hipMemcpyHostToDevice);
    hipMemcpy(bsr->dllevelsu_per_subdomain, bsr->llevelsu_per_subdomain.data(), sizeof(int) * bsr->llevelsu_per_subdomain.size(), hipMemcpyHostToDevice);
   
    hipMemcpy(bsr->dval, bsr->hdval.data(), sizeof(double) * bsr->hdval.size(), hipMemcpyHostToDevice);
    
    
}



void init_gpu(solver_data *bsr) 
{

    int N = bsr->hbsr_bval.size(); //size of b
    int Nb = bsr->mb;
    int nnzb = bsr->nnzb;
    int nnz = bsr->nnz;
    int block_size = bsr->block_size;
    
    
    int nnzb_prec = bsr->hbsr_pcol_ind.size();
    //std::cout << " pcolid size is "<<nnzb_prec<<std::endl;  
    //exit(0);
    HIP_CHECK(hipMalloc((void**)&bsr->d_r, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_rw, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_p, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_pw, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_s, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_t, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_v, sizeof(double) * N));

    HIP_CHECK(hipMalloc((void**)&bsr->d_Arows, sizeof(rocsparse_int) * (Nb + 1)));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Acols, sizeof(rocsparse_int) * nnzb));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Avals, sizeof(double) * nnz));
    HIP_CHECK(hipMalloc((void**)&bsr->d_x, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_b, sizeof(double) * N));
    //structural symmetry allows the same row pointer and column index to be used for ILU0 output
    //statement above is longer truetrue anymore as the preconditioner and A are no longer same - 
    //HIP_CHECK(hipMalloc((void**)&bsr->d_Mvals, sizeof(double) * nnzb_prec * block_size * block_size));
    
    //Create the preconditioner matrix to device for generating ILU0 decomposition 
    //bsr->d_Mcols = bsr->d_Acols;
    //bsr->d_Mrows = bsr->d_Arows;
    
    HIP_CHECK(hipMalloc((void**)&bsr->d_Mrows, sizeof(rocsparse_int) * (Nb + 1)));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Mcols, sizeof(rocsparse_int) * nnzb_prec));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Mvals, sizeof(double) * nnzb_prec*block_size*block_size));
    
    
    
    //Init handles
    ROCSPARSE_CHECK(rocsparse_create_handle(&bsr->handle));
    ROCBLAS_CHECK(rocblas_create_handle(&bsr->blas_handle));

    ROCSPARSE_CHECK(rocsparse_get_version(bsr->handle, &bsr->ver));
    ROCSPARSE_CHECK(rocsparse_get_git_rev(bsr->handle, bsr->rev));

    std::ostringstream out;
    out << "rocSPARSE version: " << bsr->ver / 100000 << "." << bsr->ver / 100 % 1000 << "."
        << bsr->ver % 100 << "-" << bsr->rev << "\n";
    std::cout<<out.str();
    
}



void init_gpu_g(solver_data2 *bsr) 
{

    int N = bsr->hbsr_bval.size(); //size of b
    int Nb = bsr->mb;
    int nnzb = bsr->nnzb;
    int nnz = bsr->nnz;
    int block_size = bsr->block_size;
    
      //std::cout<< bsr->subdomain_offsets.size()<<" \n\n\n";
    
    //for(int kk=0;kk<bsr->subdomain_offsets.size();kk++)std::cout<< bsr->subdomain_offsets[kk]<<" ";
    
    dependency_lgraph_coo_offsets(
     bsr->mb, 
     bsr->nnzb, 
     bsr->num_subdomains, 
     bsr->subdomain_offsets,   //input   //used to get start row and end of of subdomain. Size = #subdomains + 1 
     bsr->hmap2, 				//input     
     bsr->hbsr_plrow_ptr, 		//input     
     bsr->hbsr_plcol_ind, 		//input     
     bsr->llevels_per_subdomain,	//output #levels per each subdomain
     bsr->hdl2graph_rows,            //output 
     bsr->hdl2graph_cols, 			//output
     bsr->hdl2graph_levptr, 			//output #point to rows and cols for each level
     bsr->hdl2graph_data_offsets,    //output #point to the vector values
     bsr->level_start_subdomain);
  
  
     //std::cout<<"lalalalala num subds"<<  bsr->num_subdomains<<std::endl;
     //for (int i=0;i<  bsr->llevels_per_subdomain.size();i++)std::cout<<"levels in subdomain "<<i<<" "<<  bsr->llevels_per_subdomain[i]<<std::endl;
  
    dependency_ugraph_coo_offsets(
     bsr->mb, 
     bsr->nnzb, 
     bsr->num_subdomains, 
     bsr->subdomain_offsets,   //input   //used to get start row and end of of subdomain. Size = #subdomains + 1 
     bsr->hmap3, 				//input     
     bsr->hbsr_purow_ptr, 		//input     
     bsr->hbsr_pucol_ind, 		//input     
     bsr->llevelsu_per_subdomain,	//output #levels per each subdomain
     bsr->hdu2graph_rows,            //output 
     bsr->hdu2graph_cols, 			//output
     bsr->hdu2graph_levptr, 			//output #point to rows and cols for each level
     bsr->hdu2graph_data_offsets,    //output #point to the vector values
     bsr->levelu_start_subdomain);
    
    
    
    //std::cout<<"size size \n\n\n\n"<<bsr->hdu2graph_rows.size()<<std::endl;
    
    
    
    //int nnzb_prec = bsr->hbsr_col_ind.size();
    //std::cout << " pcolid size is "<<nnzb_prec<<std::endl;  
    //exit(0);
    HIP_CHECK(hipMalloc((void**)&bsr->d_r, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_rw, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_p, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_pw, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_s, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_t, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_v, sizeof(double) * N));

    HIP_CHECK(hipMalloc((void**)&bsr->d_Arows, sizeof(rocsparse_int) * (Nb + 1)));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Acols, sizeof(rocsparse_int) * nnzb));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Avals, sizeof(double) * nnz));
    HIP_CHECK(hipMalloc((void**)&bsr->d_x, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_b, sizeof(double) * N));
    //structural symmetry allows the same row pointer and column index to be used for ILU0 output
    //statement above is longer truetrue anymore as the preconditioner and A are no longer same - 
    //HIP_CHECK(hipMalloc((void**)&bsr->d_Mvals, sizeof(double) * nnzb_prec * block_size * block_size));
    
    //Create the preconditioner matrix to device for generating ILU0 decomposition 
    //bsr->d_Mcols = bsr->d_Acols;
    //bsr->d_Mrows = bsr->d_Arows;
    
    HIP_CHECK(hipMalloc((void**)&bsr->d_Mlrows, sizeof(rocsparse_int) * bsr->hbsr_plrow_ptr.size()));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Mlcols, sizeof(rocsparse_int) * bsr->hbsr_plcol_ind.size()));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Mlvals, sizeof(double) * bsr->hbsr_plval.size()));
    
    HIP_CHECK(hipMalloc((void**)&bsr->d_Murows, sizeof(rocsparse_int) * bsr->hbsr_purow_ptr.size()));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Mucols, sizeof(rocsparse_int) * bsr->hbsr_pucol_ind.size()));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Muvals, sizeof(double) * bsr->hbsr_puval.size()));
    
    hipMalloc((void**)&bsr->dval, sizeof(double) * bsr->hdval.size());
    
    hipMalloc((void**)&bsr->dlevel_start_subdomain, sizeof(int) * bsr->level_start_subdomain.size());
    hipMalloc((void**)&bsr->dhdl2graph_rows, sizeof(int) * bsr->hdl2graph_rows.size());
    hipMalloc((void**)&bsr->dhdl2graph_cols, sizeof(int) * bsr->hdl2graph_cols.size());
    hipMalloc((void**)&bsr->dhdl2graph_levptr, sizeof(int) * bsr->hdl2graph_levptr.size());
    hipMalloc((void**)&bsr->dhdl2graph_data_offsets, sizeof(int) * bsr->hdl2graph_data_offsets.size());
    hipMalloc((void**)&bsr->dllevels_per_subdomain, sizeof(int) * bsr->llevels_per_subdomain.size());
    hipMalloc((void**)&bsr->dsubdomain_offsets, sizeof(int) * bsr->subdomain_offsets.size());
    
    hipMalloc((void**)&bsr->dlevelu_start_subdomain, sizeof(int) * bsr->levelu_start_subdomain.size());
    hipMalloc((void**)&bsr->dhdu2graph_rows, sizeof(int) * bsr->hdu2graph_rows.size());
    hipMalloc((void**)&bsr->dhdu2graph_cols, sizeof(int) * bsr->hdu2graph_cols.size());
    hipMalloc((void**)&bsr->dhdu2graph_levptr, sizeof(int) * bsr->hdu2graph_levptr.size());
    hipMalloc((void**)&bsr->dhdu2graph_data_offsets, sizeof(int) * bsr->hdu2graph_data_offsets.size());
    hipMalloc((void**)&bsr->dllevelsu_per_subdomain, sizeof(int) * bsr->llevelsu_per_subdomain.size());
   
    
    
    
    //Init handles
    ROCSPARSE_CHECK(rocsparse_create_handle(&bsr->handle));
    ROCBLAS_CHECK(rocblas_create_handle(&bsr->blas_handle));

    ROCSPARSE_CHECK(rocsparse_get_version(bsr->handle, &bsr->ver));
    ROCSPARSE_CHECK(rocsparse_get_git_rev(bsr->handle, bsr->rev));

    std::ostringstream out;
    out << "rocSPARSE version: " << bsr->ver / 100000 << "." << bsr->ver / 100 % 1000 << "."
        << bsr->ver % 100 << "-" << bsr->rev << "\n";
    std::cout<<out.str();
    
}


void init_gpu_vc(solver_data_vc *bsr) 
{

    int N = bsr->hbsr_bval.size(); //size of b
    int Nb = bsr->mb;
    int nnzb = bsr->nnzb;
    int nnz = bsr->nnz;
    int block_size = bsr->block_size;
    
      //std::cout<< bsr->subdomain_offsets.size()<<" \n\n\n";
    
    //for(int kk=0;kk<bsr->subdomain_offsets.size();kk++)std::cout<< bsr->subdomain_offsets[kk]<<" ";
    
    dependency_lgraph_vc_offsets(
     bsr->mb, 
     bsr->nnzb, 
     bsr->num_subdomains, 
     bsr->subdomain_offsets,   //input   //used to get start row and end of of subdomain. Size = #subdomains + 1 
     bsr->hmap2, 				//input     
     bsr->hbsr_plrow_ptr, 		//input     
     bsr->hbsr_plcol_ind, 		//input     
     bsr->llevels_per_subdomain,	//output #levels per each subdomain
     bsr->hdl2graph_rows,            //output 
     bsr->hdl2graph_levptr, 			//output #point to rows and cols for each level
     bsr->level_start_subdomain);
  
  
     //std::cout<<"lalalalala num subds"<<  bsr->num_subdomains<<std::endl;
     //for (int i=0;i<  bsr->llevels_per_subdomain.size();i++)std::cout<<"levels in subdomain "<<i<<" "<<  bsr->llevels_per_subdomain[i]<<std::endl;
  
    dependency_ugraph_vc_offsets(
     bsr->mb, 
     bsr->nnzb, 
     bsr->num_subdomains, 
     bsr->subdomain_offsets,   //input   //used to get start row and end of of subdomain. Size = #subdomains + 1 
     bsr->hmap3, 				//input     
     bsr->hbsr_purow_ptr, 		//input     
     bsr->hbsr_pucol_ind, 		//input     
     bsr->llevelsu_per_subdomain,	//output #levels per each subdomain
     bsr->hdu2graph_rows,            //output 
     bsr->hdu2graph_levptr, 			//output #point to rows and cols for each level
     bsr->levelu_start_subdomain);
    
    
    
    //std::cout<<"size size \n\n\n\n"<<bsr->hdu2graph_rows.size()<<std::endl;
    
    
    
    //int nnzb_prec = bsr->hbsr_col_ind.size();
    //std::cout << " pcolid size is "<<nnzb_prec<<std::endl;  
    //exit(0);
    HIP_CHECK(hipMalloc((void**)&bsr->d_r, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_rw, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_p, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_pw, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_s, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_t, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_v, sizeof(double) * N));

    HIP_CHECK(hipMalloc((void**)&bsr->d_Arows, sizeof(rocsparse_int) * (Nb + 1)));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Acols, sizeof(rocsparse_int) * nnzb));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Avals, sizeof(double) * nnz));
    HIP_CHECK(hipMalloc((void**)&bsr->d_x, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&bsr->d_b, sizeof(double) * N));
    //structural symmetry allows the same row pointer and column index to be used for ILU0 output
    //statement above is longer truetrue anymore as the preconditioner and A are no longer same - 
    //HIP_CHECK(hipMalloc((void**)&bsr->d_Mvals, sizeof(double) * nnzb_prec * block_size * block_size));
    
    //Create the preconditioner matrix to device for generating ILU0 decomposition 
    //bsr->d_Mcols = bsr->d_Acols;
    //bsr->d_Mrows = bsr->d_Arows;
    
    HIP_CHECK(hipMalloc((void**)&bsr->d_Mlrows, sizeof(rocsparse_int) * bsr->hbsr_plrow_ptr.size()));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Mlcols, sizeof(rocsparse_int) * bsr->hbsr_plcol_ind.size()));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Mlvals, sizeof(double) * bsr->hbsr_plval.size()));
    
    HIP_CHECK(hipMalloc((void**)&bsr->d_Murows, sizeof(rocsparse_int) * bsr->hbsr_purow_ptr.size()));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Mucols, sizeof(rocsparse_int) * bsr->hbsr_pucol_ind.size()));
    HIP_CHECK(hipMalloc((void**)&bsr->d_Muvals, sizeof(double) * bsr->hbsr_puval.size()));
    
    hipMalloc((void**)&bsr->dval, sizeof(double) * bsr->hdval.size());
    
    hipMalloc((void**)&bsr->dlevel_start_subdomain, sizeof(int) * bsr->level_start_subdomain.size());
    hipMalloc((void**)&bsr->dhdl2graph_rows, sizeof(int) * bsr->hdl2graph_rows.size());
   // hipMalloc((void**)&bsr->dhdl2graph_cols, sizeof(int) * bsr->hdl2graph_cols.size());
    hipMalloc((void**)&bsr->dhdl2graph_levptr, sizeof(int) * bsr->hdl2graph_levptr.size());
    //hipMalloc((void**)&bsr->dhdl2graph_data_offsets, sizeof(int) * bsr->hdl2graph_data_offsets.size());
    hipMalloc((void**)&bsr->dllevels_per_subdomain, sizeof(int) * bsr->llevels_per_subdomain.size());
    hipMalloc((void**)&bsr->dsubdomain_offsets, sizeof(int) * bsr->subdomain_offsets.size());
    
    hipMalloc((void**)&bsr->dlevelu_start_subdomain, sizeof(int) * bsr->levelu_start_subdomain.size());
    hipMalloc((void**)&bsr->dhdu2graph_rows, sizeof(int) * bsr->hdu2graph_rows.size());
   // hipMalloc((void**)&bsr->dhdu2graph_cols, sizeof(int) * bsr->hdu2graph_cols.size());
    hipMalloc((void**)&bsr->dhdu2graph_levptr, sizeof(int) * bsr->hdu2graph_levptr.size());
    //hipMalloc((void**)&bsr->dhdu2graph_data_offsets, sizeof(int) * bsr->hdu2graph_data_offsets.size());
    hipMalloc((void**)&bsr->dllevelsu_per_subdomain, sizeof(int) * bsr->llevelsu_per_subdomain.size());
   
    
    
    
    //Init handles
    ROCSPARSE_CHECK(rocsparse_create_handle(&bsr->handle));
    ROCBLAS_CHECK(rocblas_create_handle(&bsr->blas_handle));

    ROCSPARSE_CHECK(rocsparse_get_version(bsr->handle, &bsr->ver));
    ROCSPARSE_CHECK(rocsparse_get_git_rev(bsr->handle, bsr->rev));

    std::ostringstream out;
    out << "rocSPARSE version: " << bsr->ver / 100000 << "." << bsr->ver / 100 % 1000 << "."
        << bsr->ver % 100 << "-" << bsr->rev << "\n";
    std::cout<<out.str();
    
}



