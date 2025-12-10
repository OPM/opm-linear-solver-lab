//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "baseline_standalone.hpp"

void preconditioner(solver_data *bsr)
{


    //ILU0 decomposition
    ROCSPARSE_CHECK(rocsparse_dbsrilu0(bsr->handle, bsr->dir, bsr->Nb, bsr->nnzbs_prec, bsr->descr_M,
                    bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, rocsparse_solve_policy_auto, bsr->d_buffer));

    // Check for zero pivot
    int zero_position = 0;
    rocsparse_status status = rocsparse_bsrilu0_zero_pivot(bsr->handle, bsr->ilu_info, &zero_position);
    if(rocsparse_status_success != status)
    {
        printf("L has structural and/or numerical zero at L(%d,%d)\n", zero_position, zero_position);
        //return false;
    }
    // 
    ///*
    std::vector<double> hilu0(bsr->nnz, 0.0f);
    hipMemcpy(hilu0.data(), bsr->d_Mvals, sizeof(double) * bsr->nnz, hipMemcpyDeviceToHost);
  
  
    std::cout << "ilu0 bsr " << std::endl;
    //for(size_t i = 0; i < hilu0.size(); ++i)
    //{
    //    std::cout << hilu0[i] << ", ";
   // }
    std::cout << std::endl;
    //*/
    /*
    if (verbosity >= 3) {
        //HIP_CHECK(hipStreamSynchronize(stream));
        hipDeviceSynchronize();
        std::ostringstream out;
        out << "rocsparseSolver::create_preconditioner(): " << t.stop() << " s";
        //OpmLog::info(out.str());
    }
    */
}

void preconditioner_g(solver_data2 *bsr)
{


    //ILU0 decomposition
    //ROCSPARSE_CHECK(rocsparse_dbsrilu0(bsr->handle, bsr->dir, bsr->Nb, bsr->nnzbs_prec, bsr->descr_M,
    //                bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, rocsparse_solve_policy_auto, bsr->d_buffer));

    // Check for zero pivot
    int zero_position = 0;
    //rocsparse_status status = rocsparse_bsrilu0_zero_pivot(bsr->handle, bsr->ilu_info, &zero_position);
    //if(rocsparse_status_success != status)
    {
    //    printf("L has structural and/or numerical zero at L(%d,%d)\n", zero_position, zero_position);
        //return false;
    }
    // 
    ///*
    std::vector<double> hilu0(bsr->nnz, 0.0f);
    //hipMemcpy(hilu0.data(), bsr->d_Mvals, sizeof(double) * bsr->nnz, hipMemcpyDeviceToHost);
  
  
    std::cout << "ilu0 bsr " << std::endl;
    //for(size_t i = 0; i < hilu0.size(); ++i)
    //{
    //    std::cout << hilu0[i] << ", ";
   // }
    std::cout << std::endl;
    //*/
    /*
    if (verbosity >= 3) {
        //HIP_CHECK(hipStreamSynchronize(stream));
        hipDeviceSynchronize();
        std::ostringstream out;
        out << "rocsparseSolver::create_preconditioner(): " << t.stop() << " s";
        //OpmLog::info(out.str());
    }
    */
}