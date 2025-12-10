//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "baseline_standalone.hpp"

void gpu_pbicgstab(solver_data *bsr, double tolerance,float maxit)
{

    float it = 0.5;
    double rho, rhop, beta, alpha, nalpha, omega, nomega, tmp1, tmp2;
    double norm, norm_0;
    double zero = 0.0;
    double one  = 1.0;
    double mone = -1.0;
    std::vector<double> norms;
   
    // HIP_VERSION is defined as (HIP_VERSION_MAJOR * 10000000 + HIP_VERSION_MINOR * 100000 + HIP_VERSION_PATCH)
#if HIP_VERSION >= 50400000
    ROCSPARSE_CHECK(rocsparse_dbsrmv_ex(bsr->handle, bsr->dir, bsr->operation,
                                        bsr->Nb, bsr->Nb, bsr->nnzb, &one, bsr->descr_M,
                                        bsr->d_Avals, bsr->d_Arows, bsr->d_Acols, bsr->block_size,
                                        bsr->spmv_info, bsr->d_x, &zero, bsr->d_r));
                                        
#else
    ROCSPARSE_CHECK(rocsparse_dbsrmv(bsr->handle, bsr->dir, bsr->operation,
                                        bsr->Nb, bsr->Nb, bsr->nnzb, &one, bsr->descr_M,
                                        bsr->d_Avals, bsr->d_Arows, bsr->d_Acols, bsr->block_size,
                                        bsr->d_x, &zero, bsr->d_r));
#endif
    
    ROCBLAS_CHECK(rocblas_dscal(bsr->blas_handle, bsr->N, &mone, bsr->d_r, 1));
    ROCBLAS_CHECK(rocblas_daxpy(bsr->blas_handle, bsr->N, &one, bsr->d_b, 1, bsr->d_r, 1));
    ROCBLAS_CHECK(rocblas_dcopy(bsr->blas_handle, bsr->N, bsr->d_r, 1, bsr->d_rw, 1));
    ROCBLAS_CHECK(rocblas_dcopy(bsr->blas_handle, bsr->N, bsr->d_r, 1, bsr->d_p, 1));
    ROCBLAS_CHECK(rocblas_dnrm2(bsr->blas_handle, bsr->N, bsr->d_r, 1, &norm_0));

	  /*
    if (verbosity >= 2) {
        std::ostringstream out;
        out << std::scientific << "rocsparseSolver initial norm: " << norm_0;
        OpmLog::info(out.str());
    }

    if (verbosity >= 3) {
        t_rest.start();
    }
	
	  */
    for (it = 0.5; it < maxit; it += 0.5) {
        rhop = rho;
        ROCBLAS_CHECK(rocblas_ddot(bsr->blas_handle, bsr->N, bsr->d_rw, 1, bsr->d_r, 1, &rho));

        if (it > 1) {
            beta = (rho / rhop) * (alpha / omega);
            nomega = -omega;
            ROCBLAS_CHECK(rocblas_daxpy(bsr->blas_handle, bsr->N, &nomega,  bsr->d_v, 1, bsr->d_p, 1));
            ROCBLAS_CHECK(rocblas_dscal(bsr->blas_handle, bsr->N, &beta, bsr->d_p, 1));
            ROCBLAS_CHECK(rocblas_daxpy(bsr->blas_handle, bsr->N, &one, bsr->d_r, 1, bsr->d_p, 1));
        }
		/*
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_rest.stop();
            t_prec.start();
        }
		*/
        // apply ilu0
        ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(bsr->handle, bsr->dir, \
                              bsr->operation, bsr->Nb, bsr->nnzbs_prec, &one, \
                              bsr->descr_L, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, bsr->d_p, bsr->d_t, rocsparse_solve_policy_auto, bsr->d_buffer));
        ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(bsr->handle, bsr->dir, \
                              bsr->operation, bsr->Nb, bsr->nnzbs_prec, &one, \
                              bsr->descr_U, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, bsr->d_t, bsr->d_pw, rocsparse_solve_policy_auto, bsr->d_buffer));
		/*
		if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_prec.stop();
            t_spmv.start();
        }
		*/
        // spmv
#if HIP_VERSION >= 50400000
        ROCSPARSE_CHECK(rocsparse_dbsrmv_ex(bsr->handle, bsr->dir, bsr->operation,
                                            bsr->Nb, bsr->Nb, bsr->nnzb, &one, bsr->descr_M,
                                            bsr->d_Avals, bsr->d_Arows, bsr->d_Acols, bsr->block_size,
                                            bsr->spmv_info, bsr->d_pw, &zero,  bsr->d_v));
#else
        ROCSPARSE_CHECK(rocsparse_dbsrmv(bsr->handle, bsr->dir, bsr->operation,
                                            bsr->Nb, bsr->Nb, bsr->nnzb, &one, bsr->descr_M,
                                            bsr->d_Avals, bsr->d_Arows, bsr->d_Acols, bsr->block_size,
                                            bsr->d_pw, &zero,  bsr->d_v));
#endif
		/*
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_spmv.stop();
            t_well.start();
        }

        // apply wellContributions
        if(wellContribs.getNumWells() > 0){
            static_cast<WellContributionsRocsparse&>(wellContribs).apply(bsr->bsr->d_pw,  bsr->d_v);
        }
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_well.stop();
            t_rest.start();
        }
		*/
        ROCBLAS_CHECK(rocblas_ddot(bsr->blas_handle, bsr->N, bsr->d_rw, 1,  bsr->d_v, 1, &tmp1));
        alpha = rho / tmp1;
        nalpha = -alpha;
        ROCBLAS_CHECK(rocblas_daxpy(bsr->blas_handle, bsr->N, &nalpha,  bsr->d_v, 1, bsr->d_r, 1));
        ROCBLAS_CHECK(rocblas_daxpy(bsr->blas_handle, bsr->N, &alpha, bsr->d_pw, 1, bsr->d_x, 1));
        ROCBLAS_CHECK(rocblas_dnrm2(bsr->blas_handle, bsr->N, bsr->d_r, 1, &norm));
        /*
		if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_rest.stop();
        }
		*/
        if (norm < tolerance * norm_0) {
            break;
        }

        it += 0.5;

        // apply ilu0
		/*
        if (verbosity >= 3) {
            t_prec.start();
        }
		*/
        ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(bsr->handle, bsr->dir, \
                              bsr->operation, bsr->Nb, bsr->nnzbs_prec, &one, \
                              bsr->descr_L, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, bsr->d_r, bsr->d_t, rocsparse_solve_policy_auto, bsr->d_buffer));
        ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(bsr->handle, bsr->dir, \
                              bsr->operation, bsr->Nb, bsr->nnzbs_prec, &one, \
                              bsr->descr_U, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, bsr->d_t, bsr->d_s, rocsparse_solve_policy_auto, bsr->d_buffer));
        /*
		if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_prec.stop();
            t_spmv.start();
        }
		*/

        // spmv
#if HIP_VERSION >= 50400000
        ROCSPARSE_CHECK(rocsparse_dbsrmv_ex(bsr->handle, bsr->dir, bsr->operation,
                                            bsr->Nb, bsr->Nb, bsr->nnzb, &one, bsr->descr_M,
                                            bsr->d_Avals, bsr->d_Arows, bsr->d_Acols, bsr->block_size,
                                            bsr->spmv_info, bsr->d_s, &zero, bsr->d_t));
#else
        ROCSPARSE_CHECK(rocsparse_dbsrmv(bsr->handle, bsr->dir, bsr->operation,
                                            bsr->Nb, bsr->Nb, bsr->nnzb, &one, bsr->descr_M,
                                            bsr->d_Avals, bsr->d_Arows, bsr->d_Acols, bsr->block_size,
                                            bsr->d_s, &zero, bsr->d_t));
#endif
		/*
        if(verbosity >= 3){
            HIP_CHECK(hipStreamSynchronize(stream));
            t_spmv.stop();
            t_well.start();
        }
		
        // apply wellContributions
        if(wellContribs.getNumWells() > 0){
            static_cast<WellContributionsRocsparse&>(wellContribs).apply(bsr->d_s, bsr->d_t);
        }
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_well.stop();
            t_rest.start();
        }
		*/
        ROCBLAS_CHECK(rocblas_ddot(bsr->blas_handle, bsr->N, bsr->d_t, 1, bsr->d_r, 1, &tmp1));
        ROCBLAS_CHECK(rocblas_ddot(bsr->blas_handle, bsr->N, bsr->d_t, 1, bsr->d_t, 1, &tmp2));
        omega = tmp1 / tmp2;
        nomega = -omega;
        ROCBLAS_CHECK(rocblas_daxpy(bsr->blas_handle, bsr->N, &omega, bsr->d_s, 1, bsr->d_x, 1));
        ROCBLAS_CHECK(rocblas_daxpy(bsr->blas_handle, bsr->N, &nomega, bsr->d_t, 1, bsr->d_r, 1));

        ROCBLAS_CHECK(rocblas_dnrm2(bsr->blas_handle, bsr->N, bsr->d_r, 1, &norm));
        /*
		if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_rest.stop();
        }
		*/
        if (norm < tolerance * norm_0) {
            break;
        }
		/*
        if (verbosity > 1) {
            std::ostringstream out;
            out << "it: " << it << std::scientific << ", norm: " << norm;
            OpmLog::info(out.str());
        }
		*/
    }

    //res.iterations = std::min(it, (float)maxit);
    //res.reduction = norm / norm_0;
    //res.conv_rate  = static_cast<double>(pow(res.reduction, 1.0 / it));
    //res.elapsed = t_total.stop();
    //res.converged = (it != (maxit + 0.5));

	  float iterations = std::min(it, (float)maxit);
    double reduction = norm / norm_0;
    double conv_rate  = static_cast<double>(pow(reduction, 1.0 / it));
    //res.elapsed = t_total.stop();
    int converged = (it != (maxit + 0.5));

    //if (verbosity >= 1) {
    //    std::ostringstream out;
    //    out << "=== converged: " << res.converged << ", conv_rate: " << res.conv_rate << ", time: " << res.elapsed << \
    //        ", time per iteration: " << res.elapsed / it << ", iterations: " << it;
    //    OpmLog::info(out.str());
    //}
	
	  //if (verbosity >= 1) {
    //    std::ostringstream out;
    //    out << "=== converged: " << res.converged << ", conv_rate: " << res.conv_rate << ", time: " << res.elapsed << \
    //        ", time per iteration: " << res.elapsed / it << ", iterations: " << it;
    //    std::cout<<out.str();
    //}
	
	  if (bsr->verbosity >= 1) {
        std::ostringstream out;
        out << "=== converged: " << converged << ", conv_rate: " << conv_rate << ", iterations: " << it << std::endl;
        std::cout<<out.str();
    }
	
	/*
    if (verbosity >= 3) {
        std::ostringstream out;
        out << "rocsparseSolver::prec_apply:  " << t_prec.elapsed() << " s\n";
        out << "rocsparseSolver::spmv:        " << t_spmv.elapsed() << " s\n";
        out << "rocsparseSolver::well:        " << t_well.elapsed() << " s\n";
        out << "rocsparseSolver::rest:        " << t_rest.elapsed() << " s\n";
        out << "rocsparseSolver::total_solve: " << res.elapsed << " s\n";
        OpmLog::info(out.str());
    }	
	*/
}
