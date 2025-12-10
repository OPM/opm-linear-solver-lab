//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "baseline_standalone.hpp"


#define MAX_THREADS_PER_BLOCK 1024
#define MIN_BLOCKS_PER_MULTIPROCESSOR 1
#define block_dim 3
#define num_wf 16 
//#define num_wf 8
#define WFSIZE 64
#define NINEWFSIZE 576

//#define lpt1024 3 //for subdomain size 1024/#define lpt1024 6 // for 2048 //rows_per_subdomain*block_dim/block_size
//#define lpt1024 6
#define lpt1024 12 //for fp32 subdomains of size 4096
#define lpt512 12
#define lpt256 24


__global__ void bsr_lts_nf_subdomain_per_block_block_levels_lds_xy_2k_fp32(
double* __restrict__ bsrl_val, 
double* __restrict__ bsru_val, 
double* __restrict__ dval,
double* __restrict__ x, 
double* __restrict__ y,
int num_subdomains,
int* __restrict__ dlevel_start_subdomain, 
int* __restrict__ dhdl2graph_rows,
int* __restrict__ dhdl2graph_cols,
int* __restrict__ dhdl2graph_levptr,
int* __restrict__ dhdl2graph_data_offsets,
int* __restrict__ dllevels_per_subdomain,
int* __restrict__ dlevelu_start_subdomain, 
int* __restrict__ dhdu2graph_rows,
int* __restrict__ dhdu2graph_cols,
int* __restrict__ dhdu2graph_levptr,
int* __restrict__ dhdu2graph_data_offsets,
int* __restrict__ dllevelsu_per_subdomain,
int* __restrict__ dsubdomain_offsets
)
{
 const int lid = hipThreadIdx_x & (WFSIZE - 1);
 const int wfid = hipThreadIdx_x / WFSIZE;
 const int bid = hipBlockIdx_x;
 const int tid = hipThreadIdx_x;
 //const int WFSIZE=64;
 //const int num_wf = rows_per_subdomain/WFSIZE; //logical number of warps per block. Doesn't matter if the actual number of active warps per block is lesser than this number. The kernel is immune to differences in physical #warps and logical#warps
 bool SLEEP =0;
 if (bid>= num_subdomains)
   return;

 //const int num_strides = rows_per_subdomain/hipBlockDim_x;
 /*
 if(tid==0)
 {
  printf("start: %d end: %d bid: %d\n",start,end,bid);
 }
 */
 //__shared__ double local_x[3072];
 //__shared__ double sum[6144];
 __shared__ float sum[12288];
 //__shared__ double sum[3072];
 //__shared__ int done[512];
 //__shared__ int done[1024];
 //__shared__ double sum[192];
 //__syncthreads();
 int rows_per_subdomain = dsubdomain_offsets[bid+1]-dsubdomain_offsets[bid];
 //if (tid ==0)
 //printf("\n#rows in subdomain %d is %d\n", bid, rows_per_subdomain);
 const int num_strides =rows_per_subdomain/hipBlockDim_x;
 int loc = tid*block_dim;
 int xloc = bid*rows_per_subdomain; //******TODO FOR NON UNIFORM SUBDOMAINS: MAKE THIS XLOC START AND XLOC END
 //done[tid]=0;
 /*
 for (int i=tid;i<rows_per_subdomain;i+=hipBlockDim_x)
 {
   //done[i]=0;
   int fetch_start = i*block_dim;
   #pragma unroll 3
   for (int j = 0; j<block_dim;j++)
   {
     sum[fetch_start + j]=x[xloc*block_dim + fetch_start + j];
     //local_x[fetch_start + j]=x[xloc*block_dim + fetch_start + j];
   }
 }
 */
 

 int fetch_start = tid*lpt1024;
 int fetch_end = fetch_start+ lpt1024;
 
 //#pragma unroll lpt1024
 
 for (int i=fetch_start;i<fetch_end;i+=1)
 {
 
   //done[i]=0;
  
   sum[i]=(float)x[xloc*block_dim + i];
   //#pragma unroll 3
   //for (int j = 0; j<block_dim;j++)
   //{
   //  sum[i + j]=x[xloc*block_dim + i + j];;
   //local_x[fetch_start + j]=x[xloc*block_dim + fetch_start + j];
   //}
 }
 
 
 /*
  //#pragma unroll lpt1024
 for (int i=bid*rows_per_subdomain*3;i<(bid+1)*rows_per_subdomain*3;i+=1)
 {
 
   //done[i]=0;
  
   sum[i-xloc*block_dim]=x[i];
   //#pragma unroll 3
   //for (int j = 0; j<block_dim;j++)
   //{
   //  sum[i + j]=x[xloc*block_dim + i + j];;
   //local_x[fetch_start + j]=x[xloc*block_dim + fetch_start + j];
   //}
 }
 */
 
 __syncthreads();
 int lev_offset= dlevel_start_subdomain[bid];
 //if(tid==0)
 //printf("lev offset %d\n ",lev_offset); 
 //if (tid ==0)
 //  printf("came here from wid %d, levels per subdomain: %d \n",bid, dllevels_per_subdomain[bid]);
 for (int i = 1;i<dllevels_per_subdomain[bid];i+=1) // i is the level
 {
   //if (tid ==0 && i==1)
   //printf("came here from wid %d\n",bid);
   int start = dhdl2graph_levptr[lev_offset + i];
   int end = dhdl2graph_levptr[lev_offset + i+1];
   
   //if (tid==0 && i==18)
   //printf ("start %d end %d diff %d\n",start,end,end-start);
   
 
   for (int j = start+tid; j<end;j+=hipBlockDim_x) //each WF processes a unique level
   {
     //printf("wavefront %d is processing row %d start %d end %d\n", wfid, i,start,end);
     //printf("wavefront %d is processing row %d\n", wfid, i);
	
     int elrow = dhdl2graph_rows[j] ;
     int elcol = dhdl2graph_cols[j] ;
	
     //int eldat = dgraph_data_offsets[j]; //Todo make this generic update this as well
     //int eldat = dgraph_data_offsets[j] + bid*nnzs_per_subdomain;  
     int eldat = dhdl2graph_data_offsets[j] ;
	 //if ((tid==0 && j==start+tid) && i==18)
	 //   printf ("elrow %d elcol %d eldat %d bid %d\n",elrow, elcol, eldat, bid);
     int start_val= eldat*block_dim*block_dim;
     double block[9] = {0,0,0,0,0,0,0,0,0};
     for (int ii=0;ii<9;ii++)
     {
      block[ii]=bsrl_val[start_val+ii];
     } 
     
     
    
     double col_els[3]={0,0,0};
     for (int ii=0;ii<3;ii++)
       col_els[ii]=sum[(block_dim*elcol- xloc*block_dim)+ii];
      //col_els[ii]=y[(block_dim*col_id)+ii];
    
     for (int k =0; k<block_dim;k++)
     {
      int offset = elrow*block_dim + k -xloc*block_dim;
      double local_update=0;
      for (int inner =0; inner<block_dim;inner++)
      {
          local_update+=  -1*block[(k*block_dim)+inner]*col_els[inner];
      }
      atomicAdd(&sum[offset],local_update); //warpsum later.
     }
   }
   __syncthreads();
   /*
   if(lid==0)
   {
     int done_row=i;
     __hip_atomic_store(&done[done_row], 1, __ATOMIC_RELEASE, __HIP_MEMORY_SCOPE_AGENT);
     //__hip_atomic_store(&done[done_row], 1, __ATOMIC_RELEASE, __HIP_MEMORY_SCOPE_WORKGROUP);
   }
   */
 }  
 __syncthreads();
 
 //int fetch_start = tid*lpt1024;
 //int fetch_end = fetch_start+ lpt1024;
 //#pragma unroll lpt1024
 //tid=0 --> sum[0,1,2,4,5,6]
 //tid=1 --> sum[7,8,9,10,11,12]
 
 
 for (int i=tid;i<rows_per_subdomain;i+=hipBlockDim_x)
 //for (int i=fetch_start;i<fetch_end;i+=1)
 {
 
   double blockd[9] = {0,0,0,0,0,0,0,0,0};
   for (int ii=0;ii<block_dim*block_dim;ii++)
   {
    blockd[ii]=dval[(xloc+i)*block_dim*block_dim+ii];
   } 
   
   
   double sumv[3] ={0,0,0};
   double ans[3] ={0,0,0};
   for (int ii=0;ii<block_dim;ii++)
   {
    sumv[ii]=sum[i*block_dim+ii];
   }
   //dense matvec
   
   
   for (int ii=0;ii<block_dim;ii++)
   {
     double tempans=0;
     for (int jj=0;jj<block_dim;jj++)
     {
         tempans+=blockd[ii*block_dim+jj]*sumv[jj];
     }
     ans[ii]=tempans;
   }
   
   
   //write back to sum
   
   for (int ii=0;ii<block_dim;ii++)
   {
      sum[i*block_dim+ii]=ans[ii];
   }

 }
 
 
 __syncthreads();
 int lev_offsetu= dlevelu_start_subdomain[bid];
 for (int i = 1;i<dllevelsu_per_subdomain[bid];i+=1) // i is the level
 {
   int start = dhdu2graph_levptr[lev_offsetu+i];
   int end = dhdu2graph_levptr[lev_offsetu+i+1];
 
   for (int j = start+tid; j<end;j+=hipBlockDim_x) //each WF processes a unique level
   {
     
     //printf("wavefront %d is processing row %d start %d end %d\n", wfid, i,start,end);
     int elrow = dhdu2graph_rows[j] ;
     int elcol = dhdu2graph_cols[j] ;
     //int eldat = dgraph_data_offsets[j]; //Todo make this generic update this as well
     //int eldat = dgraph_data_offsets[j] + bid*nnzs_per_subdomain;  
     int eldat = dhdu2graph_data_offsets[j] ;  
     int start_val= eldat*block_dim*block_dim;
     double block[9] = {0,0,0,0,0,0,0,0,0};
     for (int ii=0;ii<9;ii++)
     {
      block[ii]=bsru_val[start_val+ii];
     } 
     
    
     double col_els[3]={0,0,0};
     for (int ii=0;ii<3;ii++)
       col_els[ii]=sum[(block_dim*elcol)+ii-xloc*block_dim];
      //col_els[ii]=y[(block_dim*col_id)+ii];
    
     for (int k =0; k<block_dim;k++)
     {
      int offset = elrow*block_dim + k - (xloc*block_dim);
      double local_update=0;
      for (int inner =0; inner<block_dim;inner++)
      {
          local_update+=  -1*block[(k*block_dim)+inner]*col_els[inner];
      }
      atomicAdd(&sum[offset],local_update); //warpsum later.
     }
   }
   __syncthreads();
   /*
   if(lid==0)
   {
     int done_row=i;
     __hip_atomic_store(&done[done_row], 1, __ATOMIC_RELEASE, __HIP_MEMORY_SCOPE_AGENT);
     //__hip_atomic_store(&done[done_row], 1, __ATOMIC_RELEASE, __HIP_MEMORY_SCOPE_WORKGROUP);
   }
   */
 }  
 __syncthreads();
 
 //int fetch_start = tid*lpt1024;
 //int fetch_end = fetch_start+ lpt1024;
 //#pragma unroll lpt1024
 //tid=0 --> sum[0,1,2,4,5,6]
 //tid=1 --> sum[7,8,9,10,11,12]
 
 
 
  
 
 for (int i=fetch_start;i<fetch_end;i+=1)
 {
   //done[i]=0;
  
   y[xloc*block_dim + i] = (double)sum[i];
   //#pragma unroll 3
   //for (int j = 0; j<block_dim;j++)
   //{
   //  sum[i + j]=x[xloc*block_dim + i + j];;
   //local_x[fetch_start + j]=x[xloc*block_dim + fetch_start + j];
   //}
 }
 
 
 
 /*
 for (int i=bid*rows_per_subdomain*3;i<(bid+1)*rows_per_subdomain*3;i+=1)
 {
 
   //done[i]=0;
   y[i] = sum[i-xloc*block_dim];
  // sum[i]=x[i];
   //#pragma unroll 3
   //for (int j = 0; j<block_dim;j++)
   //{
   //  sum[i + j]=x[xloc*block_dim + i + j];;
   //local_x[fetch_start + j]=x[xloc*block_dim + fetch_start + j];
   //}
 }
 */
	
	
}

void gpu_pbicgstab_g_ec_fp32(solver_data2 *bsr, double tolerance,float maxit)
{

    float it = 0.5;
    double rho, rhop, beta, alpha, nalpha, omega, nomega, tmp1, tmp2;
    double norm, norm_0;
    double zero = 0.0;
    double one  = 1.0;
    double mone = -1.0;
    std::vector<double> norms;
    int grid_size_b_sub=bsr->num_subdomains;
    int block_size = 1024;
    //int block_size = 1;
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
			//if (rho <1e-6)
			//{
			//	ROCBLAS_CHECK(rocblas_dcopy(bsr->blas_handle, bsr->N, bsr->d_r, 1, bsr->d_rw, 1));
			//	ROCBLAS_CHECK(rocblas_dcopy(bsr->blas_handle, bsr->N, bsr->d_r, 1, bsr->d_p, 1));
			//	ROCBLAS_CHECK(rocblas_ddot(bsr->blas_handle, bsr->N, bsr->d_rw, 1,  bsr->d_r, 1, &rho));
			//	//rhop=rho;
			//}
			//else
			{
				beta = (rho / rhop) * (alpha / omega);

				nomega = -omega;
				ROCBLAS_CHECK(rocblas_daxpy(bsr->blas_handle, bsr->N, &nomega,  bsr->d_v, 1, bsr->d_p, 1));
				ROCBLAS_CHECK(rocblas_dscal(bsr->blas_handle, bsr->N, &beta, bsr->d_p, 1));
				ROCBLAS_CHECK(rocblas_daxpy(bsr->blas_handle, bsr->N, &one, bsr->d_r, 1, bsr->d_p, 1));
			}
			//std::cout<<" rho "<<rho<<std::endl;
		}
		/*
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_rest.stop();
            t_prec.start();
        }
		*/
        // apply ilu0
        /*
        ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(bsr->handle, bsr->dir, \
                              bsr->operation, bsr->Nb, bsr->nnzbs_prec, &one, \
                              bsr->descr_L, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, bsr->d_p, bsr->d_t, rocsparse_solve_policy_auto, bsr->d_buffer));
        ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(bsr->handle, bsr->dir, \
                              bsr->operation, bsr->Nb, bsr->nnzbs_prec, &one, \
                              bsr->descr_U, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, bsr->d_t, bsr->d_pw, rocsparse_solve_policy_auto, bsr->d_buffer));
                              */
                              
        bsr_lts_nf_subdomain_per_block_block_levels_lds_xy_2k_fp32<<<grid_size_b_sub,block_size>>>(
        bsr->d_Mlvals, 
        bsr->d_Muvals,
        bsr->dval,
        bsr->d_p, 
        bsr->d_pw,
        bsr->num_subdomains,
        bsr->dlevel_start_subdomain, 
        bsr->dhdl2graph_rows,
        bsr->dhdl2graph_cols,
        bsr->dhdl2graph_levptr,
        bsr->dhdl2graph_data_offsets,
        bsr->dllevels_per_subdomain,
        bsr->dlevelu_start_subdomain, 
        bsr->dhdu2graph_rows,
        bsr->dhdu2graph_cols,
        bsr->dhdu2graph_levptr,
        bsr->dhdu2graph_data_offsets,
        bsr->dllevelsu_per_subdomain,
        bsr->dsubdomain_offsets
        );
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
        //std::cout<<"norm "<<norm <<" upper bound " <<tolerance * norm_0<<std::endl;
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
        /*
        ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(bsr->handle, bsr->dir, \
                              bsr->operation, bsr->Nb, bsr->nnzbs_prec, &one, \
                              bsr->descr_L, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, bsr->d_r, bsr->d_t, rocsparse_solve_policy_auto, bsr->d_buffer));
        ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(bsr->handle, bsr->dir, \
                              bsr->operation, bsr->Nb, bsr->nnzbs_prec, &one, \
                              bsr->descr_U, bsr->d_Mvals, bsr->d_Mrows, bsr->d_Mcols, bsr->block_size, bsr->ilu_info, bsr->d_t, bsr->d_s, rocsparse_solve_policy_auto, bsr->d_buffer));
        */
        
        bsr_lts_nf_subdomain_per_block_block_levels_lds_xy_2k_fp32<<<grid_size_b_sub,block_size>>>(
        bsr->d_Mlvals, 
        bsr->d_Muvals,
        bsr->dval,
        bsr->d_r, 
        bsr->d_s,
        bsr->num_subdomains,
        bsr->dlevel_start_subdomain, 
        bsr->dhdl2graph_rows,
        bsr->dhdl2graph_cols,
        bsr->dhdl2graph_levptr,
        bsr->dhdl2graph_data_offsets,
        bsr->dllevels_per_subdomain,
        bsr->dlevelu_start_subdomain, 
        bsr->dhdu2graph_rows,
        bsr->dhdu2graph_cols,
        bsr->dhdu2graph_levptr,
        bsr->dhdu2graph_data_offsets,
        bsr->dllevelsu_per_subdomain,
        bsr->dsubdomain_offsets
        );
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
