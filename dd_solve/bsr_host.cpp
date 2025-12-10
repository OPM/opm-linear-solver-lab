//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "baseline_standalone.hpp"
#include "bsr_reader.hpp"
void init_host(solver_data *bsr, std::ostringstream &indir, int mb, int nd, int rpd, int nnzb, int block_dim) 
{
  

    int nnz;
    rocsparse_handle handle;
    
    bsr->handle =handle;
    
    //bsr->hbsr_val = {1.0, 5.0, 1.0, 4.0, 1.0, 0.0, 0.0, 0.0, 3.0, 4.0, 0.0, 1.0, 1.0, 28.0, 0.0, 1.0};
  	//bsr->hbsr_row_ptr = {0, 2, 4};
  	//bsr->hbsr_col_ind = {0, 1, 0, 1};
    //bsr->hbsr_bval = {1,1,1,1};
    
   	// This example solves
   	// | 1  5  0  0  1  2|   |y0|   |1|
   	// | 1  4  0  0  3  6|   |y1|   |1|
  	// | 0  0  1  0  6  3| * |y2| = |1|     
  	// | 0  0 28  1  8  9|   |y3|   |1|
    // | 2  3  3  5  4  5|   |y4|   |1|
    // | 1  7  6  2  9  7|   |y5|   |1| 
    
    
    read_bsr_from_txt(indir, nd, rpd, nnzb, block_dim, bsr->hbsr_val, bsr->hbsr_row_ptr, bsr->hbsr_col_ind);
    std::vector<double> hx(mb*block_dim, 1.0f);
    bsr->hbsr_bval=hx;
    
    
    //bsr->hbsr_val = {1.0, 5.0, 1.0, 4.0, 1.0, 2.0, 3.0, 6.0, 1.0, 0.0, 28.0, 1.0, 6.0, 3.0, 8.0, 9.0, 2.0, 3.0, 1.0, 7.0, 3.0, 5.0, 6.0, 2.0, 4.0, 5.0, 9.0, 7.0};
  	//bsr->hbsr_row_ptr = {0, 2, 4, 7};
  	//bsr->hbsr_col_ind = {0, 2, 1, 2, 0, 1, 2};
    //bsr->hbsr_bval = {1,1,1,1,1,1};
    
    
    
    
    
    
    //the above will be replaced by the matrices in rocky/data
    bsr->mb = bsr->hbsr_row_ptr.size()-1;
    bsr->nnzb = bsr->hbsr_col_ind.size();
    bsr->nnz =  bsr->hbsr_val.size();
    bsr->block_size = bsr->hbsr_bval.size()/bsr->mb;
    bsr->Nb = bsr->mb;
    bsr->nnzbs_prec= bsr->nnzb;
    bsr->N = bsr->hbsr_bval.size();
    
    std::cout << "created bsr with #rows "<< bsr->mb <<" #cols "<<bsr->nnzb << " block size " <<  bsr->block_size<< std::endl;
  
}





void init_host_g(solver_data2 *bsr, std::ostringstream &indir, int mb, int nd, int rpd, int nnzb, int block_dim) 
{
  

    int nnz;
    rocsparse_handle handle;
    
    bsr->handle =handle;
    
    //bsr->hbsr_val = {1.0, 5.0, 1.0, 4.0, 1.0, 0.0, 0.0, 0.0, 3.0, 4.0, 0.0, 1.0, 1.0, 28.0, 0.0, 1.0};
  	//bsr->hbsr_row_ptr = {0, 2, 4};
  	//bsr->hbsr_col_ind = {0, 1, 0, 1};
    //bsr->hbsr_bval = {1,1,1,1};
    
   	// This example solves
   	// | 1  5  0  0  1  2|   |y0|   |1|
   	// | 1  4  0  0  3  6|   |y1|   |1|
  	// | 0  0  1  0  6  3| * |y2| = |1|     
  	// | 0  0 28  1  8  9|   |y3|   |1|
    // | 2  3  3  5  4  5|   |y4|   |1|
    // | 1  7  6  2  9  7|   |y5|   |1| 
    
    
    read_bsr_from_txt(indir, nd, rpd, nnzb, block_dim, bsr->hbsr_val, bsr->hbsr_row_ptr, bsr->hbsr_col_ind);
    std::vector<double> hx(mb*block_dim, 1.0f);
    bsr->hbsr_bval=hx;
    
    
    //bsr->hbsr_val = {1.0, 5.0, 1.0, 4.0, 1.0, 2.0, 3.0, 6.0, 1.0, 0.0, 28.0, 1.0, 6.0, 3.0, 8.0, 9.0, 2.0, 3.0, 1.0, 7.0, 3.0, 5.0, 6.0, 2.0, 4.0, 5.0, 9.0, 7.0};
  	//bsr->hbsr_row_ptr = {0, 2, 4, 7};
  	//bsr->hbsr_col_ind = {0, 2, 1, 2, 0, 1, 2};
    //bsr->hbsr_bval = {1,1,1,1,1,1};
    
    
    
    
    
    
    //the above will be replaced by the matrices in rocky/data
    bsr->mb = bsr->hbsr_row_ptr.size()-1;
    bsr->nnzb = bsr->hbsr_col_ind.size();
    bsr->nnz =  bsr->hbsr_val.size();
    bsr->block_size = bsr->hbsr_bval.size()/bsr->mb;
    bsr->Nb = bsr->mb;
    bsr->nnzbs_prec= bsr->nnzb;
    bsr->N = bsr->hbsr_bval.size();
    
    std::cout << "created bsr with #rows "<< bsr->mb <<" #cols "<<bsr->nnzb << " block size " <<  bsr->block_size<< std::endl;
    
   
    
    
  
}
