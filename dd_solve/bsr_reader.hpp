//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#pragma once

void read_bsr(std::ostringstream& indir, int mb, int nnzb, int block_dim, std::vector <double> &hbsr_val, std::vector <int> &hbsr_row_ptr, std::vector <int> &hbsr_col_ind, std::vector <double> &hbsr_bval);

void read_bsr_from_txt(std::ostringstream& indir, int numblocks, int rows_per_block, int nnzb, int block_dim, std::vector <double> &hbsr_val, std::vector <int> &hbsr_row_ptr, std::vector <int> &hbsr_col_ind);