//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#pragma once
#include "baseline_standalone.hpp"
void init_host(solver_data *bsr, std::ostringstream &indir, int mb, int nd, int rpd, int nnzb, int block_dim);
void init_host_g(solver_data2 *bsr, std::ostringstream &indir, int mb, int nd, int rpd, int nnzb, int block_dim);
