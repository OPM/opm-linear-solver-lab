//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#pragma once


void rocsparse_boilerplate(solver_data *bsr) ;
void data_xfers(solver_data *bsr);
void init_gpu(solver_data *bsr);

void rocsparse_boilerplate_g(solver_data2 *bsr) ;
void data_xfers_g(solver_data2 *bsr);
void init_gpu_g(solver_data2 *bsr);

void rocsparse_boilerplate_vc(solver_data_vc *bsr) ;
void data_xfers_vc(solver_data_vc *bsr);
void init_gpu_vc(solver_data_vc *bsr);