//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "bsr.hpp"
typedef
struct bilu_prec
{
    bsr_matrix L;
    bsr_matrix U;
}
bilu_prec;

void bilu_create(bsr_matrix *M, bilu_prec *P);
void bilu_apply(bilu_prec *P, double *x);


typedef
struct bildu_prec
{
    bsr_matrix L;
    bsr_matrix D;
    bsr_matrix U;
}
bildu_prec;


void bildu_create(bsr_matrix *M, bildu_prec *P);
void bildu_apply(bildu_prec *P, double *x);
