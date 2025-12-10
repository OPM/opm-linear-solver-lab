//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#ifndef bsr_hpp
#define bsr_hpp
#include <stdbool.h>
#include "graph_bsr.hpp"
typedef
struct bsr_matrix
{
    int bsize;
    int nrows;
    int ncols;
    int nnz;
    int* colidx;
    int* rowptr;
    double *data;    
}
bsr_matrix;

void bsr_dirichlet(bsr_matrix *M, int i);
void bsr_to_graph(bsr_matrix *M, bsr_graph *G);
void bsr_reorder(bsr_matrix *M, bsr_graph *G, bsr_matrix *T);
void bsr_show(bsr_matrix *M, const char *name);
void bsr_create(bsr_matrix *M, int nrows, int ncols, int nnz, int bsize);
void bsr_identity(bsr_matrix *M, int n, int bsize);
void bsr_dot(bsr_matrix *M, const double *x, double *y);
void bsr_laplacian(bsr_matrix *M, int nx, int ny, int nz, int bsize);
void bsr_upper(bsr_matrix *M, bsr_matrix *U, bool strict);
void bsr_lower(bsr_matrix *M, bsr_matrix *L, bool strict);
void bsr_diagonal(bsr_matrix *M, bsr_matrix *D);
void bsr_transpose(bsr_matrix *M, bsr_matrix *T);
void bsr_cartesian(bsr_matrix *M, int nx, int ny, int nz, int bsize);
void bsr_flux_balance(bsr_matrix *M);
/*
void csr_copy(csr_matrix *M, csr_matrix *C);
void csr_reverse(csr_matrix *M, csr_matrix *R);
void csr_diagonal(csr_matrix *M, csr_matrix *D);
*/

#endif
