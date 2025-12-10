//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#ifndef bsr_graph_hpp
#define bsr_graph_hpp
#include <stdbool.h>
typedef
struct bsr_graph
{
    int nrows;
    int ncols;
    int nnz;
    int* colidx;
    int* rowptr;
    int *weight;    
    int post_dd_edges;
}
bsr_graph;

void bsr_graph_create(bsr_graph *M, int nrows, int ncols, int nnz);
void bsr_graph_show(bsr_graph *M, const char *name);
void bsr_graph_transpose(bsr_graph *M, bsr_graph *T);
void bsr_graph_rowcount(bsr_graph *M, int **count);
void bsr_graph_tags(bsr_graph *M, int *tag, int n);
void bsr_graph_reorder(bsr_graph *M, int **pmap,bsr_graph *P);
void bsr_graph_decompose(bsr_graph *M, bsr_graph *P, bsr_graph *D);
void bsr_graph_part(int nx, int ny, int nz,int nblkx, int nblky, int nblkz,int* part_idxs);

/*
void bsr_dirichlet(bsr_matrix *M, int i);

void bsr_copy(bsr_matrix *M, bsr_matrix *C);
void bsr_dot(bsr_matrix *M, const double *x, double *b);
void bsr_transpose(bsr_matrix *M, bsr_matrix *T);
void bsr_reverse(bsr_matrix *M, bsr_matrix *R);
void bsr_upper(bsr_matrix *M, bsr_matrix *U, bool strict);
void bsr_lower(bsr_matrix *M, bsr_matrix *L, bool strict);
void bsr_diagonal(bsr_matrix *M, bsr_matrix *D);
void bsr_identity(bsr_matrix *M, int n);
void bsr_laplacian(bsr_matrix *M, int nx, int ny, int nz);
void bsr_cartesian(bsr_matrix *M, int nx, int ny, int nz);
void bsr_flux_balance(bsr_matrix *M);
void bsr_level(bsr_matrix *M, int *level);
*/
#endif
