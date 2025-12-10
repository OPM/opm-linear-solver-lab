//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

void mat_show(const double*, int m, int n, const char *name);

void mat_identity(double *A, int m, int n);
void mat_fma(const double *A, const double *x, double *y,int n);
void mat_dia(double *A, double c, int n);
void mat_inv(const double *A, double *B, int n);
void mat_matmul(const double *A, const double *B, double *C, int n);
void mat_matfms(const double *A, const double *B, double *C, int n);

void mat_vecmul(const double *A, double *b, int n);
void mat_vecfms(const double *A, const double *b, double *c, int n);
