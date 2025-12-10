//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

void vec_show(const double*, int n, const char *name);
void vec_ishow(const int*, int n, const char *name);
void vec_fill(double *v, int n, double a);
void vec_ifill(int *v, int n, int a);
void vec_copy(const double *x, double *y,int n);
double vec_inner(const double *x, const double *y, int n);
