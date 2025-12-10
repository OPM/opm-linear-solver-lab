//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

//#ifndef gen_graph_h
//#define gen_graph_h
#pragma once

#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include "hip/hip_runtime.h"
#include <rocsparse/rocsparse.h>
#include <vector>
#include <bits/stdc++.h> 
#include <fstream>
#include <iterator>
#include <string>
#include <bits/stdc++.h> 
#include <sstream>
#include <limits>
#include <chrono>

void dependency_lgraph_coo(int mb, int nnzb, int rows_per_subdomain,  std::vector <int> &hmap, std::vector <int> &hbsr_row_ptr, std::vector <int> &hbsr_col_ind, std::vector <int> &hdgraph_rows, std::vector <int> &hdgraph_cols, std::vector <int> &hdgraph_levptr, std::vector <int> &hdgraph_data_offset);;

void dependency_ugraph_coo(int mb, int nnzb, int rows_per_subdomain,  std::vector <int> &hmap, std::vector <int> &hbsr_row_ptr, std::vector <int> &hbsr_col_ind, std::vector <int> &hdgraph_rows, std::vector <int> &hdgraph_cols, std::vector <int> &hdgraph_levptr, std::vector <int> &hdgraph_data_offsets);

void dependency_merged_coo(int mb, int nnzb, int rows_per_subdomain,  std::vector <int> &hmap, std::vector <int> &hbsr_row_ptr, std::vector <int> &hbsr_col_ind, std::vector <int> &hdgraph_rows, std::vector <int> &hdgraph_cols, std::vector <int> &hdgraph_levptr, std::vector <int> &hdgraph_data_offsets, std::vector <int> &hbsru_row_ptr, std::vector <int> &hbsru_col_ind);

void dependency_lgraph_coo_offsets(
int mb, 
int nnzb, 
int num_subdomains, 
std::vector <int> &subdomain_offsets,   //input   //used to get start row and end of of subdomain. Size = #subdomains + 1 
std::vector <int> &hmap, 				//input     
std::vector <int> &hbsr_row_ptr, 		//input     
std::vector <int> &hbsr_col_ind, 		//input     
std::vector <int> &llevels_per_subdomain,	//output #levels per each subdomain
std::vector <int> &hdgraph_rows,            //output 
std::vector <int> &hdgraph_cols, 			//output
std::vector <int> &hdgraph_levptr, 			//output #point to rows and cols for each level
std::vector <int> &hdgraph_data_offsets,    //output #point to the vector values
std::vector <int> &level_start_subdomain); //output point to the start id of levptr for each subdomain

void dependency_ugraph_coo_offsets(
int mb, 
int nnzb, 
int num_subdomains, 
std::vector <int> &subdomain_offsets,   //input   //used to get start row and end of of subdomain. Size = #subdomains + 1 
std::vector <int> &hmap, 				//input     
std::vector <int> &hbsr_row_ptr, 		//input     
std::vector <int> &hbsr_col_ind, 		//input     
std::vector <int> &llevels_per_subdomain,	//output #levels per each subdomain
std::vector <int> &hdgraph_rows,            //output 
std::vector <int> &hdgraph_cols, 			//output
std::vector <int> &hdgraph_levptr, 			//output #point to rows and cols for each level
std::vector <int> &hdgraph_data_offsets,    //output #point to the vector values
std::vector <int> &level_start_subdomain);  //output point to the start id of levptr for each subdomain


void dependency_lgraph_vc_offsets(
int mb, 
int nnzb, 
int num_subdomains, 
std::vector <int> &subdomain_offsets,   //input   //used to get start row and end of of subdomain. Size = #subdomains + 1 
std::vector <int> &hmap, 				//input     
std::vector <int> &hbsr_row_ptr, 		//input     
std::vector <int> &hbsr_col_ind, 		//input     
std::vector <int> &llevels_per_subdomain,	//output #levels per each subdomain
std::vector <int> &hdgraph_rows,            //output 
std::vector <int> &hdgraph_levptr, 			//output #point to rows and cols for each level
std::vector <int> &level_start_subdomain);	//output  #point to the start of hdgraph_rows (and cols, as both have same index) for a given subdomain

void dependency_ugraph_vc_offsets(
int mb, 
int nnzb, 
int num_subdomains, 
std::vector <int> &subdomain_offsets,   //input   //used to get start row and end of of subdomain. Size = #subdomains + 1 
std::vector <int> &hmap, 				//input     
std::vector <int> &hbsr_row_ptr, 		//input     
std::vector <int> &hbsr_col_ind, 		//input     
std::vector <int> &llevels_per_subdomain,	//output #levels per each subdomain
std::vector <int> &hdgraph_rows,            //output 
std::vector <int> &hdgraph_levptr, 			//output #point to rows and cols for each level
std::vector <int> &level_start_subdomain);	//output  #point to the start of hdgraph_rows (and cols, as both have same index) for a given subdomain


