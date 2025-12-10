//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

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
#include "gen_graph.hpp"
//#define MAX_THREADS_PER_BLOCK 512
#define MAX_THREADS_PER_BLOCK 1024
#define MIN_BLOCKS_PER_MULTIPROCESSOR 1
#define block_dim 3
#define num_wf 16 
//#define num_wf 8
#define WFSIZE 64
#define NINEWFSIZE 576
#define lpt1024 3 //rows_per_subdomain*block_dim/block_size
#define lpt512 12
#define lpt256 24


void append_bfs(std::vector<int> &hbsru_row_ptr,std::vector<int> &hbsru_col_ind, std::vector<int> & bfs_array, int colu)
{

 for(int i=0;i<bfs_array.size();i++)
 {
   if (bfs_array[i]==colu)
     return;
 } 
 bfs_array.push_back(colu);
 
 for (int k=hbsru_row_ptr[colu];k<hbsru_row_ptr[colu+1];k++)
 {
     int cold=hbsru_col_ind[k]; 
     append_bfs(hbsru_row_ptr, hbsru_col_ind, bfs_array, cold);
 }

}


bool is_subgraph_added(std::vector<int> &hbsru_row_ptr,std::vector<int> &hbsru_col_ind, std::vector<int> &added, int colu)
{

  //recursion until its exit point
  std::vector<int> bfs_array;
  
  append_bfs(hbsru_row_ptr, hbsru_col_ind, bfs_array, colu);
  
  if (colu==24)
  {
    for (int i=0;i<bfs_array.size();i++)
    {
      std::cout<<bfs_array[i]<< " haha ";
    }
    std::cout<<std::endl;
  }
  for (int i=0;i<bfs_array.size();i++)
  {
    if (added[i]==0)
      return false;
  }
  return true;

}



void dependency_merged_coo(int mb, int nnzb, int rows_per_subdomain,  std::vector <int> &hmap, std::vector <int> &hbsr_row_ptr, std::vector <int> &hbsr_col_ind, std::vector <int> &hdgraph_rows, std::vector <int> &hdgraph_cols, std::vector <int> &hdgraph_levptr, std::vector <int> &hdgraph_data_offsets, std::vector <int> &hbsru_row_ptr, std::vector <int> &hbsru_col_ind)

{

  std::vector<int> available_for_upper(rows_per_subdomain, 0);
  std::vector<int> added(rows_per_subdomain, 0);
  std::vector<int> last_level_nnzs;
  std::vector<int> hmapping(rows_per_subdomain, mb+1);
  hdgraph_levptr.push_back(0);
  int lzcount=0;
  int colid = 0;
  for (int i=0;i<rows_per_subdomain;i++)
  {
  //check for all level 0 rows fisrt and uppend them to hmap //bildu doesn't store diags
    if((hbsr_row_ptr[i+1]-hbsr_row_ptr[i])==1)
    { 
      //printf ("came here for %d\n",i);
      hmap.push_back(i);
      hmapping[i]=0;
      colid = hbsr_col_ind[hbsr_row_ptr[i]];
      hdgraph_data_offsets.push_back(hbsr_row_ptr[i]);
      hdgraph_cols.push_back(colid);
      hdgraph_rows.push_back(i);
      lzcount++;
      available_for_upper[i]=1;
    }
  }
  int lastel=hdgraph_levptr.back();
  hdgraph_levptr.push_back(lastel + lzcount);
  
  
  //for (int i=0;i<rows_per_subdomain;i++)
  // std::cout<<"mapping "<< i<<" " << hmapping[i]<<std::endl;
  
  //now traverse and keep explanding the tree. This is not just the BFS, we also need to make sure that rows added in a new  level only depend on the rows in previous levels, i.e, all paths from level 0 to a row in a given level do not contain levels after that row.
  int level = 0;
  int prev=level;
  while(true)
  {
    bool added_levels=false;
    lzcount =0 ;
    // at least one row must be added to a new level. If not, fix pivoting
    for (int i=0;i<rows_per_subdomain;i++)
    {
      //printf("came here 1\n");
      bool append=false;
      
      if(hmapping[i]>level)
      {
        for (int j=hbsr_row_ptr[i];j<hbsr_row_ptr[i+1]-1;j++)
        {
          //printf("came here 2\n");
          int col = hbsr_col_ind[j];
          //if (i==1)
          //  std::cout<<"checking if row "<<col<<  "with level " << hmapping[col]<< " has level less than "<< level <<std::endl;
          if (hmapping[col]>level)
            {
              append=false;
              //this row has dependencies beyond the rows in previous levels, --> break and wait for further expansion
              break;
            }
          append = true;
        }
        if (append)
        {
          added_levels=true;
          hmap.push_back(i);
          hmapping[i]=level+1; 
          last_level_nnzs.push_back(i);
          //append entire row in COO
          for (int j=hbsr_row_ptr[i];j<hbsr_row_ptr[i+1]-1;j++)
          {
            lzcount++;
            colid = hbsr_col_ind[j];
            hdgraph_cols.push_back(colid);
            hdgraph_rows.push_back(i);
            hdgraph_data_offsets.push_back(j);
          }
        }
      //update lvlptr; 
      }
    }
    
    //
    
    
    //note: the nnzs processed in current level do not become available to U until the execution of current level
    //now check with the current set of available unknowns for U, which rows of U can be added to existing levels;
   
    for (int i =0; i<rows_per_subdomain;i++)
    {
      bool available=true;
      for (int j=hbsru_row_ptr[i];j<hbsru_row_ptr[i+1];j++)
      {
        int colu=hbsru_col_ind[j];
        if(colu==i && available_for_upper[colu]!=1)
        {
            available=false;
            break;
        }
        
        //ALL nnz elements of i must be available for i to be added to the current level
        //not only the cols but all the vertices from the induced subgraph with col as the root node must be added
        if(colu!=i && is_subgraph_added(hbsru_row_ptr,hbsru_col_ind, added, colu))
        {
          available=false;
          break;
        }
      }
      if (available && added[i]==0)
      {
          added[i]=1;
          added_levels=true;
         //i can be inserted in the graph. We multiplu the COO pair and offset values to indicate upper
          for (int j=hbsru_row_ptr[i]+1;j<hbsru_row_ptr[i+1];j++) //no need to insert diag
          {
            lzcount++;
            int colud = hbsru_col_ind[j];
            hdgraph_cols.push_back(-1 * colud);
            hdgraph_rows.push_back(-1 * i);
            hdgraph_data_offsets.push_back(-1 * j);
          }
          //also append i to last_level_nnz to be marked as available for next level, possible that duplicates are added but doesn't matter.
          last_level_nnzs.push_back(i);
      }
    }
      
    //now mark last_level_nnzs of L as available for U for the next level  
    for (int i=0;i<last_level_nnzs.size();i++)
    {
      available_for_upper[last_level_nnzs[i]]=1;
    }

    //clear last level nnzs
    last_level_nnzs.clear();
    lastel=hdgraph_levptr.back();
    hdgraph_levptr.push_back(lastel + lzcount);
    
      
    level++;
    if(!added_levels) // no new levels added --> break
      break;
  }
  //
  //
  /*
  for (int i=0;i<hmap.size();i++)
  {
   std::cout<<"row "<< i << " level " << hmapping[i] << " hmap["<< i << "] " <<  hmap[i]<<std::endl;
  }
  
  for (int i=hbsr_row_ptr[176];i<hbsr_row_ptr[177];i++)
  {
   std::cout<< hbsr_col_ind[i]<<" ";
  }
  std::cout<<std::endl;
  */




}
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
std::vector <int> &level_start_subdomain)	//output  #point to the start of hdgraph_rows (and cols, as both have same index) for a given subdomain
{
  //std::vector<int> hmapping(rows_per_subdomain, mb+1);
hdgraph_levptr.push_back(0);
std::vector<int> hmapping;
int lzcount;
int colid = 0;
for (int ii=0;ii<num_subdomains;ii++)
	{
    
		lzcount=0;
		hmapping.resize(subdomain_offsets[ii+1]-subdomain_offsets[ii]);
		for (int i=0;i<hmapping.size();i++)hmapping[i]=mb;
		//std::vector<int> hmapping(subdomain_offsets[ii+1]-subdomain_offsets[ii], mb+1);
		//std::cout<<"size of hmapping "<<subdomain_offsets.size()<<" \n";
		bool level_added=false;
		for (int i=subdomain_offsets[ii];i<subdomain_offsets[ii+1];i++)
		{
		//check for all level 0 rows fisrt and uppend them to hmap //bildu doesn't store diags
			if((hbsr_row_ptr[i+1]-hbsr_row_ptr[i])==1)
			{ 
			  //printf ("came here for %d\n",i);
				hmap.push_back(i);
				hmapping[i-subdomain_offsets[ii]]=0;
				colid = hbsr_col_ind[hbsr_row_ptr[i]];
				//hdgraph_data_offsets.push_back(hbsr_row_ptr[i]);
				//hdgraph_cols.push_back(colid);
				hdgraph_rows.push_back(i);
				lzcount++;
			}
		}
		
		llevels_per_subdomain[ii]+=1;
		level_added=false;
		int lastel=hdgraph_levptr.back();
		hdgraph_levptr.push_back(lastel + lzcount);
		
		  
		  //for (int i=0;i<rows_per_subdomain;i++)
		  // std::cout<<"mapping "<< i<<" " << hmapping[i]<<std::endl;
		  
		  //now traverse and keep explanding the tree. This is not just the BFS, we also need to make sure that rows added in a new  level only depend on the rows in previous levels, i.e, all paths from level 0 to a row in a given level do not contain levels after that row.
		int level = 0;
		int prev=level;
		//printf("came here for subdomain %d\n",ii+1000);
		while(true)
		{
			bool added_levels=false;
			lzcount =0 ;
			// at least one row must be added to a new level. If not, fix pivoting
			for (int i=subdomain_offsets[ii];i<subdomain_offsets[ii+1];i++)
			//for (int i=0;i<rows_per_subdomain;i++)
			{
			    //if (ii==1) 
				//	printf("incremented lcount for subdomain %d\n",ii);
				bool append=false;
			  
				if(hmapping[i-subdomain_offsets[ii]]>level)
				{
					for (int j=hbsr_row_ptr[i];j<hbsr_row_ptr[i+1]-1;j++)
					{
					  //printf("came here 2\n");
					  int col = hbsr_col_ind[j];
					  //if (i==1)
					  //  std::cout<<"checking if row "<<col<<  "with level " << hmapping[col]<< " has level less than "<< level <<std::endl;
					  if (hmapping[col-subdomain_offsets[ii]]>level)
						{
						  append=false;
						  //this row has dependencies beyond the rows in previous levels, --> break and wait for further expansion
						  break;
						}
					  append = true;
					}
					if (append)
					{
					  added_levels=true;
					  hmap.push_back(i);
					  hmapping[i-subdomain_offsets[ii]]=level+1; 
					  hdgraph_rows.push_back(i);
      	    lzcount++;
            /*
					  //append entire row in COO
					  for (int j=hbsr_row_ptr[i];j<hbsr_row_ptr[i+1]-1;j++)
					  {
						lzcount++;
						colid = hbsr_col_ind[j];
						hdgraph_cols.push_back(colid);
						hdgraph_rows.push_back(i);
						hdgraph_data_offsets.push_back(j);
					  }
            */
					}
						//update lvlptr;
				}
			}
			if(!added_levels)
			{
			// no new levels added --> break
			  break;
			}
			llevels_per_subdomain[ii]+=1;
			
			lastel=hdgraph_levptr.back();
			//printf("writing %d to the tail at %d\n", lzcount, hdgraph_levptr.size());
			hdgraph_levptr.push_back(lastel + lzcount);
			//total_added +=lzcount;
			level++;	
		}
		if (ii==0)
		{
			level_start_subdomain[ii]=0;
		}
		level_start_subdomain[ii+1]=level_start_subdomain[ii]+llevels_per_subdomain[ii];
		  //level_start_subdomain.push_back(level_start_subdomain.back() + total_added);
		  //total_added=0;
	}


}

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
std::vector <int> &level_start_subdomain)	//output  #point to the start of hdgraph_rows (and cols, as both have same index) for a given subdomain
{
  //std::vector<int> hmapping(rows_per_subdomain, mb+1);
hdgraph_levptr.push_back(0);
std::vector<int> hmapping;
int lzcount;
int colid = 0;
for (int ii=0;ii<num_subdomains;ii++)
	{
    
		lzcount=0;
		hmapping.resize(subdomain_offsets[ii+1]-subdomain_offsets[ii]);
		for (int i=0;i<hmapping.size();i++)hmapping[i]=mb;
		//std::vector<int> hmapping(subdomain_offsets[ii+1]-subdomain_offsets[ii], mb+1);
		//std::cout<<"size of hmapping "<<subdomain_offsets.size()<<" \n";
		bool level_added=false;
		for (int i=subdomain_offsets[ii];i<subdomain_offsets[ii+1];i++)
		{
		//check for all level 0 rows fisrt and uppend them to hmap //bildu doesn't store diags
			if((hbsr_row_ptr[i+1]-hbsr_row_ptr[i])==1)
			{ 
			  //printf ("came here for %d\n",i);
				hmap.push_back(i);
				hmapping[i-subdomain_offsets[ii]]=0;
				//colid = hbsr_col_ind[hbsr_row_ptr[i]];
				//hdgraph_data_offsets.push_back(hbsr_row_ptr[i]);
				//hdgraph_cols.push_back(colid);
				hdgraph_rows.push_back(i);
				lzcount++;
			}
		}
		
		llevels_per_subdomain[ii]+=1;
		level_added=false;
		int lastel=hdgraph_levptr.back();
		hdgraph_levptr.push_back(lastel + lzcount);
		
		  
		  //for (int i=0;i<rows_per_subdomain;i++)
		  // std::cout<<"mapping "<< i<<" " << hmapping[i]<<std::endl;
		  
		  //now traverse and keep explanding the tree. This is not just the BFS, we also need to make sure that rows added in a new  level only depend on the rows in previous levels, i.e, all paths from level 0 to a row in a given level do not contain levels after that row.
		int level = 0;
		int prev=level;
		//printf("came here for subdomain %d\n",ii+1000);
		while(true)
		{
			bool added_levels=false;
			lzcount =0 ;
			// at least one row must be added to a new level. If not, fix pivoting
			for (int i=subdomain_offsets[ii];i<subdomain_offsets[ii+1];i++)
			//for (int i=0;i<rows_per_subdomain;i++)
			{
			    //if (ii==1) 
				//	printf("incremented lcount for subdomain %d\n",ii);
				bool append=false;
			  
				if(hmapping[i-subdomain_offsets[ii]]>level)
				{
					for (int j=hbsr_row_ptr[i]+1;j<hbsr_row_ptr[i+1];j++)
					{
					  //printf("came here 2\n");
					  int col = hbsr_col_ind[j];
					  //if (i==1)
					  //  std::cout<<"checking if row "<<col<<  "with level " << hmapping[col]<< " has level less than "<< level <<std::endl;
					  if (hmapping[col-subdomain_offsets[ii]]>level)
						{
						  append=false;
						  //this row has dependencies beyond the rows in previous levels, --> break and wait for further expansion
						  break;
						}
					  append = true;
					}
					if (append)
					{
					  added_levels=true;
					  hmap.push_back(i);
					  hmapping[i-subdomain_offsets[ii]]=level+1; 
					  hdgraph_rows.push_back(i);
      	    lzcount++;
					  //append entire row in COO
            /*
					  for (int j=hbsr_row_ptr[i]+1;j<hbsr_row_ptr[i+1];j++)
					  {
						lzcount++;
						colid = hbsr_col_ind[j];
						hdgraph_cols.push_back(colid);
						hdgraph_rows.push_back(i);
						hdgraph_data_offsets.push_back(j);
					  }
            */
					}
						//update lvlptr;
				}
			}
			if(!added_levels)
			{
			// no new levels added --> break
			  break;
			}
			llevels_per_subdomain[ii]+=1;
			
			lastel=hdgraph_levptr.back();
			//printf("writing %d to the tail at %d\n", lzcount, hdgraph_levptr.size());
			hdgraph_levptr.push_back(lastel + lzcount);
			//total_added +=lzcount;
			level++;	
		}
		if (ii==0)
		{
			level_start_subdomain[ii]=0;
		}
		level_start_subdomain[ii+1]=level_start_subdomain[ii]+llevels_per_subdomain[ii];
		  //level_start_subdomain.push_back(level_start_subdomain.back() + total_added);
		  //total_added=0;
	}
  //
  //
  /*
  for (int i=0;i<hmap.size();i++)
  {
   std::cout<<"row "<< i << " level " << hmapping[i] << " hmap["<< i << "] " <<  hmap[i]<<std::endl;
  }
  
  for (int i=hbsr_row_ptr[176];i<hbsr_row_ptr[177];i++)
  {
   std::cout<< hbsr_col_ind[i]<<" ";
  }
  std::cout<<std::endl;
  */
    
}
void dependency_lgraph_coo(int mb, int nnzb, int rows_per_subdomain,  std::vector <int> &hmap, std::vector <int> &hbsr_row_ptr, std::vector <int> &hbsr_col_ind, std::vector <int> &hdgraph_rows, std::vector <int> &hdgraph_cols, std::vector <int> &hdgraph_levptr, std::vector <int> &hdgraph_data_offsets)
{
 //dependency tracking only on the fisrt subdomain. //The remaining subdomains are identical so no need to track dependency for all subs. Update this appropriately for the generic case later on
  std::vector<int> hmapping(rows_per_subdomain, mb+1);
  hdgraph_levptr.push_back(0);
  int lzcount=0;
  int colid = 0;
  for (int i=0;i<rows_per_subdomain;i++)
  {
  //check for all level 0 rows fisrt and uppend them to hmap //bildu doesn't store diags
    if((hbsr_row_ptr[i+1]-hbsr_row_ptr[i])==1)
    { 
      //printf ("came here for %d\n",i);
      hmap.push_back(i);
      hmapping[i]=0;
      colid = hbsr_col_ind[hbsr_row_ptr[i]];
      hdgraph_data_offsets.push_back(hbsr_row_ptr[i]);
      hdgraph_cols.push_back(colid);
      hdgraph_rows.push_back(i);
      lzcount++;
    }
  }
  int lastel=hdgraph_levptr.back();
  hdgraph_levptr.push_back(lastel + lzcount);
  
  
  //for (int i=0;i<rows_per_subdomain;i++)
  // std::cout<<"mapping "<< i<<" " << hmapping[i]<<std::endl;
  
  //now traverse and keep explanding the tree. This is not just the BFS, we also need to make sure that rows added in a new  level only depend on the rows in previous levels, i.e, all paths from level 0 to a row in a given level do not contain levels after that row.
  int level = 0;
  int prev=level;
  while(true)
  {
    bool added_levels=false;
    lzcount =0 ;
    // at least one row must be added to a new level. If not, fix pivoting
    for (int i=0;i<rows_per_subdomain;i++)
    {
      //printf("came here 1\n");
      bool append=false;
      
      if(hmapping[i]>level)
      {
        for (int j=hbsr_row_ptr[i];j<hbsr_row_ptr[i+1]-1;j++)
        {
          //printf("came here 2\n");
          int col = hbsr_col_ind[j];
          //if (i==1)
          //  std::cout<<"checking if row "<<col<<  "with level " << hmapping[col]<< " has level less than "<< level <<std::endl;
          if (hmapping[col]>level)
            {
              append=false;
              //this row has dependencies beyond the rows in previous levels, --> break and wait for further expansion
              break;
            }
          append = true;
        }
        if (append)
        {
          added_levels=true;
          hmap.push_back(i);
          hmapping[i]=level+1; 
          
          //append entire row in COO
          for (int j=hbsr_row_ptr[i];j<hbsr_row_ptr[i+1]-1;j++)
          {
            lzcount++;
            colid = hbsr_col_ind[j];
            hdgraph_cols.push_back(colid);
            hdgraph_rows.push_back(i);
            hdgraph_data_offsets.push_back(j);
          }
        }
      //update lvlptr;
       
      }
    }
    
    lastel=hdgraph_levptr.back();
    hdgraph_levptr.push_back(lastel + lzcount);
    
    level++;
    if(!added_levels) // no new levels added --> break
      break;
  }
  //
  //
  /*
  for (int i=0;i<hmap.size();i++)
  {
   std::cout<<"row "<< i << " level " << hmapping[i] << " hmap["<< i << "] " <<  hmap[i]<<std::endl;
  }
  
  for (int i=hbsr_row_ptr[176];i<hbsr_row_ptr[177];i++)
  {
   std::cout<< hbsr_col_ind[i]<<" ";
  }
  std::cout<<std::endl;
  */
    
}

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
std::vector <int> &level_start_subdomain)	//output  #point to the start of hdgraph_rows (and cols, as both have same index) for a given subdomain
{
  //std::vector<int> hmapping(rows_per_subdomain, mb+1);
hdgraph_levptr.push_back(0);
std::vector<int> hmapping;
int lzcount;
int colid = 0;
for (int ii=0;ii<num_subdomains;ii++)
	{
    
		lzcount=0;
		hmapping.resize(subdomain_offsets[ii+1]-subdomain_offsets[ii]);
		for (int i=0;i<hmapping.size();i++)hmapping[i]=mb;
		//std::vector<int> hmapping(subdomain_offsets[ii+1]-subdomain_offsets[ii], mb+1);
		//std::cout<<"size of hmapping "<<subdomain_offsets.size()<<" \n";
		bool level_added=false;
		for (int i=subdomain_offsets[ii];i<subdomain_offsets[ii+1];i++)
		{
		//check for all level 0 rows fisrt and uppend them to hmap //bildu doesn't store diags
			if((hbsr_row_ptr[i+1]-hbsr_row_ptr[i])==1)
			{ 
			  //printf ("came here for %d\n",i);
				hmap.push_back(i);
				hmapping[i-subdomain_offsets[ii]]=0;
				colid = hbsr_col_ind[hbsr_row_ptr[i]];
				hdgraph_data_offsets.push_back(hbsr_row_ptr[i]);
				hdgraph_cols.push_back(colid);
				hdgraph_rows.push_back(i);
				lzcount++;
			}
		}
		
		llevels_per_subdomain[ii]+=1;
		level_added=false;
		int lastel=hdgraph_levptr.back();
		hdgraph_levptr.push_back(lastel + lzcount);
		
		  
		  //for (int i=0;i<rows_per_subdomain;i++)
		  // std::cout<<"mapping "<< i<<" " << hmapping[i]<<std::endl;
		  
		  //now traverse and keep explanding the tree. This is not just the BFS, we also need to make sure that rows added in a new  level only depend on the rows in previous levels, i.e, all paths from level 0 to a row in a given level do not contain levels after that row.
		int level = 0;
		int prev=level;
		//printf("came here for subdomain %d\n",ii+1000);
		while(true)
		{
			bool added_levels=false;
			lzcount =0 ;
			// at least one row must be added to a new level. If not, fix pivoting
			for (int i=subdomain_offsets[ii];i<subdomain_offsets[ii+1];i++)
			//for (int i=0;i<rows_per_subdomain;i++)
			{
			    //if (ii==1) 
				//	printf("incremented lcount for subdomain %d\n",ii);
				bool append=false;
			  
				if(hmapping[i-subdomain_offsets[ii]]>level)
				{
					for (int j=hbsr_row_ptr[i];j<hbsr_row_ptr[i+1]-1;j++)
					{
					  //printf("came here 2\n");
					  int col = hbsr_col_ind[j];
					  //if (i==1)
					  //  std::cout<<"checking if row "<<col<<  "with level " << hmapping[col]<< " has level less than "<< level <<std::endl;
					  if (hmapping[col-subdomain_offsets[ii]]>level)
						{
						  append=false;
						  //this row has dependencies beyond the rows in previous levels, --> break and wait for further expansion
						  break;
						}
					  append = true;
					}
					if (append)
					{
					  added_levels=true;
					  hmap.push_back(i);
					  hmapping[i-subdomain_offsets[ii]]=level+1; 
					  
					  //append entire row in COO
					  for (int j=hbsr_row_ptr[i];j<hbsr_row_ptr[i+1]-1;j++)
					  {
						lzcount++;
						colid = hbsr_col_ind[j];
						hdgraph_cols.push_back(colid);
						hdgraph_rows.push_back(i);
						hdgraph_data_offsets.push_back(j);
					  }
					}
						//update lvlptr;
				}
			}
			if(!added_levels)
			{
			// no new levels added --> break
			  break;
			}
			llevels_per_subdomain[ii]+=1;
			
			lastel=hdgraph_levptr.back();
			//printf("writing %d to the tail at %d\n", lzcount, hdgraph_levptr.size());
			hdgraph_levptr.push_back(lastel + lzcount);
			//total_added +=lzcount;
			level++;	
		}
		if (ii==0)
		{
			level_start_subdomain[ii]=0;
		}
		level_start_subdomain[ii+1]=level_start_subdomain[ii]+llevels_per_subdomain[ii];
		  //level_start_subdomain.push_back(level_start_subdomain.back() + total_added);
		  //total_added=0;
	}
  //
  //
  /*
  for (int i=0;i<hmap.size();i++)
  {
   std::cout<<"row "<< i << " level " << hmapping[i] << " hmap["<< i << "] " <<  hmap[i]<<std::endl;
  }
  
  for (int i=hbsr_row_ptr[176];i<hbsr_row_ptr[177];i++)
  {
   std::cout<< hbsr_col_ind[i]<<" ";
  }
  std::cout<<std::endl;
  */
    
}


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
std::vector <int> &level_start_subdomain)	//output  #point to the start of hdgraph_rows (and cols, as both have same index) for a given subdomain
{
  //std::vector<int> hmapping(rows_per_subdomain, mb+1);
hdgraph_levptr.push_back(0);
std::vector<int> hmapping;
int lzcount;
int colid = 0;
for (int ii=0;ii<num_subdomains;ii++)
	{
    
		lzcount=0;
		hmapping.resize(subdomain_offsets[ii+1]-subdomain_offsets[ii]);
		for (int i=0;i<hmapping.size();i++)hmapping[i]=mb;
		//std::vector<int> hmapping(subdomain_offsets[ii+1]-subdomain_offsets[ii], mb+1);
		//std::cout<<"size of hmapping "<<subdomain_offsets.size()<<" \n";
		bool level_added=false;
		for (int i=subdomain_offsets[ii];i<subdomain_offsets[ii+1];i++)
		{
		//check for all level 0 rows fisrt and uppend them to hmap //bildu doesn't store diags
			if((hbsr_row_ptr[i+1]-hbsr_row_ptr[i])==1)
			{ 
			  //printf ("came here for %d\n",i);
				hmap.push_back(i);
				hmapping[i-subdomain_offsets[ii]]=0;
				colid = hbsr_col_ind[hbsr_row_ptr[i]];
				hdgraph_data_offsets.push_back(hbsr_row_ptr[i]);
				hdgraph_cols.push_back(colid);
				hdgraph_rows.push_back(i);
				lzcount++;
			}
		}
		
		llevels_per_subdomain[ii]+=1;
		level_added=false;
		int lastel=hdgraph_levptr.back();
		hdgraph_levptr.push_back(lastel + lzcount);
		
		  
		  //for (int i=0;i<rows_per_subdomain;i++)
		  // std::cout<<"mapping "<< i<<" " << hmapping[i]<<std::endl;
		  
		  //now traverse and keep explanding the tree. This is not just the BFS, we also need to make sure that rows added in a new  level only depend on the rows in previous levels, i.e, all paths from level 0 to a row in a given level do not contain levels after that row.
		int level = 0;
		int prev=level;
		//printf("came here for subdomain %d\n",ii+1000);
		while(true)
		{
			bool added_levels=false;
			lzcount =0 ;
			// at least one row must be added to a new level. If not, fix pivoting
			for (int i=subdomain_offsets[ii];i<subdomain_offsets[ii+1];i++)
			//for (int i=0;i<rows_per_subdomain;i++)
			{
			    //if (ii==1) 
				//	printf("incremented lcount for subdomain %d\n",ii);
				bool append=false;
			  
				if(hmapping[i-subdomain_offsets[ii]]>level)
				{
					for (int j=hbsr_row_ptr[i]+1;j<hbsr_row_ptr[i+1];j++)
					{
					  //printf("came here 2\n");
					  int col = hbsr_col_ind[j];
					  //if (i==1)
					  //  std::cout<<"checking if row "<<col<<  "with level " << hmapping[col]<< " has level less than "<< level <<std::endl;
					  if (hmapping[col-subdomain_offsets[ii]]>level)
						{
						  append=false;
						  //this row has dependencies beyond the rows in previous levels, --> break and wait for further expansion
						  break;
						}
					  append = true;
					}
					if (append)
					{
					  added_levels=true;
					  hmap.push_back(i);
					  hmapping[i-subdomain_offsets[ii]]=level+1; 
					  
					  //append entire row in COO
					  for (int j=hbsr_row_ptr[i]+1;j<hbsr_row_ptr[i+1];j++)
					  {
						lzcount++;
						colid = hbsr_col_ind[j];
						hdgraph_cols.push_back(colid);
						hdgraph_rows.push_back(i);
						hdgraph_data_offsets.push_back(j);
					  }
					}
						//update lvlptr;
				}
			}
			if(!added_levels)
			{
			// no new levels added --> break
			  break;
			}
			llevels_per_subdomain[ii]+=1;
			
			lastel=hdgraph_levptr.back();
			//printf("writing %d to the tail at %d\n", lzcount, hdgraph_levptr.size());
			hdgraph_levptr.push_back(lastel + lzcount);
			//total_added +=lzcount;
			level++;	
		}
		if (ii==0)
		{
			level_start_subdomain[ii]=0;
		}
		level_start_subdomain[ii+1]=level_start_subdomain[ii]+llevels_per_subdomain[ii];
		  //level_start_subdomain.push_back(level_start_subdomain.back() + total_added);
		  //total_added=0;
	}
  //
  //
  /*
  for (int i=0;i<hmap.size();i++)
  {
   std::cout<<"row "<< i << " level " << hmapping[i] << " hmap["<< i << "] " <<  hmap[i]<<std::endl;
  }
  
  for (int i=hbsr_row_ptr[176];i<hbsr_row_ptr[177];i++)
  {
   std::cout<< hbsr_col_ind[i]<<" ";
  }
  std::cout<<std::endl;
  */
    
}


void dependency_ugraph_coo(int mb, int nnzb, int rows_per_subdomain,  std::vector <int> &hmap, std::vector <int> &hbsr_row_ptr, std::vector <int> &hbsr_col_ind, std::vector <int> &hdgraph_rows, std::vector <int> &hdgraph_cols, std::vector <int> &hdgraph_levptr, std::vector <int> &hdgraph_data_offsets)
{
 //dependency tracking only on the fisrt subdomain. //The remaining subdomains are identical so no need to track dependency for all subs. Update this appropriately for the generic case later on
  std::vector<int> hmapping(rows_per_subdomain, mb+1);
  hdgraph_levptr.push_back(0);
  int lzcount=0;
  int colid = 0;
  for (int i=0;i<rows_per_subdomain;i++)
  {
  //check for level 1 rows fisrt and uppend them to hmap
    if((hbsr_row_ptr[i+1]-hbsr_row_ptr[i])==1)
    { 
      //printf ("came here for %d\n",i);
      hmap.push_back(i);
      hmapping[i]=0;
      colid = hbsr_col_ind[hbsr_row_ptr[i]];
      hdgraph_data_offsets.push_back(hbsr_row_ptr[i]);
      hdgraph_cols.push_back(colid);
      hdgraph_rows.push_back(i);
      lzcount++;
    }
  }
  int lastel=hdgraph_levptr.back();
  hdgraph_levptr.push_back(lastel + lzcount);
  
  
  //for (int i=0;i<rows_per_subdomain;i++)
  // std::cout<<"mapping "<< i<<" " << hmapping[i]<<std::endl;
  
  //now traverse and keep explanding the tree. This is not just the BFS, we also need to make sure that rows added in a new  level only depend on the rows in previous levels, i.e, all paths from level 0 to a row in a given level do not contain levels after that row.
  int level = 0;
  int prev=level;
  while(true)
  {
    bool added_levels=false;
    lzcount =0 ;
    // at least one row must be added to a new level. If not, fix pivoting
    for (int i=0;i<rows_per_subdomain;i++)
    {
      //printf("came here 1\n");
      bool append=false;
      
      if(hmapping[i]>level)
      {
        for (int j=hbsr_row_ptr[i]+1;j<hbsr_row_ptr[i+1];j++)
        {
          //printf("came here 2\n");
          int col = hbsr_col_ind[j];
          //if (i==1)
          //  std::cout<<"checking if row "<<col<<  "with level " << hmapping[col]<< " has level less than "<< level <<std::endl;
          if (hmapping[col]>level)
            {
              append=false;
              //this row has dependencies beyond the rows in previous levels, --> break and wait for further expansion
              break;
            }
          append = true;
        }
        if (append)
        {
          added_levels=true;
          hmap.push_back(i);
          hmapping[i]=level+1; 
          
          //append entire row in COO
          for (int j=hbsr_row_ptr[i]+1;j<hbsr_row_ptr[i+1];j++) //start from rowptr[i]+1 as rowptr[i] is the diagonal
          {
            lzcount++;
            colid = hbsr_col_ind[j];
            hdgraph_cols.push_back(colid);
            hdgraph_rows.push_back(i);
            hdgraph_data_offsets.push_back(j);
          }
        }
      //update lvlptr;
       
      }
    }
    
    lastel=hdgraph_levptr.back();
    hdgraph_levptr.push_back(lastel + lzcount);
    
    level++;
    if(!added_levels) // no new levels added --> break
      break;
  }
  //
  //
  /*
  for (int i=0;i<hmap.size();i++)
  {
   std::cout<<"row "<< i << " level " << hmapping[i] << " hmap["<< i << "] " <<  hmap[i]<<std::endl;
  }
  
  for (int i=hbsr_row_ptr[176];i<hbsr_row_ptr[177];i++)
  {
   std::cout<< hbsr_col_ind[i]<<" ";
  }
  std::cout<<std::endl;
  */
    
}

