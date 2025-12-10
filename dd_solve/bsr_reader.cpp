//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "baseline_standalone.hpp"


//Compile using hipcc -o gpu_pbicgstab baseline_standalone.cpp -I/opt/rocm-5.7.1/include -L/opt/rocm-5.7.1/lib -lrocsparse -lrocblas
//The above assumes rocpsparse is installed in /opt/rocm-5.7.1
void read_bsr_from_txt(std::ostringstream& indir, int numblocks, int rows_per_block, int nnzb, int block_dim, std::vector <double> &hbsr_val, std::vector <int> &hbsr_row_ptr, std::vector <int> &hbsr_col_ind)

{

  int mb = numblocks*rows_per_block;
  std::ostringstream oss_val;
  oss_val << indir.str() <<"Avals_" << mb << "_" << nnzb<<"_"<<nnzb*block_dim*block_dim<<".txt";
  
  std::string fname_val = oss_val.str();    
  std::ostringstream oss_row;
  oss_row << indir.str() << "Arows_" << numblocks<<"_"<<rows_per_block<<".txt";
  std::string fname_row = oss_row.str();    
  
  std::ostringstream oss_col;
  oss_col << indir.str() << "Acols_" << nnzb<<".txt";
  std::string fname_col = oss_col.str();    
  
  std::cout << "Vals read from "<<fname_val << "\n"; 
  std::cout << "Rows read from "<<fname_row << "\n"; 
  std::cout << "Cols read from "<<fname_col << "\n"; 

  std::ifstream inv(fname_val.c_str());

  // Check if object is valid
  if (!inv)
  {
      std::cout << "Cannot open the File : " << fname_val << std::endl;
      exit(0);
  }
  std::string strv;
  // Read the next line from File untill it reaches the end.
  while (std::getline(inv, strv))
  {
      // Line contains string of length > 0 then save it in vector
      if (strv.size() > 0)
          hbsr_val.push_back(std::stod(strv));
  }
  inv.close();
  
  
  std::ifstream inc(fname_col.c_str());

  // Check if object is valid
  if (!inc)
  {
      std::cout << "Cannot open the File : " << fname_col << std::endl;
      exit(0);
  }
  std::string strc;
  // Read the next line from File untill it reaches the end.
  while (std::getline(inc, strc))
  {
      // Line contains string of length > 0 then save it in vector
      if (strc.size() > 0)
          hbsr_col_ind.push_back(std::stod(strc));
  }
  inc.close();
  
  
  std::ifstream inr(fname_row.c_str());

  // Check if object is valid
  if (!inr)
  {
      std::cout << "Cannot open the File : " << fname_row << std::endl;
      exit(0);
  }
  std::string strr;
  // Read the next line from File untill it reaches the end.
  while (std::getline(inr, strr))
  {
      // Line contains string of length > 0 then save it in vector
      if (strr.size() > 0)
          hbsr_row_ptr.push_back(std::stod(strr));
  }
  inr.close();
  
  
  
  
  //std::ifstream is("numbers.txt");
  //std::istream_iterator<double> start(is), end;
  //std::vector<double> numbers(start, end);
  //std::cout << "Read " << numbers.size() << " numbers" << std::endl;

  // print the numbers to stdout
  //std::cout << "numbers read in:\n";
  //std::copy(numbers.begin(), numbers.end(), 
  //          std::ostream_iterator<double>(std::cout, " "));
  //std::cout << std::endl;

  
  /*
  //printing rowptrs 
  std::cout<<"printing rowptrs\n";
  for (int i =0;i<hbsr_row_ptr.size();i++)
    std::cout<<hbsr_row_ptr[i]<<" ";
  std::cout<<"\n";
  
  
  std::cout<<"printing cols\n";
  for (int i =0;i<hbsr_col_ind.size();i++)
    std::cout<<hbsr_col_ind[i]<<" ";
  std::cout<<"\n";
  
  
  std::cout<<"printing vals\n";
  for (int i =0;i<hbsr_val.size();i++)
    std::cout<<hbsr_val[i]<<" ";
  std::cout<<"\n";
  */
}


void read_bsr(std::ostringstream& indir, int mb, int nnzb, int block_dim, std::vector <double> &hbsr_val, std::vector <int> &hbsr_row_ptr, std::vector <int> &hbsr_col_ind, std::vector <double> &hbsr_bval)

{

  std::ostringstream oss_val;
  oss_val << indir.str() <<"Avals.fp64_"<<nnzb*block_dim*block_dim;
  std::string fname_val = oss_val.str();    
  
  std::ostringstream oss_row;
  oss_row << indir.str() << "Arows.int32_" << mb+1;
  std::string fname_row = oss_row.str();    
  
  std::ostringstream oss_col;
  oss_col << indir.str() << "Acols.int32_" << nnzb;
  std::string fname_col = oss_col.str();    
  
  std::ostringstream oss_bval;
  oss_bval << indir.str() << "bvals.fp64_" <<mb*block_dim;
  std::string fname_bval = oss_bval.str();  
  
  std::cout << "BSR Vals read from "<<fname_val << "\n"; 
  std::cout << "Rows read from "<<fname_row << "\n"; 
  std::cout << "Cols read from "<<fname_col << "\n"; 
  std::cout << "B Vals read from "<<fname_bval << "\n";
 
  
  //Avals
  std::ifstream ifs(fname_val, std::ios::in | std::ios::binary);
  if (ifs) {
    std::streampos fileSize;

    size_t sizeOfBuffer;
  
    ifs.seekg(0, std::ios::end);
    fileSize = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
  
    sizeOfBuffer = fileSize / sizeof(double);
    
    hbsr_val.resize(sizeOfBuffer);
    ifs.read(reinterpret_cast<char*>(hbsr_val.data()), fileSize);
  }
  else
  {
    std::cout << fname_val << " not found, exiting "<<std::endl;
    exit(0);
  }
  
  //Arows
  std::ifstream ifs_row(fname_row, std::ios::in | std::ios::binary);
  if (ifs_row) {
    // Temporary Variables
    std::streampos fileSize;

    size_t sizeOfBuffer;
  
    ifs_row.seekg(0, std::ios::end);
    fileSize = ifs_row.tellg();
    ifs_row.seekg(0, std::ios::beg);
  
    sizeOfBuffer = fileSize / sizeof(int);
    
    hbsr_row_ptr.resize(sizeOfBuffer);
    ifs_row.read(reinterpret_cast<char*>(hbsr_row_ptr.data()), fileSize);
  }
  else
  {
    std::cout << fname_row << " not found, exiting "<<std::endl;
    exit(0);
  }
  
  //Acols
  std::ifstream ifs_cols(fname_col, std::ios::in | std::ios::binary);
  if (ifs_cols) {
    std::streampos fileSize;
 
    size_t sizeOfBuffer;
  
    ifs_cols.seekg(0, std::ios::end);
    fileSize = ifs_cols.tellg();
    ifs_cols.seekg(0, std::ios::beg);
  
    sizeOfBuffer = fileSize / sizeof(int);
    
    hbsr_col_ind.resize(sizeOfBuffer);
    ifs_cols.read(reinterpret_cast<char*>(hbsr_col_ind.data()), fileSize);
  }
  else
  {
    std::cout << fname_col << " not found, exiting "<<std::endl;
    exit(0);
  }
  
  
  //Bvals
  std::ifstream ifs_bval(fname_bval, std::ios::in | std::ios::binary);
  if (ifs_bval) {
    // Temporary Variables
    std::streampos fileSize;

    size_t sizeOfBuffer;
  
    ifs_bval.seekg(0, std::ios::end);
    fileSize = ifs_bval.tellg();
    ifs_bval.seekg(0, std::ios::beg);
  
    sizeOfBuffer = fileSize / sizeof(double);
    
    hbsr_bval.resize(sizeOfBuffer);
    ifs_bval.read(reinterpret_cast<char*>(hbsr_bval.data()), fileSize);
  }
  else
  {
    std::cout << fname_bval << " not found, exiting "<<std::endl;
    exit(0);
  }
  
  

  //printing rowptrs 
  std::cout<<"first 40 rowptrs\n";
  //for (int i =0;i<hbsr_row_ptr.size();i++)
  for (int i =0;i<40;i++)
    std::cout<<hbsr_row_ptr[i]<<" ";
  std::cout<<"\n\n";
  
  
  std::cout<<"first 40 cols\n";
  //for (int i =0;i<hbsr_col_ind.size();i++)
  for (int i =0;i<40;i++)
    std::cout<<hbsr_col_ind[i]<<" ";
  std::cout<<"\n\n";
  
  std::cout<<"first 40 vals\n";
  //for (int i =0;i<hbsr_val.size();i++)
  for (int i =0;i<40;i++)
    std::cout<<hbsr_val[i]<<" ";
  std::cout<<"\n\n";
  
  
  std::cout<<"first 40 bvals\n";
  //for (int i =0;i<hbsr_bval.size();i++)
  for (int i =0;i<40;i++)
    std::cout<<hbsr_bval[i]<<" ";
  std::cout<<"\n\n";
  
  
}

