//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "baseline_standalone.hpp"
#include "bsr_reader.hpp"
#include "bsr_ilu0_generate.hpp"
#include "bsr_host.hpp"
#include "bsr_gpu_boilerplate.hpp"
#include "solver.hpp"
#include "solver_g.hpp"
#include "solver_vc.hpp"
#include "bsr.hpp"
#include "vec.hpp"
#include "graph_bsr.hpp"
#include "chrono"
#include "bilu.hpp"
#include "solver_nolds_vc.hpp"
#include "solver_vc_fp32.hpp"
#include "solver_g_fp32.hpp"
void decompress_u(std::vector <double> &ud,std::vector <int> &ur,std::vector <int> &uc, std::vector <double> &uud,std::vector <int> &uur,std::vector <int> &uuc, int rows_per_subdomain)
{

std::vector<double> identity = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};

  uur.push_back(0);
  //uur.push_back(1);
  //uuc.push_back(0);
  //uud.insert(std::end(lld), std::begin(identity), std::end(identity));
  //std::cout<<" rows_per_subdomain "<<rows_per_subdomain << " ur size "<<ur.size()<<" urc size "<<uc.size()<<" urd size "<<ud.size()<<std::endl;
  //exit(0);
  for (int i=0;i<ur.size()-1;i++)
  {
    //if(i%rows_per_subdomain==(rows_per_subdomain-1))//TODO: for non uniform subdomain sizes, the start row of the subdomain 
    if(ur[i]==ur[i+1])
    {
      
      uud.insert(std::end(uud), std::begin(identity), std::end(identity));  
      //lct++;
      uuc.push_back(i);
      int back2=uur.back();
      uur.push_back(back2+1); 
    }
    else
    {
      int start = ur[i];
      int end = ur[i+1];
     
      bool diag_added=false;
      int uct=0;
      for (int j= start;j<end;j++)
      {
        int colid = uc[j];
        if (colid==i)
        {
         diag_added=true;
        }
      }
       if(!diag_added)
      {
        uud.insert(std::end(uud), std::begin(identity), std::end(identity));  
        uct++;
        uuc.push_back(i);
        diag_added=true;
      }
      
      for (int j= start;j<end;j++)
      {
      
        int colid2 = uc[j];
        if (colid2!=i)
        {
         //diag_added=true;
        
          uct++;
          for (int k=0;k<9;k++){
            uud.push_back(ud[j*9+k]);
          }
          uuc.push_back(colid2);  
        }   
      }
      int back=uur.back();
      uur.push_back(back+uct);
    }  
  }
  //std::cout<<" rows_per_subdomain "<<rows_per_subdomain << " uur size "<<uur.size()<<" uurc size "<<uuc.size()<<" uurd size "<<uud.size()<<std::endl;
  //exit(0);
  //int backl=uur.back();
  //uur.push_back(backl+1);
  //uuc.push_back(ur.size()-1);
  //uud.insert(std::end(uud), std::begin(identity), std::end(identity));

}


void decompress_l(std::vector <double> &ld,std::vector <int> &lr,std::vector <int> &lc, std::vector <double> &lld,std::vector <int> &llr,std::vector <int> &llc, int rows_per_subdomain)
{

  std::vector<double> identity = {1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
  llr.push_back(0);
  llr.push_back(1);
  llc.push_back(0);
  lld.insert(std::end(lld), std::begin(identity), std::end(identity));
  //std::cout<<"lr size"<< lr.size()<<std::endl;
  //exit(0);
  for (int i=1;i<lr.size()-1;i++)
  {
    int lct=0;
    //if(i%rows_per_subdomain==0)//TODO: for non uniform subdomain sizes, the start row of the subdomain 
    if(lr[i]==lr[i+1])
    {
      
      lld.insert(std::end(lld), std::begin(identity), std::end(identity));  
      //lct++;
      llc.push_back(i);
      int back2=llr.back();
      llr.push_back(back2+1); 
    }
    else
    {
      int start = lr[i];
      int end = lr[i+1];
     
      bool diag_added=false;
      
      for (int j= start;j<end;j++)
      {
        int colid = lc[j];
        //if (colid==(i+1))
        if (colid==(i))
        {
         diag_added=true;
        }
        lct++;
        for (int k=0;k<9;k++){
          lld.push_back(ld[j*9+k]);
        }
        llc.push_back(colid);     
      }
      if(!diag_added)
      {
        lld.insert(std::end(lld), std::begin(identity), std::end(identity));  
        lct++;
        llc.push_back(i);
        diag_added=true;
      }
      int back=llr.back();
      llr.push_back(back+lct); 
    } 
  }
   //std::cout<<" rows_per_subdomain "<<rows_per_subdomain << " llr size "<<llr.size()<<" llrc size "<<llc.size()<<" llrd size "<<lld.size()<<std::endl;
  // exit(0);

}

void run_solver_graph_ildu_vc_nolds(bsr_matrix *C, solver_data_vc *bsr, int nx, int ny, int nz, int bsize, int nbx, int nby, int nbz, float maxit, double tolerance, std::vector<double> &hy_og, std::vector<double> &ans, int num_subdomains)
{
   
    std::vector<int> hbsr_row_ptr;
    std::vector<int> hbsr_col_ind;
    std::vector<double> hbsr_val;
    
    std::vector<int> hbsr_plrow_ptr;
    std::vector<int> hbsr_plcol_ind;
    std::vector<double> hbsr_plval;
    
    std::vector<int> hbsr_purow_ptr;
    std::vector<int> hbsr_pucol_ind;
    std::vector<double> hbsr_puval;
    
    

    int n=nx*ny*nz;
    int nnz= nx*ny*nz +2*(nx-1)*ny*nz + 2*nx*(ny-1)*nz + 2*nx*ny*(nz-1);
    int b=bsize;
    int bb=b*b;
    int total_nnz = bb*nnz;  
    int* colidx;
    int* rowptr;
    double* data;
    int mb =n;
    //int rows_per_subdomain = nbx*nby*nbz;
    //colidx=(int*)malloc(nnz*sizeof(int));
    //rowptr=(int*)malloc((n+1)*sizeof(int));
    //data=(double*)malloc(bb*nnz*sizeof(double));
    
    
    //bsr_matrix C, D;
    bsr_cartesian(C, nx, ny, nz, bsize);

    //flux balance
    bsr_flux_balance(C);

    //boundary condititions
    bsr_dirichlet(C,0);
    bsr_dirichlet(C,n-1);
    //bsr_show(&C,"C");
    bsr_matrix  D1;
    //int nbx=16;
    //int nby=16;
    //int nbz=4;
    
   
    
    
    int *pid = (int*)malloc(n*sizeof(int));
    
    vec_ifill(pid,n,-1);
    auto start1 = std::chrono::system_clock::now();
    bsr_graph_part(nx, ny, nz, nbx, nby, nbz, pid);
    //vec_ishow(pid,n,"pid");
    
	  //reinterpret partition ids as a graph
    bsr_graph PT, PN, PC, A;
    bsr_graph_tags(&PT, pid, n);
    //bsr_graph_show(&PT,"PT");
	
	  //apply transpose to generate node permutation graph
    bsr_graph_transpose(&PT,&PN);
    //bsr_graph_show(&PN,"PN");
	
	  //extract graph of sparse matrix C and apply node permutation
    // --> the end results is a connection permutation graph
	  bsr_to_graph(C,&A);
    bsr_graph_reorder(&A, &(PN.weight), &PC);
    
    

    
	  //bsr_graph_show(&PC,"PC");
	  //    vec_ishow(PC.weight,A.nnz,"PC");

    // apply connection permutation to sparse matrix C
    // --> to get the reordered sparse matrix D
    bsr_reorder(C, &PC, &D1);
    bsr_graph B, DD;
    bsr_to_graph(&D1,&B);
    bsr_graph_decompose(&B, &PN, &DD);
    //csr_graph_show(&B,"B");
    //bsr_graph_show(&DD,"DD");

    // apply permutation
    // label conn to drop
    // apply projection
    
    DD.nnz = DD.post_dd_edges;
    bsr_matrix F;
    bsr_reorder(&D1, &DD, &F);
    //bsr_show(&F,"F");
    
    
    F.nnz = DD.post_dd_edges;
    auto end1 = std::chrono::system_clock::now();
    double elapsed_ms1= std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1).count();
    std::cout << "partition + decomposition time  "<<elapsed_ms1 << '\n';
    
    
    for(int i=0;i<(D1.nrows+1);i++)
      hbsr_row_ptr.push_back(D1.rowptr[i]);
    for(int i=0;i<D1.nnz;i++)
      hbsr_col_ind.push_back(D1.colidx[i]);
    for(int i=0;i<D1.nnz*bsize*bsize;i++)
      hbsr_val.push_back(D1.data[i]);
     
    //std::cout<<"Size "<<F.bsize<<" "<<F.nrows<<" " <<F.nrows<<" "<<F.ncols<<" "<<F.nnz<<"\n";
    bildu_prec P;
    bildu_create(&F, &P);
     
    //std::cout<<"Size D"<<D1.bsize<<" "<<D1.nrows<<" " <<D1.nrows<<" "<<D1.ncols<<" "<<D1.nnz<<"\n";
    //std::cout<<"Size P"<<P.L.bsize<<" "<<P.L.nrows<<" " <<P.L.nrows<<" "<<P.L.ncols<<" "<<P.L.nnz<<"\n";
    bsr_matrix *L= &P.L;
    bsr_matrix lt;           
    bsr_transpose(L, &lt);    
    //std::cout<<"Size lt"<<lt.bsize<<" "<<lt.nrows<<" " <<lt.nrows<<" "<<lt.ncols<<" "<<lt.nnz<<"\n";
    //&P->L = &lt;
    //bsr_show(&P.D, "L POST ildu");
    //exit(0);
  
   // int mb   = P->U.nrows;
    //int nnzb = P->U.nnz;
    //
    //std::vector <double> ld(P->L.data, P->L.data + P->L.nnz * block_dim*block_dim);
    //std::vector <int> lr(P->L.rowptr, P->L.rowptr + P->L.nrows);
    //std::vector <int> lc(P->L.colidx, P->L.colidx + P->L.nnz);
    
    int block_dim =3;
    std::vector <double> ld(lt.data, lt.data + lt.nnz * block_dim*block_dim);
    std::vector <int> lr(lt.rowptr, lt.rowptr + (lt.nrows+1));  //ideally it should be lt.nrows +1 but we are not including the first non zero in 
    std::vector <int> lc(lt.colidx, lt.colidx + lt.nnz);
    std::vector <double> hdval (P.D.data, P.D.data + P.D.nnz * block_dim*block_dim);
    //std::cout<<"size of hdval is "<<hdval.size()<<std::endl;
    
    std::vector <double> ud(P.U.data, P.U.data + P.U.nnz * block_dim*block_dim);
    std::vector <int> ur(P.U.rowptr, P.U.rowptr + (P.U.nrows+1));
    //std::cout<<"size of rows in PU is "<<P.U.nrows<<std::endl;
    std::vector <int> uc(P.U.colidx, P.U.colidx + P.U.nnz);
    
    int rows_per_subdomain = nbx*nby*nbz;
    decompress_l(ld,lr,lc, hbsr_plval,hbsr_plrow_ptr, hbsr_plcol_ind,rows_per_subdomain);
    decompress_u(ud,ur,uc, hbsr_puval,hbsr_purow_ptr, hbsr_pucol_ind,rows_per_subdomain);  
      
    /*
    for(int i=0;i<(F.nrows+1);i++)
      hbsr_prow_ptr.push_back(F.rowptr[i]);
    for(int i=0;i<F.nnz;i++)
      hbsr_pcol_ind.push_back(F.colidx[i]);
    for(int i=0;i<F.nnz*bsize*bsize;i++)
      hbsr_pval.push_back(F.data[i]);  
      */
    
   // int bsize;
   // int nrows;
   // int ncols;
   // int nnz;
   
    
    //exit(0);
    
   // std::vector<double> hy_og(mb*bsize, 1.0f);
    std::vector<double> hy(mb*bsize, 0.0f); 

    
    
    //apply row permutation to vector
    for (int i=0;i<mb;i++){
      for (int j=0;j<bsize;j++)
      {
        hy[i*bsize+j]=hy_og[PN.weight[i]*bsize+j];
      }  
    }
    
    
    std::vector <int> subdomain_offsets(num_subdomains+1, 0);
    std::vector <int> llevels_per_subdomain(num_subdomains, 0);
    std::vector <int> level_start_subdomain(num_subdomains+1, 0);
     
    std::vector <int> llevelsu_per_subdomain(num_subdomains, 0);
    std::vector <int> levelu_start_subdomain(num_subdomains+1, 0);
    
  
    //update this loop for non-uniform subdomain size
    for(int i=0;i<num_subdomains;i++)
    {
      subdomain_offsets[i+1]=subdomain_offsets[i]+rows_per_subdomain;
    }
      
    bsr->num_subdomains = num_subdomains;
    bsr->rows_per_subdomain = rows_per_subdomain;
    bsr->subdomain_offsets = subdomain_offsets;
    
    //subdomain_offsets=subdomain_offsets;
    bsr->llevels_per_subdomain=llevels_per_subdomain;
    bsr->level_start_subdomain=level_start_subdomain;
    bsr->llevelsu_per_subdomain=llevelsu_per_subdomain;
    bsr->levelu_start_subdomain=levelu_start_subdomain;
    
    bsr->hbsr_bval=hy;
    bsr->hdval=hdval;
    bsr->hbsr_val=hbsr_val; bsr->hbsr_row_ptr=hbsr_row_ptr; bsr->hbsr_col_ind=hbsr_col_ind;
    bsr->hbsr_plval=hbsr_plval; bsr->hbsr_plrow_ptr=hbsr_plrow_ptr; bsr->hbsr_plcol_ind=hbsr_plcol_ind;
    bsr->hbsr_puval=hbsr_puval; bsr->hbsr_purow_ptr=hbsr_purow_ptr; bsr->hbsr_pucol_ind=hbsr_pucol_ind;
    //std::cout<<"hval bsr size "<<hbsr_val.size()<<std::endl;
    
    bsr->mb = bsr->hbsr_row_ptr.size()-1;
    bsr->nnzb = bsr->hbsr_col_ind.size();
    bsr->nnz =  bsr->hbsr_val.size();
    bsr->block_size = bsr->hbsr_bval.size()/bsr->mb;
    bsr->Nb = bsr->mb;
    bsr->nnzbs_prec= F.nnz;
    bsr->N = bsr->hbsr_bval.size();
    
    //std::cout << "created bsr with #rows "<< bsr->mb <<" #cols "<<bsr->nnzb << " block size " <<  bsr->block_size<< std::endl;
    
    
      
    init_gpu_vc(bsr);
    data_xfers_vc(bsr);
    rocsparse_boilerplate_vc(bsr);
    //preconditioner(bsr);
    auto start = std::chrono::system_clock::now();
    gpu_pbicgstab_nolds_vc(bsr,tolerance, maxit);
    auto end = std::chrono::system_clock::now();
    double elapsed_ms= std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "solver time "<<elapsed_ms << '\n';
   
    //print result
    
    std::vector<double> hx(bsr->N, 0.0f);
    hipMemcpy(hx.data(), bsr->d_x, sizeof(double) * bsr->N, hipMemcpyDeviceToHost);
    
    
    //std::cout<<"writing to "<<s2<<s<<std::endl;
    //indir << "xvals.fp64_"<<rows <<"_"<<nnzb<<"_"<<std::scientific << tolerance<<".txt";
    //std::cout<<"writing to "<<indir.str()<<std::endl;
 
    //std::ofstream output_file(indir.str());

    //std::ostream_iterator<double> output_iterator(output_file, "\n");
    //std::copy(std::begin(hx), std::end(hx), output_iterator);
    
    //apply inverse permutation answers
    int *imap = (int*)malloc(bsr->mb*sizeof(int));
    for(int i=0;i<bsr->mb;i++) imap[PN.weight[i]]=i;
    //std::vector<double> hy(mb*bsize, 0.0f);
    
    std::cout << "actual edges "<< C->nnz<<" post partitionning  edges " << D1.nnz << " post DD edges "<<DD.post_dd_edges<<std::endl;
    //std::cout << "applying inverse permutation to answer vector " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
      int ipos = imap[i];
      for (int j=0;j<bsize;j++)
      {
        ans[i*bsize+j]=hx[ipos*bsize+j];
      }  
    //std::cout<<hx[i]<<" "; 
    }
    
    
    /*
    std::cout << "hx " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
  	std::cout<<hx[i]<<" "; 
    }
    */
    std::cout<<"\n";

  

}


void run_solver_graph_ildu_vc_fp32(bsr_matrix *C, solver_data_vc *bsr, int nx, int ny, int nz, int bsize, int nbx, int nby, int nbz, float maxit, double tolerance, std::vector<double> &hy_og, std::vector<double> &ans, int num_subdomains)
{
   
    std::vector<int> hbsr_row_ptr;
    std::vector<int> hbsr_col_ind;
    std::vector<double> hbsr_val;
    
    std::vector<int> hbsr_plrow_ptr;
    std::vector<int> hbsr_plcol_ind;
    std::vector<double> hbsr_plval;
    
    std::vector<int> hbsr_purow_ptr;
    std::vector<int> hbsr_pucol_ind;
    std::vector<double> hbsr_puval;
    
    

    int n=nx*ny*nz;
    int nnz= nx*ny*nz +2*(nx-1)*ny*nz + 2*nx*(ny-1)*nz + 2*nx*ny*(nz-1);
    int b=bsize;
    int bb=b*b;
    int total_nnz = bb*nnz;  
    int* colidx;
    int* rowptr;
    double* data;
    int mb =n;
    //int rows_per_subdomain = nbx*nby*nbz;
    //colidx=(int*)malloc(nnz*sizeof(int));
    //rowptr=(int*)malloc((n+1)*sizeof(int));
    //data=(double*)malloc(bb*nnz*sizeof(double));
    
    
    //bsr_matrix C, D;
    bsr_cartesian(C, nx, ny, nz, bsize);

    //flux balance
    bsr_flux_balance(C);

    //boundary condititions
    bsr_dirichlet(C,0);
    bsr_dirichlet(C,n-1);
    //bsr_show(&C,"C");
    bsr_matrix  D1;
    //int nbx=16;
    //int nby=16;
    //int nbz=4;
    
   
    
    
    int *pid = (int*)malloc(n*sizeof(int));
    
    vec_ifill(pid,n,-1);
    auto start1 = std::chrono::system_clock::now();
    bsr_graph_part(nx, ny, nz, nbx, nby, nbz, pid);
    //vec_ishow(pid,n,"pid");
    
	  //reinterpret partition ids as a graph
    bsr_graph PT, PN, PC, A;
    bsr_graph_tags(&PT, pid, n);
    //bsr_graph_show(&PT,"PT");
	
	  //apply transpose to generate node permutation graph
    bsr_graph_transpose(&PT,&PN);
    //bsr_graph_show(&PN,"PN");
	
	  //extract graph of sparse matrix C and apply node permutation
    // --> the end results is a connection permutation graph
	  bsr_to_graph(C,&A);
    bsr_graph_reorder(&A, &(PN.weight), &PC);
    
    

    
	  //bsr_graph_show(&PC,"PC");
	  //    vec_ishow(PC.weight,A.nnz,"PC");

    // apply connection permutation to sparse matrix C
    // --> to get the reordered sparse matrix D
    bsr_reorder(C, &PC, &D1);
    bsr_graph B, DD;
    bsr_to_graph(&D1,&B);
    bsr_graph_decompose(&B, &PN, &DD);
    //csr_graph_show(&B,"B");
    //bsr_graph_show(&DD,"DD");

    // apply permutation
    // label conn to drop
    // apply projection
    
    DD.nnz = DD.post_dd_edges;
    bsr_matrix F;
    bsr_reorder(&D1, &DD, &F);
    //bsr_show(&F,"F");
    
    
    F.nnz = DD.post_dd_edges;
    auto end1 = std::chrono::system_clock::now();
    double elapsed_ms1= std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1).count();
    std::cout << "partition + decomposition time  "<<elapsed_ms1 << '\n';
    
    
    for(int i=0;i<(D1.nrows+1);i++)
      hbsr_row_ptr.push_back(D1.rowptr[i]);
    for(int i=0;i<D1.nnz;i++)
      hbsr_col_ind.push_back(D1.colidx[i]);
    for(int i=0;i<D1.nnz*bsize*bsize;i++)
      hbsr_val.push_back(D1.data[i]);
     
    //std::cout<<"Size "<<F.bsize<<" "<<F.nrows<<" " <<F.nrows<<" "<<F.ncols<<" "<<F.nnz<<"\n";
    bildu_prec P;
    bildu_create(&F, &P);
     
    //std::cout<<"Size D"<<D1.bsize<<" "<<D1.nrows<<" " <<D1.nrows<<" "<<D1.ncols<<" "<<D1.nnz<<"\n";
    //std::cout<<"Size P"<<P.L.bsize<<" "<<P.L.nrows<<" " <<P.L.nrows<<" "<<P.L.ncols<<" "<<P.L.nnz<<"\n";
    bsr_matrix *L= &P.L;
    bsr_matrix lt;           
    bsr_transpose(L, &lt);    
    //std::cout<<"Size lt"<<lt.bsize<<" "<<lt.nrows<<" " <<lt.nrows<<" "<<lt.ncols<<" "<<lt.nnz<<"\n";
    //&P->L = &lt;
    //bsr_show(&P.D, "L POST ildu");
    //exit(0);
  
   // int mb   = P->U.nrows;
    //int nnzb = P->U.nnz;
    //
    //std::vector <double> ld(P->L.data, P->L.data + P->L.nnz * block_dim*block_dim);
    //std::vector <int> lr(P->L.rowptr, P->L.rowptr + P->L.nrows);
    //std::vector <int> lc(P->L.colidx, P->L.colidx + P->L.nnz);
    
    int block_dim =3;
    std::vector <double> ld(lt.data, lt.data + lt.nnz * block_dim*block_dim);
    std::vector <int> lr(lt.rowptr, lt.rowptr + (lt.nrows+1));  //ideally it should be lt.nrows +1 but we are not including the first non zero in 
    std::vector <int> lc(lt.colidx, lt.colidx + lt.nnz);
    std::vector <double> hdval (P.D.data, P.D.data + P.D.nnz * block_dim*block_dim);
    //std::cout<<"size of hdval is "<<hdval.size()<<std::endl;
    
    std::vector <double> ud(P.U.data, P.U.data + P.U.nnz * block_dim*block_dim);
    std::vector <int> ur(P.U.rowptr, P.U.rowptr + (P.U.nrows+1));
    //std::cout<<"size of rows in PU is "<<P.U.nrows<<std::endl;
    std::vector <int> uc(P.U.colidx, P.U.colidx + P.U.nnz);
    
    int rows_per_subdomain = nbx*nby*nbz;
    decompress_l(ld,lr,lc, hbsr_plval,hbsr_plrow_ptr, hbsr_plcol_ind,rows_per_subdomain);
    decompress_u(ud,ur,uc, hbsr_puval,hbsr_purow_ptr, hbsr_pucol_ind,rows_per_subdomain);  
      
    /*
    for(int i=0;i<(F.nrows+1);i++)
      hbsr_prow_ptr.push_back(F.rowptr[i]);
    for(int i=0;i<F.nnz;i++)
      hbsr_pcol_ind.push_back(F.colidx[i]);
    for(int i=0;i<F.nnz*bsize*bsize;i++)
      hbsr_pval.push_back(F.data[i]);  
      */
    
   // int bsize;
   // int nrows;
   // int ncols;
   // int nnz;
   
    
    //exit(0);
    
   // std::vector<double> hy_og(mb*bsize, 1.0f);
    std::vector<double> hy(mb*bsize, 0.0f); 

    
    
    //apply row permutation to vector
    for (int i=0;i<mb;i++){
      for (int j=0;j<bsize;j++)
      {
        hy[i*bsize+j]=hy_og[PN.weight[i]*bsize+j];
      }  
    }
    
    
    std::vector <int> subdomain_offsets(num_subdomains+1, 0);
    std::vector <int> llevels_per_subdomain(num_subdomains, 0);
    std::vector <int> level_start_subdomain(num_subdomains+1, 0);
     
    std::vector <int> llevelsu_per_subdomain(num_subdomains, 0);
    std::vector <int> levelu_start_subdomain(num_subdomains+1, 0);
    
  
    //update this loop for non-uniform subdomain size
    for(int i=0;i<num_subdomains;i++)
    {
      subdomain_offsets[i+1]=subdomain_offsets[i]+rows_per_subdomain;
    }
      
    bsr->num_subdomains = num_subdomains;
    bsr->rows_per_subdomain = rows_per_subdomain;
    bsr->subdomain_offsets = subdomain_offsets;
    
    //subdomain_offsets=subdomain_offsets;
    bsr->llevels_per_subdomain=llevels_per_subdomain;
    bsr->level_start_subdomain=level_start_subdomain;
    bsr->llevelsu_per_subdomain=llevelsu_per_subdomain;
    bsr->levelu_start_subdomain=levelu_start_subdomain;
    
    bsr->hbsr_bval=hy;
    bsr->hdval=hdval;
    bsr->hbsr_val=hbsr_val; bsr->hbsr_row_ptr=hbsr_row_ptr; bsr->hbsr_col_ind=hbsr_col_ind;
    bsr->hbsr_plval=hbsr_plval; bsr->hbsr_plrow_ptr=hbsr_plrow_ptr; bsr->hbsr_plcol_ind=hbsr_plcol_ind;
    bsr->hbsr_puval=hbsr_puval; bsr->hbsr_purow_ptr=hbsr_purow_ptr; bsr->hbsr_pucol_ind=hbsr_pucol_ind;
    //std::cout<<"hval bsr size "<<hbsr_val.size()<<std::endl;
    
    bsr->mb = bsr->hbsr_row_ptr.size()-1;
    bsr->nnzb = bsr->hbsr_col_ind.size();
    bsr->nnz =  bsr->hbsr_val.size();
    bsr->block_size = bsr->hbsr_bval.size()/bsr->mb;
    bsr->Nb = bsr->mb;
    bsr->nnzbs_prec= F.nnz;
    bsr->N = bsr->hbsr_bval.size();
    
    //std::cout << "created bsr with #rows "<< bsr->mb <<" #cols "<<bsr->nnzb << " block size " <<  bsr->block_size<< std::endl;
    
    
      
    init_gpu_vc(bsr);
    data_xfers_vc(bsr);
    rocsparse_boilerplate_vc(bsr);
    //preconditioner(bsr);
    auto start = std::chrono::system_clock::now();
     gpu_pbicgstab_vc_fp32(bsr,tolerance, maxit);
    auto end = std::chrono::system_clock::now();
    double elapsed_ms= std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "solver time "<<elapsed_ms << '\n';
   
    //print result
    
    std::vector<double> hx(bsr->N, 0.0f);
    hipMemcpy(hx.data(), bsr->d_x, sizeof(double) * bsr->N, hipMemcpyDeviceToHost);
    
    
    //std::cout<<"writing to "<<s2<<s<<std::endl;
    //indir << "xvals.fp64_"<<rows <<"_"<<nnzb<<"_"<<std::scientific << tolerance<<".txt";
    //std::cout<<"writing to "<<indir.str()<<std::endl;
 
    //std::ofstream output_file(indir.str());

    //std::ostream_iterator<double> output_iterator(output_file, "\n");
    //std::copy(std::begin(hx), std::end(hx), output_iterator);
    
    //apply inverse permutation answers
    int *imap = (int*)malloc(bsr->mb*sizeof(int));
    for(int i=0;i<bsr->mb;i++) imap[PN.weight[i]]=i;
    //std::vector<double> hy(mb*bsize, 0.0f);
    
    std::cout << "actual edges "<< C->nnz<<" post partitionning  edges " << D1.nnz << " post DD edges "<<DD.post_dd_edges<<std::endl;
    //std::cout << "applying inverse permutation to answer vector " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
      int ipos = imap[i];
      for (int j=0;j<bsize;j++)
      {
        ans[i*bsize+j]=hx[ipos*bsize+j];
      }  
    //std::cout<<hx[i]<<" "; 
    }
    
    
    /*
    std::cout << "hx " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
  	std::cout<<hx[i]<<" "; 
    }
    */
    std::cout<<"\n";

  

}



void run_solver_graph_ildu_vc(bsr_matrix *C, solver_data_vc *bsr, int nx, int ny, int nz, int bsize, int nbx, int nby, int nbz, float maxit, double tolerance, std::vector<double> &hy_og, std::vector<double> &ans, int num_subdomains)
{
   
    std::vector<int> hbsr_row_ptr;
    std::vector<int> hbsr_col_ind;
    std::vector<double> hbsr_val;
    
    std::vector<int> hbsr_plrow_ptr;
    std::vector<int> hbsr_plcol_ind;
    std::vector<double> hbsr_plval;
    
    std::vector<int> hbsr_purow_ptr;
    std::vector<int> hbsr_pucol_ind;
    std::vector<double> hbsr_puval;
    
    

    int n=nx*ny*nz;
    int nnz= nx*ny*nz +2*(nx-1)*ny*nz + 2*nx*(ny-1)*nz + 2*nx*ny*(nz-1);
    int b=bsize;
    int bb=b*b;
    int total_nnz = bb*nnz;  
    int* colidx;
    int* rowptr;
    double* data;
    int mb =n;
    //int rows_per_subdomain = nbx*nby*nbz;
    //colidx=(int*)malloc(nnz*sizeof(int));
    //rowptr=(int*)malloc((n+1)*sizeof(int));
    //data=(double*)malloc(bb*nnz*sizeof(double));
    
    
    //bsr_matrix C, D;
    bsr_cartesian(C, nx, ny, nz, bsize);

    //flux balance
    bsr_flux_balance(C);

    //boundary condititions
    bsr_dirichlet(C,0);
    bsr_dirichlet(C,n-1);
    //bsr_show(&C,"C");
    bsr_matrix  D1;
    //int nbx=16;
    //int nby=16;
    //int nbz=4;
    
   
    
    
    int *pid = (int*)malloc(n*sizeof(int));
    
    vec_ifill(pid,n,-1);
    auto start1 = std::chrono::system_clock::now();
    bsr_graph_part(nx, ny, nz, nbx, nby, nbz, pid);
    //vec_ishow(pid,n,"pid");
    
	  //reinterpret partition ids as a graph
    bsr_graph PT, PN, PC, A;
    bsr_graph_tags(&PT, pid, n);
    //bsr_graph_show(&PT,"PT");
	
	  //apply transpose to generate node permutation graph
    bsr_graph_transpose(&PT,&PN);
    //bsr_graph_show(&PN,"PN");
	
	  //extract graph of sparse matrix C and apply node permutation
    // --> the end results is a connection permutation graph
	  bsr_to_graph(C,&A);
    bsr_graph_reorder(&A, &(PN.weight), &PC);
    
    

    
	  //bsr_graph_show(&PC,"PC");
	  //    vec_ishow(PC.weight,A.nnz,"PC");

    // apply connection permutation to sparse matrix C
    // --> to get the reordered sparse matrix D
    bsr_reorder(C, &PC, &D1);
    bsr_graph B, DD;
    bsr_to_graph(&D1,&B);
    bsr_graph_decompose(&B, &PN, &DD);
    //csr_graph_show(&B,"B");
    //bsr_graph_show(&DD,"DD");

    // apply permutation
    // label conn to drop
    // apply projection
    
    DD.nnz = DD.post_dd_edges;
    bsr_matrix F;
    bsr_reorder(&D1, &DD, &F);
    //bsr_show(&F,"F");
    
    
    F.nnz = DD.post_dd_edges;
    auto end1 = std::chrono::system_clock::now();
    double elapsed_ms1= std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1).count();
    std::cout << "partition + decomposition time  "<<elapsed_ms1 << '\n';
    
    
    for(int i=0;i<(D1.nrows+1);i++)
      hbsr_row_ptr.push_back(D1.rowptr[i]);
    for(int i=0;i<D1.nnz;i++)
      hbsr_col_ind.push_back(D1.colidx[i]);
    for(int i=0;i<D1.nnz*bsize*bsize;i++)
      hbsr_val.push_back(D1.data[i]);
     
    //std::cout<<"Size "<<F.bsize<<" "<<F.nrows<<" " <<F.nrows<<" "<<F.ncols<<" "<<F.nnz<<"\n";
    bildu_prec P;
    bildu_create(&F, &P);
     
    //std::cout<<"Size D"<<D1.bsize<<" "<<D1.nrows<<" " <<D1.nrows<<" "<<D1.ncols<<" "<<D1.nnz<<"\n";
    //std::cout<<"Size P"<<P.L.bsize<<" "<<P.L.nrows<<" " <<P.L.nrows<<" "<<P.L.ncols<<" "<<P.L.nnz<<"\n";
    bsr_matrix *L= &P.L;
    bsr_matrix lt;           
    bsr_transpose(L, &lt);    
    //std::cout<<"Size lt"<<lt.bsize<<" "<<lt.nrows<<" " <<lt.nrows<<" "<<lt.ncols<<" "<<lt.nnz<<"\n";
    //&P->L = &lt;
    //bsr_show(&P.D, "L POST ildu");
    //exit(0);
  
   // int mb   = P->U.nrows;
    //int nnzb = P->U.nnz;
    //
    //std::vector <double> ld(P->L.data, P->L.data + P->L.nnz * block_dim*block_dim);
    //std::vector <int> lr(P->L.rowptr, P->L.rowptr + P->L.nrows);
    //std::vector <int> lc(P->L.colidx, P->L.colidx + P->L.nnz);
    
    int block_dim =3;
    std::vector <double> ld(lt.data, lt.data + lt.nnz * block_dim*block_dim);
    std::vector <int> lr(lt.rowptr, lt.rowptr + (lt.nrows+1));  //ideally it should be lt.nrows +1 but we are not including the first non zero in 
    std::vector <int> lc(lt.colidx, lt.colidx + lt.nnz);
    std::vector <double> hdval (P.D.data, P.D.data + P.D.nnz * block_dim*block_dim);
    //std::cout<<"size of hdval is "<<hdval.size()<<std::endl;
    
    std::vector <double> ud(P.U.data, P.U.data + P.U.nnz * block_dim*block_dim);
    std::vector <int> ur(P.U.rowptr, P.U.rowptr + (P.U.nrows+1));
    //std::cout<<"size of rows in PU is "<<P.U.nrows<<std::endl;
    std::vector <int> uc(P.U.colidx, P.U.colidx + P.U.nnz);
    
    int rows_per_subdomain = nbx*nby*nbz;
    decompress_l(ld,lr,lc, hbsr_plval,hbsr_plrow_ptr, hbsr_plcol_ind,rows_per_subdomain);
    decompress_u(ud,ur,uc, hbsr_puval,hbsr_purow_ptr, hbsr_pucol_ind,rows_per_subdomain);  
      
    /*
    for(int i=0;i<(F.nrows+1);i++)
      hbsr_prow_ptr.push_back(F.rowptr[i]);
    for(int i=0;i<F.nnz;i++)
      hbsr_pcol_ind.push_back(F.colidx[i]);
    for(int i=0;i<F.nnz*bsize*bsize;i++)
      hbsr_pval.push_back(F.data[i]);  
      */
    
   // int bsize;
   // int nrows;
   // int ncols;
   // int nnz;
   
    
    //exit(0);
    
   // std::vector<double> hy_og(mb*bsize, 1.0f);
    std::vector<double> hy(mb*bsize, 0.0f); 

    
    
    //apply row permutation to vector
    for (int i=0;i<mb;i++){
      for (int j=0;j<bsize;j++)
      {
        hy[i*bsize+j]=hy_og[PN.weight[i]*bsize+j];
      }  
    }
    
    
    std::vector <int> subdomain_offsets(num_subdomains+1, 0);
    std::vector <int> llevels_per_subdomain(num_subdomains, 0);
    std::vector <int> level_start_subdomain(num_subdomains+1, 0);
     
    std::vector <int> llevelsu_per_subdomain(num_subdomains, 0);
    std::vector <int> levelu_start_subdomain(num_subdomains+1, 0);
    
  
    //update this loop for non-uniform subdomain size
    for(int i=0;i<num_subdomains;i++)
    {
      subdomain_offsets[i+1]=subdomain_offsets[i]+rows_per_subdomain;
    }
      
    bsr->num_subdomains = num_subdomains;
    bsr->rows_per_subdomain = rows_per_subdomain;
    bsr->subdomain_offsets = subdomain_offsets;
    
    //subdomain_offsets=subdomain_offsets;
    bsr->llevels_per_subdomain=llevels_per_subdomain;
    bsr->level_start_subdomain=level_start_subdomain;
    bsr->llevelsu_per_subdomain=llevelsu_per_subdomain;
    bsr->levelu_start_subdomain=levelu_start_subdomain;
    
    bsr->hbsr_bval=hy;
    bsr->hdval=hdval;
    bsr->hbsr_val=hbsr_val; bsr->hbsr_row_ptr=hbsr_row_ptr; bsr->hbsr_col_ind=hbsr_col_ind;
    bsr->hbsr_plval=hbsr_plval; bsr->hbsr_plrow_ptr=hbsr_plrow_ptr; bsr->hbsr_plcol_ind=hbsr_plcol_ind;
    bsr->hbsr_puval=hbsr_puval; bsr->hbsr_purow_ptr=hbsr_purow_ptr; bsr->hbsr_pucol_ind=hbsr_pucol_ind;
    //std::cout<<"hval bsr size "<<hbsr_val.size()<<std::endl;
    
    bsr->mb = bsr->hbsr_row_ptr.size()-1;
    bsr->nnzb = bsr->hbsr_col_ind.size();
    bsr->nnz =  bsr->hbsr_val.size();
    bsr->block_size = bsr->hbsr_bval.size()/bsr->mb;
    bsr->Nb = bsr->mb;
    bsr->nnzbs_prec= F.nnz;
    bsr->N = bsr->hbsr_bval.size();
    
    //std::cout << "created bsr with #rows "<< bsr->mb <<" #cols "<<bsr->nnzb << " block size " <<  bsr->block_size<< std::endl;
    
    
      
    init_gpu_vc(bsr);
    data_xfers_vc(bsr);
    rocsparse_boilerplate_vc(bsr);
    //preconditioner(bsr);
    auto start = std::chrono::system_clock::now();
     gpu_pbicgstab_vc(bsr,tolerance, maxit);
    auto end = std::chrono::system_clock::now();
    double elapsed_ms= std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "solver time "<<elapsed_ms << '\n';
   
    //print result
    
    std::vector<double> hx(bsr->N, 0.0f);
    hipMemcpy(hx.data(), bsr->d_x, sizeof(double) * bsr->N, hipMemcpyDeviceToHost);
    
    
    //std::cout<<"writing to "<<s2<<s<<std::endl;
    //indir << "xvals.fp64_"<<rows <<"_"<<nnzb<<"_"<<std::scientific << tolerance<<".txt";
    //std::cout<<"writing to "<<indir.str()<<std::endl;
 
    //std::ofstream output_file(indir.str());

    //std::ostream_iterator<double> output_iterator(output_file, "\n");
    //std::copy(std::begin(hx), std::end(hx), output_iterator);
    
    //apply inverse permutation answers
    int *imap = (int*)malloc(bsr->mb*sizeof(int));
    for(int i=0;i<bsr->mb;i++) imap[PN.weight[i]]=i;
    //std::vector<double> hy(mb*bsize, 0.0f);
    
    std::cout << "actual edges "<< C->nnz<<" post partitionning  edges " << D1.nnz << " post DD edges "<<DD.post_dd_edges<<std::endl;
    //std::cout << "applying inverse permutation to answer vector " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
      int ipos = imap[i];
      for (int j=0;j<bsize;j++)
      {
        ans[i*bsize+j]=hx[ipos*bsize+j];
      }  
    //std::cout<<hx[i]<<" "; 
    }
    
    
    /*
    std::cout << "hx " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
  	std::cout<<hx[i]<<" "; 
    }
    */
    std::cout<<"\n";

  

}




void run_solver_graph_ildu(bsr_matrix *C, solver_data2 *bsr, int nx, int ny, int nz, int bsize, int nbx, int nby, int nbz, float maxit, double tolerance, std::vector<double> &hy_og, std::vector<double> &ans, int num_subdomains)
{
   
    std::vector<int> hbsr_row_ptr;
    std::vector<int> hbsr_col_ind;
    std::vector<double> hbsr_val;
    
    std::vector<int> hbsr_plrow_ptr;
    std::vector<int> hbsr_plcol_ind;
    std::vector<double> hbsr_plval;
    
    std::vector<int> hbsr_purow_ptr;
    std::vector<int> hbsr_pucol_ind;
    std::vector<double> hbsr_puval;
    
    

    int n=nx*ny*nz;
    int nnz= nx*ny*nz +2*(nx-1)*ny*nz + 2*nx*(ny-1)*nz + 2*nx*ny*(nz-1);
    int b=bsize;
    int bb=b*b;
    int total_nnz = bb*nnz;  
    int* colidx;
    int* rowptr;
    double* data;
    int mb =n;
    //int rows_per_subdomain = nbx*nby*nbz;
    //colidx=(int*)malloc(nnz*sizeof(int));
    //rowptr=(int*)malloc((n+1)*sizeof(int));
    //data=(double*)malloc(bb*nnz*sizeof(double));
    
    
    //bsr_matrix C, D;
    bsr_cartesian(C, nx, ny, nz, bsize);

    //flux balance
    bsr_flux_balance(C);

    //boundary condititions
    bsr_dirichlet(C,0);
    bsr_dirichlet(C,n-1);
    //bsr_show(&C,"C");
    bsr_matrix  D1;
    //int nbx=16;
    //int nby=16;
    //int nbz=4;
    
   
    
    
    int *pid = (int*)malloc(n*sizeof(int));
    
    vec_ifill(pid,n,-1);
    auto start1 = std::chrono::system_clock::now();
    bsr_graph_part(nx, ny, nz, nbx, nby, nbz, pid);
    //vec_ishow(pid,n,"pid");
    
	  //reinterpret partition ids as a graph
    bsr_graph PT, PN, PC, A;
    bsr_graph_tags(&PT, pid, n);
    //bsr_graph_show(&PT,"PT");
	
	  //apply transpose to generate node permutation graph
    bsr_graph_transpose(&PT,&PN);
    //bsr_graph_show(&PN,"PN");
	
	  //extract graph of sparse matrix C and apply node permutation
    // --> the end results is a connection permutation graph
	  bsr_to_graph(C,&A);
    bsr_graph_reorder(&A, &(PN.weight), &PC);
    
    

    
	  //bsr_graph_show(&PC,"PC");
	  //    vec_ishow(PC.weight,A.nnz,"PC");

    // apply connection permutation to sparse matrix C
    // --> to get the reordered sparse matrix D
    bsr_reorder(C, &PC, &D1);
    bsr_graph B, DD;
    bsr_to_graph(&D1,&B);
    bsr_graph_decompose(&B, &PN, &DD);
    //csr_graph_show(&B,"B");
    //bsr_graph_show(&DD,"DD");

    // apply permutation
    // label conn to drop
    // apply projection
    
    DD.nnz = DD.post_dd_edges;
    bsr_matrix F;
    bsr_reorder(&D1, &DD, &F);
    //bsr_show(&F,"F");
    
    
    F.nnz = DD.post_dd_edges;
    auto end1 = std::chrono::system_clock::now();
    double elapsed_ms1= std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1).count();
    std::cout << "partition + decomposition time  "<<elapsed_ms1 << '\n';
    
    
    for(int i=0;i<(D1.nrows+1);i++)
      hbsr_row_ptr.push_back(D1.rowptr[i]);
    for(int i=0;i<D1.nnz;i++)
      hbsr_col_ind.push_back(D1.colidx[i]);
    for(int i=0;i<D1.nnz*bsize*bsize;i++)
      hbsr_val.push_back(D1.data[i]);
     
    //std::cout<<"Size "<<F.bsize<<" "<<F.nrows<<" " <<F.nrows<<" "<<F.ncols<<" "<<F.nnz<<"\n";
    bildu_prec P;
    bildu_create(&F, &P);
     
    //std::cout<<"Size D"<<D1.bsize<<" "<<D1.nrows<<" " <<D1.nrows<<" "<<D1.ncols<<" "<<D1.nnz<<"\n";
    //std::cout<<"Size P"<<P.L.bsize<<" "<<P.L.nrows<<" " <<P.L.nrows<<" "<<P.L.ncols<<" "<<P.L.nnz<<"\n";
    bsr_matrix *L= &P.L;
    bsr_matrix lt;           
    bsr_transpose(L, &lt);    
    //std::cout<<"Size lt"<<lt.bsize<<" "<<lt.nrows<<" " <<lt.nrows<<" "<<lt.ncols<<" "<<lt.nnz<<"\n";
    //&P->L = &lt;
    //bsr_show(&P.D, "L POST ildu");
    //exit(0);
  
   // int mb   = P->U.nrows;
    //int nnzb = P->U.nnz;
    //
    //std::vector <double> ld(P->L.data, P->L.data + P->L.nnz * block_dim*block_dim);
    //std::vector <int> lr(P->L.rowptr, P->L.rowptr + P->L.nrows);
    //std::vector <int> lc(P->L.colidx, P->L.colidx + P->L.nnz);
    
    int block_dim =3;
    std::vector <double> ld(lt.data, lt.data + lt.nnz * block_dim*block_dim);
    std::vector <int> lr(lt.rowptr, lt.rowptr + (lt.nrows+1));  //ideally it should be lt.nrows +1 but we are not including the first non zero in 
    std::vector <int> lc(lt.colidx, lt.colidx + lt.nnz);
    std::vector <double> hdval (P.D.data, P.D.data + P.D.nnz * block_dim*block_dim);
    //std::cout<<"size of hdval is "<<hdval.size()<<std::endl;
    
    std::vector <double> ud(P.U.data, P.U.data + P.U.nnz * block_dim*block_dim);
    std::vector <int> ur(P.U.rowptr, P.U.rowptr + (P.U.nrows+1));
    //std::cout<<"size of rows in PU is "<<P.U.nrows<<std::endl;
    std::vector <int> uc(P.U.colidx, P.U.colidx + P.U.nnz);
    
    int rows_per_subdomain = nbx*nby*nbz;
    decompress_l(ld,lr,lc, hbsr_plval,hbsr_plrow_ptr, hbsr_plcol_ind,rows_per_subdomain);
    decompress_u(ud,ur,uc, hbsr_puval,hbsr_purow_ptr, hbsr_pucol_ind,rows_per_subdomain);  
      
    /*
    for(int i=0;i<(F.nrows+1);i++)
      hbsr_prow_ptr.push_back(F.rowptr[i]);
    for(int i=0;i<F.nnz;i++)
      hbsr_pcol_ind.push_back(F.colidx[i]);
    for(int i=0;i<F.nnz*bsize*bsize;i++)
      hbsr_pval.push_back(F.data[i]);  
      */
    
   // int bsize;
   // int nrows;
   // int ncols;
   // int nnz;
   
    
    //exit(0);
    
   // std::vector<double> hy_og(mb*bsize, 1.0f);
    std::vector<double> hy(mb*bsize, 0.0f); 

    
    
    //apply row permutation to vector
    for (int i=0;i<mb;i++){
      for (int j=0;j<bsize;j++)
      {
        hy[i*bsize+j]=hy_og[PN.weight[i]*bsize+j];
      }  
    }
    
    
    std::vector <int> subdomain_offsets(num_subdomains+1, 0);
    std::vector <int> llevels_per_subdomain(num_subdomains, 0);
    std::vector <int> level_start_subdomain(num_subdomains+1, 0);
     
    std::vector <int> llevelsu_per_subdomain(num_subdomains, 0);
    std::vector <int> levelu_start_subdomain(num_subdomains+1, 0);
    
  
    //update this loop for non-uniform subdomain size
    for(int i=0;i<num_subdomains;i++)
    {
      subdomain_offsets[i+1]=subdomain_offsets[i]+rows_per_subdomain;
    }
      
    bsr->num_subdomains = num_subdomains;
    bsr->rows_per_subdomain = rows_per_subdomain;
    bsr->subdomain_offsets = subdomain_offsets;
    
    //subdomain_offsets=subdomain_offsets;
    bsr->llevels_per_subdomain=llevels_per_subdomain;
    bsr->level_start_subdomain=level_start_subdomain;
    bsr->llevelsu_per_subdomain=llevelsu_per_subdomain;
    bsr->levelu_start_subdomain=levelu_start_subdomain;
    
    bsr->hbsr_bval=hy;
    bsr->hdval=hdval;
    bsr->hbsr_val=hbsr_val; bsr->hbsr_row_ptr=hbsr_row_ptr; bsr->hbsr_col_ind=hbsr_col_ind;
    bsr->hbsr_plval=hbsr_plval; bsr->hbsr_plrow_ptr=hbsr_plrow_ptr; bsr->hbsr_plcol_ind=hbsr_plcol_ind;
    bsr->hbsr_puval=hbsr_puval; bsr->hbsr_purow_ptr=hbsr_purow_ptr; bsr->hbsr_pucol_ind=hbsr_pucol_ind;
    //std::cout<<"hval bsr size "<<hbsr_val.size()<<std::endl;
    
    bsr->mb = bsr->hbsr_row_ptr.size()-1;
    bsr->nnzb = bsr->hbsr_col_ind.size();
    bsr->nnz =  bsr->hbsr_val.size();
    bsr->block_size = bsr->hbsr_bval.size()/bsr->mb;
    bsr->Nb = bsr->mb;
    bsr->nnzbs_prec= F.nnz;
    bsr->N = bsr->hbsr_bval.size();
    
    //std::cout << "created bsr with #rows "<< bsr->mb <<" #cols "<<bsr->nnzb << " block size " <<  bsr->block_size<< std::endl;
    
    
      
    init_gpu_g(bsr);
    data_xfers_g(bsr);
    rocsparse_boilerplate_g(bsr);
    //preconditioner(bsr);
    auto start = std::chrono::system_clock::now();
     gpu_pbicgstab_g(bsr,tolerance, maxit);
    auto end = std::chrono::system_clock::now();
    double elapsed_ms= std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "solver time "<<elapsed_ms << '\n';
   
    //print result
    
    std::vector<double> hx(bsr->N, 0.0f);
    hipMemcpy(hx.data(), bsr->d_x, sizeof(double) * bsr->N, hipMemcpyDeviceToHost);
    
    
    //std::cout<<"writing to "<<s2<<s<<std::endl;
    //indir << "xvals.fp64_"<<rows <<"_"<<nnzb<<"_"<<std::scientific << tolerance<<".txt";
    //std::cout<<"writing to "<<indir.str()<<std::endl;
 
    //std::ofstream output_file(indir.str());

    //std::ostream_iterator<double> output_iterator(output_file, "\n");
    //std::copy(std::begin(hx), std::end(hx), output_iterator);
    
    //apply inverse permutation answers
    int *imap = (int*)malloc(bsr->mb*sizeof(int));
    for(int i=0;i<bsr->mb;i++) imap[PN.weight[i]]=i;
    //std::vector<double> hy(mb*bsize, 0.0f);
    
    std::cout << "actual edges "<< C->nnz<<" post partitionning  edges " << D1.nnz << " post DD edges "<<DD.post_dd_edges<<std::endl;
    //std::cout << "applying inverse permutation to answer vector " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
      int ipos = imap[i];
      for (int j=0;j<bsize;j++)
      {
        ans[i*bsize+j]=hx[ipos*bsize+j];
      }  
    //std::cout<<hx[i]<<" "; 
    }
    
    
    /*
    std::cout << "hx " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
  	std::cout<<hx[i]<<" "; 
    }
    */
    std::cout<<"\n";

  

}


void run_solver_graph_ildu_fp32(bsr_matrix *C, solver_data2 *bsr, int nx, int ny, int nz, int bsize, int nbx, int nby, int nbz, float maxit, double tolerance, std::vector<double> &hy_og, std::vector<double> &ans, int num_subdomains)
{
   
    std::vector<int> hbsr_row_ptr;
    std::vector<int> hbsr_col_ind;
    std::vector<double> hbsr_val;
    
    std::vector<int> hbsr_plrow_ptr;
    std::vector<int> hbsr_plcol_ind;
    std::vector<double> hbsr_plval;
    
    std::vector<int> hbsr_purow_ptr;
    std::vector<int> hbsr_pucol_ind;
    std::vector<double> hbsr_puval;
    
    

    int n=nx*ny*nz;
    int nnz= nx*ny*nz +2*(nx-1)*ny*nz + 2*nx*(ny-1)*nz + 2*nx*ny*(nz-1);
    int b=bsize;
    int bb=b*b;
    int total_nnz = bb*nnz;  
    int* colidx;
    int* rowptr;
    double* data;
    int mb =n;
    //int rows_per_subdomain = nbx*nby*nbz;
    //colidx=(int*)malloc(nnz*sizeof(int));
    //rowptr=(int*)malloc((n+1)*sizeof(int));
    //data=(double*)malloc(bb*nnz*sizeof(double));
    
    
    //bsr_matrix C, D;
    bsr_cartesian(C, nx, ny, nz, bsize);

    //flux balance
    bsr_flux_balance(C);

    //boundary condititions
    bsr_dirichlet(C,0);
    bsr_dirichlet(C,n-1);
    //bsr_show(&C,"C");
    bsr_matrix  D1;
    //int nbx=16;
    //int nby=16;
    //int nbz=4;
    
   
    
    
    int *pid = (int*)malloc(n*sizeof(int));
    
    vec_ifill(pid,n,-1);
    auto start1 = std::chrono::system_clock::now();
    bsr_graph_part(nx, ny, nz, nbx, nby, nbz, pid);
    //vec_ishow(pid,n,"pid");
    
	  //reinterpret partition ids as a graph
    bsr_graph PT, PN, PC, A;
    bsr_graph_tags(&PT, pid, n);
    //bsr_graph_show(&PT,"PT");
	
	  //apply transpose to generate node permutation graph
    bsr_graph_transpose(&PT,&PN);
    //bsr_graph_show(&PN,"PN");
	
	  //extract graph of sparse matrix C and apply node permutation
    // --> the end results is a connection permutation graph
	  bsr_to_graph(C,&A);
    bsr_graph_reorder(&A, &(PN.weight), &PC);
    
    

    
	  //bsr_graph_show(&PC,"PC");
	  //    vec_ishow(PC.weight,A.nnz,"PC");

    // apply connection permutation to sparse matrix C
    // --> to get the reordered sparse matrix D
    bsr_reorder(C, &PC, &D1);
    bsr_graph B, DD;
    bsr_to_graph(&D1,&B);
    bsr_graph_decompose(&B, &PN, &DD);
    //csr_graph_show(&B,"B");
    //bsr_graph_show(&DD,"DD");

    // apply permutation
    // label conn to drop
    // apply projection
    
    DD.nnz = DD.post_dd_edges;
    bsr_matrix F;
    bsr_reorder(&D1, &DD, &F);
    //bsr_show(&F,"F");
    
    
    F.nnz = DD.post_dd_edges;
    auto end1 = std::chrono::system_clock::now();
    double elapsed_ms1= std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1).count();
    std::cout << "partition + decomposition time  "<<elapsed_ms1 << '\n';
    
    
    for(int i=0;i<(D1.nrows+1);i++)
      hbsr_row_ptr.push_back(D1.rowptr[i]);
    for(int i=0;i<D1.nnz;i++)
      hbsr_col_ind.push_back(D1.colidx[i]);
    for(int i=0;i<D1.nnz*bsize*bsize;i++)
      hbsr_val.push_back(D1.data[i]);
     
    //std::cout<<"Size "<<F.bsize<<" "<<F.nrows<<" " <<F.nrows<<" "<<F.ncols<<" "<<F.nnz<<"\n";
    bildu_prec P;
    bildu_create(&F, &P);
     
    //std::cout<<"Size D"<<D1.bsize<<" "<<D1.nrows<<" " <<D1.nrows<<" "<<D1.ncols<<" "<<D1.nnz<<"\n";
    //std::cout<<"Size P"<<P.L.bsize<<" "<<P.L.nrows<<" " <<P.L.nrows<<" "<<P.L.ncols<<" "<<P.L.nnz<<"\n";
    bsr_matrix *L= &P.L;
    bsr_matrix lt;           
    bsr_transpose(L, &lt);    
    //std::cout<<"Size lt"<<lt.bsize<<" "<<lt.nrows<<" " <<lt.nrows<<" "<<lt.ncols<<" "<<lt.nnz<<"\n";
    //&P->L = &lt;
    //bsr_show(&P.D, "L POST ildu");
    //exit(0);
  
   // int mb   = P->U.nrows;
    //int nnzb = P->U.nnz;
    //
    //std::vector <double> ld(P->L.data, P->L.data + P->L.nnz * block_dim*block_dim);
    //std::vector <int> lr(P->L.rowptr, P->L.rowptr + P->L.nrows);
    //std::vector <int> lc(P->L.colidx, P->L.colidx + P->L.nnz);
    
    int block_dim =3;
    std::vector <double> ld(lt.data, lt.data + lt.nnz * block_dim*block_dim);
    std::vector <int> lr(lt.rowptr, lt.rowptr + (lt.nrows+1));  //ideally it should be lt.nrows +1 but we are not including the first non zero in 
    std::vector <int> lc(lt.colidx, lt.colidx + lt.nnz);
    std::vector <double> hdval (P.D.data, P.D.data + P.D.nnz * block_dim*block_dim);
    //std::cout<<"size of hdval is "<<hdval.size()<<std::endl;
    
    std::vector <double> ud(P.U.data, P.U.data + P.U.nnz * block_dim*block_dim);
    std::vector <int> ur(P.U.rowptr, P.U.rowptr + (P.U.nrows+1));
    //std::cout<<"size of rows in PU is "<<P.U.nrows<<std::endl;
    std::vector <int> uc(P.U.colidx, P.U.colidx + P.U.nnz);
    
    int rows_per_subdomain = nbx*nby*nbz;
    decompress_l(ld,lr,lc, hbsr_plval,hbsr_plrow_ptr, hbsr_plcol_ind,rows_per_subdomain);
    decompress_u(ud,ur,uc, hbsr_puval,hbsr_purow_ptr, hbsr_pucol_ind,rows_per_subdomain);  
      
    /*
    for(int i=0;i<(F.nrows+1);i++)
      hbsr_prow_ptr.push_back(F.rowptr[i]);
    for(int i=0;i<F.nnz;i++)
      hbsr_pcol_ind.push_back(F.colidx[i]);
    for(int i=0;i<F.nnz*bsize*bsize;i++)
      hbsr_pval.push_back(F.data[i]);  
      */
    
   // int bsize;
   // int nrows;
   // int ncols;
   // int nnz;
   
    
    //exit(0);
    
   // std::vector<double> hy_og(mb*bsize, 1.0f);
    std::vector<double> hy(mb*bsize, 0.0f); 

    
    
    //apply row permutation to vector
    for (int i=0;i<mb;i++){
      for (int j=0;j<bsize;j++)
      {
        hy[i*bsize+j]=hy_og[PN.weight[i]*bsize+j];
      }  
    }
    
    
    std::vector <int> subdomain_offsets(num_subdomains+1, 0);
    std::vector <int> llevels_per_subdomain(num_subdomains, 0);
    std::vector <int> level_start_subdomain(num_subdomains+1, 0);
     
    std::vector <int> llevelsu_per_subdomain(num_subdomains, 0);
    std::vector <int> levelu_start_subdomain(num_subdomains+1, 0);
    
  
    //update this loop for non-uniform subdomain size
    for(int i=0;i<num_subdomains;i++)
    {
      subdomain_offsets[i+1]=subdomain_offsets[i]+rows_per_subdomain;
    }
      
    bsr->num_subdomains = num_subdomains;
    bsr->rows_per_subdomain = rows_per_subdomain;
    bsr->subdomain_offsets = subdomain_offsets;
    
    //subdomain_offsets=subdomain_offsets;
    bsr->llevels_per_subdomain=llevels_per_subdomain;
    bsr->level_start_subdomain=level_start_subdomain;
    bsr->llevelsu_per_subdomain=llevelsu_per_subdomain;
    bsr->levelu_start_subdomain=levelu_start_subdomain;
    
    bsr->hbsr_bval=hy;
    bsr->hdval=hdval;
    bsr->hbsr_val=hbsr_val; bsr->hbsr_row_ptr=hbsr_row_ptr; bsr->hbsr_col_ind=hbsr_col_ind;
    bsr->hbsr_plval=hbsr_plval; bsr->hbsr_plrow_ptr=hbsr_plrow_ptr; bsr->hbsr_plcol_ind=hbsr_plcol_ind;
    bsr->hbsr_puval=hbsr_puval; bsr->hbsr_purow_ptr=hbsr_purow_ptr; bsr->hbsr_pucol_ind=hbsr_pucol_ind;
    //std::cout<<"hval bsr size "<<hbsr_val.size()<<std::endl;
    
    bsr->mb = bsr->hbsr_row_ptr.size()-1;
    bsr->nnzb = bsr->hbsr_col_ind.size();
    bsr->nnz =  bsr->hbsr_val.size();
    bsr->block_size = bsr->hbsr_bval.size()/bsr->mb;
    bsr->Nb = bsr->mb;
    bsr->nnzbs_prec= F.nnz;
    bsr->N = bsr->hbsr_bval.size();
    
    //std::cout << "created bsr with #rows "<< bsr->mb <<" #cols "<<bsr->nnzb << " block size " <<  bsr->block_size<< std::endl;
    
    
      
    init_gpu_g(bsr);
    data_xfers_g(bsr);
    rocsparse_boilerplate_g(bsr);
    //preconditioner(bsr);
    auto start = std::chrono::system_clock::now();
     gpu_pbicgstab_g_ec_fp32(bsr,tolerance, maxit);
    auto end = std::chrono::system_clock::now();
    double elapsed_ms= std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "solver time "<<elapsed_ms << '\n';
   
    //print result
    
    std::vector<double> hx(bsr->N, 0.0f);
    hipMemcpy(hx.data(), bsr->d_x, sizeof(double) * bsr->N, hipMemcpyDeviceToHost);
    
    
    //std::cout<<"writing to "<<s2<<s<<std::endl;
    //indir << "xvals.fp64_"<<rows <<"_"<<nnzb<<"_"<<std::scientific << tolerance<<".txt";
    //std::cout<<"writing to "<<indir.str()<<std::endl;
 
    //std::ofstream output_file(indir.str());

    //std::ostream_iterator<double> output_iterator(output_file, "\n");
    //std::copy(std::begin(hx), std::end(hx), output_iterator);
    
    //apply inverse permutation answers
    int *imap = (int*)malloc(bsr->mb*sizeof(int));
    for(int i=0;i<bsr->mb;i++) imap[PN.weight[i]]=i;
    //std::vector<double> hy(mb*bsize, 0.0f);
    
    std::cout << "actual edges "<< C->nnz<<" post partitionning  edges " << D1.nnz << " post DD edges "<<DD.post_dd_edges<<std::endl;
    //std::cout << "applying inverse permutation to answer vector " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
      int ipos = imap[i];
      for (int j=0;j<bsize;j++)
      {
        ans[i*bsize+j]=hx[ipos*bsize+j];
      }  
    //std::cout<<hx[i]<<" "; 
    }
    
    
    /*
    std::cout << "hx " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
  	std::cout<<hx[i]<<" "; 
    }
    */
    std::cout<<"\n";

  

}








void run_solver_unmodified(bsr_matrix *C, solver_data *bsr, int nx, int ny, int nz, int bsize, float maxit, double tolerance, std::vector<double> &hy, std::vector<double> &ans)
{
    std::vector<int> hbsr_row_ptr;
    std::vector<int> hbsr_col_ind;
    std::vector<double> hbsr_val;
    
    
    std::vector<int> hbsr_prow_ptr;
    std::vector<int> hbsr_pcol_ind;
    std::vector<double> hbsr_pval;
       
    int n=nx*ny*nz;
    int nnz= nx*ny*nz +2*(nx-1)*ny*nz + 2*nx*(ny-1)*nz + 2*nx*ny*(nz-1);
    int b=bsize;
    int bb=b*b;
    int total_nnz = bb*nnz;  
    int* colidx;
    int* rowptr;
    double* data;
    int mb =n;
    
    //colidx=(int*)malloc(nnz*sizeof(int));
    //rowptr=(int*)malloc((n+1)*sizeof(int));
    //data=(double*)malloc(bb*nnz*sizeof(double));
    
    
    //bsr_matrix C, D;
    bsr_cartesian(C, nx, ny, nz, bsize);

    //flux balance
    bsr_flux_balance(C);

    //boundary condititions
    bsr_dirichlet(C,0);
    bsr_dirichlet(C,n-1);
    //bsr_show(&C,"C");
    
    
    
    for(int i=0;i<=n;i++)
      hbsr_row_ptr.push_back(C->rowptr[i]);
    for(int i=0;i<nnz;i++)
      hbsr_col_ind.push_back(C->colidx[i]);
    for(int i=0;i<total_nnz;i++)
      hbsr_val.push_back(C->data[i]);
      
      
      
    for(int i=0;i<=n;i++)
      hbsr_prow_ptr.push_back(C->rowptr[i]);
    for(int i=0;i<nnz;i++)
      hbsr_pcol_ind.push_back(C->colidx[i]);
    for(int i=0;i<total_nnz;i++)
      hbsr_pval.push_back(C->data[i]);
    
    
  
    bsr->hbsr_bval=hy;
       
    
    bsr->hbsr_val=hbsr_val; bsr->hbsr_row_ptr=hbsr_row_ptr; bsr->hbsr_col_ind=hbsr_col_ind;
    bsr->hbsr_pval=hbsr_pval; bsr->hbsr_prow_ptr=hbsr_prow_ptr; bsr->hbsr_pcol_ind=hbsr_pcol_ind;
    
    //std::cout<<"hcol bsr size "<<hbsr_col_ind.size()<<std::endl;
    // std::cout<<"hcol bsr size "<<bsr->hbsr_pcol_ind.size()<<std::endl;
  
    bsr->mb = bsr->hbsr_row_ptr.size()-1;
    bsr->nnzb = bsr->hbsr_col_ind.size();
    bsr->nnz =  bsr->hbsr_val.size();
    bsr->block_size = bsr->hbsr_bval.size()/bsr->mb;
    bsr->Nb = bsr->mb;
    bsr->nnzbs_prec= bsr->nnzb;
    bsr->N = bsr->hbsr_bval.size();
    
    //std::cout << "created bsr with #rows "<< bsr->mb <<" #cols "<<bsr->nnzb << " block size " <<  bsr->block_size<< std::endl;
    
      
    init_gpu(bsr);
    data_xfers(bsr);
    rocsparse_boilerplate(bsr);
    preconditioner(bsr);
    
    auto start = std::chrono::system_clock::now();
    gpu_pbicgstab(bsr,tolerance, maxit);
    auto end = std::chrono::system_clock::now();
    double elapsed_ms= std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "solver time "<<elapsed_ms << '\n';
   
    //print result
    
    //std::vector<double> hx(bsr->N, 0.0f);
    hipMemcpy(ans.data(), bsr->d_x, sizeof(double) * bsr->N, hipMemcpyDeviceToHost);
    
    
    //std::cout<<"writing to "<<s2<<s<<std::endl;
    //indir << "xvals.fp64_"<<rows <<"_"<<nnzb<<"_"<<std::scientific << tolerance<<".txt";
    //std::cout<<"writing to "<<indir.str()<<std::endl;
 
    //std::ofstream output_file(indir.str());

    //std::ostream_iterator<double> output_iterator(output_file, "\n");
    //std::copy(std::begin(hx), std::end(hx), output_iterator);
    
    
    
    /*
    std::cout << "ans " <<std::endl;
    for (int i =0; i< bsr->N; i++)
    {
  	std::cout<<ans[i]<<" "; 
    }
    */
    std::cout<<"\n";
  

}

void run_solver_reordered(bsr_matrix *C, solver_data *bsr, int nx, int ny, int nz, int bsize, int nbx, int nby, int nbz, float maxit, double tolerance, std::vector<double> &hy_og, std::vector<double> &ans)
{

    std::vector<int> hbsr_row_ptr;
    std::vector<int> hbsr_col_ind;
    std::vector<double> hbsr_val;
    
    std::vector<int> hbsr_prow_ptr;
    std::vector<int> hbsr_pcol_ind;
    std::vector<double> hbsr_pval;

    int n=nx*ny*nz;
    int nnz= nx*ny*nz +2*(nx-1)*ny*nz + 2*nx*(ny-1)*nz + 2*nx*ny*(nz-1);
    int b=bsize;
    int bb=b*b;
    int total_nnz = bb*nnz;  
    int* colidx;
    int* rowptr;
    double* data;
    int mb =n;
    
    //colidx=(int*)malloc(nnz*sizeof(int));
    //rowptr=(int*)malloc((n+1)*sizeof(int));
    //data=(double*)malloc(bb*nnz*sizeof(double));
    
    
    //bsr_matrix C, D;
    bsr_cartesian(C, nx, ny, nz, bsize);

    //flux balance
    bsr_flux_balance(C);

    //boundary condititions
    bsr_dirichlet(C,0);
    bsr_dirichlet(C,n-1);
    //bsr_show(&C,"C");
    bsr_matrix  D;
    

    
     //int nbx=16;
    //int nby=16;
    //int nbz=4;
    
    
    
    int *pid = (int*)malloc(n*sizeof(int));
    vec_ifill(pid,n,-1);
    bsr_graph_part(nx, ny, nz, nbx, nby, nbz, pid);
    //metis function(matrix, max parition, )
    //vec_ishow(pid,n,"pid");
    
	  //reinterpret partition ids as a graph
    bsr_graph PT, PN, PC, A;
    bsr_graph_tags(&PT, pid, n);
    //bsr_graph_show(&PT,"PT");
	
	  //apply transpose to generate node permutation graph
    bsr_graph_transpose(&PT,&PN);
    //bsr_graph_show(&PN,"PN");
	
	  //extract graph of sparse matrix C and apply node permutation
    // --> the end results is a connection permutation graph
	  bsr_to_graph(C,&A);
    bsr_graph_reorder(&A, &(PN.weight), &PC);
    
    
	  //bsr_graph_show(&PC,"PC");
	  //    vec_ishow(PC.weight,A.nnz,"PC");

    // apply connection permutation to sparse matrix C
    // --> to get the reordered sparse matrix D
    bsr_reorder(C, &PC, &D);
    
    
    for(int i=0;i<(D.nrows+1);i++)
      hbsr_row_ptr.push_back(D.rowptr[i]);
    for(int i=0;i<D.nnz;i++)
      hbsr_col_ind.push_back(D.colidx[i]);
    for(int i=0;i<D.nnz*bsize*bsize;i++)
      hbsr_val.push_back(D.data[i]);
      
      
    for(int i=0;i<(D.nrows+1);i++)
      hbsr_prow_ptr.push_back(D.rowptr[i]);
    for(int i=0;i<D.nnz;i++)
      hbsr_pcol_ind.push_back(D.colidx[i]);
    for(int i=0;i<D.nnz*bsize*bsize;i++)
      hbsr_pval.push_back(D.data[i]);  
      
    
    //std::vector<double> hy_og(mb*bsize, 1.0f);
    std::vector<double> hy(mb*bsize, 0.0f); 

    
    
    //apply row permutation to vector
    for (int i=0;i<mb;i++){
      for (int j=0;j<bsize;j++)
      {
        hy[i*bsize+j]=hy_og[PN.weight[i]*bsize+j];
      }  
    }
    
    bsr->hbsr_bval=hy;
    
    bsr->hbsr_val=hbsr_val; bsr->hbsr_row_ptr=hbsr_row_ptr; bsr->hbsr_col_ind=hbsr_col_ind;
    bsr->hbsr_pval=hbsr_pval; bsr->hbsr_prow_ptr=hbsr_prow_ptr; bsr->hbsr_pcol_ind=hbsr_pcol_ind;
    //std::cout<<"hval bsr size "<<hbsr_val.size()<<std::endl;
    
    bsr->mb = bsr->hbsr_row_ptr.size()-1;
    bsr->nnzb = bsr->hbsr_col_ind.size();
    bsr->nnz =  bsr->hbsr_val.size();
    bsr->block_size = bsr->hbsr_bval.size()/bsr->mb;
    bsr->Nb = bsr->mb;
    bsr->nnzbs_prec= bsr->nnzb;
    bsr->N = bsr->hbsr_bval.size();
    
    //std::cout << "created bsr with #rows "<< bsr->mb <<" #cols "<<bsr->nnzb << " block size " <<  bsr->block_size<< std::endl;
    
    
      
    init_gpu(bsr);
    data_xfers(bsr);
    rocsparse_boilerplate(bsr);
    preconditioner(bsr);
    auto start = std::chrono::system_clock::now();
    gpu_pbicgstab(bsr,tolerance, maxit);
    auto end = std::chrono::system_clock::now();
    double elapsed_ms= std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "solver time "<<elapsed_ms << '\n';
   
    //print result
    
    std::vector<double> hx(bsr->N, 0.0f);
    hipMemcpy(hx.data(), bsr->d_x, sizeof(double) * bsr->N, hipMemcpyDeviceToHost);
    
    
    //std::cout<<"writing to "<<s2<<s<<std::endl;
    //indir << "xvals.fp64_"<<rows <<"_"<<nnzb<<"_"<<std::scientific << tolerance<<".txt";
    //std::cout<<"writing to "<<indir.str()<<std::endl;
 
    //std::ofstream output_file(indir.str());

    //std::ostream_iterator<double> output_iterator(output_file, "\n");
    //std::copy(std::begin(hx), std::end(hx), output_iterator);
    
    //apply inverse permutation answers
    int *imap = (int*)malloc(bsr->mb*sizeof(int));
    for(int i=0;i<bsr->mb;i++) imap[PN.weight[i]]=i;
    //std::vector<double> hy(mb*bsize, 0.0f);
    
    
    //std::cout << "applying inverse permutation to answer vector " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
      int ipos = imap[i];
      for (int j=0;j<bsize;j++)
      {
        ans[i*bsize+j]=hx[ipos*bsize+j];
      }  
    //std::cout<<hx[i]<<" "; 
    }
    
    
    /*
    std::cout << "hx " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
  	std::cout<<hx[i]<<" "; 
    }
    */
    std::cout<<"\n";
    
}

void run_solver_dd(bsr_matrix *C, solver_data *bsr, int nx, int ny, int nz, int bsize, int nbx, int nby, int nbz, float maxit, double tolerance, std::vector<double> &hy_og, std::vector<double> &ans)
{

    std::vector<int> hbsr_row_ptr;
    std::vector<int> hbsr_col_ind;
    std::vector<double> hbsr_val;
    
    std::vector<int> hbsr_prow_ptr;
    std::vector<int> hbsr_pcol_ind;
    std::vector<double> hbsr_pval;

    int n=nx*ny*nz;
    int nnz= nx*ny*nz +2*(nx-1)*ny*nz + 2*nx*(ny-1)*nz + 2*nx*ny*(nz-1);
    int b=bsize;
    int bb=b*b;
    int total_nnz = bb*nnz;  
    int* colidx;
    int* rowptr;
    double* data;
    int mb =n;
    
    //colidx=(int*)malloc(nnz*sizeof(int));
    //rowptr=(int*)malloc((n+1)*sizeof(int));
    //data=(double*)malloc(bb*nnz*sizeof(double));
    
    
    //bsr_matrix C, D;
    bsr_cartesian(C, nx, ny, nz, bsize);

    //flux balance
    bsr_flux_balance(C);

    //boundary condititions
    bsr_dirichlet(C,0);
    bsr_dirichlet(C,n-1);
    //bsr_show(&C,"C");
    bsr_matrix  D;
    //int nbx=16;
    //int nby=16;
    //int nbz=4;
    
   
    
    
    int *pid = (int*)malloc(n*sizeof(int));
    
    vec_ifill(pid,n,-1);
    auto start1 = std::chrono::system_clock::now();
    bsr_graph_part(nx, ny, nz, nbx, nby, nbz, pid);
    //vec_ishow(pid,n,"pid");
    
	  //reinterpret partition ids as a graph
    bsr_graph PT, PN, PC, A;
    bsr_graph_tags(&PT, pid, n);
    //bsr_graph_show(&PT,"PT");
	
	  //apply transpose to generate node permutation graph
    bsr_graph_transpose(&PT,&PN);
    //bsr_graph_show(&PN,"PN");
	
	  //extract graph of sparse matrix C and apply node permutation
    // --> the end results is a connection permutation graph
	  bsr_to_graph(C,&A);
    bsr_graph_reorder(&A, &(PN.weight), &PC);
    

    
	  //bsr_graph_show(&PC,"PC");
	  //    vec_ishow(PC.weight,A.nnz,"PC");

    // apply connection permutation to sparse matrix C
    // --> to get the reordered sparse matrix D
    bsr_reorder(C, &PC, &D);
    bsr_graph B, DD;
    bsr_to_graph(&D,&B);
    bsr_graph_decompose(&B, &PN, &DD);
    //csr_graph_show(&B,"B");
    //bsr_graph_show(&DD,"DD");

    // apply permutation
    // label conn to drop
    // apply projection
    bsr_matrix F;
    bsr_reorder(&D, &DD, &F);
    //bsr_show(&F,"F");
    
    
    auto end1 = std::chrono::system_clock::now();
    double elapsed_ms1= std::chrono::duration_cast<std::chrono::duration<double>>(end1 - start1).count();
    std::cout << "partition + decomposition time  "<<elapsed_ms1 << '\n';
    
    for(int i=0;i<(D.nrows+1);i++)
      hbsr_row_ptr.push_back(D.rowptr[i]);
    for(int i=0;i<D.nnz;i++)
      hbsr_col_ind.push_back(D.colidx[i]);
    for(int i=0;i<D.nnz*bsize*bsize;i++)
      hbsr_val.push_back(D.data[i]);
      
      
    for(int i=0;i<(F.nrows+1);i++)
      hbsr_prow_ptr.push_back(F.rowptr[i]);
    for(int i=0;i<F.nnz;i++)
      hbsr_pcol_ind.push_back(F.colidx[i]);
    for(int i=0;i<F.nnz*bsize*bsize;i++)
      hbsr_pval.push_back(F.data[i]);  
      
    
   // std::vector<double> hy_og(mb*bsize, 1.0f);
    std::vector<double> hy(mb*bsize, 0.0f); 

    
    
    //apply row permutation to vector
    for (int i=0;i<mb;i++){
      for (int j=0;j<bsize;j++)
      {
        hy[i*bsize+j]=hy_og[PN.weight[i]*bsize+j];
      }  
    }
    
    bsr->hbsr_bval=hy;
    
    bsr->hbsr_val=hbsr_val; bsr->hbsr_row_ptr=hbsr_row_ptr; bsr->hbsr_col_ind=hbsr_col_ind;
    bsr->hbsr_pval=hbsr_pval; bsr->hbsr_prow_ptr=hbsr_prow_ptr; bsr->hbsr_pcol_ind=hbsr_pcol_ind;
    std::cout<<"hval bsr size "<<hbsr_val.size()<<std::endl;
    
    bsr->mb = bsr->hbsr_row_ptr.size()-1;
    bsr->nnzb = bsr->hbsr_col_ind.size();
    bsr->nnz =  bsr->hbsr_val.size();
    bsr->block_size = bsr->hbsr_bval.size()/bsr->mb;
    bsr->Nb = bsr->mb;
    bsr->nnzbs_prec= F.nnz;
    bsr->N = bsr->hbsr_bval.size();
    
    std::cout << "created bsr with #rows "<< bsr->mb <<" #cols "<<bsr->nnzb << " block size " <<  bsr->block_size<< std::endl;
    
    
      
    init_gpu(bsr);
    data_xfers(bsr);
    rocsparse_boilerplate(bsr);
    preconditioner(bsr);
    auto start = std::chrono::system_clock::now();
    gpu_pbicgstab(bsr,tolerance, maxit);
    auto end = std::chrono::system_clock::now();
    double elapsed_ms= std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "solver time "<<elapsed_ms << '\n';
   
    //print result
    
    std::vector<double> hx(bsr->N, 0.0f);
    hipMemcpy(hx.data(), bsr->d_x, sizeof(double) * bsr->N, hipMemcpyDeviceToHost);
    
    
    //std::cout<<"writing to "<<s2<<s<<std::endl;
    //indir << "xvals.fp64_"<<rows <<"_"<<nnzb<<"_"<<std::scientific << tolerance<<".txt";
    //std::cout<<"writing to "<<indir.str()<<std::endl;
 
    //std::ofstream output_file(indir.str());

    //std::ostream_iterator<double> output_iterator(output_file, "\n");
    //std::copy(std::begin(hx), std::end(hx), output_iterator);
    
    //apply inverse permutation answers
    int *imap = (int*)malloc(bsr->mb*sizeof(int));
    for(int i=0;i<bsr->mb;i++) imap[PN.weight[i]]=i;
    //std::vector<double> hy(mb*bsize, 0.0f);
    
    std::cout << "actual edges "<< C->nnz<<" post partitionning  edges " << D.nnz << " post DD edges "<<DD.post_dd_edges<<std::endl;
    //std::cout << "applying inverse permutation to answer vector " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
      int ipos = imap[i];
      for (int j=0;j<bsize;j++)
      {
        ans[i*bsize+j]=hx[ipos*bsize+j];
      }  
    //std::cout<<hx[i]<<" "; 
    }
    
    
    /*
    std::cout << "hx " <<std::endl;
    for (int i =0; i< bsr->mb; i++)
    {
  	std::cout<<hx[i]<<" "; 
    }
    */
    std::cout<<"\n";


}


void compare_vectors(std::vector<double> &ans1,std::vector<double> &ans2)
{

  double thresh =1e-4;
  int j=0;
 
  for (int k=0;k<ans1.size();k++)
    if(abs(ans1[k]-ans2[k])>thresh)
    //if(abs(ans1[k]-ans2[k])>thresh*(abs(ans1[k]+ans2[k])))
    {
      if(j<10)
        std::cout<<std::setprecision(20) <<" answers don't match at element " << k << " rocsparse dd " <<ans1[k] << " ildu dd "<<ans2[k] << " relative diff "<< abs(ans1[k] - ans2[k])/abs(ans1[k] + ans2[k])  <<" "  <<std::endl;
      j++;
    }
  std::cout<<"# vector elements with abs(answer-OPM_answer) > 1e-4> = " << j <<std::endl;
}




int main(int argc, char** argv)
{
	 
   	int rows=4;
    int bsize = 3;
    //int nx=32;
    //int ny=8;
    //int nz=16;
    
    
    //Try to include recursive inertial bisection 
    //Default:  
    int nx=128;
    int ny=128;
    int nz=128;
  
    int nbx=8;;
    int nby=16;
    int nbz=16;
    
    int runvariant=8; 
    
    //runvariant
    //1: Undomidied solver with rocsparse
    //2: Rocsparse Solver with paritions but without dropping inter-parition edges(non zeros outside of the parition are not dropped)
    //3: Rosparse Solver with domain decomposition ( dropping the nonzeros outside of paritions)
    //4: Solver with our ILDU (edge-centric - non-deterministic atomic adds)
    //5: Solver with our ILDU (vertex-centric - non-deterministic atomic adds)
    
    if (argc>6)
    {
      nx = atoi(argv[1]);
      ny = atoi(argv[2]);
      nz = atoi(argv[3]);
      nbx = atoi(argv[4]);
      nby = atoi(argv[5]);
      nbz = atoi(argv[6]);
    } 
    
    if (argc>7)
      runvariant= atoi(argv[7]);
    
    
    
    if (runvariant>8)
      runvariant=8; //only 5 variants are added here. 6 --> run all. >6 --> run all
    int runcorrectness = 0;
    if(argc>8)
      runcorrectness = 1;
    
    if(nbx*nby*nbz!=2048)
    {
      std::cout<<"Subdomain size !=2048 not supported ... exiting\n";
      exit(0);
    }
    
    int mb=nx*ny*nz;
    
    if(mb%(nbx*nby*nbz)!=0)
    {
      std::cout<<"Subdomain dimentions do not divide the grid not ... exiting\n";
      exit(0);
    }
    
    int num_subdomains = mb/(nbx*nby*nbz);
    std::vector<double> ans1(mb*bsize, 0.0f);
    std::vector<double> ans2(mb*bsize, 0.0f);
    std::vector<double> ans3(mb*bsize, 0.0f);
    std::vector<double> ans4(mb*bsize, 0.0f);
    std::vector<double> ans5(mb*bsize, 0.0f);
    std::vector<double> ans6(mb*bsize, 0.0f);
    std::vector<double> ans7(mb*bsize, 0.0f);
   	//int nnzb = 0;
    //int num_jacobi_blocks;
    //int rows_per_jacobi_block;
  	/*
    if (argc >= 2)
    {
  		num_jacobi_blocks = atoi(argv[1]);    
    }
    rows_per_jacobi_block = num_jacobi_blocks;
    
    if (argc >= 3)
    {
  		rows_per_jacobi_block = atoi(argv[2]);    
    }
  
  	if (argc >= 4)
    {
  		nnzb = atoi(argv[3]);    
    }
  	
    if (argc >= 5)
    {
  		block_dim = atoi(argv[4]);    
    }
  
    rows = num_jacobi_blocks * rows_per_jacobi_block;
   */
 
    //float maxit = 2000.0;
    float maxit = 2000.0;
    //double tolerance = 1e-10;
    double tolerance = 1e-8;
    //solver_data bsr;
    std::ostringstream indir, ss; 
    //indir << "/home/agondhal/linalg_standalone/rocky/data/";
    indir << "/home/agondhal/linalg_standalone/rocky/rocsparse_tests/subdomains/";
    //init_host(&bsr, indir, rows, num_jacobi_blocks, rows_per_jacobi_block,nnzb, block_dim);
    
    bsr_matrix C, C_1, C_2, C_3, C_4, C_5, C_6;
    solver_data bsr, bsr_1, bsr_2;
    solver_data2 bsr_3, bsr_6;
    solver_data_vc bsr_4,bsr_5;
    
     //std::vector<double> hy_og(mb*bsize, 1.0f);
    std::vector<double> hy(mb*bsize, 1e-6); 

    //for(int i=0;i<mb*bsize/2;i++) 
    //  hy[i]=1.0;
    for(int i=0;i<bsize;i++) 
    {
     hy[bsize*mb-1-i]=-1;
     hy[i]=1;
    }
    
    if (runvariant==1 || runvariant==8)
    {
    std::cout<<"----------------------------------------------------------------------------------\n";
    std::cout<<"Running OPM unmodified\n";
    run_solver_unmodified(&C, &bsr, nx, ny, nz, bsize, maxit, tolerance, hy,ans1);
    std::cout<<"----------------------------------------------------------------------------------\n\n";
    }
    if (runvariant==2 || runvariant==8)
    {
    std::cout<<"----------------------------------------------------------------------------------\n";
    std::cout<<"Running OPM with partitioning\n";
    run_solver_reordered(&C_1, &bsr_1, nx, ny, nz, bsize, nbx, nby, nbz, maxit, tolerance, hy, ans2);
    std::cout<<"----------------------------------------------------------------------------------\n\n";
    }
    if (runvariant==3 || runvariant==8)
    {
    std::cout<<"----------------------------------------------------------------------------------\n";
    std::cout<<"Running OPM with partitioning + domain decomposition \n";
    run_solver_dd(&C_2, &bsr_2, nx, ny, nz, bsize, nbx, nby, nbz, maxit, tolerance, hy, ans3);
    std::cout<<"----------------------------------------------------------------------------------\n\n";
    }
    if (runvariant==4 || runvariant==8)
    {
    std::cout<<"----------------------------------------------------------------------------------\n";
    std::cout<<"Running solver with edge-centric traversal for triangular solves\n";
    run_solver_graph_ildu(&C_3, &bsr_3, nx, ny, nz, bsize, nbx, nby, nbz, maxit, tolerance, hy, ans4, num_subdomains);
    std::cout<<"----------------------------------------------------------------------------------\n\n";
    }
    if (runvariant==5 || runvariant==8)
    {
    std::cout<<"----------------------------------------------------------------------------------\n";
    std::cout<<"Running solver with edge-centric traversal for triangular solves with fp32 \n";
    int fp32_nbx = 2*nbx;
    int fp32_num_subdomains = num_subdomains/2;
    run_solver_graph_ildu_fp32(&C_6, &bsr_6, nx, ny, nz, bsize, fp32_nbx, nby, nbz, maxit, tolerance, hy, ans7, fp32_num_subdomains);
    std::cout<<"----------------------------------------------------------------------------------\n\n";
    }
     
    if (runvariant==6 || runvariant==8)
    {
    std::cout<<"----------------------------------------------------------------------------------\n";
    std::cout<<"Running solver with vertex-centric traversal for triangular solves\n";
    run_solver_graph_ildu_vc(&C_4, &bsr_4, nx, ny, nz, bsize, nbx, nby, nbz, maxit, tolerance, hy, ans5, num_subdomains);
    std::cout<<"----------------------------------------------------------------------------------\n\n";
    }
    
    if (runvariant==7 || runvariant==8)
    {
    std::cout<<"----------------------------------------------------------------------------------\n";
    std::cout<<"Running solver with vertex-centric traversal for triangular solves with fp32\n";
    int fp32_nbx = 2*nbx;
    int fp32_num_subdomains = num_subdomains/2;
    run_solver_graph_ildu_vc_fp32(&C_5, &bsr_5, nx, ny, nz, bsize, fp32_nbx, nby, nbz, maxit, tolerance, hy, ans6, fp32_num_subdomains);
    std::cout<<"----------------------------------------------------------------------------------\n\n";
    }
    
    //run_solver_graph_ildu_vc_nolds(&C_5, &bsr_5, nx, ny, nz, bsize, nbx, nby, nbz, maxit, tolerance, hy, ans6, num_subdomains);
    //compare_vectors(ans1,ans6);
    if(runcorrectness==1 && runvariant==8)
    {
    compare_vectors(ans1,ans2);
    compare_vectors(ans1,ans3);
    compare_vectors(ans1,ans4);
    compare_vectors(ans1,ans5);
    compare_vectors(ans1,ans6);
    compare_vectors(ans1,ans7);
    }
    

}  

