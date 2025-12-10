# Domain-Decomposed BiCGSTAB Solver

## Overview
This program performs **domain decomposition using geometric cuts** on 3D Laplacian matrices and uses the docomposed matrices to perform ILU0 preconditioning for a **preconditioned BiCGSTAB solver**.  
The matrix is represented using a block sparse row format where block size =3.
It implements multiple solver configurations, including baseline rocSPARSE variants, reordered/partitioned solvers, and custom ILDU-based solvers with both edge-centric and vertex-centric triangular solves.  

The implementation leverages **fine-grained domain decomposition** to increase parallelism on GPUs by mapping subdomains to thread blocks and reducing irregular memory accesses.  

---

## Compilation
Simply run:
```bash
make
```

---

## Running
The program takes the following arguments:  

```bash
./solver nx ny nz nbx nby nbz [runvariant] [check_correctness]
Currently supported subdomain size = 2048 rows for a block sparse matrix with block size =3
```

- `nx ny nz` → Grid dimensions of the 3D Laplacian (default: 128×128×128).  
- `nbx nby nbz` → Subdomain decomposition along x, y, z directions (must multiply to 2048).  
- `runvariant` → Which solver configuration to run (default: 8 = run all).  
- `check_correctness` → Optional flag (set `1` to compare solutions across variants).  

Example:  
```bash
./solver 128 128 128 8 16 16 8
```

---

## Implemented Configurations
The following solver variants are supported:  

1. **Unmodified rocSPARSE solver** (baseline).  
2. **Reordered solver with partitioning only** (no dropping of inter-partition edges).  
3. **Domain-decomposed solver** (dropping inter-partition edges).  
4. **ILDU solver with edge-centric triangular solves (FP64)**.  
5. **ILDU solver with edge-centric triangular solves (FP32)**.  
6. **ILDU solver with vertex-centric triangular solves (FP64)**.  
7. **ILDU solver with vertex-centric triangular solves (FP32)**.  
8. **Run all of the above** in sequence.  

