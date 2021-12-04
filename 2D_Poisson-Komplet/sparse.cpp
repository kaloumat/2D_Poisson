#define DEBUG
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include "sparse.h"
/* ----------------------------------------------------------------------------------- */
void Sparse::Create(Triplet *ptr){
    ptr->Quicksort(0, ptr->nz);
    ptr->Unique();

    n  = ptr->n;
    nz = ptr->nz;

    J.resize(nz);
    VAL.resize(nz);

    int i, j = 0;

    PI.push_back(0);

    for(i = 0; i < n; i++){
        while(ptr->I[j] == i + 1 && j < ptr->nz){
            J[j]   = ptr->J[j] - 1;
            VAL[j] = ptr->VAL[j];
            j++;
        }
        PI.push_back(j);
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse::GaussSeidel(std::vector<double>& b, std::vector<double>& x, int niter, double err){
    int i, iter;
    double aii, rez, rezi;

    for(iter = 0; iter < niter; iter++){
        rez = 0;
        for(i = 0; i < n; i++){
            rezi = GetRowRezidiuum(b, x, i, &aii);
            x[i] += rezi / aii;
            rez = std::max(rez, fabs(rezi));
        }

#ifdef DEBUG
        printf("rez[%d] = %e\n", iter, rez);
#endif // DEBUG

        if(rez < err){
            break;
        }
    }
}
/* ----------------------------------------------------------------------------------- */
double Sparse::GetRowRezidiuum(std::vector<double>& b, std::vector<double>& x, int i, double *aii){
    int k, j;
    double rezi, aij;

    rezi = b[i];
    *aii = 0;

    for(k = PI[i]; k < PI[i + 1]; k++){
        aij   = VAL[k];
        j     = J[k];
        rezi -= aij * x[j];

        if(i == j)
            *aii = aij;
    }
    return rezi;
}
