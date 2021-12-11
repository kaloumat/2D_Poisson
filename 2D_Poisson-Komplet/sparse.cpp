#define DEBUG
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include "sparse.h"
/* ----------------------------------------------------------------------------------- */
void Sparse::Create(Triplet *ptr){
    ptr->Unique();

    n  = ptr->n;
    nz = ptr->nz;

    int i, j = 0;

    PI.push_back(0);

    for(i = 0; i < n; i++){
        while(ptr->triplet[j].I == i + 1 && j < ptr->nz){
            J_VAL.push_back(std::make_pair(ptr->triplet[j].J - 1, ptr->triplet[j].VAL));
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
            rezi = GetRowReziduum(b, x, i, &aii);
            x[i] += rezi / aii;
            rez = std::max(rez, fabs(rezi));
        }

#ifdef DEBUG
        std::cout << "rez[" << iter << "] = " << std::scientific << rez << std::endl;
#endif // DEBUG

        if(rez < err){
            break;
        }
    }
}
/* ----------------------------------------------------------------------------------- */
double Sparse::GetRowReziduum(std::vector<double>& b, std::vector<double>& x, int i, double *aii){
    int k, j;
    double rezi, aij;

    rezi = b[i];
    *aii = 0;

    for(k = PI[i]; k < PI[i + 1]; k++){
        aij   = J_VAL[k].second;
        j     = J_VAL[k].first;
        rezi -= aij * x[j];

        if(i == j)
            *aii = aij;
    }
    return rezi;
}
