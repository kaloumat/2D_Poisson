#ifndef SPARSE_H_INCLUDED
#define SPARSE_H_INCLUDED

#include <vector>
#include "triplet.h"

class Sparse{
public:
    std::vector<int> PI, J;     // indexy radku a sloupcu matice
    std::vector<double> VAL;    // hodnoty nenulovych prvku matice

    int     n;                  // rozmer matice
    int     nz;                 // pocet nenulovych prvku

    void Create(Triplet *ptr);
    void GaussSeidel(std::vector<double>& b, std::vector<double>& x, int niter, double err);
    double GetRowRezidiuum(std::vector<double>& b, std::vector<double>& x, int i, double *aii);
};
#endif // SPARSE_H_INCLUDED
