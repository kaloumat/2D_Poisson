#ifndef SPARSE_H_INCLUDED
#define SPARSE_H_INCLUDED

#include <vector>
#include "triplet.h"

class Sparse{
public:
    int n;                                      // rozmer matice
    int nz;                                     // pocet nenulovych prvku

    std::vector<int> PI;                        // pole pro indexy radku
    std::vector<std::pair<int, double>> J_VAL;  // indexy sloupcu matice a hodnoty nenulovych prvku matice


    void Create(Triplet *ptr);
    void GaussSeidel(std::vector<double>& b, std::vector<double>& x, int niter, double err);
    double GetRowReziduum(std::vector<double>& b, std::vector<double>& x, int i, double *aii);
};
#endif // SPARSE_H_INCLUDED
