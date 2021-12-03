#ifndef TRIPLET_H_INCLUDED
#define TRIPLET_H_INCLUDED

#include <vector>

class Triplet{
public:
    std::vector<int> I, J;      // indexy radku a sloupcu matice
    std::vector<double> VAL;    // hodnoty nenulovych prvku matice

    int      n;                 // rozmer matice
    int      nz;                // pocet nenulovych prvku

    Triplet(int n, int nz){
        this->n = n;
        this->nz = nz;
    };

    void Add(int _i, int _j, double _val);
    void AddRHS(std::vector<double>& b, int i, int j, double qval);
    void Quicksort(int first, int last);
    void Unique();
};
#endif // TRIPLET_H_INCLUDED
