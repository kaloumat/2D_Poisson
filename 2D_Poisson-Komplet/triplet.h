#ifndef TRIPLET_H_INCLUDED
#define TRIPLET_H_INCLUDED

#include <vector>

class Triplet{
public:
    int n;                                              // rozmer matice
    int nz;                                             // pocet nenulovych prvku

    struct formatTriplet{
        int   I, J;                                     // indexy radku a sloupcu matice
        double VAL;                                     // hodnoty nenulovych prvku matice

        bool operator<(const formatTriplet& a) const{
            if(I < a.I)
                return I < a.I;
            if(I == a.I && J < a.J)
                return J < a.J;
            return false;
        }
    };

    std::vector<formatTriplet> triplet;

    Triplet(int n, int nz){
        this->n  = n;
        this->nz = nz;
    };

    void Add(int _i, int _j, double _val);
    void AddRHS(std::vector<double>& b, int i, int j, double qval);
    void Unique();
};
#endif // TRIPLET_H_INCLUDED
