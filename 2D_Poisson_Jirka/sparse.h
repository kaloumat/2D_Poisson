#ifndef SPARSE_H_INCLUDED
#define SPARSE_H_INCLUDED
typedef struct {
    int     *PI;    // pole pro ukazatel na radky I
        /*
            Pole PI pouzivame nasledovne. Misto toho abysme meli pole I, kde jsou serazeny
            indexy vsech radku, vytvorime pomocne pole PI. Do pole PI budeme ukladat
            pouze polohy v poli I, na kterych vzdy konci a zacina novy index.
            Priklad: Sparse matice 3x3
                I   J   VAL
                0   0   0.912354
                0   2   1.232145
                1   0   -0.15676
                2   0   0.51348
                2   1   0.4156
                2   2   0.81169

                Pole PI potom vypada nasledovne
                PI[0] = 0
                PI[1] = 2
                PI[2] = 3
                PI[3] = 5
                PS: Pole PI je tak vzdy o jedno vetsi nez velikost matice.
        */
    int     *J;     // pole pro indexy sloupcu v matici tuhosti
    double  *VAL;   // pole pro nenulove hodnoty matice tuhosti
    int     n;      // rozmer matice (v prikladu je n = 3)
    int     nz;     // pocet nenulovych prvku (pocet radku ve sparse, v prikladu nz = 6)
    int     nalloc; // cislo pro alokoci pameti
} sparse;

/*
    Definice funkci ktere pouzivame v sparse.cpp abysme je tam nemuseli mit napsany poporade
    a nestalo se, ze pri kompilaci nebyla jeste nejaka funkce jakoby definovana, je to taky
    fajn pro prehled, kdyz uz je tech funkci hodne

*/
void Sparse_Allocate(sparse *p, int n, int nalloc);
void Sparse_Free(sparse *p);
void Sparse_GaussSeidel(sparse *A, double *b, double *x, int niter, double err);
double max(double a, double b);
double Sparse_GetRowRezidiuum(sparse *A, double *b, double *x, int i, double *aii);
#endif // SPARSE_H_INCLUDED
