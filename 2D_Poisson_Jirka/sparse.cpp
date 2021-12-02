#include <cstdlib>
#include <iostream>
#include <cmath>
#include "sparse.h"
/* ----------------------------------------------------------------------------------- */
void Sparse_Allocate(sparse *p, int n, int nalloc){
    // funkce pro alokaci pameti
    p->PI = (int *)malloc((n + 1) * sizeof(int)); // malloc vycleni v pocitacove pameti misto o velikost (n+1) bytu
    p->J = (int *)malloc(nalloc * sizeof(int)); // malloc vycleni v pocitacove pameti misto o velikost nalloc bytu
    p->VAL = (double *)malloc(nalloc * sizeof(double));

    p->n = n;
    p->nz = 0;
    p->nalloc = nalloc;
}
/* ----------------------------------------------------------------------------------- */
void Sparse_Free(sparse *p){
    // funkce pro dealokaci pameti
    if (p->nalloc > 0){
        free(p->PI);
        free(p->J);
        free(p->VAL);

        p->nalloc = 0;
        p->nz = 0;
    }
}
/* ----------------------------------------------------------------------------------- */
void Sparse_GaussSeidel(sparse *A, double *b, double *x, int niter, double err){
    // funkce pro vypocet soustavy lin. rovnic pomoci GaussSeidel
    int i, iter; // incializace pomocnych indexu
    double aii, rez, rezi; // incializace pomocnych hodnot typu double

    for(iter = 0; iter < niter; iter++){ // cyklus probiha dokud neni iterace mensi nez max. pocet iteraci (zadavame jako parametr)
        rez = 0; // pro kazdy radek pocitame rezidium zvlast, takze ho vzdy pro dalsi radek vynulujeme
        for(i = 0; i < A->n; i++){ // vypoce
            /*
                Vypocet funkce na aktualnim radku pomoci funkce Sparse_GetRowRezidiuum
                Jako parametry predame matici A, vektor b, vektor x(u), radek i
                    a cislo na diagonale aii, ktere funkce Sparse_GetRowRezidiuum vyplni
                rez(i) = b(i) - suma{a(ij)*x(j)}
            */
            rezi = Sparse_GetRowRezidiuum(A, b, x, i, &aii);
            /*
                Vypocet i-te hodnoty vektoru x(u) v aktualni iteraci
                x(i) = x(i) + 1/a(ii) * [b(i) - suma{a(ij)*x(j)}]    <=>    x(i) = x(i) + 1/a(ii) * rez(i)
                -> tenhle zapis GS ti tady pres komentar asi nevysvetlim :D, takze nad tim moc nelamej hlavu, doladime na pristi schuzce :D
            */
            x[i] += rezi / aii;
            /*
                do promenne "rez" ulozime max rezidium,
                    to je bud rezidium z predchoziho radku, nebo rezidium z aktualniho radku
            */
            rez = max(rez, fabs(rezi));
        }
        printf("rez[%d] = %.16lf\n", iter, rez); // tisk rezidia do terminalu pro kazdou iteraci (pro kontrolu)
        if(rez < err){ // cyklus skonci, kdyz bude rezidium mensi nez minimalni chyba metody (zadavame jako parametr)
            break;
        }
    }
}
/* ----------------------------------------------------------------------------------- */
double Sparse_GetRowRezidiuum(sparse *A, double *b, double *x, int i, double *aii){
    // funkce pro vypocet rezidia na jednom radku matice
    int k, j; // inicializace pomocnych neznamych
    double rezi, aij;

    /* rez(i) = b(i) - suma{a(ij)*x(j)}
        a abych mohl pouzit iteracni vypocet, polozim rez(i) = b(i)
        potom rez(i) = rez(i) - suma{a(ij)*x(j)}
    */
    rezi = b[i];
    *aii = 0; // inicializace nezname, do ktere ukladam prvek na diagonale

    for(k = A->PI[i]; k < A->PI[i + 1]; k++){ // provadim cyklus pred idexy radku I, dokud se index I o jedna nezvedne
        aij   = A->VAL[k];
        j     = A->J[k];
        rezi -= aij * x[j]; // vypocet rezidua na i-tem radku matice

        if(i == j)
            *aii = aij; //ulozeni prvku na diagonale do nezname aii
    }
    return rezi; // vraci rezidium na i-tem radku
}
/* ----------------------------------------------------------------------------------- */
double max(double a, double b){
    // pomocna funkce pro nalezeni maxima
    double max_val;

    if(a > b)
        max_val = a;
    else
        max_val = b;
    return max_val;
}
