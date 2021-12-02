#ifndef VECTOR_H_INCLUDED
#define VECTOR_H_INCLUDED
/*
    Definice funkci ktere pouzivame v vector.cpp abysme je tam nemuseli mit napsany poporade
    a nestalo se, ze pri kompilaci nebyla jeste nejaka funkce jakoby definovana, je to taky
    fajn pro prehled, kdyz uz je tech funkci hodne

*/
typedef double * vector;
double *Vector_Allocate(int n);
void    Vector_Free(double *x);
void    Vector_Save(double *x, int n, char const *fname);
#endif // VECTOR_H_INCLUDED
