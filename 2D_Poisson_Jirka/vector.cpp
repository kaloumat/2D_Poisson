#include <cstdlib>
#include <iostream>
#include "vector.h"

using namespace std;
/* ----------------------------------------------------------------------------------- */
double *Vector_Allocate(int n){
    // funkce pro alokaci pameti
    double *p;

    p = (double *)calloc(n, sizeof(double));

    return p;
}
/* ----------------------------------------------------------------------------------- */
void Vector_Free(double *x){
    // funkce pro dealokaci pameti
    free(x);
}
/* ----------------------------------------------------------------------------------- */
void Vector_Save(double *x, int n, char const *fname){
    // funkce pro ulozeni vektoru x do textoveho souboru(lze potom napriklad otevrit v matlabu)
    int i;
    FILE *fid;  // syntaxe pro praci se souborem, soubor pojmenuje interne fid (muze byt jakekoliv jmeno)

    fid = fopen(fname, "w"); // otevreni souboru se jmenem fname (zadavame v main.cpp) v modu w-writing(zapisujeme nove hodnoty)
    if (fid == NULL){ // pro kontrolu, kdyby se nahodou z nejakeho duvodu nepodarilo otevrit soubor
        cout << "Error, nepodarilo se otevrit soubor pro VECTOR" << endl;
        exit(1);
    }

    for(i = 0; i < n; i++){ // zapis hodnot z pole x do sloupce
        fprintf(fid, "%lf\n", x[i]);
    }

    fclose(fid); // uzavreni souboru fid
}
