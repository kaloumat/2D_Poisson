#ifndef MESH_H_INCLUDED
#define MESH_H_INCLUDED
typedef struct{
    int nbNods;              // pocet uzlu site
    int nbBndrEdges;         // pocet hranicnich usecek site
    int nbTriangles;         // pocet trojuhelniku site

    double *x, *y;           // x, y souradnice
    int *BndrA, *BndrB;      // indexy uzlu danych hranicnich usecek
    int *BndrMark;           // speficka cisla pro okrajove podminky (cislo nastavuje v gmsh)
    int *TriA, *TriB, *TriC; // indexy uzlu danych trojuhelniku
    int *TriMark;            // specificka cisla oznacujici dane oblasti (nastavuje v gmsh, v nasem pripade
                             //     pouze jedna oblast, tedy jedno cislo pro vsechy trojuhelniky)

     int *isDirichlet;        // pole jednicek a nul, odpovida velikosti matice, neboli poctu uzlu
                             //     kdyz jednicka -> v uzlu je predepsana Dirichletova okr pod
                             //     kdyz nula -> v uzlu neni predepsana zadna okr pod
} mesh;

/*
    Definice funkci ktere pouzivame v mesh.cpp abysme je tam nemuseli mit napsany poporade
    a nestalo se, ze pri kompilaci nebyla jeste nejaka funkce jakoby definovana, je to taky
    fajn pro prehled, kdyz uz je tech funkci hodne

*/
void Mesh_Allocate(mesh *p);
void Mesh_Free(mesh *p);
void Mesh_SaveVTK(mesh *p, double *u, const char *fname);
#endif // MESH_H_INCLUDED
