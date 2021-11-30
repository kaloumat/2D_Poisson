#ifndef MESH_H_INCLUDED
#define MESH_H_INCLUDED
typedef struct{
    int nbNods;              // pocet uzlu site
    int nbBndrEdges;         // pocet hranicnich usecek site
    int nbTriangles;         // pocet trojuhelniku site

    double *x, *y;           // x, y souradnice
    int *BndrA, *BndrB;      // indexy uzlu danych hranicnich usecek
    int *BndrMark;           // speficka cisla pro okrajove podminky
    int *TriA, *TriB, *TriC; // indexy uzlu danych trojuhelniku
    int *TriMark;            // specificka cisla oznacujici dane oblasti (v nasem pripade
                             //     pouze jedna oblast, tedy jedno cislo pro vsechy trojuhelniky)

    int *isDirichlet;        // pole jednicek a nul, odpovida velikosti matice, neboli poctu uzlu 
                             //     kdyz jednicka -> v uzlu je predepsana Dirichletova okr pod
                             //     kdyz nula -> v uzlu neni predepsana zadna okr pod
} mesh; //pojmenovani struktury (muze byt jakekoliv)

// pomocna struktura pro nacitani dat ze souboru .msh
typedef struct{
    int idx;            // index elementu
    int etyp;           // typ elementu
    int markPhysical;   // specificke cislo
    int ilist[3];       // pole pro tri hodnoty -> odpovida trem vrcholum pro trojuhelnik
                        //      v pripade ctyruhelniku bysme museli zvetsit
} gmshline;

// definovani velicin, ktere se v celem skriptu nemeni
#define GMSH_SEGMENT 1  // znamena GMSH_SEGMENT = 1 a odpovida usecce
#define GMSH_TRIANGLE 2 // znamena GMSH_TRIANGLE = 2 a odpovida trojuhelniku
#define BUFFER_SIZE 500 // znamena max BUFFER_SIZE = 500, predpokladame, ze radek
                        //      v .txt souboru bude mit max 500 znaku
/*
    Definice funkci ktere pouzivame v mesh.cpp abysme je tam nemuseli mit napsany poporade
    a nestalo se, ze pri kompilaci nebyla jeste nejaka funkce jakoby definovana, je to taky
    fajn pro prehled, kdyz uz je tech funkci hodne

*/
void Mesh_Allocate(mesh *p);
void Mesh_Free(mesh *p);
void Mesh_Load1(mesh *p, const char *fname);
int GmshLine_Read(gmshline *p, const char *buf);
void Mesh_Read(mesh *p, const char *fname);
void Mesh_Load2(mesh *p, const char *fname);
#endif // MESH_H_INCLUDED
