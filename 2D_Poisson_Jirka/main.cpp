#include <iostream>
#include <cstdlib>
#include "mesh.h"
#include "vector.h"
#include "sparse.h"

using namespace std;

void Nacti_fakeMesh(mesh *p, const char *fname);
void Nacti_fakeSparseMat_a_fakeVekb(mesh *p, sparse *A, double* b, const char *fname);

int main()
{
    int i, j; // inicializace indexu, ktere budeme pouzivat pri ruznych cyklech

    /*
        inicializace structury typu mesh, ktera je definovana v hlavickovem
        souboru mesh.h, a jeji pojmenovani "sit" (muze mit jakykoliv nazev)
    */
    mesh sit;

    Nacti_fakeMesh(&sit, "fakeMesh.txt"); // nacte data, ktera by jsme jinak dostali nacteni souboru .msh

    /*
        inicializace pole se jmenem b, puvodne byla syntaxe double* b, tu jsme ale
        nahradili v hlavickovem souboru vector.h za "vector", tedy misto kompilator
        vyhodnocuje "double*" a "vector" jako stejnou vec
    */
    vector b;
    b = Vector_Allocate(sit.nbNods); // funkce Vector_Allocate alokuje pamet o velikosti sit.nbNods

    /*
        inicializace structury typu sparse, ktera je definovana v hlavickovem
        souboru sparse.h, a jeji pojmenovani "A" (muze mit jakykoliv nazev)
    */
    sparse A;
    Nacti_fakeSparseMat_a_fakeVekb(&sit, &A, b, "fakeSparseMatice_fakeVekb.txt"); // nacte data ve formatu sparse, ktere by jinak vznikly vypoctem lokalnich prispevku do matice A a vektoru b

    vector u; // stejne jak u vektoru b
    u = Vector_Allocate(sit.nbNods); // funkce Vector_Allocate alokuje pamet o velikosti sit.nbNods

    /*
        Tisk do terminalu, slouzi pro kontrolu jestli se data z textovych souboru nacetli spravne
    */
    /*-------------------------------------------------------------------------*/
    cout << "souradnice x,y:" << endl;
    for (i = 0; i < sit.nbNods; i++)
        cout << "%lf %lf\n" << sit.x[i] << sit.y[i] << endl;

    cout << "hranicni usecky:" << endl;
    for (i = 0; i < sit.nbBndrEdges; i++)
        cout << "%d %d %d\n" << sit.BndrA[i] << sit.BndrB[i] << sit.BndrMark[i] << endl;

    cout << "trojuhelniky:" << endl;
    for (i = 0; i < sit.nbTriangles; i++)
        cout << "%d %d %d %d\n" << sit.TriA[i] << sit.TriB[i] << sit.TriC[i] << sit.TriMark[i] << endl;

    cout << "sparse matice:" << endl;
    for(i = 0; i < A.n; i++){
        for(j = A.PI[i]; j < A.PI[i + 1]; j++){
            cout << "%d %d %lf\n" << i << A.J[j] << A.VAL[j] << endl;
        }
    }

    cout << "\nvektor b:" << endl;
    for(i = 0; i < sit.nbNods; i++)
        cout << "%lf\n" << b[i] << endl;
    /*-------------------------------------------------------------------------*/

    cout << "\nreziduum GaussSeidel:" << endl;
     /*
        Funkce pro vypocet soustavy lin. rovnic pomoci GaussSeidel
        Jako parametry ji predame sparse matici "A", vektor "b" a nulovy vektor(pole) "u" ktery funkce vyplni
        4. parametry odpovida maximalnimu poctu iteraci
        5. parametry odpovida minimalni chybe metody
    */
    Sparse_GaussSeidel(&A, b, u, 100, 1e-10);

    /*
        Vytisk vysledneho vektoru "u" do terminalu
    */
    cout << "\nvysledny vektor u:" << endl;
    for(int i = 0; i < A.n; i++){
        cout << "%lf\n" << u[i] << endl;
    }


    /*
        Ulozeni vektoru "u" do textove souboru se jmenem "poisson.txt" a do souboru .vtk se jmenem "poisson.vtk"
    */
    Vector_Save(u, A.n, "poisson.txt");
    Mesh_SaveVTK(&sit, u, "poisson.vtk");


    /*
        Kdyz uz s uvedenymi poli nepotrebujeme pracovat dealokujeme pamet (v kompletnim
        skriptu bude prikaz uveden az na konci), viz napr. funkce Mesh_Free v mesh.cpp
    */
    Mesh_Free(&sit);
    Vector_Free(u);
    Vector_Free(b);
    Sparse_Free(&A);

    return 0;
}

//Tuhle funkci neres, pouze pomocna pro nacteni dat, ktere predstavuji data z .msh formatu
void Nacti_fakeMesh(mesh *p, const char *fname){
    p->nbNods = 12;
    p->nbBndrEdges = 8;
    p->nbTriangles = 14;
    Mesh_Allocate(p);

    int i;
    char radek[101];
    FILE *fid;

    fid = fopen(fname, "r");

    for(i = 0; i < p->nbNods; i++){
        fgets(radek, 100, fid);
        sscanf(radek, "%lf %lf", &p->x[i], &p->y[i]);
    }

    int pom = i;
    for(i = pom; i < p->nbBndrEdges + pom; i++){
        fgets(radek, 100, fid);
        sscanf(radek, "%d %d %d", &p->BndrA[i - pom], &p->BndrB[i - pom], &p->BndrMark[i - pom]);
    }

    for(i = pom; i < p->nbTriangles + pom; i++){
        fgets(radek, 100, fid);
        sscanf(radek, "%d %d %d %d", &p->TriA[i - pom], &p->TriB[i - pom], &p->TriC[i - pom], &p->TriMark[i - pom]);
    }

    fclose(fid);
}

//Tuhle funkci neres, pouze pomocna pro nacteni dat, ktere predstavuji sparse matici a vektor b
void Nacti_fakeSparseMat_a_fakeVekb(mesh *p, sparse *A, double *b, const char *fname){
    Sparse_Allocate(A, p->nbNods, 22);

    int i;
    char radek[101];
    FILE *fid;

    fid = fopen(fname, "r");

    for(i = 0; i < p->nbNods + 1; i++){
        fgets(radek, 100, fid);
        sscanf(radek, "%d", &A->PI[i]);
    }

    int pom = i;
    for(i = pom; i < A->nalloc + pom; i++){
        fgets(radek, 100, fid);
        sscanf(radek, "%d %lf", &A->J[i - pom], &A->VAL[i - pom]);
    }

    for(i = pom; i < p->nbNods + pom; i++){
        fgets(radek, 100, fid);
        sscanf(radek, "%lf", &b[i - pom]);
    }

    fclose(fid);
}


