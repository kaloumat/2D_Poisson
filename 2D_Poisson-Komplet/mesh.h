#ifndef MESH_H_INCLUDED
#define MESH_H_INCLUDED

#define GMSH_SEGMENT 1
#define GMSH_TRIANGLE 2

#include <vector>

class Mesh;
class element;
class gradfi;
class quadrature;
class gmshline;

class Mesh{
public:
    int nbNods;                                             // pocet uzlu site
    int nbBndrEdges;                                        // pocet hranicnich usecek site
    int nbTriangles;                                        // pocet trojuhelniku site

    std::vector<double> x, y;                               // x, y souradnice
    std::vector<int>    BndrA, BndrB, BndrMark;             // indexy uzlu hranicnich usecek a jejich "physical" marker z gmsh
    std::vector<int>    TriA, TriB, TriC, TriMark;          // indexy uzlu trojuhelniku a jejich "physical" marker z gmsh
    std::vector<bool>   isDirichlet;                        // jestli je v uzlu predepsan Dirichlet

    void Load(const std::string& filename);
    void Read(const std::string& filename);
    void VectorSave(std::vector<double>& u, const std::string& filename);
    void VectorData(std::vector<double>& u);
};

class element{
public:
    int typ;                                // typ elementu
    int idx;                                // index elementu
    int mark;                               // oznaceni "physical" z gmsh
    int idxA, idxB, idxC;                   // indexy uzlu
    int idxBaseFn[3];                       // pomocne pole pro Dirichleta
    double A[2],  B[2],  C[2];              // souradnice vrcholu
    double Sa[2], Sb[2], Sc[2];             // souradnice stredu stran
    double T[2];                            // souradnice teziste
    double matB[2][2];                      // transformaci matice
    double invB[2][2];                      // inverze transformacni matice
    double detB;                            // determinant transformaci matice
    double vol;                             // obsah elementu

    void GetElement(int i, Mesh *ptr);
    void Info(int i, Mesh *ptr);
};

class gradfi{
public:
    double fiA[2];                          // gradient bazove funkce fiA
    double fiB[2];                          // gradient bazove funkce fiB
    double fiC[2];                          // gradient bazove funkce fiC

    double fiAfiA;                          // skalarni soucin gradientu bazovych funkci fiA a fiA
    double fiAfiB;                          // skalarni soucin gradientu bazovych funkci fiA a fiB
    double fiAfiC;                          // skalarni soucin gradientu bazovych funkci fiA a fiC
    double fiBfiB;                          // skalarni soucin gradientu bazovych funkci fiB a fiB
    double fiBfiC;                          // skalarni soucin gradientu bazovych funkci fiB a fiC
    double fiCfiC;                          // skalarni soucin gradientu bazovych funkci fiC a fiC

    void GetBaseFnGradient(element *K);
};

class quadrature{
public:
    double xh[3];       // x-souradnice uzlu numericke kvadratury na ref. trojuhelniku
    double yh[3];       // y-souradnice uzlu numericke kvadratury na ref. trojuhelniku
    double wh[3];       // vahy hodnot pro jednotlive uzly numericke kvadratury na ref. trojuhelniku

    double fi[3][3];    // matice hodnot bazovych funkci v uzlech numericke kvadratury

    quadrature();
};

class gmshline{
public:
    int idx;                            // index elementu
    int etyp;                           // typ elementu
    int markPhysical;                   // oznaceni "physical" z gmsh
    int ilist[3];                       // pomocne pole pro uzly elementu

    int Read(const std::string& radek);
};

#endif // MESH_H_INCLUDED
