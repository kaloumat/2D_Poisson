#define DEBUG
#include <iostream>
#include <string.h>
#include <cmath>
#include <vector>
#include "mesh.h"

/* ----------------------------- Mesh ------------------------------------------------ */
void Mesh::Load(const std::string& filename){
    Mesh grid;
    grid.Read(filename);

    nbNods      = grid.nbNods;
    nbBndrEdges = grid.nbBndrEdges;
    nbTriangles = grid.nbTriangles;

    double xcoord, ycoord;

    BndrA.resize(nbBndrEdges);
    BndrB.resize(nbBndrEdges);
    BndrMark.resize(nbBndrEdges);

    TriA.resize(nbTriangles);
    TriB.resize(nbTriangles);
    TriC.resize(nbTriangles);
    TriMark.resize(nbTriangles);

    isDirichlet.resize(nbNods);

#ifdef DEBUG
    printf("pocet uzlu = %d\n", nbNods);
    printf("pocet usecek = %d\n", nbBndrEdges);
    printf("pocet trojuhelniku = %d\n\n", nbTriangles);
#endif // DEBUG

    int      i, nbGeomElements, idx;
    int      i1, i2; i1 = i2 = 0;
    FILE     *fid;
    char     radek[500];
    gmshline current;

    fid = fopen(filename.c_str(), "r");
    if (fid == NULL){
        printf("Error, nepodarilo se otevrit soubor ve formatu .MSH");
        exit(1);
    }
    while (!feof(fid)){
        fgets(radek, sizeof(radek), fid);
        if (memcmp(radek, "$Nodes", 6) == 0){
            fgets(radek, sizeof(radek), fid); // dalsi radek
            for (i = 0; i < nbNods; i++){
                fgets(radek, sizeof(radek), fid);
                sscanf(radek, "%d %lf %lf", &idx, &xcoord, &ycoord);
                x.push_back(xcoord);
                y.push_back(ycoord);
            }
        }
        if (memcmp(radek, "$Elements", 9) == 0){
            fgets(radek, sizeof(radek), fid); // dalsi radek
            sscanf(radek, "%d", &nbGeomElements);
            for(i = 0; i < nbGeomElements; i++){
                fgets(radek, sizeof(radek), fid);
                current.Read(radek);
                switch (current.etyp){
                case GMSH_SEGMENT:
                    BndrA[i1]    = current.ilist[0];
                    BndrB[i1]    = current.ilist[1];
                    BndrMark[i1] = current.markPhysical;
                    if(current.markPhysical == 111)
                        isDirichlet[BndrA[i1]] = true;
                    i1++;
                    break;
                case GMSH_TRIANGLE:
                    TriA[i2]    = current.ilist[0];
                    TriB[i2]    = current.ilist[1];
                    TriC[i2]    = current.ilist[2];
                    TriMark[i2] = current.markPhysical;
                    i2++;
                    break;
                default:
                    printf("\n\nError[MSH]:\t Unknown element %d, type %d!\n", i + 1, current.etyp);
                    exit(1);
                    break;
                }
            }
        }
    }

    fclose(fid);
}
/* ----------------------------------------------------------------------------------- */
void Mesh::Read(const std::string& filename){
    int      i, nbGeomElements;
    FILE     *fid;
    char     radek[500];
    gmshline current;

    nbNods = 0;
    nbBndrEdges = 0;
    nbTriangles = 0;

    fid = fopen(filename.c_str(), "r");
    if (fid == NULL){
        printf("Error, nepodarilo se otevrit soubor ve formatu .MSH");
        exit(1);
    }
    while (!feof(fid)){
        fgets(radek, sizeof(radek), fid);
        if (memcmp(radek, "$Nodes", 6) == 0){
            fgets(radek, sizeof(radek), fid); //dalsi radek
            sscanf(radek, "%d", &nbNods);
        }
        if (memcmp(radek, "$Elements", 9) == 0){
            fgets(radek, sizeof(radek), fid); //dalsi radek
            sscanf(radek, "%d", &nbGeomElements);
            for(i = 0; i < nbGeomElements; i++){
                fgets(radek, sizeof(radek), fid);
                current.Read(radek);
                switch (current.etyp){
                case GMSH_SEGMENT:     nbBndrEdges++; break;
                case GMSH_TRIANGLE:    nbTriangles++; break;
                default:
                    printf("\n\nError[MSH]:\t Neznamy element %d, type %d!\n", i + 1, current.etyp);
                    exit(1);
                    break;
                }
            }
        }
    }
    fclose(fid);
}


/* ----------------------------- element --------------------------------------------- */
void element::GetElement(int i, Mesh *ptr){
    typ  = GMSH_TRIANGLE;
    idx  = i;

    idxA = ptr->TriA[i];
    idxB = ptr->TriB[i];
    idxC = ptr->TriC[i];
    mark = ptr->TriMark[i];

    idxBaseFn[0] = idxA + 1;
    idxBaseFn[1] = idxB + 1;
    idxBaseFn[2] = idxC + 1;

    if(ptr->isDirichlet[idxA] == true)
        idxBaseFn[0] = - idxBaseFn[0];

    if(ptr->isDirichlet[idxB] == true)
        idxBaseFn[1] = - idxBaseFn[1];

    if(ptr->isDirichlet[idxC] == true)
        idxBaseFn[2] = - idxBaseFn[2];

    A[0] = ptr->x[idxA];     A[1] = ptr->y[idxA];
    B[0] = ptr->x[idxB];     B[1] = ptr->y[idxB];
    C[0] = ptr->x[idxC];     C[1] = ptr->y[idxC];

    Sa[0] = 1/2. * (B[0] + C[0]);      Sa[1] = 1/2. * (B[1] + C[1]);
    Sb[0] = 1/2. * (A[0] + C[0]);      Sb[1] = 1/2. * (A[1] + C[1]);
    Sc[0] = 1/2. * (A[0] + B[0]);      Sc[1] = 1/2. * (A[1] + B[1]);

    T[0] = 2/3. * (Sa[0] - A[0]) + A[0];
    T[1] = 2/3. * (Sa[1] - A[1]) + A[1];

    matB[0][0] = B[0] - A[0];      matB[0][1] = C[0] - A[0];
    matB[1][0] = B[1] - A[1];      matB[1][1] = C[1] - A[1];

    detB = matB[0][0] * matB[1][1] - matB[0][1] * matB[1][0];

    invB[0][0] =  1 / detB * matB[1][1];    invB[0][1] = -1 / detB * matB[0][1];
    invB[1][0] = -1 / detB * matB[1][0];    invB[1][1] =  1 / detB * matB[0][0];

    vol = detB/2.;
}
/* ----------------------------------------------------------------------------------- */
void element::Info(int i, Mesh *ptr){
    GetElement(i, ptr);
    printf("\ntrojuhelnik:\n");
    printf("typ elementu: %d\n", typ);
    printf("index elementu: %d\n", idx);
    printf("oznaceni: %d\n", mark);
    printf("indexy uzlu: A = %d, B = %d, C = %d\n", idxA, idxB, idxC);
    printf("indexy uzlu v ramci Dirichleta: A = %d, B = %d, C = %d\n", idxBaseFn[0], idxBaseFn[1], idxBaseFn[2]);
    printf("souradnice vrcholu: A[%lf, %lf], B[%lf, %lf] C[%lf, %lf]\n", A[0], A[1], B[0], B[1], C[0], C[1]);
    printf("souradnice stredu: Sa[%lf, %lf], Sb[%lf, %lf] Sc[%lf, %lf]\n", Sa[0], Sa[1], Sb[0], Sb[1], Sc[0], Sc[1]);
    printf("souradnice teziste: T[%lf, %lf]\n", T[0], T[1]);
    printf("transformacni matice: [%lf   %lf]\n", matB[0][0], matB[0][1]);
    printf("transformacni matice: [%lf   %lf]\n", matB[1][0], matB[1][1]);
    printf("inverzni matice: [%lf   %lf]\n", invB[0][0], invB[0][1]);
    printf("inverzni matice: [%lf   %lf]\n", invB[1][0], invB[1][1]);
    printf("determinant matice: %lf\n", detB);
    printf("obsah elementu: %lf\n", vol);
}


/* ------------------------------ gradfi --------------------------------------------- */
void gradfi::GetBaseFnGradient(element *K){
    fiA[0] = - K->invB[0][0] - K->invB[1][0];
    fiA[1] = - K->invB[0][1] - K->invB[1][1];

    fiB[0] = K->invB[0][0];
    fiB[1] = K->invB[0][1];

    fiC[0] = K->invB[1][0];
    fiC[1] = K->invB[1][1];

    fiAfiA = fiA[0] * fiA[0] + fiA[1] * fiA[1];
    fiAfiB = fiA[0] * fiB[0] + fiA[1] * fiB[1];
    fiAfiC = fiA[0] * fiC[0] + fiA[1] * fiC[1];
    fiBfiB = fiB[0] * fiB[0] + fiB[1] * fiB[1];
    fiBfiC = fiB[0] * fiC[0] + fiB[1] * fiC[1];
    fiCfiC = fiC[0] * fiC[0] + fiC[1] * fiC[1];
}


/* ---------------------------- quadrature ------------------------------------------- */
quadrature::quadrature(){
    //uzly - stredy stran ref. trojuhelnika
        //Sa              Sb              Sc
    xh[0] = 1/2.;    xh[1] = 0;       xh[2] = 1/2.;
    yh[0] = 1/2.;    yh[1] = 1/2.;    yh[2] = 0;
    wh[0] = 1/3.;    wh[1] = 1/3.;    wh[2] = 1/3.;

    fi[0][0] = - xh[0] - yh[0] + 1;   fi[0][1] = - xh[1] - yh[1] + 1;    fi[0][2] = - xh[2] - yh[2] + 1;
    fi[1][0] = xh[0];                 fi[1][1] = xh[1];                  fi[1][2] = xh[2];
    fi[2][0] = yh[0];                 fi[2][1] = yh[1];                  fi[2][2] = yh[2];
}


/* ----------------------------- gmshline -------------------------------------------- */
int gmshline::Read(const std::string& radek){
    int idx, etp, tgs, tgs1; //index //typ elementu //pocet oznaceni v gmsh
    int val;
    int tags[10], pom[10];
    int nred;

    nred = sscanf(radek.c_str(), "%d %d %d", &idx, &etp, &tgs);
    if (nred != 3){
        printf("Error[MSH]:\t GMSH format je spatne, nebo spatny radek\n");
        exit(1);
    }
    etyp = etp;
    this->idx  = idx;
    tags[0] = tags[1] = tags[2] = tags[3] = tags[4] = tags[5] = pom[0] = pom[1] = pom[2] = pom[3] = pom[4] = 0;

    val  = tgs + etp * 10;
    nred = 0;
    switch (val){
    case 12:
        nred = sscanf(radek.c_str(), "%d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1]);
        break;
    case 22:
        nred = sscanf(radek.c_str(), "%d %d %d %d %d %d %d %d", &idx, &etp, &tgs1, &tags[0], &tags[1], &pom[0], &pom[1], &pom[2]);
        break;
    default:
        printf("Error[MSH]:\t Incorrect number of tags (NTAGS %d)\n", tgs);
        break;
    }

    markPhysical = tags[0];
    switch (etp){
    case GMSH_SEGMENT:
        if (nred != tgs + 3 + 2)
            printf("Error[MSH]:\t Spatny format site pro usecky\n");
        ilist[0] = pom[0] - 1;
        ilist[1] = pom[1] - 1;
        break;
    case GMSH_TRIANGLE:
        if (nred != tgs + 3 + 3)
            printf("Error[MSH]:\t Spatny format site pro trojuhelniky\n");
        ilist[0] = pom[0] - 1;
        ilist[1] = pom[1] - 1;
        ilist[2] = pom[2] - 1;
        break;
    }
    return 0;
}
