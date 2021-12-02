#include <cstdlib>
#include <iostream>
#include "mesh.h"
/* ----------------------------------------------------------------------------------- */
void Mesh_Allocate(mesh *p){
    // funkce pro alokaci pameti
    p->x = (double *)malloc(p->nbNods * sizeof(double)); // malloc vycleni v pocitacove pameti misto o velikost p->nbNods bytu
    p->y = (double *)malloc(p->nbNods * sizeof(double));

    p->BndrA = (int *)malloc(p->nbBndrEdges * sizeof(int));
    p->BndrB = (int *)malloc(p->nbBndrEdges * sizeof(int));
    p->BndrMark = (int *)malloc(p->nbBndrEdges * sizeof(int));

    p->TriA = (int *)malloc(p->nbTriangles * sizeof(int));
    p->TriB = (int *)malloc(p->nbTriangles * sizeof(int));
    p->TriC = (int *)malloc(p->nbTriangles * sizeof(int));
    p->TriMark = (int *)malloc(p->nbTriangles * sizeof(int));

    p->isDirichlet = (int *)calloc(p->nbNods, sizeof(int)); // calloc totez jako malloc, akorat vyclenou pamet jeste vynuluje (malloc tam necha nahodne hodnoty)
}
/* ----------------------------------------------------------------------------------- */
void Mesh_Free(mesh *p){
    // funkce pro dealokaci pameti
    free(p->x);
    free(p->y);

    free(p->BndrA);
    free(p->BndrB);
    free(p->BndrMark);

    free(p->TriA);
    free(p->TriB);
    free(p->TriC);
    free(p->TriMark);
}
/* ----------------------------------------------------------------------------------- */
void Mesh_SaveVTK(mesh *p, double *u, const char *fname){
    // funkce pro ulozeni vektoru u do formatu .vtk, ktery otevre paraview (funkce je skoro identicka s matlabovskou)
    int i;
    FILE *fid; // syntaxe pro praci se souborem, soubor pojmenuje interne fid (muze byt jakekoliv jmeno)

    fid = fopen(fname, "w"); // otevreni souboru se jmenem fname (zadavame v main.cpp) v modu w-writing(zapisujeme nove hodnoty)
    if(fid == NULL){ // pro kontrolu, kdyby se nahodou z nejakeho duvodu nepodarilo otevrit soubor
        printf("Error, nepodarilo se otevrit soubor pro format .VTK\n");
        exit(1);
    }

    fprintf(fid, "# vtk DataFile Version 2.0\n");
    fprintf(fid, "Comment: Scalar data\n");
    fprintf(fid, "ASCII\n");
    fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fid, "POINTS %d float\n", p->nbNods);
    for(i = 0; i < p->nbNods; i++){
        fprintf(fid, "%g %g %g\n", p->x[i], p->y[i], 0.);
    }

    fprintf(fid, "\nCELLS %d %d\n", p->nbTriangles, 4 * p->nbTriangles);
    for(i = 0; i < p->nbTriangles; i++){
        fprintf(fid, "3 %d %d %d\n", p->TriA[i], p->TriB[i], p->TriC[i]);
    }

    fprintf(fid, "\n\nCELL_TYPES %d\n", p->nbTriangles);
    for(i = 0; i < p->nbTriangles; i++){
        fprintf(fid, "%d\n", 5); //VTK_TRIANGLE
    }

    fprintf(fid, "\nPOINT_DATA %d\n", p->nbNods);
    fprintf(fid, "SCALARS %s float 1\n", "u");
    fprintf(fid, "LOOKUP_TABLE default\n");
    for(i = 0; i < p->nbNods; i++){
        fprintf(fid, "%g\n", u[i]);
    }

    fclose(fid); // uzavreni souboru fid
}
