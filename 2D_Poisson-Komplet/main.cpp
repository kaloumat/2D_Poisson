#define DEBUG
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include "triplet.h"
#include "mesh.h"
#include "sparse.h"

using namespace std;

double mojefce(double x, double y){
        return 2 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y);
}

double Dirichlet(double x, double y){
    return sin(M_PI * x) * sin(M_PI * y);
}

int main(){
    int i, k;
    double x, y;

    Mesh sit;
    sit.Load("ctverecnxn.msh");

    Triplet tri(sit.nbNods, 0);

    vector<double> b, u;
    b.resize(sit.nbNods);
    u.resize(sit.nbNods);

    element K;
    quadrature Q;
    for(k = 0; k < sit.nbTriangles; k++){
        K.GetElement(k, &sit);

        double qval[3] = {0, 0, 0};
        for(i = 0; i < 3; i++){
            x = K.A[0] + K.matB[0][0] * Q.xh[i] + K.matB[0][1] * Q.yh[i];
            y = K.A[1] + K.matB[1][0] * Q.xh[i] + K.matB[1][1] * Q.yh[i];

            qval[0] += Q.wh[i] * K.vol * mojefce(x, y) * Q.fi[0][i];
            qval[1] += Q.wh[i] * K.vol * mojefce(x, y) * Q.fi[1][i];
            qval[2] += Q.wh[i] * K.vol * mojefce(x, y) * Q.fi[2][i];
        }

        b[K.idxA] += qval[0];
        b[K.idxB] += qval[1];
        b[K.idxC] += qval[2];
    }

    for(i = 0; i < sit.nbBndrEdges; i++){
        if(sit.isDirichlet[sit.BndrA[i]] == true)
            b[sit.BndrA[i]] = Dirichlet(sit.x[sit.BndrA[i]], sit.y[sit.BndrA[i]]);
    }

    for(i = 0; i < sit.nbNods; i++){
        if(sit.isDirichlet[i] == true)
            tri.Add(i + 1, i + 1, 1);
    }

    gradfi grad;
    for(k = 0; k < sit.nbTriangles; k++){
        K.GetElement(k, &sit);
        grad.GetBaseFnGradient(&K);

        tri.AddRHS(b, K.idxBaseFn[0], K.idxBaseFn[0], K.vol * grad.fiAfiA);
        tri.AddRHS(b, K.idxBaseFn[0], K.idxBaseFn[1], K.vol * grad.fiAfiB);
        tri.AddRHS(b, K.idxBaseFn[0], K.idxBaseFn[2], K.vol * grad.fiAfiC);

        tri.AddRHS(b, K.idxBaseFn[1], K.idxBaseFn[0], K.vol * grad.fiAfiB);
        tri.AddRHS(b, K.idxBaseFn[1], K.idxBaseFn[1], K.vol * grad.fiBfiB);
        tri.AddRHS(b, K.idxBaseFn[1], K.idxBaseFn[2], K.vol * grad.fiBfiC);

        tri.AddRHS(b, K.idxBaseFn[2], K.idxBaseFn[0], K.vol * grad.fiAfiC);
        tri.AddRHS(b, K.idxBaseFn[2], K.idxBaseFn[1], K.vol * grad.fiBfiC);
        tri.AddRHS(b, K.idxBaseFn[2], K.idxBaseFn[2], K.vol * grad.fiCfiC);
    }

    Sparse A;
    A.Create(&tri);

#ifdef DEBUG
    cout << "\nvektor b:" << endl;
    for(i = 0; i < sit.nbNods; i++)
        cout << "b[" << i << "] = " << fixed << b[i] << endl;

    cout << "\nA.n = " << A.n << endl;
    cout << "A.nz = " << A.nz << endl;

    cout << "\nmatice tuhosti A:" << endl;
    for(i = 0; i < A.n; i++){
        for(int j = A.PI[i]; j < A.PI[i + 1]; j++){
            cout << i << " " << A.J_VAL[j].first << " " << A.J_VAL[j].second << endl;
        }
    }

    cout << "\nreziduum GaussSeidel:" << endl;
#endif // DEBUG

    A.GaussSeidel(b, u, 1e5, 1e-10);

#ifdef DEBUG
    cout << "\nvysledny vektor u:" << endl;
    for(i = 0; i < A.n; i++){
        cout << fixed << u[i] << endl;
    }
#endif // DEBUG

    sit.VectorData(u);
    sit.VectorSave(u, "poisson2D.txt");

    return 0;
}
