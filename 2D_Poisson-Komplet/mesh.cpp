#define DEBUG
#include <iostream>
#include <iterator>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>
#include "mesh.h"

using namespace std;

/* ----------------------------- Mesh ------------------------------------------------ */
void Mesh::Load(const string& filename){
    Mesh grid;
    grid.Read(filename);

    nbNods      = grid.nbNods;
    nbBndrEdges = grid.nbBndrEdges;
    nbTriangles = grid.nbTriangles;

    isDirichlet.resize(nbNods);

#ifdef DEBUG
    cout << "pocet uzlu = " << nbNods << endl;
    cout << "pocet usecek = " << nbBndrEdges << endl;
    cout << "pocet trojuhelniku = " << nbTriangles << endl;
#endif // DEBUG

    int      i, nbGeomElements, idx, ii = 0;
    char     buff[100];
    double   xcoord, ycoord;
    string   radek;

    gmshline current;

    ifstream vstup(filename);
    if(!vstup.is_open()){
        cout << "Error, nepodarilo se otevrit soubor ve formatu .MSH" << endl;
        exit(1);
    }
    while (getline(vstup, radek)){
        if (radek == "$Nodes"){
            vstup.getline(buff, sizeof(buff));
            for (i = 0; i < nbNods; i++){
                getline(vstup, radek);
                istringstream ss(radek);
                ss >> idx >> xcoord >> ycoord;
                x.push_back(xcoord);
                y.push_back(ycoord);
            }
        }
        if (radek == "$Elements"){
            getline(vstup, radek);
            istringstream ss(radek);
            ss >> nbGeomElements;
            for(i = 0; i < nbGeomElements; i++){
                getline(vstup, radek);
                current.Read(radek);
                switch (current.etyp){
                case GMSH_SEGMENT:
                    BndrA.push_back(current.ilist[0]);
                    BndrB.push_back(current.ilist[1]);
                    BndrMark.push_back(current.markPhysical);
                    if(current.markPhysical == 111)
                        isDirichlet[BndrA[ii]] = true;
                    ii++;
                    break;
                case GMSH_TRIANGLE:
                    TriA.push_back(current.ilist[0]);
                    TriB.push_back(current.ilist[1]);
                    TriC.push_back(current.ilist[2]);
                    TriMark.push_back(current.markPhysical);
                    break;
                default:
                    cout << "\nError[Mesh]: Neznamy element " << i + 1 << ", typ " << current.etyp << endl;
                    exit(1);
                    break;
                }
            }
        }
    }
    vstup.close();
}
/* ----------------------------------------------------------------------------------- */
void Mesh::Read(const string& filename){
    int      i, nbGeomElements;
    string   radek;

    gmshline current;

    nbNods      = 0;
    nbBndrEdges = 0;
    nbTriangles = 0;

    ifstream vstup(filename);
    if(!vstup.is_open()){
        cout << "Error, nepodarilo se otevrit soubor ve formatu .MSH" << endl;
        exit(1);
    }
    while (getline(vstup, radek)){
        if (radek == "$Nodes"){
            getline(vstup, radek);
            istringstream ss(radek);
            ss >> nbNods;
        }
        if (radek == "$Elements"){
            getline(vstup, radek);
            istringstream ss(radek);
            ss >> nbGeomElements;
            for(i = 0; i < nbGeomElements; i++){
                getline(vstup, radek);
                current.Read(radek);
                switch (current.etyp){
                case GMSH_SEGMENT:
                    nbBndrEdges++;
                    break;
                case GMSH_TRIANGLE:
                    nbTriangles++;
                    break;
                default:
                    cout << "\nError[Mesh]: Neznamy element " << i + 1 << ", typ " << current.etyp << endl;
                    exit(1);
                    break;
                }
            }
        }
    }
    vstup.close();
}
/* ----------------------------------------------------------------------------------- */
void Mesh::VectorSave(vector<double>& u, const string& filename){
    ofstream vystup(filename);
    vystup << fixed;
    copy(u.begin(), u.end(), ostream_iterator<double>(vystup, "\n"));
    vystup.close();
}
/* ----------------------------------------------------------------------------------- */
void Mesh::VectorData(vector<double>& u){
    int i;

    struct Data{
        double Ydata;
        double Udata;
    };

    vector<Data> DataX0_2;
    for(i = 0; i < nbNods; i++){
        if(x[i] == 0.2){
            DataX0_2.push_back({y[i], u[i]});
        }
    }

    vector<Data> DataX0_5;
    for(i = 0; i < nbNods; i++){
        if(x[i] == 0.5){
            DataX0_5.push_back({y[i], u[i]});
        }
    }

    vector<Data> DataX0_9;
    for(i = 0; i < nbNods; i++){
        if(x[i] == 0.9){
            DataX0_9.push_back({y[i], u[i]});
        }
    }

    sort(DataX0_2.begin(), DataX0_2.end(), [](const auto& j, const auto& k) { return j.Ydata < k.Ydata; });
    sort(DataX0_5.begin(), DataX0_5.end(), [](const auto& j, const auto& k) { return j.Ydata < k.Ydata; });
    sort(DataX0_9.begin(), DataX0_9.end(), [](const auto& j, const auto& k) { return j.Ydata < k.Ydata; });

    ofstream vystupX0_2("X0_2.txt");
    for(size_t i = 0; i < DataX0_2.size(); i++)
        vystupX0_2 << fixed << DataX0_2[i].Ydata << " " << DataX0_2[i].Udata << endl;
    vystupX0_2.close();

    ofstream vystupX0_5("X0_5.txt");
    for(size_t i = 0; i < DataX0_5.size(); i++)
        vystupX0_5 << fixed << DataX0_5[i].Ydata << " " << DataX0_5[i].Udata << endl;
    vystupX0_5.close();

    ofstream vystupX0_9("X0_9.txt");
    for(size_t i = 0; i < DataX0_9.size(); i++)
        vystupX0_9 << fixed << DataX0_9[i].Ydata << " " << DataX0_9[i].Udata << endl;
    vystupX0_9.close();
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
    cout << "\ntrojuhelnik "<< idx << ": " << endl;
    cout << "oznaceni: " << mark << endl;
    cout << "indexy uzlu: A = " << idxA << ", B = " << idxB << ", C = " << idxC << endl;
    cout << "indexy uzlu v ramci Dirichleta: A = " << idxBaseFn[0] << ", B = " << idxBaseFn[1] << ", C = " << idxBaseFn[2] << endl;
    cout << "souradnice vrcholu: A[" << A[0] << ", " << A[1] << "], B[" << B[0] << ", " << B[1] << "], C[" << C[0] << ", " << C[1] << "]" << endl;
    cout << "souradnice stredu: Sa[" << Sa[0] << ", " << Sa[1] << "], Sb[" << Sb[0] << ", " << Sb[1] << "], Sc[" << Sc[0] << ", " << Sc[1] << "]" << endl;
    cout << "souradnice teziste: T[" << T[0] << ", " << T[1] << "]" << endl;
    cout << "transformacni matice: [" << matB[0][0] << " " << matB[0][1] << "]" << endl;
    cout << "transformacni matice: [" << matB[1][0] << " " << matB[1][1] << "]" << endl;
    cout << "inverzni matice: [" << invB[0][0] << ", " << invB[0][1] << "]" << endl;
    cout << "inverzni matice: [" << invB[1][0] << ", " << invB[1][1] << "]" << endl;
    cout << "determinant transformacni matice: " << detB << endl;
    cout << "obsah trojuhelniku: " << vol << endl;
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
int gmshline::Read(const string& radek){
    int idx, etyp, tgs;    // index // typ elementu // pocet oznaceni v gmsh
    int val;               // ridici cislo pro rozliseni elementu
    int tags[2], pom[3];
    int nred;

    nred = sscanf(radek.c_str(), "%d %d %d", &idx, &etyp, &tgs);
    if (nred != 3){
        cout << "Error[Mesh]: GMSH format je spatne, nebo spatny radek" << endl;
        exit(1);
    }

    this->idx  = idx;
    this->etyp = etyp;

    tags[0] = tags[1] = pom[0] = pom[1] = pom[2] = 0;
    val     = tgs + etyp * 10;
    nred    = 0;

    switch (val){
    case 12:
        nred = sscanf(radek.c_str(), "%d %d %d %d %d %d %d", &idx, &etyp, &tgs, &tags[0], &tags[1], &pom[0], &pom[1]);
        break;
    case 22:
        nred = sscanf(radek.c_str(), "%d %d %d %d %d %d %d %d", &idx, &etyp, &tgs, &tags[0], &tags[1], &pom[0], &pom[1], &pom[2]);
        break;
    default:
        cout << "Error[Mesh]: Spatny pocet tagu: " << tgs << endl;
        break;
    }

    markPhysical = tags[0];

    switch (etyp){
    case GMSH_SEGMENT:
        if (nred != tgs + 3 + 2)
            cout << "Error[Mesh]: Spatny format site pro usecky" << endl;
        ilist[0] = pom[0] - 1;
        ilist[1] = pom[1] - 1;
        break;
    case GMSH_TRIANGLE:
        if (nred != tgs + 3 + 3)
            cout << "Error[Mesh]: Spatny format site pro trojuhelniky" << endl;
        ilist[0] = pom[0] - 1;
        ilist[1] = pom[1] - 1;
        ilist[2] = pom[2] - 1;
        break;
    }
    return 0;
}
