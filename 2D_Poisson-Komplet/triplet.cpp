#include <iostream>
#include <vector>
#include <algorithm>
#include "triplet.h"
/* ----------------------------------------------------------------------------------- */
void Triplet::Add(int _i, int _j, double _val){
    triplet.push_back({_i, _j, _val});

    nz++;
}
/* ----------------------------------------------------------------------------------- */
void Triplet::AddRHS(std::vector<double>& b, int i, int j, double qval){
    if(i > 0){
        if(j > 0){
            Add(i, j, qval);
        }
        else if(j < 0){
            b[i - 1] -= qval * b[-j - 1];
        }
    }
}
/* ----------------------------------------------------------------------------------- */
void Triplet::Unique(){
    sort(triplet.begin(), triplet.end());

    int i, j = 1;

    for(i = 0; i < nz; i++){
        while(triplet[i].I == triplet[j].I && triplet[i].J == triplet[j].J){
            triplet[i].VAL += triplet[j].VAL;
            j++;
        }
        if(j < nz){
            triplet[i + 1].I   = triplet[j].I;
            triplet[i + 1].J   = triplet[j].J;
            triplet[i + 1].VAL = triplet[j].VAL;
            j++;
        }
        else
            break;
    }
    nz = i + 1;
}
