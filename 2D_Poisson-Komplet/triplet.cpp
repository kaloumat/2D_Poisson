#include <vector>
#include "triplet.h"
/* ----------------------------------------------------------------------------------- */
void Triplet::Add(int _i, int _j, double _val){
    I.push_back(_i);
    J.push_back(_j);
    VAL.push_back(_val);

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
void Triplet::Quicksort(int first, int last){
    int i, pivot;
    double temp;

    if(first < last){
        pivot = first;
        for(i = first + 1; i < last; i++){
            if(I[i] < I[first]){
                pivot++;

                temp = I[pivot];
                I[pivot] = I[i];
                I[i] = temp;

                temp = J[pivot];
                J[pivot] = J[i];
                J[i] = temp;

                temp = VAL[pivot];
                VAL[pivot] = VAL[i];
                VAL[i] = temp;
            }
            if(I[i] == I[first] && J[i] < J[first]){
                pivot++;

                temp = I[pivot];
                I[pivot] = I[i];
                I[i] = temp;

                temp = J[pivot];
                J[pivot] = J[i];
                J[i] = temp;

                temp = VAL[pivot];
                VAL[pivot] = VAL[i];
                VAL[i] = temp;
            }
        }

        temp = I[pivot];
        I[pivot] = I[first];
        I[first] = temp;

        temp = J[pivot];
        J[pivot] = J[first];
        J[first] = temp;

        temp = VAL[pivot];
        VAL[pivot] = VAL[first];
        VAL[first] = temp;

        Quicksort(first, pivot);
        Quicksort(pivot + 1, last);
    }
}
/* ----------------------------------------------------------------------------------- */
void Triplet::Unique(){
    int i, j = 1;

    for(i = 0; i < nz; i++){
        while(I[i] == I[j] && J[i] == J[j]){
            VAL[i] += VAL[j];
            j++;
        }
        if(j < nz){
            I[i + 1] = I[j];
            J[i + 1] = J[j];
            VAL[i + 1] = VAL[j];
            j++;
        }
        else
            break;
    }
    nz = i + 1;
}
