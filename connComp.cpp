/**
 * USAGE: [group_idx group_size]= affinity_connected_comp(A)
 * This function finds out the connected components from an affinity matrix
 * idx(i) stores the id of the connected component that contains element i.
 * group_size(i) stores the number of elements in this group
 * We use Union-Find Algorithm
 *
 * Shugao Ma
 * Last modify: 9-24-2012
 **/
#include "mex.h"
#include "math.h"
#include <iostream>
#include <cstdio>
using namespace std;

void group(double* I, int x, int y);
int find(double* I, int x);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i = 0, j = 0;
    //check for proper number of input and output arguments
    if(nrhs<1){
        mexErrMsgTxt("Please input the adjacency matrix. \n USAGE: [group_idx]= connComp(A)");
        return;
    }
    if(nlhs!=1){
        mexErrMsgTxt("Only two output arguments are allowed. \n USAGE: [group_idx] = connComp(A)");
        return;
    }
    //check the dimensions of the input trajectory matrix
    int nDim = mxGetNumberOfDimensions(prhs[0]);
    if(nDim!=2){
        mexErrMsgTxt("The adjacency matrix should be a two dimensional square matrix.");
        return;
    }
    int rows = mxGetM(prhs[0]);
    int cols = mxGetN(prhs[0]);
    if(rows != cols){
        mexErrMsgTxt("The adjacency matrix should be a two dimensional square matrix.");
        return;
    }
    //find the connected components
    mxArray* idx = mxCreateNumericMatrix(cols, 1, mxDOUBLE_CLASS, mxREAL);
    double* A = (double*)mxGetData(prhs[0]);
    double* I = (double*)mxGetData(idx);
    for(i = 0; i < cols; i++){
        I[i] = i;
    }
    for(i = 0; i < cols-1; i++){
        for(j = i+1; j < cols; j++){
            if(A[i * cols + j] > 0){
                group(I, i, j);
            }   
        }
    }
    for(i = 0; i < cols; i++){
        int xp = find(I,i);
    }
    for(i = 0; i < cols; i++){
        I[i]+=1;
    }   
    plhs[0] = idx;
}

void group(double* I, int x, int y)
{
    int xp = find(I, x);
    int yp = find(I, y);
    if(xp == yp){
        return;
    }
    if(xp > yp){
        I[xp] = yp;
        find(I, x);
        return;
    }
    if(xp < yp){
        I[yp] = xp;
        find(I, y);
        return;
    }
}

int find(double* I, int x)
{
    if(I[x] == x){
        return (int)x;
    }
    else{
        int xp = find(I, I[x]);
        I[x] = xp;
        return (int)xp;
    }
}