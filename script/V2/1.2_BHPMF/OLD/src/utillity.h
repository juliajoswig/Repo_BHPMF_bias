//
//  utillity.h
//  
//
//  Created by Farideh Fazayeli on 7/18/13.
//
//

#ifndef ____utillity__
#define ____utillity__

#include <iostream>
//#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <iterator>
// #include <iomanip>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <algorithm>

#include <R.h>      // R functions
#include <Rmath.h>  // R math
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>


//#include <gsl/gsl_math.h>
//#include <gsl/gsl_cblas.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>


//#include <cblas.h>
//#include "source_axpy_r.h"

#include <map>

#include <sys/time.h>
#include <sys/resource.h>

using namespace std;


typedef vector<vector<float> > matVecFloat;
typedef vector<vector<int> > matVecInt;
typedef vector<set<int> >  vecSet;


double getcputime();
void PrintArr(int *myArr, int nRow);
void PrintArr(double *myArr, int nRow);
void printMat(float *myArr, int nRow, int nCol);
void printMat(float **myMat, int nRow, int nCol);
void printMat(int **myMat, int nRow, int nCol);

float **dlmreadRFloat(SEXP env, char *fileName, int &nRow, int &nCol, char *delim);
int **dlmreadRInt(SEXP env, char *fileName, int &nRow, int &nCol, char *delim);

float **dlmread1(char *filePath, int nRow, int nCol);
float **dlmread2(char *filePath, int &nRow, int &nCol, char delim);
float **dlmread2(char *filePath, int &nRow, int &nCol, char delim, int flag);
void *dlmread3(char *filePath, char delim, int &nRow, int &nCol, int flag );
float **dlmread3(char *filePath, char delim, int &nRow, int &nCol);
void init_R();

int **readFile(char *filePath, int &nRow, char delim);
SEXP getListElement (SEXP list, char *str);
void mvrnorm(double *des, double *mu, double *cholCov, int dim, bool upper);


void dlmreadVec(char *filePath, matVecFloat &myMat, int &nRow, int nCol);
void dlmreadVec(char *filePath, matVecInt &myMat, int &nRow, int nCol);
void dlmreadVec(char *filePath, matVecInt &myMat, int nRow, vector<int> nCols);
void dlmreadVec(char *filePath, vector<int> &myArr);
void dlmwriteVec(const char *inFilePath, const char *outFilePath, int nRow, int nCol, int gap, int burn, int nSam);

template <class T>
inline void SafeDelete(T x)
{
    assert(x != NULL);
    delete[] x;
    x = NULL;
}

#endif /* defined(____utillity__) */
