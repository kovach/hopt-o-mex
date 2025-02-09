// sxMatrixInterface.cpp

#include "math.h"
#include <cstdio>
#include <sstream>
#include <string.h>

#include <casadi/sx/sx.hpp>
#include <casadi/sx/sx_tools.hpp>
#include <casadi/matrix/matrix_tools.hpp>

#include "sxMatrixInterface.hpp"

//#define COUT_MEMORY_MANAGEMENT

using namespace std;
using namespace CasADi;


/******************** memory management *******************/
int sxMatrixSizeOfAddress(){
  return sizeof(SXMatrix*);
}

SXMatrix * sxMatrixCreateSymbolic(char * charPrefix, int n, int m){
  string prefix;
  prefix.assign(charPrefix);
  SXMatrix * out = new SXMatrix(create_symbolic(prefix, n, m));
  #ifdef COUT_MEMORY_MANAGEMENT
  cout << "(cpp) new sx matrix at " << out << ", val: " << *out << endl;
  #endif
  return out;
}

SXMatrix * sxMatrixDuplicate(const SXMatrix & old){
  SXMatrix * out = new SXMatrix(old);
  #ifdef COUT_MEMORY_MANAGEMENT
  cout << "(cpp) duplicate " << out << ", val: " << *out << endl;
  #endif
  return out;
}

void sxMatrixDuplicateAt(const SXMatrix & source, int idx, SXMatrix * destination[]){
  destination[idx] = new SXMatrix(source);

  #ifdef COUT_MEMORY_MANAGEMENT
  cout << "(cpp) duplicateAt destination {head: " << destination << ", address: " << destination[idx] << ", val: " << *(destination[idx]) << "}\n";
  #endif
}

void sxMatrixFreeArray(SXMatrix * matArray[], int n){
  for (int idx=0; idx<n; idx++){
    #ifdef COUT_MEMORY_MANAGEMENT
    cout << "(cpp) freeArray deleting {head: " << matArray << ", address: " << matArray[idx] << ", val: " << *(matArray[idx]) << "}\n";
    #endif
    delete matArray[idx];
  }
}


void sxMatrixDelete(SXMatrix * const sx){
  #ifdef COUT_MEMORY_MANAGEMENT
  cout << "(cpp) deleting sx matrix at " << sx << endl;
  #endif
  delete sx;
}

SXMatrix * sxMatrixZeros(int n, int m){
  SXMatrix * out = new SXMatrix( zerosSX(n,m) );
  #ifdef COUT_MEMORY_MANAGEMENT
  cout << "(cpp) new sx zeros at " << out << ", val: " << *out << endl;
  #endif
  return out;
}


/******************** show *******************/
void sxMatrixShow(char * stringOut, int strLen, const SXMatrix & sx){
  ostringstream sxOutStream;
  sxOutStream << sx;
  strncpy(stringOut, sxOutStream.str().c_str(), strLen);

  if (sxOutStream.str().length() > strLen)
    cerr << "(cpp) ERROR - sxMatrixShow trying to write " << sxOutStream.str().length() << " characters to output string with capacity of " << strLen << " characters\n";
}


/******************** accessors *******************/
void sxMatrixAt(const SXMatrix & mat, int n, int m, SX & out){
  out = SX(mat.indexed(n,m));
}

void sxMatrixSet(const SX & sx, int n, int m, SXMatrix & mat){
  mat.indexed_assignment(n, m, sx);
}

int sxMatrixSize1(const SXMatrix & mat){
  return mat.size1();
}

int sxMatrixSize2(const SXMatrix & mat){
  return mat.size2();
}


/******************** math *******************/
void sxMatrixPlus(const SXMatrix & m0, const SXMatrix & m1, SXMatrix & mOut){
  mOut = m0 + m1;
}

void sxMatrixMinus(const SXMatrix & m0, const SXMatrix & m1, SXMatrix & mOut){
  mOut = m0 - m1;
}

void sxMM(const SXMatrix & m0, const SXMatrix & m1, SXMatrix & mOut){
  mOut = prod(m0, m1);
}

void sxMatrixTranspose(const SXMatrix & mIn, SXMatrix & mOut){
  mOut = SXMatrix(mIn.trans());
}

int sxMatrixIsEqual(const SXMatrix & m0, const SXMatrix & m1){
  if (isEqual(m0, m1))
    return 1;
  return 0;
}
