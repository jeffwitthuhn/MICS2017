/**
  * @file matrix.h
  * Matrix im class implementation
  * used to make a nx1 vector 
  * @authors Jeff Witthuhn
*/

#ifndef MATRIX
#define MATRIX
#include <iostream>
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/random/variate_generator.hpp"
#include "boost/random/uniform_real_distribution.hpp"
#include <ctime>
#include <chrono> 
#include <vector>
#include <cmath>
#include "newbitmatrix.h"
using namespace std;

unsigned seed2=std::chrono::system_clock::now().time_since_epoch().count();
boost::random::mt19937 gen2 (seed2);

template <class T> 
double randomNum(T min,T max) {
    boost::random::uniform_int_distribution<> dist(min, max);
    return dist(gen2);
}


/**
  * \breif utility class to implement the algorithms Matrix
  *
  * Matrix can be used for arbitrary matrices however it is used
  * as a nx1 vector for purposes of the algorithms. 
  */
class Matrix {

  private:
  unsigned long long int rows;/**< number of rows in the matrix*/
  unsigned long long int columns; /**< number of columns in the matrix*/
  std::vector<std::vector<double> > v; /**< representation of the matrix*/
  public:


  /// Constructors

  /**
    default constructor
  */
  Matrix() {
    rows=0;columns=0;
  }

  /**
    itialize size constructor
    \param rows the number of rows in the matrix
  */
  Matrix(int rows) {
    this->rows=rows;
    this->columns=1;
    v.resize(rows+1);
    for(int i=0;i<rows;i++) {
      v[i].resize(1);
      v[i][0]=0;
    }
  }

  /**
    itialize size constructor
    \param rows the number of rows in the matrix
    \param columns the number of columns in the matrix
  */
  Matrix(int rows, int columns) {
    this->rows=rows;
    this->columns=columns;
    v.resize(rows);
    for(int i=0;i<rows;i++) {
      v[i].resize(columns);
      v[i][0]=0;
    }
  }

  /// Setters


  /**
    sets a specified element to an integer value

    \param row the row index of the desired element in the matrix
    \param column the column index of the desired element in the matrix 
    \param value the value to set the desired element
  */
  void set(int r, int c, double input) {
    v[r][c]=input;
  }

  /**
    Assignment operation for Matrix
  */
  void operator =(Matrix in) {
    for(int i=0;i<rows;i++)
      for(int j=0;j<columns;j++)
        v[i][j]=in.getElement(i,j);
  }

  /**
    reinitialize the matrix
    \param rows the number of rows in the matrix
    \param columns the number of columns in the matrix
  */
  void init(int rows, int columns) {
    this->rows=rows;
    this->columns=columns;
    v.resize(rows);
    for(int i=0;i<rows;i++) {
      v[i].resize(columns);
      v[i][0]=0;
    }
  }

  /**
    sets all elements in the matrix to zero. 
  */
  void setzero() {
    for(int i=0;i<rows;i++)
      for(int j=0;j<columns;j++)
        v[i][j]=0;
  }

  /**
    sets all elements in the matrix to one. 
  */
  void setone() {
  for(int i=0;i<rows;i++)
      for(int j=0;j<columns;j++)
        v[i][j]=1;
  }

  /**
    sets all elements in the matrix to a "random" integer. 
  */
  void randomize(int lower, int upper) {
    for(int i=0;i<rows;i++)
      for(int j=0;j<columns;j++)
        v[i][j]=randomNum<int>(lower,upper);
  }

  /// Getters

  /**
    returns the number of rows in the matrix
  */
  int getrows() {
    return rows;
  }

  /**
    returns the specified element assuming it is a vector
    \param r row of the desired element
  */
  double operator()(int r) {
    return v[r][0];
  }

  /**
    returns the specified element
    \param r row of the desired element
    \param c column of the desired element
  */
  double operator()(int r, int c) {
    return v[r][c];
  }

  /**
    returns the specified element
    \param r row of the desired element
    \param c column of the desired element
  */
  double getElement(int r,int c) {
    return v[r][c];
  }

  /// Testing Methods

  /**
    prints a representation of the matrix to std out
  */
  void rprint(ostream& outfile) {
    for(int i=0;i<rows;i++) {
      for(int j=0;j<columns;j++)
        outfile<<v[i][j];
      outfile<<std::endl;
    }
  }

  /// Operations

  /**
    For the matrix-vector product algorithms,used to combine
    the resultant vectors of the partitioned products

    \param d vector of Matrices (nx1 vectors) to combine
    \param amount the amount of Matrices to combine
  */
  int addto(vector<Matrix> d, int amount) {
    setzero();
    int counter=0;
    for(int m=0;m<amount;m++)
      for(int i=0;i<rows;i++)
        for(int j=0;j<columns;j++) {
          if(d[m].getElement(i,j)) {
                  counter++;
                  v[i][j]+=d[m].getElement(i,j);
                }
        }
      return counter;
  }


  /**
    checks for equality of two Matrix

    \param in the other matrix
  */
  bool eq(Matrix in) {
    for(int i=0;i<rows;i++)
      for(int j=0;j<columns;j++) {
        if(v[i][j]!=in.getElement(i,j))
          return false;
      }
      return true;
  }

};

#endif