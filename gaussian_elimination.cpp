/**
 * File Name: gaussian_elemination.cpp
 * 
 * Author: Justin "Tanner" Broaddus 
 * 
 * Description:
 * 
 * (For Math 5600 Numerical Analysis I Final Project) 
 * (Problem 1)
 * 
 * Representing Gaussian Elemination on the following matrix A:
 *  1 -1  2  1
 *  3  2  1  4
 *  5  8  6  3
 *  4  2  5  3  
 * 
 * With the following column vector b:
 *  10.3041
 *  29.8941
 *  59.1774
 *  40.4211 
 * 
 * ... to find the solution to the equation of the form Ax=b.
 * 
 * NOTICE: This program is using partial pivoting, hence the 
 * row_swap() function defined below.
 *   
 **/

#include <iostream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <future>
#include <functional>


// Helper function to print 4x4 matrix A
void print_matrix(double (A)[4][4]) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      std::cout << A[i][j] << ' ';
    }
    std::cout << '\n';
  }
  std::cout << std::flush;
}

// Helper function to print 4x4 matrix A with column vector b
void print_matrix(double (A)[4][4], double (x)[4]) {
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      std::cout << A[i][j] << ' ';
    }
    std::cout << x[i] << '\n';
  }
  std::cout << std::flush;
}


// Helper function for row swap to perform partial pivoting
void row_swap(double (&A)[4][4], const int row1, const int row2) {
  std::swap(A[row1], A[row2]);
}

// Helper function to find largest value in a column
const int largest_index(double (&A)[4][4], const int column) {
  int index = 0;
  double temp;
  for (int i = 0 + column; i <= 3; i++) {
    if (i == 0) {
      temp = A[i][column];
    } else {
      if (temp < std::abs(A[i][column])) {
        temp = A[i][column];
        index = i;
      } 
    }
  }
  return index;
}

// Function for gaussian elimination with partial pivoting
// Passing matrix A and column vectors x and b by reference
void gaussian_elimination (double (&A)[4][4], double (&x)[4], double (&b)[4]) {
  // For loop to iterate through each column
  for (int j = 0; j < 4; j++) {
    // Finding the index with the largest value in order to perform a row swap for partial pivoting       
    int largest = largest_index(A, j);
    if (largest != j) {   
      // Using an std::async with a deferred launch to make this a synchronous operation.
      // For some reason, a simple call to row_swap is not providing expected behavior
      // and places what should be the top row of the resulting matrix A after 
      // gaussian elimination into the 2nd row. Weird. 
      // This leaves me with strange results such as -inf for some of the values in x. 
      // I believe this is due to an anomaly in the std::swap function 
      // when called on the rows of a multidimensional array!
      // TODO: Find out why!
      std::future<void> f = std::async(std::launch::deferred, row_swap, std::ref(A), j, largest);
      f.wait();
      std::swap(b[j],b[largest]);
    }
    for (int i = 1 + j; i < 4; i++) {
      double multiplier = -(A[i][j]/A[0+j][j]);
      for (int k = j; k < 4; k++) {
        A[i][k] += multiplier*A[0+j][k];
      }
      b[i] = b[i] + multiplier*b[0+j];
    }
  }
  std::cout << "Augmented Matrix A (with Column vector b) after gaussian elimination:\n";
  print_matrix(A,b);
  std::cout << std::endl;
  // At this point, we are done with gaussian elimination and are ready
  // for backwards substitution.
  // For loop to iterate from the bottom up of matrix A for backwards substitution.
  for (int j = 3; j>=0; j--) {
    double numerator = 0;
    for (int i = 3; i >= 0 + j; i--) {
      if (j==i) {
        // calculating the entries in the solution vector x
        x[j] = static_cast<double>((b[j] - numerator)/A[j][i]);
      } else {
        numerator += (x[i] * A[j][i]);
      }
    }
  }
}


// Driver function
int main() {

  // Setting precision to 5 decimal places
  std::cout << std::fixed << std::setprecision(5);

  // Matrix A; represented as a 4x4 2D array of type double 
  double A[4][4] = {{1,-1,2,1}, {3,2,1,4}, {5,8,6,3}, {4,2,5,3}};

  // Column vector b; represented as an 1D array of type double
  double b[4] = {10.3041, 29.8941, 59.1774, 40.4211};

  // Column vector x; represented as an 
  double x[4];  

  std::cout << "\nPerforming gaussian elimination and backwards substitution\n"
    << "to find the solution to the following equation Ax=b.\n\n";
  // Printing Matrix A  
  std::cout << "Matrix A:\n";
  print_matrix(A, b);
  std::cout << "\n";

  // Printing column vector b
  std::cout << "Col Vector b:\n";
  for (int i = 0; i < 4; i++) {
    std::cout << b[i] << '\n';
  }
  std::cout << "\n\n" << std::flush;


  // Performing gaussian elemination with partial pivoting
  gaussian_elimination(A,x,b); 

  // Printing solution vector x
  std::cout << "Values to solution vector x after backwards substitution:\n";
  for (int i = 0; i < 4; i++) {
    std::cout << "x" << i+1 << " = " << x[i] << '\n';
  }

  std::cout << std::endl;

  return 0;
}