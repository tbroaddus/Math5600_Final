/**
 * File Name: mod_newton.cpp
 * 
 * Author: Justin "Tanner" Broaddus
 * 
 * Description:
 * 
 * (For Math 5600 Numerical Analysis I Final Project)
 * (Problem 2)
 * 
 * Representing modified newton's method for approximating the 
 * root of a function f(x) that has a multiplicity greater than 1.
 * 
 * f(x) = x^3 - 3 * x^2 + 4
 * f'(x) = 3x^2 - 6x
 * 
 * x_n = x_n-1 - (m * (f(x_n-1)/f'(x_n-1)))
 * 
 **/

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>

// Helper function to represent function f(x)
// f(x) = x^3 - 3 * x^2 + 4
double f (const double x) {
  return (pow(x,3) - (3.0 * pow(x,2)) + 4.0);
}

// Helper function to represent function f'(x)
//f'(x) = 3x^2 - 6x
double f_prime (const double x) {
  return ((3.0 * pow(x,2)) - (6.0 * x));
}


// Function to perform modified newton method for x sub i
// x_n = x_n-1 - (m * (f(x_n-1)/f'(x_n-1)))
double mod_newton(const double m, const double x) {
  return x - (m * (f(x)/f_prime(x)));
}

// Driver function
// Params: m,n
int main(int argc, char** argv) {

  if (argc != 3) {
    std::cerr << "Invalid number of arguments, please pass in value of m and n (number of iterations) respectively!" << std::endl;
    return -1;
  }
  // Value of m from command line args
  const double m = static_cast<double>(atof(argv[1]));
  // Value of n (max iter) from command line args
  const int n = static_cast<int>(atoi(argv[2]));
  // Value of x; x_0 = 3 (initial guess)
  double x = 3.0;

  std::cout << "\nParameters:\n";
  std::cout << "m = " << static_cast<int>(m) << '\n';
  std::cout << "n = " << n << '\n';
  std::cout << "x_b 0 = " << static_cast<int>(x) << "\n\n\n";
  if (m == 1)
    std::cout << "Step n\t\tx_n\t\t\t|f(x_c)|\t\t|f(x_n)/f(x_n-1)|\n\n";
  else 
    std::cout << "Step n\t\tx_n\t\t\t|f(x_c)|\t\t|f(x_n)/f(x_n-1)^2|\n\n";

  std::cout << std::fixed << std::setprecision(8);
  // Printing initial step 0
  std::cout << 0 << ":\t\t" << std::fixed << std::setprecision(8)
      << x << "\t\t" << std::fixed << std::setprecision(12) << std::abs(f(x)) 
       << '\n';

  // Performing n iterations of modified newtons method for quadratic convergence
  for (int i = 0; i < n; i++) {
    double x_prev = x;
    x = mod_newton(m,x);
    // backwards error = |f(x sub c)|
    // where x sub c = approximation of true value of root
    double backwards_error = std::abs(f(x));
    std::cout << i+1 << ":\t\t" << std::fixed << std::setprecision(8)
      << x << "\t\t" << std::fixed << std::setprecision(12) << backwards_error 
      << "\t\t" << std::fixed << std::setprecision(8) 
      << ( m == 1 ? abs(f(x)/f(x_prev)) : abs(f(x)/(pow(f(x_prev),2)))) << '\n';
    if (backwards_error < 1.0e-12) 
      break;
  }

  std::cout << std::endl;

  return 0;
}