/**
 * File Name: third_derivative.cpp
 * 
 * Author: Justin "Tanner" Broaddus
 * 
 * Description:
 *
 * (For Math 5600 Numerical Analysis I Final Project)
 * (Problem 3)
 * 
 * Using a second-order-accurate centered difference formula to 
 * approximate the third derivative of a function f(x)
 * where f(x) = sin(2x) and the actual value of f'''(x)
 * is acquired from f'''(x) = -8cos(2x).
 * 
 * M_PI is a macro that is an approximate value of pi.
 * 
 **/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>

// Representing function f(x) = sin(2x)
double f(const double x) {
  return sin(2*x);
}

// Representing the approximate f'''(x)
double third_dir_f(const double x, const double h) {
  return (f(x + 2 * h) - (2 * f(x + h)) + (2 * f(x - h)) - f(x - 2 * h))/(2 * pow(h,3));
}

int main() {

  // f'''(x) = -8cos(2x) when f(x) = sin(2x)
  // Calculating actual value of f'''(x) when x = pi/8
  double actual = -8 * cos(2 * (M_PI/8));

  std::cout << std::fixed << std::setprecision(8);
  std::cout << "\n\nApproximating f'''(x) for f(x) = (sin(2x)) when x = pi/8\nActual value = " 
    << actual << "\n\n\n";
  std::cout << "h\t\tf'''(x)\t\terror\t\terror/h^2\n\n";

  // For loop to iterate through each h value
  for (int i = 1; i <= 5; i++) {
    // Calculating the approximate value
    double approximate =  third_dir_f(M_PI/8, 1/pow(4,i));
    // Calculating the absolute error
    double error = std::abs(actual - approximate);
    std::cout <<  "1/4^" << i << "\t\t" 
      << approximate 
      << "\t" << error 
      << "\t" << error/pow(1/pow(4,i),2) << '\n'; 
  }

  std::cout << std::endl;

  return 0;
}