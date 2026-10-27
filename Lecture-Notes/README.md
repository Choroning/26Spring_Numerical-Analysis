# Lecture Notes

[Repository Guide](../README.md) | [English](README.md) | [한국어](README.ko.md)

> **Last Updated:** 2026-10-27

## Part 1 Modeling, Computers, and Error Analysis

| Unit | Topics | Documents |
|:-----|:-----|:-----|
| [Chapter 1 Mathematical Modeling and Engineering Problem Solving](Ch01_Mathematical-Modeling/) | Mathematical models, Newton's second law and the bungee jumper example, Euler's method, types of engineering problems solved numerically, NumPy arrays and matrix operations, Matplotlib visualization, views and copies in NumPy | [Lecture notes](Ch01_Mathematical-Modeling/Concepts_Lecture.md), [Lab notes](Ch01_Mathematical-Modeling/Concepts_Lab.md) |
| [Chapter 2 Python Fundamentals and Programming](Ch02_Python-Fundamentals-and-Programming/) | Programming environment, Python basics and NumPy and SciPy functions, control structures and functions, the bungee jumper example in Python, variable scope and the LEGB rule, structured programming, exception handling, decorators, object-oriented basics, type hints, generators | [Lecture notes](Ch02_Python-Fundamentals-and-Programming/Concepts_Lecture.md), [Lab notes](Ch02_Python-Fundamentals-and-Programming/Concepts_Lab.md) |
| [Chapter 3 Roundoff Errors](Ch03_Roundoff-Errors/) | Significant figures, accuracy and precision, error definitions and the Scarborough criterion, round-off errors and problems with floating-point arithmetic, error metrics, truncation error and the Taylor series | [Lecture notes](Ch03_Roundoff-Errors/Concepts_Lecture.md), [Lab notes](Ch03_Roundoff-Errors/Concepts_Lab.md) |
| [Chapter 4 Truncation Errors](Ch04_Truncation-Errors/) | Taylor series expansion, order of approximation, numerical differentiation using the Taylor series, total numerical error, stability and condition number, finite differences in Python | [Lecture notes](Ch04_Truncation-Errors/Concepts_Lecture.md) |

## Part 2 Roots and Optimization

| Unit | Topics | Documents |
|:-----|:-----|:-----|
| [Chapter 5 Roots: Bracketing Methods](Ch05_Roots-Bracketing-Methods/) | Roots problems and graphical methods, bracketing versus open methods, incremental search, the bisection method and its convergence, the false position and modified false position methods | [Lecture notes](Ch05_Roots-Bracketing-Methods/Concepts_Lecture.md) |
| [Chapter 6 Roots: Open Methods](Ch06_Roots-Open-Methods/) | Fixed-point iteration, the Wegstein method, the Newton-Raphson method, the secant method, inverse quadratic interpolation, Brent's method | [Lecture notes](Ch06_Roots-Open-Methods/Concepts_Lecture.md) |
| [Chapter 7 Optimization](Ch07_Optimization/) | Local and global optimality, existence of an optimal solution, one-dimensional optimization, golden-section search, parabolic interpolation, level sets and the gradient, convex problems, the steepest descent method | [Lecture notes](Ch07_Optimization/Concepts_Lecture.md) |

## Part 3 Linear Systems

| Unit | Topics | Documents |
|:-----|:-----|:-----|
| [Chapter 8 Solving Linear Algebraic Equations](Ch08_Linear-Algebraic-Equations/) | Graphical method for small systems, determinants and Cramer's rule, naive Gauss elimination and its computational cost, pitfalls of naive Gauss elimination, a spring-mass system application | [Lecture notes](Ch08_Linear-Algebraic-Equations/Concepts_Lecture.md), [Lab notes](Ch08_Linear-Algebraic-Equations/Concepts_Lab.md) |
| [Chapter 9 Gauss Elimination: Improved Methods](Ch09_Gauss-Elimination/) | Partial pivoting and scaled partial pivoting, the Gauss-Jordan method, tridiagonal systems and the Thomas algorithm, condition number, heat conduction | [Lecture notes](Ch09_Gauss-Elimination/Concepts_Lecture.md), [Lab notes](Ch09_Gauss-Elimination/Concepts_Lab.md) |
| [Chapter 10 LU Factorization](Ch10_LU-Factorization/) | Gauss elimination as LU factorization, forward and back substitution, LU factorization with pivoting, Cholesky factorization, Doolittle's algorithm, multiple right-hand sides and a thermal resistor network | [Lecture notes](Ch10_LU-Factorization/Concepts_Lecture.md), [Lab notes](Ch10_LU-Factorization/Concepts_Lab.md) |
| [Chapter 11 Matrix Inverse and Condition Number](Ch11_Matrix-Inverse-and-Condition/) | The matrix inverse and its computation via LU factorization, physical interpretation with a mass-spring system, detecting ill-conditioning, vector and matrix norms, the matrix condition number and the error bound theorem, the Hilbert matrix | [Lecture notes](Ch11_Matrix-Inverse-and-Condition/Concepts_Lecture.md), [Lab notes](Ch11_Matrix-Inverse-and-Condition/Concepts_Lab.md) |
| [Chapter 12 Iterative Methods for Linear Systems](Ch12_Iterative-Methods/) | The Gauss-Seidel method, diagonal dominance and convergence, Jacobi iteration, relaxation, preconditioning, the power method for eigenvalues, the inverse and shifted power methods | [Lecture notes](Ch12_Iterative-Methods/Concepts_Lecture.md), [Lab notes](Ch12_Iterative-Methods/Concepts_Lab.md) |
| [Chapter 13 Eigenvalue Methods: Symmetric Matrices](Ch13_Eigenvalues-Symmetric-Matrices/) | Computing eigenvalues and eigenvectors, eigenvalues and the stability of ordinary differential equations, coupled second-order ODE systems, mass-spring systems, the power method, spectral decomposition of symmetric matrices, vibration analysis | [Lecture notes](Ch13_Eigenvalues-Symmetric-Matrices/Concepts_Lecture.md), [Lab notes](Ch13_Eigenvalues-Symmetric-Matrices/Concepts_Lab.md) |

## Part 4 Curve Fitting

| Unit | Topics | Documents |
|:-----|:-----|:-----|
| [Chapter 14 Linear Regression](Ch14_Linear-Regression/) | Statistics review and the normal distribution, straight-line least-squares regression, quantifying the quality of the fit | [Lecture notes](Ch14_Linear-Regression/Concepts_Lecture.md) |
| [Chapter 15 Nonlinear Regression](Ch15_Nonlinear-Regression/) | Standard error and coefficient of determination, linearization of nonlinear relationships, polynomial regression, multiple linear regression, general linear least-squares regression, Gauss-Newton nonlinear regression | [Lecture notes](Ch15_Nonlinear-Regression/Concepts_Lecture.md) |
| [Chapter 16 Fourier Analysis](Ch16_Fourier-Analysis/) | Even and odd functions, Fourier series and the Fourier integral, the Fourier transform, least-squares fit of a sinusoid, frequency and time domains, the discrete Fourier transform, the Nyquist frequency, the fast Fourier transform, the power spectrum | [Lecture notes](Ch16_Fourier-Analysis/Concepts_Lecture.md) |
| [Chapter 17 Polynomial Interpolation](Ch17_Interpolation/) | Polynomial interpolation, the Newton interpolating polynomial, the Lagrange interpolating polynomial, inverse interpolation, extrapolation and oscillations | [Lecture notes](Ch17_Interpolation/Concepts_Lecture.md) |
| [Chapter 18 Splines and Piecewise Interpolation](Ch18_Splines/) | Continuity conditions, linear splines, quadratic splines, cubic splines, bilinear interpolation | [Lecture notes](Ch18_Splines/Concepts_Lecture.md) |

## Part 5 Integration and Differentiation

| Unit | Topics | Documents |
|:-----|:-----|:-----|
| [Chapter 19 Numerical Integration Formulas](Ch19_Numerical-Integration-Formulas/) | The trapezoidal rule, Simpson's 1/3 and 3/8 rules, composite integration formulas, combining integration rules | [Lecture notes](Ch19_Numerical-Integration-Formulas/Concepts_Lecture.md) |
| [Chapter 20 Numerical Integration of Equations](Ch20_Numerical-Integration-of-Equations/) | Romberg integration, Gauss quadrature | [Lecture notes](Ch20_Numerical-Integration-of-Equations/Concepts_Lecture.md) |
| [Chapter 21 Numerical Differentiation](Ch21_Numerical-Differentiation/) | Finite-difference approximations of first and second derivatives, higher-order forward differences, Richardson extrapolation, derivatives of unequally spaced data, numerical partial derivatives, the method of undetermined coefficients, differentiable programming | [Lecture notes](Ch21_Numerical-Differentiation/Concepts_Lecture.md) |

## Part 6 Ordinary Differential Equations

| Unit | Topics | Documents |
|:-----|:-----|:-----|
| [Chapter 22 Initial-Value Problems for ODEs](Ch22_Initial-Value-Problems/) | Euler's method, local and global truncation errors, second-order Runge-Kutta methods, the classical fourth-order Runge-Kutta method, systems of ODEs, numerical stability | [Lecture notes](Ch22_Initial-Value-Problems/Concepts_Lecture.md) |
| [Chapter 23 Stiff ODEs and Adaptive Methods](Ch23_Stiff-ODEs/) | Adaptive Runge-Kutta methods, the embedded RK method RKF45 and step size control, multistep methods, stiff systems | [Lecture notes](Ch23_Stiff-ODEs/Concepts_Lecture.md) |
| [Chapter 24 Boundary-Value Problems](Ch24_Boundary-Value-Problems/) | Boundary-value problems and the analytical solution, the shooting method, finite-difference methods | [Lecture notes](Ch24_Boundary-Value-Problems/Concepts_Lecture.md) |
