# HW4 Assignment Explanation — Iterative Solvers

[Repository Guide](../../README.md) | [English](Assignment-Explanation.md) | [한국어](Assignment-Explanation.ko.md)

> **Last Updated:** 2026-11-16
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 12

> **Learning Objectives:**
> 1. Split a matrix into its diagonal, strictly lower triangular, and strictly upper triangular parts and relate them to the matrix forms of iterative methods
> 2. Implement Jacobi iteration as a single vector operation and verify it with an iteration table
> 3. Implement the Gauss-Seidel method with in-place updates and compare its convergence rate with Jacobi iteration
> 4. Explain the difference in iteration counts between the two methods with the spectral radius of the iteration matrix
> 5. Implement the diagonal dominance check and successive over-relaxation (SOR) and explain how the iteration count changes with the weighting factor and the grid size

---

## Table of Contents

- [1. Assignment Overview](#1-assignment-overview)
- [2. Requirements](#2-requirements)
- [3. Problem 1: Matrix Splitting](#3-problem-1-matrix-splitting)
  - [3.1 Approach](#31-approach)
  - [3.2 Code Written](#32-code-written)
  - [3.3 Results and Verification](#33-results-and-verification)
- [4. Problem 2: Jacobi Iteration](#4-problem-2-jacobi-iteration)
  - [4.1 Approach](#41-approach)
  - [4.2 Code Written](#42-code-written)
  - [4.3 Results and Verification](#43-results-and-verification)
- [5. Problem 3: The Gauss-Seidel Method](#5-problem-3-the-gauss-seidel-method)
  - [5.1 Approach](#51-approach)
  - [5.2 Code Written](#52-code-written)
  - [5.3 Results and Verification](#53-results-and-verification)
- [6. Problem 4: Diagonal Dominance Check](#6-problem-4-diagonal-dominance-check)
  - [6.1 Approach](#61-approach)
  - [6.2 Code Written](#62-code-written)
  - [6.3 Results and Verification](#63-results-and-verification)
- [7. Problem 5: Successive Over-Relaxation](#7-problem-5-successive-over-relaxation)
  - [7.1 Approach](#71-approach)
  - [7.2 Code Written](#72-code-written)
  - [7.3 Results and Verification](#73-results-and-verification)
- [8. Problem 6: Grid SOR for 2D Heat Conduction](#8-problem-6-grid-sor-for-2d-heat-conduction)
  - [8.1 Approach](#81-approach)
  - [8.2 Code Written](#82-code-written)
  - [8.3 Results and Verification](#83-results-and-verification)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Assignment Overview

> **Try It First:** This document contains the answers and solutions. Solve the problems in the [assignment handout](Original/Assignment.md) yourself before reading on.

The assignment is to fill in the six blanks of the lab notebook `Lab05-IterativeSolver.ipynb` distributed by the instructor and run it. The notebook has six parts: matrix splitting and the fixed-point framework, Jacobi iteration, the Gauss-Seidel method, convergence analysis, successive over-relaxation (SOR), and an application to 2D steady-state heat conduction, and each part contains one blank. The problem statements are in the [assignment handout](Original/Assignment.md), the distributed blank notebook is [Original/Lab05-IterativeSolver.ipynb](Original/Lab05-IterativeSolver.ipynb), and the notebook with the blanks filled in and run is [Solution/Lab05-IterativeSolver_Completed.ipynb](Solution/Lab05-IterativeSolver_Completed.ipynb).

| Problem | Notebook section | Blanks to fill |
|:-----|:-----|:-----|
| Problem 1 | 1.3 Fixed-point framework | `L` and `U` in the function `split_DLU` |
| Problem 2 | 2.2 Jacobi iteration | `x_new` in the function `jacobi` |
| Problem 3 | 3.1 Gauss-Seidel method | `sigma` and `x[i]` in the function `gauss_seidel` |
| Problem 4 | 4.1 Diagonal dominance | `offdiag_sum` and `margins` in the function `is_diag_dominant` |
| Problem 5 | 5.1 Successive over-relaxation (SOR) | `x[i]` in the function `sor` |
| Problem 6 | 6.2 Gauss-Seidel and SOR on the grid | `T[i, j]` in the function `solve_laplace_2d` |

*Table 1. The six blanks that make up the assignment*

Jacobi iteration, the Gauss-Seidel method, matrix splitting, diagonal dominance, relaxation, and the spectral radius are covered in the [Chapter 12 lecture note](../../Lecture-Notes/Ch12_Iterative-Methods/Concepts_Lecture.md). This document repeats only what the solutions need.

---

<br>

## 2. Requirements

The completed notebook must satisfy the following requirements. The rest of each cell stays as distributed, and only the blanks are filled in.

1. **Problem 1:** `L` must keep only the entries of $`A`$ below the diagonal and `U` only those above it, and the test in the same cell must give the zero matrix for $`D + L + U - A`$.
2. **Problem 2:** `x_new` must be the new vector of a Jacobi iteration computed from the previous iterate `x` only. The error list returned by the function is used as the relative approximate error column of the iteration table that the test cell computes separately.
3. **Problem 3:** `sigma` must be the sum of the terms of equation $`i`$ other than the diagonal term, with the components for $`j < i`$ taken from the values updated in this iteration. `x[i]` is updated in place right away.
4. **Problem 4:** `offdiag_sum` must be the sum of the absolute off-diagonal entries of each row and `margins` the absolute diagonal entry minus that sum, with all rows computed at once by array operations.
5. **Problem 5:** `x[i]` must blend the Gauss-Seidel value `x_gs_i` and the old value with the weighting factor `omega`, and with `omega = 1` it must equal the Gauss-Seidel method.
6. **Problem 6:** For `method='sor'`, `T[i, j]` must blend the grid's Gauss-Seidel value `T_gs` and the node's current value with `omega`.

---

<br>

## 3. Problem 1: Matrix Splitting

### 3.1 Approach

`np.tril(A, k)` keeps the $`k`$-th diagonal and everything below it, and `np.triu(A, k)` keeps the $`k`$-th diagonal and everything above it. The main diagonal is $`k = 0`$, so `np.tril(A, -1)`, which keeps entries from just below the diagonal, is the strictly lower triangular part $`L`$, and `np.triu(A, 1)`, which keeps entries from just above the diagonal, is the strictly upper triangular part $`U`$. The distributed code already obtains the diagonal part with `np.diag(np.diag(A))`: the inner `np.diag` extracts the vector of diagonal entries, and the outer `np.diag` builds a diagonal matrix from that vector. Matrix splitting is covered in Section 4.1 of the [Chapter 12 lecture note](../../Lecture-Notes/Ch12_Iterative-Methods/Concepts_Lecture.md).

> **Caution:** The $`L`$ and $`U`$ here are different objects from the $`L`$ and $`U`$ of the LU factorization in HW3. The LU factorization is the product $`A = LU`$ with 1s on the diagonal of $`L`$, whereas matrix splitting is the sum $`A = D + L + U`$ with zeros on the diagonals of $`L`$ and $`U`$.

### 3.2 Code Written

```python
    L = np.tril(A, -1)                               # (strictly lower)
    U = np.triu(A,  1)                               # (strictly upper)
```

### 3.3 Results and Verification

```text
D =
 [[ 3.  0.  0.]
 [ 0.  7.  0.]
 [ 0.  0. 10.]]
L =
 [[ 0.   0.   0. ]
 [ 0.1  0.   0. ]
 [ 0.3 -0.2  0. ]]
U =
 [[ 0.  -0.1 -0.2]
 [ 0.   0.  -0.3]
 [ 0.   0.   0. ]]

Sanity check  D + L + U - A =
 [[0. 0. 0.]
 [0. 0. 0.]
 [0. 0. 0.]]
```

The output matches the output in the distributed notebook. Only $`a_{21} = 0.1`$, $`a_{31} = 0.3`$, and $`a_{32} = -0.2`$ remain in $`L`$, only $`a_{12} = -0.1`$, $`a_{13} = -0.2`$, and $`a_{23} = -0.3`$ remain in $`U`$, and the three parts add up exactly to $`A`$. The function `iteration_matrices` of Section 4.2 of the notebook uses this splitting directly to build the iteration matrices $`T_J = -D^{-1}(L+U)`$ and $`T_{GS} = -(D+L)^{-1}U`$.

---

<br>

## 4. Problem 2: Jacobi Iteration

### 4.1 Approach

The function prepares `d = np.diag(A)` (the vector of diagonal entries) and `R = A - np.diag(d)` (the part without the diagonal, $`L + U`$) in advance. In the matrix form of Jacobi iteration, $`\mathbf{x}^{(k+1)} = D^{-1}(\mathbf{b} - (L+U)\mathbf{x}^{(k)})`$, multiplying by $`D^{-1}`$ is the same as dividing each component by the corresponding diagonal entry, so one iteration can be written as `(b - R @ x) / d`. The $`i`$-th component of `R @ x` is $`\sum_{j \neq i} a_{ij} x_j`$, and `/ d` is componentwise division.

Every component is computed from the previous values `x` only, so the result goes into a new array `x_new`. The code that follows computes the relative approximate error from the difference between `x_new` and `x` before replacing `x` with `x = x_new`. Jacobi iteration is covered in Section 5 of the [Chapter 12 lecture note](../../Lecture-Notes/Ch12_Iterative-Methods/Concepts_Lecture.md).

### 4.2 Code Written

```python
        x_new = (b - R @ x) / d
```

### 4.3 Results and Verification

```text
Jacobi converged in 8 iterations
Solution : [ 3.  -2.5  7. ]
Reference: [ 3.  -2.5  7. ]
Max error: 2.24e-10

 Iter            x1            x2            x3        ea (%)
------------------------------------------------------------
    1      2.616667     -2.757143      7.140000    100.000000
    2      3.000762     -2.488524      7.006357      6.078749
    3      3.000806     -2.499738      7.000207      0.159563
    4      3.000022     -2.500003      6.999981      0.010697
    5      2.999999     -2.500001      6.999999      0.000375
    6      3.000000     -2.500000      7.000000      0.000023
    7      3.000000     -2.500000      7.000000      0.000001
    8      3.000000     -2.500000      7.000000      0.000000
```

The first iteration starts from $`\mathbf{x}^{(0)} = \mathbf{0}`$, so $`x_1 = 7.85/3 = 2.616667`$, $`x_2 = -19.3/7 = -2.757143`$, and $`x_3 = 71.4/10 = 7.14`$, matching the output. The $`x`$ values in the iteration table are computed separately by the test cell with the same formula, while the error column is the list returned by the function, so the agreement of the two columns means the function performs the same update. The error of the first iteration is 100% because the previous iterate is the zero vector, so $`\lVert \mathbf{x}^{(1)} - \mathbf{0} \rVert / \lVert \mathbf{x}^{(1)} \rVert = 1`$.

The tolerance `tol=1e-6` is interpreted as a percentage, so the loop stops at the 8th iteration, where $`\varepsilon_a < 10^{-6}\%`$. Section 4.2 of the notebook computes the spectral radius of the Jacobi iteration matrix of this matrix as $`\rho(T_J) = 0.0482`$, which means each iteration gains about $`-\log_{10} 0.0482 \approx 1.3`$ digits. This matches the error column shrinking by a factor of 0.03-0.07 per iteration.

---

<br>

## 5. Problem 3: The Gauss-Seidel Method

### 5.1 Approach

When component $`i`$ is updated, the sum of the terms other than the diagonal term splits into two inner products, one to the left of the diagonal and one to the right.

$$
\sigma_i = \sum_{j < i} a_{ij} x_j + \sum_{j > i} a_{ij} x_j
$$

Because the array `x` is updated in place, when component $`i`$ is computed, `x[:i]` holds the values already updated in this iteration and `x[i+1:]` holds the not-yet-updated values of the previous iteration. Therefore `A[i, :i] @ x[:i] + A[i, i+1:] @ x[i+1:]` implements the two sums of the Gauss-Seidel formula exactly, and the right part could equally be written with `x_old[i+1:]`. Then `x[i] = (b[i] - sigma) / A[i, i]` replaces the entry in place, so the new value is used for the next component. The Gauss-Seidel method is covered in Sections 2 and 4 of the [Chapter 12 lecture note](../../Lecture-Notes/Ch12_Iterative-Methods/Concepts_Lecture.md).

> **Caution:** Using `x_old[:i]` instead of `x[:i]` computes every component from old values, which turns the method into Jacobi iteration. The difference between the two methods lies not in the formula but in which values are used and when.

### 5.2 Code Written

```python
            sigma = A[i, :i] @ x[:i] + A[i, i+1:] @ x[i+1:]
            x[i]  = (b[i] - sigma) / A[i, i]
```

### 5.3 Results and Verification

```text
Jacobi       converged in   8 iterations
Gauss-Seidel converged in   6 iterations
GS speedup (iterations): 1.3x

GS solution: [ 3.  -2.5  7. ]
True value : [ 3.  -2.5  7. ]
```

The test cell does not print an iteration table, so the iterates of the same function, printed separately, are as follows.

| Iteration | $`x_1`$ | $`x_2`$ | $`x_3`$ | $`\varepsilon_a`$ (%) |
|:-----|:-----|:-----|:-----|:-----|
| 1 | 2.616667 | -2.794524 | 7.005610 | 100.000000 |
| 2 | 2.990557 | -2.499625 | 7.000291 | 5.943698 |
| 3 | 3.000032 | -2.499988 | 6.999999 | 0.118354 |
| 4 | 3.000000 | -2.500000 | 7.000000 | 0.000421 |
| 5 | 3.000000 | -2.500000 | 7.000000 | 0.000004 |
| 6 | 3.000000 | -2.500000 | 7.000000 | 0.000000 |

*Table 2. Iterates and relative approximate errors of the Gauss-Seidel method*

The first-iteration values $`2.6167`$, $`-2.7945`$, and $`7.0056`$ equal those computed by hand in Section 3.1 of the notebook. From the second component on, the updated $`x_1`$ is used, so they are closer to the true values than the first Jacobi iterate ($`-2.757143`$, $`7.14`$).

The spectral radius explains the difference in convergence rate. According to the output of Section 4.2 of the notebook, this matrix has $`\rho(T_J) = 0.0482`$ and $`\rho(T_{GS}) = 0.0065`$. The number of iterations needed for a given accuracy is roughly proportional to $`1/(-\log\rho)`$, so Gauss-Seidel needs about $`\log 0.0482 / \log 0.0065 \approx 0.6`$ times as many as Jacobi, and the actual counts were 8 and 6.

For the heat network of Section 3.2 of the notebook, Gauss-Seidel also converges in 13 iterations to the tolerance $`10^{-8}\%`$, and the node temperatures 54.464286°C, 17.857143°C, and 16.964286°C equal the exact values $`1525/28`$, $`125/7`$, and $`475/28`$. In the comparison of Section 5 of the notebook, Jacobi iteration took 23 iterations on the same system. This matrix is tridiagonal, so $`\rho(T_{GS}) = \rho(T_J)^2`$ $`(0.125 = 0.3536^2)`$ holds, and Gauss-Seidel converges in about half as many iterations as Jacobi.

---

<br>

## 6. Problem 4: Diagonal Dominance Check

### 6.1 Approach

The row sums of `np.abs(A)`, `np.sum(np.abs(A), axis=1)`, equal $`\sum_j |a_{ij}|`$, so subtracting the absolute diagonal entries `diag_abs` gives the sum of the absolute off-diagonal entries $`\sum_{j \neq i} |a_{ij}|`$. The margin is $`|a_{ii}| - \sum_{j \neq i} |a_{ij}|`$, and the function decides by whether all margins are positive when `strict=True` or all are nonnegative when `strict=False`. Because these are array operations, the values for all rows are computed at once. Diagonal dominance is covered in Sections 3 and 8.4 of the [Chapter 12 lecture note](../../Lecture-Notes/Ch12_Iterative-Methods/Concepts_Lecture.md).

### 6.2 Code Written

```python
    offdiag_sum = np.sum(np.abs(A), axis=1) - diag_abs  # row sums of |A| minus the diagonal
    margins = diag_abs - offdiag_sum
```

### 6.3 Results and Verification

```text
A_good — strictly diagonally dominant: True
  Row margins: [2.5 3.8 5.2]

A_bad  — strictly diagonally dominant: False
  Row margins: [-4.  0. -1.]

A_ex   — strictly diagonally dominant: True,  margins: [2.7 6.6 9.5]
```

The output matches the output in the distributed notebook. By hand, the margins of `A_ex` are $`3 - 0.3 = 2.7`$, $`7 - 0.4 = 6.6`$, and $`10 - 0.5 = 9.5`$, those of `A_good` are $`4 - 1.5 = 2.5`$, $`5 - 1.2 = 3.8`$, and $`6 - 0.8 = 5.2`$, and those of `A_bad` are $`1 - 5 = -4`$, $`2 - 2 = 0`$, and $`1 - 2 = -1`$. Row 1 of `A_bad` has margin 0, so it is not dominant under the strict test.

Comparing with the spectral radii of Section 4.2 of the notebook shows what the test means. `A_good` has $`\rho(T_J) = 0.2`$ and $`\rho(T_{GS}) = 0.0817`$, so both methods converge, while `A_bad` has $`\rho(T_J) = 2.4256`$ and $`\rho(T_{GS}) = 7.0707`$, so both diverge, and in the following cell the norm of the Jacobi iterates grows to 4404.97 in 10 iterations. Diagonal dominance, however, is only a sufficient condition. The 2D Laplace matrix of Problem 6 has $`|-4| = 1 + 1 + 1 + 1`$ in its interior rows and is therefore not strictly diagonally dominant, yet the Gauss-Seidel method converges.

---

<br>

## 7. Problem 5: Successive Over-Relaxation

### 7.1 Approach

The SOR formula can be written as a weighted combination of the Gauss-Seidel value $`\tilde{x}_i`$ and the old value.

$$
x_i^{(k+1)} = (1 - \omega)\, x_i^{(k)} + \omega\, \tilde{x}_i
$$

The distributed code already computes the Gauss-Seidel value `x_gs_i`, so the blank is this formula written out, `(1 - omega) * x_old[i] + omega * x_gs_i`. When component $`i`$ is updated, `x_old[i]` and `x[i]` hold the same value (not yet changed in this iteration), so either can be used. With $`\omega = 1`$ the update becomes `x_gs_i` itself, which is the Gauss-Seidel method. Relaxation is covered in Section 6 of the [Chapter 12 lecture note](../../Lecture-Notes/Ch12_Iterative-Methods/Concepts_Lecture.md).

### 7.2 Code Written

```python
            x[i] = (1 - omega) * x_old[i] + omega * x_gs_i
```

### 7.3 Results and Verification

```text
Method            Iterations
----------------------------
Jacobi                    23
GS (ω=1.0)                13
SOR ω=1.1                 11
SOR ω=1.3                 21
SOR ω=1.5                 35
SOR ω=1.7                 67
```

```text
Empirical optimal omega : 1.045  (10 iters)
Theoretical optimal omega: 1.033
GS iterations            : 13
Best SOR iterations      : 10
SOR speedup over GS      : 1.3x
```

Running the same function with $`\omega = 1`$ converges in 13 iterations, with the same iteration count and solution as the Gauss-Seidel method. The theoretical optimum is $`2/(1 + \sqrt{1 - 0.125}) = 1.033`$, obtained by putting $`\rho(T_J) = 0.3536`$ into the notebook's formula, and it is close to the empirically found optimum 1.045. Solving with $`\omega = 1.033`$ also takes 10 iterations. This system has a small $`\rho(T_J)`$ and converges quickly with the Gauss-Seidel method alone, so the optimal $`\omega`$ is close to 1, and values of $`\omega`$ of 1.3 and above overshoot and increase the iteration count (21, 35, and 67 iterations).

---

<br>

## 8. Problem 6: Grid SOR for 2D Heat Conduction

### 8.1 Approach

The grid is swept in place, so when `T_gs` is computed at node $`(i, j)`$, the upper and left neighbors `T[i-1, j]` and `T[i, j-1]` already hold values updated in this iteration, while the lower and right neighbors `T[i+1, j]` and `T[i, j+1]` and the node's own `T[i, j]` hold values of the previous iteration. `T_gs` is therefore the Gauss-Seidel value, and `T[i, j]` plays the role of the old value, so the same weighted combination as in Problem 5, `(1 - omega) * T[i, j] + omega * T_gs`, is the SOR update. Without assembling a matrix, each equation of the 5-point stencil is itself the update of one node.

### 8.2 Code Written

```python
                    T[i, j] = (1 - omega) * T[i, j] + omega * T_gs
```

### 8.3 Results and Verification

```text
2D Laplace (n=30, 900 unknowns):
  Gauss-Seidel : 1267 iterations
  SOR (ω=1.7)  :  243 iterations
  SOR speedup  : 5.2x
  Max |GS - SOR|: 8.22e-05
```

The next cell in Section 6 of the notebook compares the iteration counts for different grid sizes with the tolerance $`10^{-5}`$.

```text
   n    unknowns    GS iters   SOR iters   SOR speedup
------------------------------------------------------
  10         100         157          49           3.2x
  20         400         514          91           5.6x
  30         900        1043         207           5.0x
  50        2500        2553         530           4.8x
```

For this model problem the spectral radius of the Jacobi iteration matrix is $`\rho(T_J) = \cos(\pi h)`$ $`(h = 1/(n+1))`$. Putting it into the notebook's optimal formula gives $`\omega_{\text{opt}} = 2/(1 + \sin(\pi h))`$, about 1.816 for $`n = 30`$. Running the same function with this value converges in 110 iterations, fewer than the 243 iterations at $`\omega = 1.7`$ and the 1267 iterations of Gauss-Seidel.

- **Growth with grid size:** $`\rho(T_{GS}) = \cos^2(\pi h)`$ approaches 1 as $`0.921, 0.978, 0.990, 0.996`$ for $`n = 10, 20, 30, 50`$, so the Gauss-Seidel iteration count grows quickly from 157 to 2553.
- **Fixed weighting factor:** The optimum $`2/(1 + \sin(\pi h))`$ is $`1.560, 1.741, 1.816, 1.884`$ for the same $`n`$. The fixed $`\omega = 1.7`$ is closest to the optimum for $`n = 20`$, where the speedup is largest at 5.6 times.
- **Difference between the GS and SOR results:** Both methods stop when a single update becomes smaller than $`10^{-6}`$, but when convergence is slow the remaining error is much larger than the update. The error is roughly $`\rho/(1-\rho)`$ times the update, about 97 times for $`\rho(T_{GS}) \approx 0.990`$. Compared with a solution computed to a very small tolerance, the error of the Gauss-Seidel result is $`9.6 \times 10^{-5}`$ and that of the SOR ($`\omega = 1.7`$) result is $`1.4 \times 10^{-5}`$, so the printed difference of $`8.22 \times 10^{-5}`$ is mostly the error remaining on the Gauss-Seidel side.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Matrix splitting | `np.tril(A, -1)` and `np.triu(A, 1)` give the strictly lower and strictly upper triangular parts, and the $`L`$ and $`U`$ of $`A = D + L + U`$ are different objects from the factors of the LU factorization. |
| Jacobi iteration | One iteration divides $`(\mathbf{b} - (L+U)\mathbf{x})`$ componentwise by the diagonal entries, and because every component uses old values the result goes into a new array. |
| Gauss-Seidel method | In-place updates make the components for $`j < i`$ use the new values of this iteration, and using `x_old[:i]` turns it into Jacobi iteration. |
| Convergence test | The relative approximate error is computed as a percentage and the loop stops when it falls below the tolerance, and the first iteration always gives 100% because the previous iterate is the zero vector. |
| Spectral radius and iteration count | The number of iterations needed is roughly proportional to $`1/(-\log\rho(T))`$, and for a tridiagonal matrix $`\rho(T_{GS}) = \rho(T_J)^2`$, so Gauss-Seidel converges in about half as many iterations. |
| Diagonal dominance check | For each row the margin is the absolute diagonal entry minus the sum of the absolute off-diagonal entries, and strict diagonal dominance is only a sufficient condition for convergence. |
| Successive over-relaxation (SOR) | The Gauss-Seidel value and the old value are blended as $`(1-\omega)x^{(k)} + \omega\tilde{x}`$, $`\omega = 1`$ is the Gauss-Seidel method, and the optimal $`\omega`$ is $`2/(1 + \sqrt{1 - \rho(T_J)^2})`$. |
| Grid SOR | An in-place sweep over the grid performs SOR without a matrix, and the finer the grid, the closer $`\rho(T_{GS})`$ gets to 1, so the optimal $`\omega`$ and the advantage of SOR both grow. |

---

<br>

## Self-Check Questions

1. **Two kinds of L and U:** How does the $`L`$ of the matrix splitting $`A = D + L + U`$ differ from the $`L`$ of the LU factorization $`A = LU`$?

   <details>
   <summary>Answer</summary>

   The $`L`$ of the matrix splitting copies the entries of $`A`$ below the diagonal and has zeros on the diagonal. The $`L`$ of the LU factorization is a unit lower triangular matrix made of elimination factors, with 1s on the diagonal, and in general its entries are not those of $`A`$.

   </details>

2. **Meaning of the division:** In `(b - R @ x) / d` of Problem 2, what matrix operation does `/ d` correspond to?

   <details>
   <summary>Answer</summary>

   `d` is the vector of diagonal entries and `/` is componentwise division, so it is the same as multiplying the vector by $`D^{-1}`$. The inverse of a diagonal matrix is the diagonal matrix of the reciprocals of its diagonal entries.

   </details>

3. **x and x_old:** In Problem 3, does the result change if the right part of `sigma` is replaced with `A[i, i+1:] @ x_old[i+1:]`? What if the left part is replaced with `x_old[:i]`?

   <details>
   <summary>Answer</summary>

   The right part does not change the result, because when component $`i`$ is computed, `x[i+1:]` has not yet been updated and equals `x_old[i+1:]`. Replacing the left part with `x_old[:i]` stops using the new values of this iteration, which turns the method into Jacobi iteration.

   </details>

4. **Ratio of iteration counts:** For the heat network, Jacobi iteration took 23 iterations and the Gauss-Seidel method 13. Explain this ratio with the spectral radius.

   <details>
   <summary>Answer</summary>

   The number of iterations is roughly proportional to $`1/(-\log\rho)`$. For this tridiagonal matrix $`\rho(T_J) = 0.3536`$ and $`\rho(T_{GS}) = 0.125 = \rho(T_J)^2`$, so $`\log\rho(T_{GS}) = 2\log\rho(T_J)`$ and Gauss-Seidel needs about half as many iterations. The ratio $`13/23 \approx 0.57`$ agrees with this.

   </details>

5. **Effect of the weighting factor:** On the heat network $`\omega = 1.7`$ was slower than Gauss-Seidel (67 iterations against 13), but on the 2D grid it was much faster (243 iterations against 1267). Explain why.

   <details>
   <summary>Answer</summary>

   The optimal $`\omega`$ is $`2/(1 + \sqrt{1 - \rho(T_J)^2})`$ and depends on $`\rho(T_J)`$. The heat network has a small $`\rho(T_J) = 0.354`$, so its optimum is 1.033 and 1.7 is far too large. The grid with $`n = 30`$ has $`\rho(T_J) = \cos(\pi/31) \approx 0.995`$, close to 1, so its optimum is about 1.816, and 1.7, being close to it, reduces the iterations greatly.

   </details>

6. **Tolerance and actual error:** In Problem 6 both methods stopped when the update became smaller than $`10^{-6}`$, so why do their results differ by $`8.22 \times 10^{-5}`$?

   <details>
   <summary>Answer</summary>

   A small update does not mean that convergence is complete. The error remains at roughly $`\rho/(1-\rho)`$ times the update, about 97 times for the Gauss-Seidel value $`\rho \approx 0.990`$. In fact the error of the Gauss-Seidel result is $`9.6 \times 10^{-5}`$, which accounts for most of the printed difference.

   </details>
