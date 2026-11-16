# HW3 Assignment Explanation — LU Factorization and Matrix Inverse

[Repository Guide](../../README.md) | [English](Assignment-Explanation.md) | [한국어](Assignment-Explanation.ko.md)

> **Last Updated:** 2026-11-16
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 10, 11

> **Learning Objectives:**
> 1. Compute $`L`$ and $`U`$ with Doolittle's algorithm and compare them with SciPy's pivoted LU factorization
> 2. Solve an LU-factored system by forward and backward substitution and solve several right-hand sides with one factorization
> 3. Solve a symmetric positive definite system with the Cholesky factorization and verify the result against a direct solve
> 4. Compute the columns of a matrix inverse one at a time with the LU factorization
> 5. Diagnose an ill-conditioned system with the condition number and estimate the number of significant digits lost

---

## Table of Contents

- [1. Assignment Overview](#1-assignment-overview)
- [2. Requirements](#2-requirements)
- [3. Problem 1: Doolittle's Algorithm](#3-problem-1-doolittles-algorithm)
  - [3.1 Approach](#31-approach)
  - [3.2 Code Written](#32-code-written)
  - [3.3 Results and Verification](#33-results-and-verification)
- [4. Problem 2: Solving with L and U](#4-problem-2-solving-with-l-and-u)
  - [4.1 Approach](#41-approach)
  - [4.2 Code Written](#42-code-written)
  - [4.3 Results and Verification](#43-results-and-verification)
- [5. Problem 3: Spring-Mass System under Multiple Load Cases](#5-problem-3-spring-mass-system-under-multiple-load-cases)
  - [5.1 Approach](#51-approach)
  - [5.2 Code Written](#52-code-written)
  - [5.3 Results and Verification](#53-results-and-verification)
- [6. Problem 4: Solving with the Cholesky Factorization](#6-problem-4-solving-with-the-cholesky-factorization)
  - [6.1 Approach](#61-approach)
  - [6.2 Code Written](#62-code-written)
  - [6.3 Results and Verification](#63-results-and-verification)
- [7. Problem 5: Matrix Inverse via LU](#7-problem-5-matrix-inverse-via-lu)
  - [7.1 Approach](#71-approach)
  - [7.2 Code Written](#72-code-written)
  - [7.3 Results and Verification](#73-results-and-verification)
- [8. Problem 6: Diagnosing an Ill-Conditioned System](#8-problem-6-diagnosing-an-ill-conditioned-system)
  - [8.1 Approach](#81-approach)
  - [8.2 Code Written](#82-code-written)
  - [8.3 Results and Verification](#83-results-and-verification)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Assignment Overview

> **Try It First:** This document contains the answers and solutions. Solve the problems in the [assignment handout](Original/Assignment.md) yourself before reading on.

The assignment is to fill in the six blanks of the lab notebook `Lab04-LU.ipynb` distributed by the instructor and run it. The notebook has five parts: LU factorization, multiple right-hand sides, Cholesky factorization, the matrix inverse, and the condition number with ill-conditioning. The problem statements are in the [assignment handout](Original/Assignment.md), the distributed blank notebook is [Original/Lab04-LU.ipynb](Original/Lab04-LU.ipynb), and the notebook with the blanks filled in and run is [Solution/Lab04-LU_Completed.ipynb](Solution/Lab04-LU_Completed.ipynb).

| Problem | Notebook section | Blanks to fill |
|:-----|:-----|:-----|
| Problem 1 | 1.2 Doolittle's algorithm | `U[k, j]` and `L[i, k]` in the function `lu_doolittle` |
| Problem 2 | 1.4 Solving with $`L`$ and $`U`$ | `d[i]` in the function `forward_sub` and `x[i]` in the function `backward_sub` |
| Problem 3 | 2.2 Multiple load cases | `x_gravity`, `x_wind`, `x_combined` |
| Problem 4 | 3.2 Cholesky factorization | `d` and `x` in the function `cholesky_solve` |
| Problem 5 | 4.2 Matrix inverse via LU | `A_inv[:, j]` in the function `lu_inverse` |
| Problem 6 | 5.4 Diagnosing an ill-conditioned system | `kappa`, `digits_lost`, `x_susp` |

*Table 1. The six blanks that make up the assignment*

LU factorization, forward and backward substitution, LU with pivoting, and the Cholesky factorization are covered in the [Chapter 10 lecture note](../../Lecture-Notes/Ch10_LU-Factorization/Concepts_Lecture.md) and the [Chapter 10 lab note](../../Lecture-Notes/Ch10_LU-Factorization/Concepts_Lab.md), and the matrix inverse and the condition number are covered in the [Chapter 11 lecture note](../../Lecture-Notes/Ch11_Matrix-Inverse-and-Condition/Concepts_Lecture.md) and the [Chapter 11 lab note](../../Lecture-Notes/Ch11_Matrix-Inverse-and-Condition/Concepts_Lab.md). This document repeats only what the solutions need.

---

<br>

## 2. Requirements

The completed notebook must satisfy the following requirements. The rest of each cell stays as distributed, and only the blanks are filled in.

1. **Problem 1:** `U[k, j]` and `L[i, k]` must implement the sums of Doolittle's formulas with the entries already computed, and `L @ U` must reproduce the original matrix in the test cell.
2. **Problem 2:** `forward_sub` must be forward substitution for a unit lower triangular matrix and `backward_sub` backward substitution for an upper triangular matrix. The solution from `lu_solve` must equal that of `np.linalg.solve`.
3. **Problem 3:** The three load cases must be solved with the single factorization `lu_K`, and the solution for the combined load must equal the sum of the solutions for the two loads.
4. **Problem 4:** As the comments instruct, the two triangular systems must be solved with `forward_sub` and `backward_sub`, and the cell that follows compares the result with `np.linalg.solve`.
5. **Problem 5:** Column $`j`$ of `A_inv` must be the solution of $`A\mathbf{x}_j = \mathbf{e}_j`$, and the factorization is done only once, outside the loop.
6. **Problem 6:** The condition number, the number of significant digits lost, and the solution must be computed as the comments instruct, and the results are used to diagnose the system and think of a remedy.

---

<br>

## 3. Problem 1: Doolittle's Algorithm

### 3.1 Approach

The sum $`\sum_{m=0}^{k-1} \ell_{km} u_{mj}`$ in the formula is the inner product of the leading part of row $`k`$ of $`L`$ and the upper part of column $`j`$ of $`U`$, so it is written as `L[k, :k] @ U[:k, j]`. Likewise, $`\sum_{m=0}^{k-1} \ell_{im} u_{mk}`$ is `L[i, :k] @ U[:k, k]`. When $`k = 0`$ both slices are empty, and the inner product of arrays of length 0 is 0, so $`u_{0j} = a_{0j}`$ and $`\ell_{i0} = a_{i0}/u_{00}`$ hold without special handling.

The order of computation also matters. Row $`k`$ of $`U`$ uses only the leading part of row $`k`$ of $`L`$ (columns found in earlier steps), and column $`k`$ of $`L`$ divides by the newly computed $`u_{kk}`$, so each step computes the row of $`U`$ first and the column of $`L`$ second. The structure of the LU factorization and its relation to the elimination factors are covered in Sections 1 and 2 of the [Chapter 10 lecture note](../../Lecture-Notes/Ch10_LU-Factorization/Concepts_Lecture.md) and Section 1 of the [Chapter 10 lab note](../../Lecture-Notes/Ch10_LU-Factorization/Concepts_Lab.md).

### 3.2 Code Written

```python
            U[k, j] = A[k, j] - L[k, :k] @ U[:k, j]
```

```python
            L[i, k] = (A[i, k] - L[i, :k] @ U[:k, k]) / U[k, k]
```

### 3.3 Results and Verification

```text
L =
 [[1.         0.         0.        ]
 [0.25       1.         0.        ]
 [0.5        0.66666667 1.        ]]
U =
 [[ 8.          7.          9.        ]
 [ 0.         -0.75       -1.25      ]
 [ 0.          0.         -0.66666667]]
LU =
 [[8. 7. 9.]
 [2. 1. 1.]
 [4. 3. 3.]]
L_sc =
 [[1.         0.         0.        ]
 [0.25       1.         0.        ]
 [0.5        0.66666667 1.        ]]
U_sc =
 [[ 8.          7.          9.        ]
 [ 0.         -0.75       -1.25      ]
 [ 0.          0.         -0.66666667]]
L_scU_sc =
 [[8. 7. 9.]
 [2. 1. 1.]
 [4. 3. 3.]]
```

By hand, row 0 of $`U`$ is $`[8, 7, 9]`$, $`\ell_{10} = 2/8 = 0.25`$, $`\ell_{20} = 4/8 = 0.5`$, and then $`u_{11} = 1 - 0.25 \cdot 7 = -0.75`$, $`u_{12} = 1 - 0.25 \cdot 9 = -1.25`$, $`\ell_{21} = (3 - 0.5 \cdot 7)/(-0.75) = 2/3`$, and $`u_{22} = 3 - 0.5 \cdot 9 - (2/3)(-1.25) = -2/3`$, matching the output, and `L @ U` reproduces the test matrix.

Contrary to the comment in the cell, the test matrix `A_test` is the example $`A`$ of Section 1.2 of the notebook with its rows reordered. In this order the first entry 8 of column 0 is the largest, and after elimination the upper of the lower values $`-0.75`$ and $`-0.5`$ in column 1 also has the larger absolute value, so `sc.linalg.lu`, which pivots partially, does not swap rows either, $`P = I`$, and the two results coincide. Running the function on the original order of Section 1.2 gives the $`L`$ and $`U`$ found by hand in the notebook, while applying `sc.linalg.lu` to that matrix gives the $`L`$ and $`U`$ above together with a $`P \neq I`$ that moves the row containing 8 to the top.

---

<br>

## 4. Problem 2: Solving with L and U

### 4.1 Approach

The sum $`\sum_{j=0}^{i-1} \ell_{ij} d_j`$ of forward substitution is `L[i, :i] @ d[:i]`. As its docstring states, `forward_sub` assumes a unit lower triangular matrix with 1s on the diagonal, so it does not divide by $`\ell_{ii}`$. The sum $`\sum_{j=i+1}^{n-1} u_{ij} x_j`$ of backward substitution is `U[i, i+1:] @ x[i+1:]`, and the diagonal entries of $`U`$ are not 1, so it divides by `U[i, i]`.

The function `sc.linalg.lu` used by `lu_solve` returns a $`P`$ that satisfies $`A = PLU`$. Then $`A\mathbf{x} = \mathbf{b}`$ becomes $`LU\mathbf{x} = P^{-1}\mathbf{b} = P^T\mathbf{b}`$, so the right-hand side is replaced with `P.T @ b` and the two functions are applied in turn. Forward and backward substitution are covered in Sections 4, 5, and 8 of the [Chapter 10 lecture note](../../Lecture-Notes/Ch10_LU-Factorization/Concepts_Lecture.md) and Section 2 of the [Chapter 10 lab note](../../Lecture-Notes/Ch10_LU-Factorization/Concepts_Lab.md).

### 4.2 Code Written

```python
        d[i] = b[i] - L[i, :i] @ d[:i]
```

```python
        x[i] = (d[i] - U[i, i+1:] @ x[i+1:]) / U[i, i]
```

### 4.3 Results and Verification

```text
Solution via LU solve : [ 3.  -2.5  7. ]
Solution via np.linalg: [ 3.  -2.5  7. ]
Max difference        : 8.881784197001252e-16
Residual ||Ax - b||   : 1.4210854715202004e-14
```

Substituting the solution $`(3, -2.5, 7)`$ into the first equation gives $`3 \cdot 3 - 0.1(-2.5) - 0.2 \cdot 7 = 7.85`$, equal to the right-hand side. The difference of $`8.9 \times 10^{-16}`$ from `np.linalg.solve` is at the level of machine epsilon relative to the size of the solution (at most 7).

The thermal resistor network of Section 1.5 of the notebook is also solved with `lu_solve`. The completed notebook prints $`T_2 = 100`$°C and $`T_3 = 50`$°C, and $`1.75 \cdot 100 - 0.5 \cdot 50 = 150`$ and $`-0.5 \cdot 100 + 2 \cdot 50 = 50`$, so both equations hold.

> **Note:** The same cell ends by printing `Energy balance at node 2 (should = 0): 100.000000 W`. The expression is `Q2 + q12_out + q23_out + q20_out`, and the heat leaving node 2 is $`0 + 25 + 25 = 50`$ W, equal to the supplied heat $`Q_2 = 50`$ W, so the energy balance holds. The quantity that becomes 0 is the outgoing heat minus $`Q_2`$; because the printed expression adds $`Q_2`$, it prints $`2Q_2 = 100`$ W at balance.

---

<br>

## 5. Problem 3: Spring-Mass System under Multiple Load Cases

### 5.1 Approach

`sc.linalg.lu_factor` packs $`L`$ and $`U`$ into one matrix and returns the pivot information with it. Passing this result `lu_K` to `sc.linalg.lu_solve` gives the solution by forward and backward substitution alone, without factoring again. The three load cases are therefore solved three times with the same `lu_K`, changing only the right-hand side. The advantage of the LU factorization for multiple right-hand sides is covered in Section 4 of the [Chapter 10 lab note](../../Lecture-Notes/Ch10_LU-Factorization/Concepts_Lab.md).

### 5.2 Code Written

```python
x_gravity  = sc.linalg.lu_solve(lu_K, F_gravity)
x_wind     = sc.linalg.lu_solve(lu_K, F_wind)
x_combined = sc.linalg.lu_solve(lu_K, F_combined)
```

### 5.3 Results and Verification

```text
Load Case         x1 (m)      x2 (m)      x3 (m)
------------------------------------------------
Gravity          0.52000     0.88000     0.92000
Wind             1.46000     1.74000     1.16000
Combined         1.98000     2.62000     2.08000

Superposition check ||x_comb - (x_grav + x_wind)||: 4.97e-16
```

The output matches the output in the distributed notebook. Substituting the displacements under gravity gives $`70(0.52) - 30(0.88) = 10`$, $`-30(0.52) + 50(0.88) - 20(0.92) = 10`$, and $`-20(0.88) + 30(0.92) = 10`$, equal to the right-hand side $`[10, 10, 10]`$. Because the system is linear, $`K(\mathbf{x}_g + \mathbf{x}_w) = \mathbf{F}_g + \mathbf{F}_w`$, and the solution for the combined load must be the sum of the two solutions (the principle of superposition). The difference of $`4.97 \times 10^{-16}`$ is at the level of roundoff error.

The cell of Section 2.1 of the notebook compares, as the number of right-hand sides $`m`$ grows, the time of calling `np.linalg.solve` for each system with that of factoring once and calling `lu_solve`. In the completed notebook the speedup grows from 0.9 times at $`m = 1`$ to 6.0 times at $`m = 100`$. The more right-hand sides there are, the more solves share the cost of the single factorization.

---

<br>

## 6. Problem 4: Solving with the Cholesky Factorization

### 6.1 Approach

With $`A = LL^T`$, the system $`A\mathbf{x} = \mathbf{b}`$ becomes the two triangular systems $`L\mathbf{d} = \mathbf{b}`$ and $`L^T\mathbf{x} = \mathbf{d}`$. The comments above the blanks instruct using `forward_sub` defined in Section 1 for the first system and noting that $`L^T`$ is upper triangular for the second, and the code written accordingly uses `forward_sub(L, b)` and `backward_sub(L.T, d)`. The Cholesky factorization is covered in Section 9 of the [Chapter 10 lecture note](../../Lecture-Notes/Ch10_LU-Factorization/Concepts_Lecture.md) and Section 5 of the [Chapter 10 lab note](../../Lecture-Notes/Ch10_LU-Factorization/Concepts_Lab.md).

### 6.2 Code Written

```python
    d = forward_sub(L, b)
```

```python
    x = backward_sub(L.T, d)
```

### 6.3 Results and Verification

```text
Is K_fem SPD? True

Displacements (m):
  x1 = 222949.7266 mm
  x2 = 445899.4532 mm
  x3 = 668849.1799 mm
  x4 = 886025.4038 mm

Max error vs np.linalg.solve: 885.4754037844386

Solve time (n=500):  Cholesky = 0.0063s,  LU = 0.0012s,  speedup = 0.2x
```

The displacements reach hundreds of meters and the largest difference from `np.linalg.solve` is 885 m, so this result is wrong. The correct displacements are $`[150, 300, 450, 550]`$ mm, which can be checked by substituting them into each equation, for example $`10^4 (2 \cdot 0.15 - 0.3) = 0`$ and $`10^4(-0.3 + 0.9 - 0.55) = 500`$.

The cause lies in the assumption of `forward_sub`. As its docstring states, the `forward_sub` of Problem 2 assumes a unit lower triangular matrix with 1s on the diagonal, so it stops at `d[i] = b[i] - L[i, :i] @ d[:i]` and does not divide by $`\ell_{ii}`$. The $`L`$ of the LU factorization satisfies this assumption, but the diagonal entries of the Cholesky factor $`L`$ are $`141.42, 122.47, 115.47, 50`$ for this matrix, not 1. The $`\mathbf{d}`$ of the first step is therefore wrong, and the result is wrong even though the backward substitution that follows works correctly. `backward_sub` divides by `U[i, i]`, so it can be used as is for $`L^T`$, whose diagonal entries are not 1.

> **Note:** Solving with `forward_sub(L, b)`, using the Section 1 `forward_sub` as the comment above the blank instructs, does not hold for a Cholesky factor. Forward substitution that divides by the diagonal entries gives the correct displacements $`[150, 300, 450, 550]`$ mm, with a difference of $`3.3 \times 10^{-16}`$ m from `np.linalg.solve`. The implementation below is corrected that way and is not the code of the submitted notebook.

```python
def forward_sub_general(L, b):
    """Solve L d = b for a lower triangular L with any nonzero diagonal."""
    n = len(b)
    d = np.zeros(n)
    for i in range(n):
        d[i] = (b[i] - L[i, :i] @ d[:i]) / L[i, i]
    return d
```

If only the one blank line is to change, a SciPy function that solves a general lower triangular system, as in `d = sc.linalg.solve_triangular(L, b, lower=True)`, gives the same result.

In the timing comparison on the last line, the SPD solver was measured as slower. A single measurement is strongly affected by function-call overhead and by the machine, so this result cannot be used to judge the theoretical cost advantage of the Cholesky factorization (about $`n^3/6`$ against $`n^3/3`$).

---

<br>

## 7. Problem 5: Matrix Inverse via LU

### 7.1 Approach

`I[:, j]` is the $`j`$-th column $`\mathbf{e}_j`$ of the identity matrix, so solving `sc.linalg.lu_solve(lu_obj, I[:, j])` with the `lu_obj` created once outside the loop gives the solution of $`A\mathbf{x}_j = \mathbf{e}_j`$, which is column $`j`$ of $`A^{-1}`$. This is stored in `A_inv[:, j]`. The factorization costs $`O(n^3)`$ and each column solve $`O(n^2)`$. The matrix inverse via LU is covered in Section 2 of the [Chapter 11 lecture note](../../Lecture-Notes/Ch11_Matrix-Inverse-and-Condition/Concepts_Lecture.md) and Section 1 of the [Chapter 11 lab note](../../Lecture-Notes/Ch11_Matrix-Inverse-and-Condition/Concepts_Lab.md).

### 7.2 Code Written

```python
        A_inv[:, j] = sc.linalg.lu_solve(lu_obj, I[:, j])
```

### 7.3 Results and Verification

```text
A^{-1} via lu_inverse =
 [[ 0.332489  0.004944  0.006798]
 [-0.005182  0.142903  0.004183]
 [-0.010078  0.00271   0.09988 ]]

A^{-1} via np.linalg.inv =
 [[ 0.332489  0.004944  0.006798]
 [-0.005182  0.142903  0.004183]
 [-0.010078  0.00271   0.09988 ]]

Max error: 0.0
Sanity check  A @ A_inv - I  max = 2.220446049250313e-16

Determinant of singular matrix: 0.0
Error: Singular matrix
```

The largest entry of $`AA^{-1} - I`$ is $`2.2 \times 10^{-16}`$, at the level of machine epsilon, so the inverse is correct. The difference from `np.linalg.inv` being exactly 0 can be attributed to `np.linalg.inv` also using a LAPACK routine that solves $`AX = I`$ with a partially pivoted LU factorization. The singular matrix has a second row equal to twice the first, so its determinant is 0, and requesting its inverse raises a `LinAlgError`.

Section 4.3 of the notebook demonstrates the rule of solving the system directly instead of forming the inverse explicitly. In the completed notebook, multiplying by the inverse took 0.0011 s and `np.linalg.solve` took 0.0004 s, and the residuals were similar at $`1.59 \times 10^{-14}`$ and $`1.58 \times 10^{-14}`$.

---

<br>

## 8. Problem 6: Diagnosing an Ill-Conditioned System

### 8.1 Approach

Following the order of the numbered comments between the blanks, the 2-norm condition number is computed with `np.linalg.cond`, the integer part of $`\log_{10}\kappa`$ is taken as the number of significant digits lost following the rule of thumb "if $`\kappa \approx 10^k`$, about $`k`$ digits are lost," and the solution is obtained with `np.linalg.solve`. The condition number and the error bound are covered in Sections 7-9 of the [Chapter 11 lecture note](../../Lecture-Notes/Ch11_Matrix-Inverse-and-Condition/Concepts_Lecture.md) and Section 2 of the [Chapter 11 lab note](../../Lecture-Notes/Ch11_Matrix-Inverse-and-Condition/Concepts_Lab.md).

### 8.2 Code Written

```python
kappa = np.linalg.cond(K_susp)
```

```python
digits_lost = int(np.log10(kappa))
```

```python
x_susp = np.linalg.solve(K_susp, f_susp)
```

### 8.3 Results and Verification

```text
Condition number kappa = 2.00e+08
Significant digits lost ≈ 8
Solution x = [0.0005 0.0005 0.1   ]
Residual ||Kx - f|| = 0.00e+00

Diagnosis: The upper-left 2x2 block has nearly identical rows
(rows differ by only 1 N/m out of 1e6), making the system near-singular.
The solution in the first two DOFs is unreliable.
```

The output matches the output in the distributed notebook. $`K`$ is symmetric, so its 2-norm condition number is the ratio of the largest to the smallest absolute eigenvalue. The upper-left $`2 \times 2`$ block has the eigenvalue $`2 \times 10^6 - 1`$ in the direction $`(1, 1)`$ and $`1`$ in the direction $`(1, -1)`$, and the third degree of freedom has the eigenvalue $`0.01`$, so $`\kappa_2 = (2 \times 10^6 - 1)/0.01 \approx 2.0 \times 10^8`$, and $`\log_{10}\kappa_2 \approx 8.3`$ gives 8 digits. The solution is $`x_1 = x_2 = 1000/(2 \times 10^6 - 1) \approx 5.0 \times 10^{-4}`$ m and $`x_3 = 0.001/0.01 = 0.1`$ m.

To check the diagnosis concretely, solving again with $`f_2`$ changed by only 0.001 N (a relative change of $`10^{-6}`$) gives the solution $`[5.0 \times 10^{-10}, 1.0 \times 10^{-3}, 0.1]`$. A change of one part in a million in the right-hand side changes $`x_1`$ and $`x_2`$ by nearly 100%. Adding and subtracting the two equations reveals why. The sum $`s = x_1 + x_2`$ is determined with stiffness $`2 \times 10^6 - 1`$ by $`(2 \times 10^6 - 1)s = f_1 + f_2`$, but the difference $`d = x_1 - x_2`$ is determined by $`d = f_1 - f_2`$ with a stiffness of only 1 N/m, so small errors in $`f_1`$ and $`f_2`$ pass directly into $`d`$.

The remedies can be summarized as follows.

- **Checking the model and data:** First confirm whether the stiffness difference of 1 N/m between the two rows reflects the actual structure or comes from a unit or input error. If the two degrees of freedom are effectively connected rigidly, rebuild the model with them merged into one.
- **Change of variables:** Using $`s`$ and $`d`$ as unknowns turns each into a well-conditioned scalar equation. The sensitivity of $`d`$ to $`f_1 - f_2`$, however, is a property of the problem itself, so if $`d`$ is needed, the difference of the loads must be known with matching precision.
- **Scaling:** Multiplying the third equation by $`10^8`$ to balance its size reduces the condition number to the block's value of about $`2 \times 10^6`$. This removes only the part that comes from the size difference with the third degree of freedom, and the ill-conditioning of the block itself remains.

> **Note:** Of the condition number $`2 \times 10^8`$, a factor of 100 comes from the small stiffness 0.01 of the unconnected third degree of freedom. The solution $`x_3 = 0.1`$ of the third degree of freedom comes from an independent scalar equation and is accurate; what is actually unreliable is, as the notebook's diagnosis says, the first two degrees of freedom, and in particular their difference $`x_1 - x_2`$.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Doolittle's algorithm | Each step computes the row of $`U`$ first and the column of $`L`$ second, and the sums in the formulas are computed as inner products of entries already found. |
| Comparison with pivoting | If each pivot is already the largest in its column, partial pivoting swaps no rows, so Doolittle's factorization without pivoting and SciPy's factorization coincide. |
| Forward and backward substitution | $`L\mathbf{d} = P^T\mathbf{b}`$ is solved from the top and $`U\mathbf{x} = \mathbf{d}`$ from the bottom, and a triangular matrix whose diagonal entries are not 1 must be divided by them. |
| Multiple right-hand sides | One factorization with `lu_factor` is followed by an $`O(n^2)`$ `lu_solve` for each right-hand side, and because the system is linear the solutions satisfy the principle of superposition. |
| Solving with the Cholesky factorization | With $`A = LL^T`$ two triangular systems are solved, but the diagonal entries of the Cholesky factor are not 1, so forward substitution that assumes a unit diagonal gives a wrong solution. |
| Matrix inverse via LU | After one factorization, each column of the identity matrix is solved as a right-hand side to give a column of the inverse, and when solving a system it is better not to form the inverse. |
| Condition number diagnosis | If $`\kappa(A) \approx 10^k`$, about $`k`$ digits are lost, and two nearly identical rows leave the difference of two unknowns almost undetermined. |

---

<br>

## Self-Check Questions

1. **Empty slices:** In Problem 1, what value does `L[k, :k] @ U[:k, j]` take when $`k = 0`$, and why does it agree with the formula?

   <details>
   <summary>Answer</summary>

   Both slices become arrays of length 0, and their inner product is 0. The sum $`\sum_{m=0}^{-1}`$ in the formula is also an empty sum equal to 0, so $`u_{0j} = a_{0j}`$ and the two agree.

   </details>

2. **Transpose of P:** Why does `lu_solve` in Problem 2 replace the right-hand side with `P.T @ b` rather than `P @ b`?

   <details>
   <summary>Answer</summary>

   `sc.linalg.lu` returns a $`P`$ that satisfies $`A = PLU`$. Then $`LU\mathbf{x} = P^{-1}\mathbf{b}`$, and a permutation matrix is orthogonal, so $`P^{-1} = P^T`$ and the right-hand side is $`P^T\mathbf{b}`$.

   </details>

3. **Error in the Cholesky solve:** Explain the cause of the wrong result in Problem 4 and how to fix it.

   <details>
   <summary>Answer</summary>

   `forward_sub` assumes a unit lower triangular matrix and does not divide by the diagonal entries, but the diagonal entries of the Cholesky factor $`L`$ are $`141.42, 122.47, 115.47, 50`$. Dividing by `L[i, i]` in forward substitution or using `sc.linalg.solve_triangular(L, b, lower=True)` gives the correct displacements $`[150, 300, 450, 550]`$ mm.

   </details>

4. **Where to factor:** In Problem 5, how would the result and the cost change if `sc.linalg.lu_factor(A)` were placed inside the loop?

   <details>
   <summary>Answer</summary>

   The result would be the same, but the $`O(n^3)`$ factorization would be repeated for every column, so the total cost would grow from $`O(n^3) + nO(n^2)`$ to $`nO(n^3) = O(n^4)`$.

   </details>

5. **Composition of the condition number:** The condition number $`2 \times 10^8`$ in Problem 6 is the ratio of which two values, and why does scaling the third equation reduce it to $`2 \times 10^6`$?

   <details>
   <summary>Answer</summary>

   The 2-norm condition number of a symmetric matrix is the ratio of the largest to the smallest absolute eigenvalue, here the block's $`2 \times 10^6 - 1`$ to the third degree of freedom's $`0.01`$. Multiplying the third equation by $`10^8`$ makes that eigenvalue $`10^6`$, so the smallest eigenvalue becomes the block's 1 and the condition number becomes about $`2 \times 10^6`$.

   </details>
