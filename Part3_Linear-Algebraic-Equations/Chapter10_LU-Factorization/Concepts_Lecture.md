# Chapter 10 Lecture — LU Factorization

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 10

> **Prerequisites**: [Linear Algebra] Gaussian elimination (Ch 8-9).
>
> **Learning Objectives**:
> 1. Decompose matrices into L and U factors
> 2. Apply LU factorization for efficient multiple solves
> 3. Implement Cholesky decomposition for symmetric matrices

---

<br>

## Table of Contents

- [1. Overview of LU Factorization](#1-overview-of-lu-factorization)
  - [1.1 The Core Idea](#11-the-core-idea)
  - [1.2 Two-Phase Solution Process](#12-two-phase-solution-process)
  - [1.3 Computational Advantage](#13-computational-advantage)
- [2. Gauss Elimination as LU Factorization](#2-gauss-elimination-as-lu-factorization)
  - [2.1 Elimination Factors Form L](#21-elimination-factors-form-l)
  - [2.2 Structure of L and U for a 3x3 System](#22-structure-of-l-and-u-for-a-3x3-system)
  - [2.3 Verification](#23-verification)
- [3. Storage Efficiency](#3-storage-efficiency)
- [4. Forward Substitution](#4-forward-substitution)
- [5. Back Substitution](#5-back-substitution)
- [6. Example 10.1: LU Decomposition (Complete Step-by-Step)](#6-example-101-lu-decomposition-complete-step-by-step)
  - [6.1 Problem Setup](#61-problem-setup)
  - [6.2 Step 1 — First Column of L and First Row of U](#62-step-1--first-column-of-l-and-first-row-of-u)
  - [6.3 Step 2 — Second Column of L and Second Row of U](#63-step-2--second-column-of-l-and-second-row-of-u)
  - [6.4 Result](#64-result)
- [7. Example 10.2: Solving with LU](#7-example-102-solving-with-lu)
  - [7.1 Forward Substitution](#71-forward-substitution)
  - [7.2 Back Substitution](#72-back-substitution)
- [8. LU with Pivoting: PA = LU](#8-lu-with-pivoting-pa--lu)
  - [8.1 Why Pivoting is Needed](#81-why-pivoting-is-needed)
  - [8.2 The Permutation Matrix P](#82-the-permutation-matrix-p)
  - [8.3 Procedure](#83-procedure)
  - [8.4 Example 10.3: 2x2 with Small Pivot](#84-example-103-2x2-with-small-pivot)
- [9. Cholesky Factorization](#9-cholesky-factorization)
  - [9.1 Symmetric Positive Definite Matrices](#91-symmetric-positive-definite-matrices)
  - [9.2 Cholesky Decomposition Formulas](#92-cholesky-decomposition-formulas)
  - [9.3 Example 10.5: Cholesky Decomposition](#93-example-105-cholesky-decomposition)
- [10. Python Code](#10-python-code)
- [Summary](#summary)

---

<br>

## 1. Overview of LU Factorization

### 1.1 The Core Idea

The goal of LU factorization is to decompose a square matrix $A$ into the product of two triangular matrices:

$$A = LU$$

where $L$ is a **lower triangular** matrix (all entries above the diagonal are zero) and $U$ is an **upper triangular** matrix (all entries below the diagonal are zero).

> **[Linear Algebra]** LU factorization is mathematically equivalent to Gaussian elimination. The elimination factors are exactly the entries of L. This connection is deep: any algorithm that performs Gaussian elimination implicitly computes an LU factorization.

### 1.2 Two-Phase Solution Process

Once $A = LU$ is known, solving $Ax = b$ is decomposed into two simpler triangular solves:

1. **Forward substitution:** Solve $Ld = b$ for the intermediate vector $d$
2. **Back substitution:** Solve $Ux = d$ for the solution vector $x$

Each of these triangular solves is straightforward because the triangular structure allows us to compute each unknown sequentially without any further elimination.

To verify correctness: $Ax = (LU)x = L(Ux) = Ld = b$.

### 1.3 Computational Advantage

The key advantage of LU factorization is the separation of the expensive decomposition step from the relatively cheap solution step:

| Phase | Cost | Depends on |
|:------|:-----|:-----------|
| Decomposition: $A = LU$ | $O(n^3)$ | Only $A$ |
| Forward substitution: $Ld = b$ | $O(n^2)$ | $L$ and $b$ |
| Back substitution: $Ux = d$ | $O(n^2)$ | $U$ and $d$ |

When the same coefficient matrix $A$ is used with **multiple right-hand side vectors** $b_1, b_2, \ldots, b_m$, we decompose $A$ only once (at $O(n^3)$ cost) and then solve each new right-hand side with just $O(n^2)$ work. This is a dramatic saving compared to performing full Gaussian elimination ($O(n^3)$) for each right-hand side separately.

This situation arises frequently in practice — for example, in structural analysis where the stiffness matrix remains the same but different load cases are applied, or in circuit analysis where the network topology is fixed but different input signals are considered.

---

<br>

## 2. Gauss Elimination as LU Factorization

### 2.1 Elimination Factors Form L

Recall from Chapter 8 that during Gaussian elimination, we compute factors:

$$f_{ij} = \frac{a_{ij}}{a_{jj}}$$

to eliminate element $a_{ij}$ by subtracting $f_{ij}$ times row $j$ from row $i$. The remarkable insight is that these same factors, collected into a lower triangular matrix with 1s on the diagonal, give us the matrix $L$.

### 2.2 Structure of L and U for a 3x3 System

For a $3 \times 3$ system, the factorization produces:

$$L = \begin{bmatrix} 1 & 0 & 0 \\ f_{21} & 1 & 0 \\ f_{31} & f_{32} & 1 \end{bmatrix}, \quad U = \begin{bmatrix} a_{11} & a_{12} & a_{13} \\ 0 & a'_{22} & a'_{23} \\ 0 & 0 & a''_{33} \end{bmatrix}$$

where:

- $f_{21}, f_{31}$ are the factors computed in the first elimination step (eliminating the first column)
- $f_{32}$ is the factor computed in the second elimination step (eliminating the second column)
- $a'_{22}, a'_{23}$ are the modified entries after the first elimination step
- $a''_{33}$ is the modified entry after both elimination steps

$L$ is **unit lower triangular** (diagonal entries are all 1), and $U$ is the upper triangular matrix that results from forward elimination — it is exactly the row-echelon form of $A$.

### 2.3 Verification

The factorization can be verified by performing the matrix multiplication $LU$ and confirming that it equals $A$. For the $3 \times 3$ case:

$$(LU)_{11} = 1 \cdot a_{11} = a_{11}$$

$$(LU)_{21} = f_{21} \cdot a_{11} = \frac{a_{21}}{a_{11}} \cdot a_{11} = a_{21}$$

$$(LU)_{22} = f_{21} \cdot a_{12} + 1 \cdot a'_{22} = f_{21} \cdot a_{12} + (a_{22} - f_{21} \cdot a_{12}) = a_{22}$$

and so on for all entries. Each multiplication reconstructs the original entry of $A$ because the elimination factors encode exactly the row operations performed.

---

<br>

## 3. Storage Efficiency

Since the diagonal of $L$ is always 1 (by construction of the Doolittle method), there is no need to store these values explicitly. This means $L$ and $U$ can be stored together in a single $n \times n$ matrix:

$$\begin{bmatrix} u_{11} & u_{12} & u_{13} \\ l_{21} & u_{22} & u_{23} \\ l_{31} & l_{32} & u_{33} \end{bmatrix}$$

The upper triangular part (including the diagonal) stores $U$, and the strictly lower triangular part stores $L$ (with the implicit understanding that the diagonal of $L$ is all 1s).

This compact storage format is used by most numerical libraries (including LAPACK and SciPy). It is memory-efficient and avoids any redundancy. In practice, the decomposition can even be performed **in place**, overwriting the original matrix $A$ with the combined $[L \setminus U]$ storage.

---

<br>

## 4. Forward Substitution

Given $L$ and $b$, we solve $Ld = b$ for $d$ by working from top to bottom:

$$d_1 = b_1$$

$$d_i = b_i - \sum_{j=1}^{i-1} L_{ij} \, d_j, \quad i = 2, 3, \ldots, n$$

Since $L$ is unit lower triangular (diagonal entries are 1), there is no division by a diagonal element. Each $d_i$ depends only on previously computed values $d_1, d_2, \ldots, d_{i-1}$, so the computation proceeds sequentially from the first equation to the last.

The cost of forward substitution is:

$$\sum_{i=2}^{n}(i-1) = \frac{n(n-1)}{2} \approx \frac{n^2}{2}$$

multiplications and the same number of subtractions, giving a total cost of $O(n^2)$.

---

<br>

## 5. Back Substitution

Given $U$ and $d$, we solve $Ux = d$ for $x$ by working from bottom to top:

$$x_n = \frac{d_n}{U_{nn}}$$

$$x_i = \frac{d_i - \sum_{j=i+1}^{n} U_{ij} \, x_j}{U_{ii}}, \quad i = n-1, n-2, \ldots, 1$$

Each $x_i$ depends only on previously computed values $x_{i+1}, x_{i+2}, \ldots, x_n$. The division by $U_{ii}$ requires that all diagonal entries of $U$ be nonzero — equivalently, that $A$ is nonsingular.

The cost of back substitution is also $O(n^2)$, similar to forward substitution.

---

<br>

## 6. Example 10.1: LU Decomposition (Complete Step-by-Step)

### 6.1 Problem Setup

Decompose the following matrix into $A = LU$:

$$A = \begin{bmatrix} 3 & -0.1 & -0.2 \\ 0.1 & 7 & -0.3 \\ 0.3 & -0.2 & 10 \end{bmatrix}$$

### 6.2 Step 1 — First Column of L and First Row of U

The first row of $U$ is identical to the first row of $A$:

$$U_{1,:} = \begin{bmatrix} 3 & -0.1 & -0.2 \end{bmatrix}$$

Compute the elimination factors for the first column:

$$f_{21} = \frac{a_{21}}{a_{11}} = \frac{0.1}{3} = 0.0333\overline{3}$$

$$f_{31} = \frac{a_{31}}{a_{11}} = \frac{0.3}{3} = 0.1$$

Eliminate the first column (compute the remaining entries of U's first "pass"):

- $a'_{22} = 7 - (0.0333\overline{3})(-0.1) = 7.00333\overline{3}$
- $a'_{23} = -0.3 - (0.0333\overline{3})(-0.2) = -0.29333\overline{3}$
- $a'_{32} = -0.2 - (0.1)(-0.1) = -0.19$
- $a'_{33} = 10 - (0.1)(-0.2) = 10.02$

### 6.3 Step 2 — Second Column of L and Second Row of U

The second row of $U$ is:

$$U_{2,:} = \begin{bmatrix} 0 & 7.00333 & -0.29333 \end{bmatrix}$$

Compute the elimination factor for the second column:

$$f_{32} = \frac{a'_{32}}{a'_{22}} = \frac{-0.19}{7.00333} = -0.0271$$

Eliminate the second column:

- $a''_{33} = 10.02 - (-0.0271)(-0.29333) = 10.02 - 0.00795 = 10.0120$

### 6.4 Result

$$L = \begin{bmatrix} 1 & 0 & 0 \\ 0.0333 & 1 & 0 \\ 0.1 & -0.0271 & 1 \end{bmatrix}, \quad U = \begin{bmatrix} 3 & -0.1 & -0.2 \\ 0 & 7.00333 & -0.29333 \\ 0 & 0 & 10.0120 \end{bmatrix}$$

Verification: multiplying $LU$ reproduces the original matrix $A$.

---

<br>

## 7. Example 10.2: Solving with LU

Using the $L$ and $U$ from Example 10.1, solve $Ax = b$ where:

$$b = \begin{bmatrix} 7.85 \\ -19.3 \\ 71.4 \end{bmatrix}$$

### 7.1 Forward Substitution

Solve $Ld = b$:

$$d_1 = b_1 = 7.85$$

$$d_2 = b_2 - L_{21} \, d_1 = -19.3 - (0.0333)(7.85) = -19.3 - 0.2617 = -19.5617$$

$$d_3 = b_3 - L_{31} \, d_1 - L_{32} \, d_2 = 71.4 - (0.1)(7.85) - (-0.0271)(-19.5617)$$

$$= 71.4 - 0.785 - 0.5303 = 70.0847$$

> **Note:** The slight differences in the last digits (e.g., $70.0847$ vs. $70.0868$) depend on how many decimal places are carried through the intermediate calculations. With full precision, $d_3 \approx 70.0843$.

### 7.2 Back Substitution

Solve $Ux = d$:

$$x_3 = \frac{d_3}{U_{33}} = \frac{70.0843}{10.0120} = 7.0003$$

$$x_2 = \frac{d_2 - U_{23} \, x_3}{U_{22}} = \frac{-19.5617 - (-0.29333)(7.0003)}{7.00333} = \frac{-19.5617 + 2.0534}{7.00333} = \frac{-17.5083}{7.00333} = -2.5001$$

$$x_1 = \frac{d_1 - U_{12} \, x_2 - U_{13} \, x_3}{U_{11}} = \frac{7.85 - (-0.1)(-2.5001) - (-0.2)(7.0003)}{3}$$

$$= \frac{7.85 - 0.2500 + 1.4001}{3} = \frac{9.0001}{3} = 3.0000$$

**Final solution:**

$$x_1 = 3.0000, \quad x_2 = -2.5001, \quad x_3 = 7.0003$$

This matches the solution obtained by direct Gaussian elimination in Chapter 8, confirming the equivalence of the two approaches. The small deviations from the exact values $(3, -2.5, 7)$ are due to round-off in the intermediate calculations.

---

<br>

## 8. LU with Pivoting: $PA = LU$

### 8.1 Why Pivoting is Needed

Just as with Gaussian elimination, LU factorization without pivoting can fail or produce inaccurate results when:

- A pivot element is zero (causing division by zero)
- A pivot element is very small (amplifying round-off errors)

Partial pivoting — selecting the largest absolute value in the current column as the pivot — resolves both issues.

### 8.2 The Permutation Matrix P

When row swaps are performed during the factorization, they are recorded in a **permutation matrix** $P$. The factorization becomes:

$$PA = LU$$

where $P$ is an $n \times n$ matrix obtained by reordering the rows of the identity matrix according to the row swaps performed during elimination. The permutation matrix has exactly one entry of 1 in each row and each column, with all other entries being 0.

Key properties of permutation matrices:

- $P^{-1} = P^T$ (orthogonal matrix)
- $\det(P) = \pm 1$
- $PP^T = I$

### 8.3 Procedure

The solution procedure with pivoting is:

1. **Factorize:** Decompose $A$ with partial pivoting to obtain $P$, $L$, and $U$ such that $PA = LU$
2. **Permute the RHS:** Compute $Pb$ (reorder the entries of $b$ according to the same row swaps)
3. **Forward substitution:** Solve $Ld = Pb$
4. **Back substitution:** Solve $Ux = d$

### 8.4 Example 10.3: 2x2 with Small Pivot

Consider the system:

$$\begin{bmatrix} 0.0003 & 3.0000 \\ 1.0000 & 1.0000 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 2.0001 \\ 1.0000 \end{bmatrix}$$

The exact solution is $x_1 = 1/3$ and $x_2 = 2/3$.

**Without pivoting:**

The factor is $f_{21} = 1.0000 / 0.0003 = 3333.3$, an extremely large number. When this factor multiplies the first row and is subtracted from the second row:

$$a'_{22} = 1.0000 - 3333.3 \times 3.0000 = 1 - 9999.9 = -9998.9$$

$$b'_2 = 1.0000 - 3333.3 \times 2.0001 = 1 - 6667.0 = -6666.0$$

With limited precision (e.g., 4 significant digits), the subtraction of large nearly-equal numbers causes significant loss of accuracy.

**With pivoting:**

Swap rows to put the largest element in the pivot position:

$$P = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}, \quad PA = \begin{bmatrix} 1.0000 & 1.0000 \\ 0.0003 & 3.0000 \end{bmatrix}, \quad Pb = \begin{bmatrix} 1.0000 \\ 2.0001 \end{bmatrix}$$

Now the factor is $f_{21} = 0.0003 / 1.0000 = 0.0003$, a small number that does not amplify errors:

$$a'_{22} = 3.0000 - 0.0003 \times 1.0000 = 2.9997$$

$$b'_2 = 2.0001 - 0.0003 \times 1.0000 = 1.9998$$

Back substitution gives:

$$x_2 = \frac{1.9998}{2.9997} = 0.6667 \approx \frac{2}{3}$$

$$x_1 = 1.0000 - 1.0000 \times 0.6667 = 0.3333 \approx \frac{1}{3}$$

The pivoted version produces accurate results, while the unpivoted version suffers from catastrophic error amplification.

---

<br>

## 9. Cholesky Factorization

### 9.1 Symmetric Positive Definite Matrices

A special and important class of matrices are **symmetric positive definite (SPD)** matrices. A matrix $A$ is SPD if:

1. **Symmetric:** $A = A^T$ (i.e., $a_{ij} = a_{ji}$ for all $i, j$)
2. **Positive definite:** $x^T A x > 0$ for all nonzero vectors $x$

Equivalently, all eigenvalues of an SPD matrix are strictly positive.

> **[Linear Algebra]** SPD matrices arise naturally in physics (mass/stiffness matrices), statistics (covariance matrices), and optimization (Hessians of convex functions). Cholesky is about 2x faster and more numerically stable than general LU.

For SPD matrices, LU factorization can be simplified to the **Cholesky factorization**:

$$A = U^T U$$

where $U$ is an upper triangular matrix. Only the upper triangle of $A$ needs to be stored and processed, resulting in roughly half the computation and half the storage compared to general LU factorization.

### 9.2 Cholesky Decomposition Formulas

The entries of $U$ are computed using the following formulas:

**Diagonal entries:**

$$U_{ii} = \sqrt{a_{ii} - \sum_{k=1}^{i-1} U_{ki}^2}$$

**Off-diagonal entries** (for $j > i$):

$$U_{ij} = \frac{a_{ij} - \sum_{k=1}^{i-1} U_{ki} \, U_{kj}}{U_{ii}}$$

The computation proceeds row by row, from top to bottom. For each row $i$, the diagonal entry $U_{ii}$ is computed first (requiring a square root), and then the off-diagonal entries $U_{ij}$ for $j > i$ are computed.

The square root in the diagonal formula guarantees that $U_{ii} > 0$ for all $i$ if and only if $A$ is positive definite. If $A$ is not positive definite, the algorithm will encounter a negative argument under the square root, which serves as a built-in test for positive definiteness.

### 9.3 Example 10.5: Cholesky Decomposition

Decompose the following SPD matrix:

$$A = \begin{bmatrix} 6 & 15 & 55 \\ 15 & 55 & 225 \\ 55 & 225 & 979 \end{bmatrix}$$

**Row 1:**

$$U_{11} = \sqrt{a_{11}} = \sqrt{6} = 2.4495$$

$$U_{12} = \frac{a_{12}}{U_{11}} = \frac{15}{2.4495} = 6.1237$$

$$U_{13} = \frac{a_{13}}{U_{11}} = \frac{55}{2.4495} = 22.4537$$

**Row 2:**

$$U_{22} = \sqrt{a_{22} - U_{12}^2} = \sqrt{55 - (6.1237)^2} = \sqrt{55 - 37.5} = \sqrt{17.5} = 4.1833$$

$$U_{23} = \frac{a_{23} - U_{12} \, U_{13}}{U_{22}} = \frac{225 - (6.1237)(22.4537)}{4.1833} = \frac{225 - 137.5}{4.1833} = \frac{87.5}{4.1833} = 20.9165$$

**Row 3:**

$$U_{33} = \sqrt{a_{33} - U_{13}^2 - U_{23}^2} = \sqrt{979 - (22.4537)^2 - (20.9165)^2}$$

$$= \sqrt{979 - 504.17 - 437.50} = \sqrt{37.33} = 6.1101$$

**Result:**

$$U = \begin{bmatrix} 2.4495 & 6.1237 & 22.4537 \\ 0 & 4.1833 & 20.9165 \\ 0 & 0 & 6.1101 \end{bmatrix}$$

Verification:

$$U^T U = \begin{bmatrix} 2.4495 & 0 & 0 \\ 6.1237 & 4.1833 & 0 \\ 22.4537 & 20.9165 & 6.1101 \end{bmatrix} \begin{bmatrix} 2.4495 & 6.1237 & 22.4537 \\ 0 & 4.1833 & 20.9165 \\ 0 & 0 & 6.1101 \end{bmatrix} = A$$

---

<br>

## 10. Python Code

SciPy provides efficient implementations of both LU and Cholesky factorizations:

```python
import numpy as np
from scipy.linalg import lu, cholesky

# LU factorization with partial pivoting
A = np.array([[3, -0.1, -0.2],
              [0.1, 7, -0.3],
              [0.3, -0.2, 10]])

P, L, U = lu(A)
print(f'P = \n{P}')
print(f'L = \n{L}')
print(f'U = \n{U}')
print(f'||PA - LU|| = {np.linalg.norm(P @ A - L @ U):.2e}')

# Cholesky factorization (for SPD matrices)
A_spd = np.array([[6, 15, 55],
                   [15, 55, 225],
                   [55, 225, 979]], dtype=float)

U_chol = cholesky(A_spd)          # upper triangular (default)
print(f'U = \n{U_chol}')
print(f'||A - U^T U|| = {np.linalg.norm(A_spd - U_chol.T @ U_chol):.2e}')
```

> **Note:** SciPy's `lu` returns the permutation matrix `P` such that `P @ A = L @ U` (some references define it as `A = P @ L @ U`; always check the documentation for your library). The `cholesky` function returns the upper triangular factor by default; pass `lower=True` to get the lower triangular factor.

---

<br>

## Summary

| Method | Condition | Cost | Key Property |
|:-------|:---------|:-----|:------------|
| LU (Doolittle) | General | $n^3/3$ | $A = LU$ |
| LU with pivoting | General | $n^3/3$ | $PA = LU$ |
| Cholesky | SPD only | $n^3/6$ | $A = U^T U$, unconditionally stable |

| Topic | Key Point |
|:------|:----------|
| Core idea | Decompose $A = LU$, then solve two triangular systems |
| Key advantage | Factor once ($O(n^3)$), solve many RHS ($O(n^2)$ each) |
| Gauss elimination link | Elimination factors $f_{ij}$ become the entries of $L$ |
| Storage | $L$ and $U$ stored in one matrix (diagonal 1s of $L$ are implicit) |
| Forward substitution | $d_i = b_i - \sum_{j=1}^{i-1} L_{ij} d_j$ |
| Back substitution | $x_i = (d_i - \sum_{j=i+1}^{n} U_{ij} x_j) / U_{ii}$ |
| Pivoting | $PA = LU$ with permutation matrix $P$ for numerical stability |
| Cholesky | $A = U^T U$ for SPD matrices — half the cost, unconditionally stable |

---
