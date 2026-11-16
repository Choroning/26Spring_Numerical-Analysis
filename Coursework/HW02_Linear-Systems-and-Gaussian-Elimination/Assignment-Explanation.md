# HW2 Assignment Explanation — Linear Systems and Gaussian Elimination

[Repository Guide](../../README.md) | [English](Assignment-Explanation.md) | [한국어](Assignment-Explanation.ko.md)

> **Last Updated:** 2026-11-16
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 8, 9

> **Learning Objectives:**
> 1. Assemble the right-hand side vector of a linear system from force-equilibrium equations in the order of the rows and verify the solution through the residual
> 2. Implement forward elimination and back substitution on the augmented matrix and compare the result with a hand trace
> 3. Explain how naive elimination fails on a system with a zero pivot and implement the row swap
> 4. Implement the forward sweep of the Thomas algorithm for tridiagonal systems and verify it on a heat conduction problem

---

## Table of Contents

- [1. Assignment Overview](#1-assignment-overview)
- [2. Requirements](#2-requirements)
- [3. Problem 1: Right-Hand Side of the Spring-Mass System](#3-problem-1-right-hand-side-of-the-spring-mass-system)
  - [3.1 Approach](#31-approach)
  - [3.2 Code Written](#32-code-written)
  - [3.3 Results and Verification](#33-results-and-verification)
- [4. Problem 2: Naive Gaussian Elimination](#4-problem-2-naive-gaussian-elimination)
  - [4.1 Approach](#41-approach)
  - [4.2 Code Written](#42-code-written)
  - [4.3 Results and Verification](#43-results-and-verification)
- [5. Problem 3: Partial Pivoting](#5-problem-3-partial-pivoting)
  - [5.1 Approach](#51-approach)
  - [5.2 Code Written](#52-code-written)
  - [5.3 Results and Verification](#53-results-and-verification)
- [6. Problem 4: The Thomas Algorithm](#6-problem-4-the-thomas-algorithm)
  - [6.1 Approach](#61-approach)
  - [6.2 Code Written](#62-code-written)
  - [6.3 Results and Verification](#63-results-and-verification)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Assignment Overview

> **Try It First:** This document contains the answers and solutions. Solve the problems in the [assignment handout](Original/Assignment.md) yourself before reading on.

The assignment is to fill in the four blanks of the lab notebook `Lab03-LinearSystem.ipynb` distributed by the instructor and run it. The notebook has four parts: setting up and solving $`A\mathbf{x} = \mathbf{b}`$, naive Gaussian elimination, partial pivoting, and tridiagonal systems with the Thomas algorithm, and each part contains one blank. The problem statements are in the [assignment handout](Original/Assignment.md), the distributed blank notebook is [Original/Lab03-LinearSystem.ipynb](Original/Lab03-LinearSystem.ipynb), and the notebook with the blanks filled in and run is [Solution/Lab03-LinearSystem_Completed.ipynb](Solution/Lab03-LinearSystem_Completed.ipynb).

| Problem | Notebook section | Blanks to fill |
|:-----|:-----|:-----|
| Problem 1 | 1.3 Spring-mass system | The right-hand side vector `b` |
| Problem 2 | 2.4 Implementing naive Gaussian elimination | `factor`, the row update, `x[n-1]`, and `x[i]` in the function `gauss_naive` |
| Problem 3 | 3.2 Partial pivoting | The row swap in the function `gauss_pivot` |
| Problem 4 | 4.2 Thomas algorithm | `factor`, `f[k]`, and `r[k]` in the function `thomas` |

*Table 1. The four blanks that make up the assignment*

The augmented matrix, the stiffness matrix of the spring-mass system, and naive Gaussian elimination are covered in the [Chapter 8 lecture note](../../Lecture-Notes/Ch08_Linear-Algebraic-Equations/Concepts_Lecture.md) and the [Chapter 8 lab note](../../Lecture-Notes/Ch08_Linear-Algebraic-Equations/Concepts_Lab.md), and partial pivoting and the Thomas algorithm are covered in the [Chapter 9 lecture note](../../Lecture-Notes/Ch09_Gauss-Elimination/Concepts_Lecture.md) and the [Chapter 9 lab note](../../Lecture-Notes/Ch09_Gauss-Elimination/Concepts_Lab.md). This document repeats only what the solutions need.

---

<br>

## 2. Requirements

The completed notebook must satisfy the following requirements. The rest of each cell (test cells, output format, plots) stays as distributed, and only the blanks are filled in.

1. **Problem 1:** `b` must be an array of length 3 holding the forces on the masses in the same order as the rows of the matrix `A`. The next cell computes the displacements with `np.linalg.solve(A, b)` and prints the residual.
2. **Problem 2:** Forward elimination must be a row update that zeros the entries below pivot column `k`, and back substitution must use the transformed right-hand side left in the last column (`nb-1`) of the augmented matrix to obtain the unknowns from the last one upward. The function must reproduce the solution $`(2, 3, -1)`$ of the $`3 \times 3`$ example traced by hand in the notebook and the `np.linalg.solve` result for the spring-mass system.
3. **Problem 3:** Rows `k` and `imax` must be swapped across the entire augmented matrix (coefficients and right-hand side). The function must give the same solution as `np.linalg.solve` on the example whose leading pivot is zero.
4. **Problem 4:** The three forward-sweep expressions must follow the formulas of Section 4.2 of the notebook, and in the heat conduction application that follows, the numerical solution must match the linear temperature profile.

---

<br>

## 3. Problem 1: Right-Hand Side of the Spring-Mass System

### 3.1 Approach

Row $`i`$ of the matrix `A` is the force-equilibrium equation of mass $`m_i`$, and the right-hand side of that equation is the external force $`F_i`$ acting on mass $`m_i`$. The right-hand side vector is therefore the external forces listed in the order of the masses, $`[F_1, F_2, F_3]`$. If the order were mixed up, a force would enter the equilibrium equation of another mass and a different system would be solved. The derivation of the stiffness matrix is covered in Section 2.11 of the [Chapter 8 lecture note](../../Lecture-Notes/Ch08_Linear-Algebraic-Equations/Concepts_Lecture.md) and Section 2 of the [Chapter 8 lab note](../../Lecture-Notes/Ch08_Linear-Algebraic-Equations/Concepts_Lab.md).

### 3.2 Code Written

```python
b = np.array([F1, F2, F3])
```

### 3.3 Results and Verification

```text
Coefficient matrix A:
[[ 70. -30.   0.]
 [-30.  50. -20.]
 [  0. -20.  30.]]

Right-hand side b: [20. 40. 30.]
```

```text
Displacements:
  x1 = 1.520000 m
  x2 = 2.880000 m
  x3 = 2.920000 m

Residual A*x - b = [ 0.00000000e+00 -2.84217094e-14 -1.42108547e-14]
Max |residual|   = 2.84e-14  (machine precision: ~2e-16)
```

Substituting the displacements gives $`70(1.52) - 30(2.88) = 20`$, $`-30(1.52) + 50(2.88) - 20(2.92) = 40`$, and $`-20(2.88) + 30(2.92) = 30`$, so all three equations hold. The largest residual, $`2.84 \times 10^{-14}`$, looks much larger than the machine epsilon $`2.2 \times 10^{-16}`$, but the right-hand side is of size about 40, so relative to it the residual is $`2.84 \times 10^{-14} / 40 \approx 7 \times 10^{-16}`$, only a few machine epsilons. A residual must be judged against the size of the right-hand side.

---

<br>

## 4. Problem 2: Naive Gaussian Elimination

### 4.1 Approach

**Forward elimination:** In pivot column `k`, the elimination factor for row `i` is computed from the entries of the augmented matrix already transformed by the earlier steps.

$$
f_{ik} = \frac{a_{ik}^{(k)}}{a_{kk}^{(k)}}
$$

Then $`f_{ik}`$ times row `k` is subtracted from row `i`, which zeros entry $`(i, k)`$. Columns `0` through `k-1` of row `i` already became zero in earlier steps, so the update only needs the slice `k:nb`, from column `k` to the last column `nb-1`. The last column holds the right-hand side, so the same update is applied to the right-hand side as well. It also matters that the factor is stored in the variable `factor` first, because the row update overwrites `Aug[i, k]`.

**Back substitution:** After forward elimination, the transformed right-hand side remains in the last column (`nb-1`) of the augmented matrix. The last row has only one unknown, so `x[n-1]` is obtained directly, and each row above it is solved by substituting the unknowns already found. The sum $`\sum_{j=i+1}^{n-1} a_{ij} x_j`$ of Section 2.3 of the notebook is computed in one step as the inner product `@` of the row vector `Aug[i, i+1:n]` and `x[i+1:n]`. Naive Gaussian elimination is covered in Section 6 of the [Chapter 8 lecture note](../../Lecture-Notes/Ch08_Linear-Algebraic-Equations/Concepts_Lecture.md) and Section 3 of the [Chapter 8 lab note](../../Lecture-Notes/Ch08_Linear-Algebraic-Equations/Concepts_Lab.md).

### 4.2 Code Written

```python
            factor = Aug[i, k] / Aug[k, k]
            Aug[i, k:nb] = Aug[i, k:nb] - factor * Aug[k, k:nb]
```

```python
    x[n-1] = Aug[n-1, nb-1] / Aug[n-1, n-1]
    for i in range(n - 2, -1, -1):
        x[i] = (Aug[i, nb-1] - Aug[i, i+1:n] @ x[i+1:n]) / Aug[i, i]
```

### 4.3 Results and Verification

```text
Upper triangular augmented matrix after elimination:
[[ 2.   1.  -1.   8. ]
 [ 0.   0.5  0.5  1. ]
 [ 0.   0.  -1.   1. ]]

Solution x  = [ 2.  3. -1.]
Expected    = [2, 3, -1]
Residual    = [0. 0. 0.]
```

```text
Spring-mass displacements:
  Gauss Elimination : [1.52 2.88 2.92]
  np.linalg.solve   : [1.52 2.88 2.92]
  Max difference    : 0.00e+00
```

The augmented matrix after elimination equals the hand trace of Section 2.4 of the notebook, $`[0, 1/2, 1/2 \mid 1]`$ and $`[0, 0, -1 \mid 1]`$, and the solution matches the true solution $`(2, 3, -1)`$. For the spring-mass system, the difference from `np.linalg.solve` is also printed as 0. In this system the pivot at each step (70 and $`50 - 30 \cdot 30/70 \approx 37.1`$) is already the largest in absolute value in its column, so naive elimination without row swaps solves it safely.

---

<br>

## 5. Problem 3: Partial Pivoting

### 5.1 Approach

First, consider what happens without pivoting. The cell of Section 3.1 of the notebook applies `gauss_naive` to the system with a zero leading pivot, and the completed notebook prints the following.

```text
True solution (via numpy):
  x = [ 0.5942029  -0.01449275  0.56521739]

Naive Gaussian Elimination result:
  x = [nan nan nan]
Augmented matrix after elimination:
[[  0.  -3.   7.   4.]
 [ nan  inf -inf -inf]
 [ nan  nan  nan  nan]]
```

The factor `Aug[1, 0] / Aug[0, 0] = 1/0` becomes `inf`, and the following operations produce `nan` from $`0 \cdot \infty`$ and $`\infty - \infty`$, so the whole solution becomes `nan`. The cell also prints a `RuntimeWarning` about the division by zero.

Partial pivoting swaps rows `k` and `imax` before elimination. In NumPy, selecting rows with an integer array such as `Aug[[k, imax]]` creates a copy of those rows, so `Aug[[k, imax]] = Aug[[imax, k]]` first builds the copy on the right and then assigns it to the two rows in reverse order, which swaps them safely. Because entire rows are swapped, the right-hand side (the last column) is swapped along with them. Partial pivoting is covered in Section 1 of the [Chapter 9 lecture note](../../Lecture-Notes/Ch09_Gauss-Elimination/Concepts_Lecture.md) and Section 1 of the [Chapter 9 lab note](../../Lecture-Notes/Ch09_Gauss-Elimination/Concepts_Lab.md).

> **Caution:** Selecting rows with single integers, as in `Aug[k], Aug[imax] = Aug[imax], Aug[k]`, puts views of the original array on the right instead of copies. The first assignment overwrites row `k` with the values of row `imax`, and the second writes the already changed row `k` back into row `imax`, so both rows end up as the original row `imax` and no swap takes place.

### 5.2 Code Written

```python
            Aug[[k, imax]] = Aug[[imax, k]]
```

### 5.3 Results and Verification

```text
Gauss with partial pivoting:
  x       = [ 0.5942029  -0.01449275  0.56521739]
  x (ref) = [ 0.5942029  -0.01449275  0.56521739]
  Max error: 0.00e+00

Upper triangular augmented matrix after pivoting:
[[ 5.  -2.   0.   3. ]
 [ 0.  -3.   7.   4. ]
 [ 0.   0.   4.6  2.6]]
```

Following the elimination by hand produces the upper triangular matrix in the output.

- **$`k = 0`$:** Among the values $`0, 1, 5`$ in column 0, 5 is the largest, so rows 0 and 2 are swapped. With the factor $`1/5 = 0.2`$, row 1 becomes $`[0, 2.4, -1 \mid -0.6]`$, and the original row 0, now row 2, already has 0 in column 0 and stays $`[0, -3, 7 \mid 4]`$.
- **$`k = 1`$:** Of the lower values $`2.4`$ and $`-3`$ in column 1, $`-3`$ has the larger absolute value and becomes the pivot, so the two rows are swapped again. With the factor $`2.4/(-3) = -0.8`$, row 2 becomes $`[0, 0, 4.6 \mid 2.6]`$.
- **Back substitution:** $`x_3 = 2.6/4.6 = 0.565217`$, $`x_2 = (4 - 7 x_3)/(-3) = -0.014493`$, and $`x_1 = (3 + 2 x_2)/5 = 0.594203`$, the same as the `np.linalg.solve` solution.

At $`k = 1`$ the swap happens even though the pivot 2.4 is not zero. Partial pivoting asks not whether the pivot is zero but whether it is the largest in absolute value in its column, so the elimination factors always satisfy $`|f_{ik}| \leq 1`$.

---

<br>

## 6. Problem 4: The Thomas Algorithm

### 6.1 Approach

In a tridiagonal matrix, the only nonzero entries of pivot row $`k-1`$ are the diagonal entry $`f_{k-1}`$ and the superdiagonal entry $`g_{k-1}`$. Gaussian elimination of the subdiagonal entry $`e_k`$ of row $`k`$ therefore only needs the factor $`e_k / f_{k-1}`$ and two updates, one to the diagonal entry and one to the right-hand side.

$$
f_k \leftarrow f_k - \frac{e_k}{f_{k-1}} g_{k-1}, \qquad r_k \leftarrow r_k - \frac{e_k}{f_{k-1}} r_{k-1}
$$

Because the loop runs in order of $`k`$, `f[k-1]` and `r[k-1]` are the values already updated in the previous iteration, which correspond to $`a_{kk}^{(k)}`$ in Gaussian elimination. The function first copies `f` and `r` so that it does not change the caller's arrays. The Thomas algorithm is covered in Section 5 of the [Chapter 9 lecture note](../../Lecture-Notes/Ch09_Gauss-Elimination/Concepts_Lecture.md) and Section 3 of the [Chapter 9 lab note](../../Lecture-Notes/Ch09_Gauss-Elimination/Concepts_Lab.md).

### 6.2 Code Written

```python
        factor = e[k] / f[k-1]
        f[k] = f[k] - factor * g[k-1]
        r[k] = r[k] - factor * r[k-1]
```

### 6.3 Results and Verification

The heat conduction cell only draws a plot, so the same function was evaluated for three interior nodes, with the following result.

| Node | $`x`$ (m) | Numerical $`T`$ (°C) | Linear profile $`50 + 100x`$ (°C) |
|:-----|:-----|:-----|:-----|
| 1 | 0.25 | 75.0 | 75.0 |
| 2 | 0.50 | 100.0 | 100.0 |
| 3 | 0.75 | 125.0 | 125.0 |

*Table 2. Temperatures of steady-state 1D heat conduction solved with three interior nodes*

With six and fifteen interior nodes, the difference from the linear profile is also at most $`3 \times 10^{-14}`$. A linear function makes the second difference $`T_{i-1} - 2T_i + T_{i+1}`$ exactly zero, so in this problem the nodal values must equal the exact solution with no discretization error, and they do. In addition, this matrix is not strictly diagonally dominant because $`|f_i| = |e_i| + |g_i|`$ in the interior rows, but the diagonal entries after the forward sweep are $`f_k = -(k+2)/(k+1)`$, that is, $`-2, -1.5, -1.333, \ldots`$, which never become zero, so the system can be solved without pivoting.

The cell after Section 4.3 of the notebook solves the same tridiagonal system with the Thomas algorithm and with `np.linalg.solve` and compares the times. The completed notebook prints the following.

```text
     n    Thomas (s)    Full Gauss (s)     Speedup
----------------------------------------------------
   100      0.000071          0.000094         1.3x
   500      0.000403          0.001319         3.3x
  2000      0.001310          0.032362        24.7x
  5000      0.003320          0.553740       166.8x
```

The times depend on the machine, but the gap widens as $`n`$ grows. The time of the Thomas algorithm grows roughly in proportion to $`n`$, while the dense solve is $`O(n^3)`$ and grows much faster.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Right-hand side and residual | The right-hand side vector holds the external force of each row's equation in the same order, and a residual must be compared with the size of the right-hand side to judge whether it is at the level of machine epsilon. |
| Forward elimination | The factor $`f_{ik} = a_{ik}^{(k)} / a_{kk}^{(k)}`$ is stored first, and then $`f_{ik}`$ times row `k` is subtracted from row `i` from column `k` to the end, with the right-hand side receiving the same operation. |
| Back substitution | The unknowns are found from the last one upward, and the sum with the unknowns already found is computed as the inner product of a row slice and the solution vector. |
| System with a zero pivot | Naive elimination divides by zero, `inf` and `nan` spread, and no solution is obtained. |
| Partial pivoting | Each step moves the row with the largest absolute value in the column up to the pivot row, and NumPy swaps the rows safely through the copy made by integer-array indexing. |
| Thomas algorithm | In the tridiagonal structure each step only updates the diagonal entry and the right-hand side, so the system is solved in $`O(n)`$, and the linear temperature profile is reproduced exactly. |

---

<br>

## Self-Check Questions

1. **Why the factor is stored first:** In Problem 2, would the result change if the row update were written in one line without storing the factor, as in `Aug[i, k:nb] = Aug[i, k:nb] - Aug[i, k] / Aug[k, k] * Aug[k, k:nb]`?

   <details>
   <summary>Answer</summary>

   Not in this one line. NumPy evaluates the entire right-hand side before assigning to the left, so the value of `Aug[i, k]` before it is overwritten is used. In an implementation that loops over the columns and updates `Aug[i, j]` one at a time, however, `Aug[i, k]` becomes 0 at `j = k`, the factor for the remaining columns becomes 0, and the elimination goes wrong. Storing the factor first is safe regardless of the implementation.

   </details>

2. **System with a zero pivot:** Explain how naive elimination produces `nan` on the example of Problem 3.

   <details>
   <summary>Answer</summary>

   At $`k = 0`$ the factor for row 1 becomes $`1/0 = \infty`$, and the update of row 1 produces `nan` from operations such as $`1 - \infty \cdot 0`$ and `inf` in the other columns. These values spread through the later elimination and back substitution, so every component of the solution becomes `nan`.

   </details>

3. **The row-swap trap:** Why does `Aug[k], Aug[imax] = Aug[imax], Aug[k]` fail to swap the two rows?

   <details>
   <summary>Answer</summary>

   `Aug[imax]` and `Aug[k]` selected with single integers are views of the original array. After the first assignment makes row `k` equal to row `imax`, the second assignment writes the changed row `k` into row `imax`, so both rows become the original row `imax`. `Aug[[k, imax]] = Aug[[imax, k]]` does not have this problem because the right-hand side is a copy.

   </details>

4. **Pivots of the Thomas algorithm:** The heat conduction matrix of Problem 4 is not strictly diagonally dominant, yet it is solved without pivoting. Explain why using the diagonal entries after the forward sweep.

   <details>
   <summary>Answer</summary>

   Starting from $`f_0 = -2`$ and repeating $`f_k = -2 - 1/f_{k-1}`$ gives $`f_k = -(k+2)/(k+1)`$, that is, $`-2, -1.5, -1.333, \ldots`$, which approaches $`-1`$ but never becomes zero. The pivots used in the divisions therefore never vanish.

   </details>

5. **Why there is no discretization error:** Why does the numerical solution of the heat conduction problem in Problem 4 equal the exact solution regardless of the number of interior nodes?

   <details>
   <summary>Answer</summary>

   The exact solution $`T = 50 + 100x`$ is linear, so it satisfies $`T_{i-1} - 2T_i + T_{i+1} = 0`$ exactly at any three equally spaced points. The discretized equations therefore have the nodal values of the exact solution as their solution, and only roundoff-level differences remain.

   </details>
