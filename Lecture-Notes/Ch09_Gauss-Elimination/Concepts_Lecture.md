# Chapter 9 — Gauss Elimination: Improved Methods

[Repository Guide](../../README.md) | [English](Concepts_Lecture.md) | [한국어](Concepts_Lecture.ko.md)

> **Last Updated:** 2026-11-19
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 9

> **Learning Objectives:**
> 1. Perform Gaussian elimination with partial pivoting
> 2. Analyze computational complexity of elimination
> 3. Explain the principle and use of scaled partial pivoting and the Gauss-Jordan method
> 4. Solve tridiagonal systems with the Thomas algorithm
> 5. Compute a determinant by Gauss elimination and account for the sign of row interchanges
> 6. Derive the operation count of Gauss elimination with summation formulas and assess the extra cost of pivoting

---

<br>

## Table of Contents

- [1. Partial Pivoting](#1-partial-pivoting)
  - [1.1 The Problem with Naive Gauss Elimination](#11-the-problem-with-naive-gauss-elimination)
  - [1.2 Partial Pivoting Strategy](#12-partial-pivoting-strategy)
  - [1.3 Modified Algorithm Pseudocode](#13-modified-algorithm-pseudocode)
  - [1.4 Step-by-Step Trace of Partial Pivoting](#14-step-by-step-trace-of-partial-pivoting)
  - [1.5 Pivoting and Error Amplification](#15-pivoting-and-error-amplification)
- [2. Scaled Partial Pivoting](#2-scaled-partial-pivoting)
- [3. Example Demonstrating Pivoting Necessity](#3-example-demonstrating-pivoting-necessity)
  - [3.1 Without Pivoting (4-Digit Arithmetic)](#31-without-pivoting-4-digit-arithmetic)
  - [3.2 With Pivoting](#32-with-pivoting)
  - [3.3 Results for Different Numbers of Significant Digits](#33-results-for-different-numbers-of-significant-digits)
- [4. Gauss-Jordan Method](#4-gauss-jordan-method)
- [5. Tridiagonal Systems: Thomas Algorithm](#5-tridiagonal-systems-thomas-algorithm)
  - [5.1 Tridiagonal Matrix Structure](#51-tridiagonal-matrix-structure)
  - [5.2 Thomas Algorithm](#52-thomas-algorithm)
  - [5.3 Computational Cost](#53-computational-cost)
  - [5.4 Normalized Form and Stability Condition](#54-normalized-form-and-stability-condition)
  - [5.5 Thomas Algorithm Example](#55-thomas-algorithm-example)
  - [5.6 Tridiagonal Systems from a Heat Conduction Problem](#56-tridiagonal-systems-from-a-heat-conduction-problem)
- [6. Python Implementation](#6-python-implementation)
  - [6.1 Implementing Gauss Elimination with Partial Pivoting](#61-implementing-gauss-elimination-with-partial-pivoting)
  - [6.2 Implementing the Thomas Algorithm and the Banded Solver](#62-implementing-the-thomas-algorithm-and-the-banded-solver)
- [7. Determinant Evaluation with Gauss Elimination](#7-determinant-evaluation-with-gauss-elimination)
  - [7.1 Determinant of a Triangular Matrix and Row Operations](#71-determinant-of-a-triangular-matrix-and-row-operations)
  - [7.2 Determinant Examples](#72-determinant-examples)
  - [7.3 Cost Compared with Cofactor Expansion](#73-cost-compared-with-cofactor-expansion)
- [8. Deriving the Operation Count](#8-deriving-the-operation-count)
  - [8.1 Operation Count of Forward Elimination](#81-operation-count-of-forward-elimination)
  - [8.2 Operation Count of Back Substitution](#82-operation-count-of-back-substitution)
  - [8.3 Total Cost and the Extra Cost of Pivoting](#83-total-cost-and-the-extra-cost-of-pivoting)
- [9. Common Mistakes](#9-common-mistakes)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)
- [References](#references)

---

<br>

## 1. Partial Pivoting

### 1.1 The Problem with Naive Gauss Elimination

In naive Gauss elimination, the pivot element $`a_{ii}`$ is used directly as the divisor during each elimination step. This causes two critical problems:

1. **Division by zero**: If $`a_{ii} = 0`$, the algorithm fails outright.
2. **Large roundoff errors**: If $`a_{ii}`$ is very small (close to zero), the multiplier $`a_{ki}/a_{ii}`$ becomes very large, amplifying every roundoff error in subsequent operations by that factor.

Both problems stem from the same root cause: the pivot element is too small relative to other elements in the column.

### 1.2 Partial Pivoting Strategy

The solution is straightforward: before performing elimination on column $`i`$, search the rows from $`i`$ to $`n`$ to find the row with the **largest absolute value** in column $`i`$, and swap that row with row $`i`$. This ensures the pivot is always the largest available element in the column, which minimizes the magnitude of the multipliers.

Since the multiplier is computed as $`\text{factor} = a_{ki}/a_{ii}`$, making $`|a_{ii}|`$ as large as possible guarantees that all factors satisfy $`|\text{factor}| \leq 1`$. This bounds the error amplification at each step.

> **Definition:** **Partial pivoting** is the strategy of finding, at each step $`k`$ of Gauss elimination, the row among rows $`k`$ through $`n`$ whose entry in column $`k`$ has the largest absolute value, interchanging it with row $`k`$, and then eliminating. Only rows are interchanged and the column order is unchanged, so the order of the unknowns is preserved.

### 1.3 Modified Algorithm Pseudocode

The only modification to naive Gauss elimination is the addition of a pivot search and row swap before each elimination step:

```text
for i = 1 to n-1:
    // Find max |a[k][i]| for k = i, ..., n
    imax = argmax(|a[k][i]|) for k = i..n
    swap row i with row imax
    // Proceed with standard elimination
    for k = i+1 to n:
        factor = a[k][i] / a[i][i]
        for j = i to n+1:
            a[k][j] = a[k][j] - factor * a[i][j]
```

> **Note:** Row swapping is an $`O(n)`$ operation per elimination step, so the total overhead of pivoting is $`O(n^2)`$, which is negligible compared to the $`O(n^3)`$ cost of the elimination itself. There is essentially no reason to skip pivoting.

### 1.4 Step-by-Step Trace of Partial Pivoting

The following system has no $`x_1`$ in its first equation, so $`a_{11} = 0`$. The system itself has a unique solution, yet naive Gauss elimination divides by zero the moment it computes the first multiplier $`a_{21}/a_{11}`$, and it stops.

```math
\left[\begin{array}{ccc|c} 0 & 2 & 3 & 8 \\ 4 & 6 & 7 & -3 \\ 2 & -3 & 6 & 5 \end{array}\right]
```

**Step 1 (column 1).** The candidates in column 1 have absolute values $`|0|`$, $`|4|`$, and $`|2|`$, so row 2, which holds the largest value 4, is the pivot row. After interchanging rows 1 and 2, the multipliers are $`l_{21} = 0/4 = 0`$ and $`l_{31} = 2/4 = 0.5`$. The new row 2 (the original row 1) already has a zero in column 1 and does not change. Row 3 becomes $`R_3 - 0.5R_1 = (2 - 2,\ -3 - 3,\ 6 - 3.5,\ 5 + 1.5) = (0,\ -6,\ 2.5,\ 6.5)`$.

```math
\left[\begin{array}{ccc|c} 4 & 6 & 7 & -3 \\ 0 & 2 & 3 & 8 \\ 0 & -6 & 2.5 & 6.5 \end{array}\right]
```

**Step 2 (column 2).** The candidates in column 2 from row 2 down are $`|2|`$ and $`|-6|`$, so row 3 becomes the pivot row and rows 2 and 3 are interchanged. The multiplier is $`l_{32} = 2/(-6) = -1/3`$, and computing $`R_3 - (-1/3)R_2`$ gives $`2 + (1/3)(-6) = 0`$ in column 2, $`3 + (1/3)(2.5) = 23/6 \approx 3.8333`$ in column 3, and $`8 + (1/3)(6.5) = 61/6 \approx 10.1667`$ on the right-hand side.

```math
\left[\begin{array}{ccc|c} 4 & 6 & 7 & -3 \\ 0 & -6 & 2.5 & 6.5 \\ 0 & 0 & 23/6 & 61/6 \end{array}\right]
```

| Step | Absolute values of candidates | Pivot row | Row interchange | Multipliers |
|:-----|:-----|:-----|:-----|:-----|
| 1 | 0, 4, 2 | Row 2 | Rows 1 and 2 | $`l_{21} = 0`$, $`l_{31} = 0.5`$ |
| 2 | 2, 6 | Row 3 | Rows 2 and 3 | $`l_{32} = -1/3`$ |

*Table 1. Pivots and multipliers chosen by partial pivoting in the 3x3 example*

**Back substitution.** The upper triangular system is solved from the bottom up.

```math
x_3 = \frac{61/6}{23/6} = \frac{61}{23} \approx 2.652174
```

```math
x_2 = \frac{6.5 - 2.5x_3}{-6} = \frac{6.5 - 6.630435}{-6} = \frac{1}{46} \approx 0.021739
```

```math
x_1 = \frac{-3 - 6x_2 - 7x_3}{4} = \frac{-3 - 0.130435 - 18.565217}{4} = -\frac{499}{92} \approx -5.423913
```

It is worth noting that every multiplier has absolute value at most 1. Because the pivot has the largest absolute value among the candidates in its column, $`|l_{ik}| = |a_{ik}|/|a_{kk}| \leq 1`$ always holds. Two row interchanges occurred in this example; this count is used in Section 7 to fix the sign of the determinant.

### 1.5 Pivoting and Error Amplification

The $`k`$th step of elimination changes each entry to $`a_{ij}^{(k)} = a_{ij}^{(k-1)} - l_{ik}a_{kj}^{(k-1)}`$. When partial pivoting guarantees $`|l_{ik}| \leq 1`$, the following inequality holds.

```math
|a_{ij}^{(k)}| \leq |a_{ij}^{(k-1)}| + |l_{ik}|\,|a_{kj}^{(k-1)}| \leq 2\max_{p,q}|a_{pq}^{(k-1)}|
```

That is, each step can at most double the largest entry of the matrix. How much the entries grow during elimination is measured by the following ratio.

> **Definition:** The **growth factor** is the largest absolute value of any entry that appears in any step of elimination divided by the largest absolute value of the entries of the original matrix, $`\rho = \max_{i,j,k}|a_{ij}^{(k)}| / \max_{i,j}|a_{ij}|`$. With partial pivoting, the inequality above gives $`\rho \leq 2^{n-1}`$.

The growth factor matters because the size of roundoff errors is proportional to the size of the numbers handled during elimination. One floating-point operation has a relative error of about the machine epsilon $`\varepsilon_{\text{mach}}`$, so when the entries grow by a factor of $`\rho`$, the absolute error left by that operation grows by the same factor. In terms of backward error analysis, the solution $`\hat{x}`$ computed by Gauss elimination is the exact solution of a slightly changed system $`(A + E)\hat{x} = b`$, and the perturbation is bounded by $`\lVert E \rVert \leq c(n)\,\rho\,\varepsilon_{\text{mach}}\lVert A \rVert`$ (with $`c(n)`$ a constant that grows like a low-degree polynomial in $`n`$). When $`\rho`$ is small, the computed solution is the exact solution of a problem almost identical to the original one, and the algorithm is called stable.

Applying this measure to the example of Section 3 makes the difference clear. Without pivoting, $`a_{22}`$ becomes $`1 - 3333.3 \times 3 \approx -9999`$, about 3333 times the original largest entry 3 ($`\rho \approx 3333`$). After the row interchange, $`a_{22} = 3 - 0.0003 \times 1 = 2.9997`$ and the entries do not grow ($`\rho = 1`$). The following code computes the growth factor in both cases.

```python
import numpy as np

def growth_factor(A, pivot=True):
    U = np.array(A, dtype=float)
    n = U.shape[0]
    amax0 = np.abs(U).max()
    amax = amax0
    for k in range(n - 1):
        if pivot:
            p = k + np.argmax(np.abs(U[k:, k]))
            U[[k, p]] = U[[p, k]]
        for i in range(k + 1, n):
            U[i, k:] -= (U[i, k] / U[k, k]) * U[k, k:]
        amax = max(amax, np.abs(U).max())
    return amax / amax0

B = [[0.0003, 3.0], [1.0, 1.0]]
print(growth_factor(B, pivot=False), growth_factor(B, pivot=True))
```

```text
3333.0 1.0
```

The function remembers the largest absolute value `amax0` of the original matrix and updates `amax` with the largest absolute value of the whole matrix after each elimination step. With `pivot=True` it finds the pivot row with `np.argmax` and swaps it in; with `False` it eliminates without swaps. At the end it returns the ratio of the two maxima.

> **Note:** **Worst-case growth and complete pivoting.** For the $`n \times n`$ matrix with ones on the diagonal, $`-1`$ everywhere below the diagonal, and ones in the last column, partial pivoting performs no interchanges and the last column doubles at every step, reaching $`\rho = 2^{n-1}`$ (16 for $`n = 5`$). Such matrices almost never arise in applications, which is why partial pivoting is the standard. Complete pivoting, which interchanges both rows and columns to choose the largest entry of the whole remaining submatrix as the pivot, restrains growth more strongly, but it needs $`(n-k+1)^2`$ comparisons per step and must track the order of the unknowns, so it is rarely used.

> **Caution:** Pivoting only prevents the amplification of errors made by the algorithm; it cannot change the sensitivity of the problem itself. If the coefficient matrix is ill-conditioned, the solution can be inaccurate even with pivoting, and this sensitivity is assessed with the condition number of Chapter 11.

---

<br>

## 2. Scaled Partial Pivoting

Sometimes partial pivoting alone is insufficient. This happens when the rows of the matrix have **very different magnitudes**. A row with all large entries might dominate the pivot selection even when it is not the best choice relative to its own scale.

**Scaled partial pivoting** addresses this by normalizing each row before comparing pivot candidates:

1. Compute the **scale factor** for each row: $`s_i = \max_j |a_{ij}|`$.
2. When selecting the pivot for column $`i`$, choose the row $`k`$ that maximizes $`|a_{ki}| / s_k`$ (the ratio of the candidate pivot to the row's scale factor).
3. Swap and proceed with standard elimination.

The scale factors $`s_i`$ are computed once at the beginning; when rows are swapped they move with their rows and are not recomputed. This normalization ensures that the pivot selection accounts for the overall magnitude of each row, preventing a row with uniformly large entries from always being chosen.

> **Note:** **Implicit scaling.** Scaled partial pivoting is equivalent to implicitly scaling the rows of the matrix before applying partial pivoting. The scale factors are used only to choose the pivot, and the elimination itself uses the original coefficients. The key insight is that the best pivot is not necessarily the largest in absolute value, but the largest relative to the other entries in its row.

**Example.** The following system is solved with 3-significant-digit arithmetic. The exact solution is $`x_1 = 1.00002`$, $`x_2 = 0.99998`$.

```math
\left[\begin{array}{cc|c} 2 & 100000 & 100000 \\ 1 & 1 & 2 \end{array}\right]
```

Partial pivoting keeps row 1 as the pivot row because $`|2| > |1|`$ in column 1. The multiplier is $`1/2 = 0.5`$, $`a_{22} = 1 - 0.5 \times 100000 = -49999`$ rounds to $`-50000`$ with three digits, and $`b_2 = 2 - 50000 = -49998`$ also becomes $`-50000`$, so $`x_2 = 1.00`$. Back substitution then gives $`x_1 = (100000 - 100000 \times 1.00)/2 = 0`$, which is completely wrong.

Scaled partial pivoting compares the ratios $`2/100000 = 0.00002`$ and $`1/1 = 1`$ with $`s_1 = 100000`$ and $`s_2 = 1`$ and selects row 2 as the pivot row. After the interchange the multiplier is $`2/1 = 2`$, $`a_{22} = 100000 - 2 = 99998`$ rounds to $`100000`$ with three digits, and $`b_2 = 100000 - 4 = 99996`$ also becomes $`100000`$, so $`x_2 = 1.00`$. Back substitution gives $`x_1 = (2 - 1.00)/1 = 1.00`$, which is accurate to three digits.

The first case fails for the same reason as the example of Section 3. Dividing row 1 by 100000 gives $`0.00002x_1 + x_2 = 1`$, so within its own row the pivot 2 is even smaller than the pivot 0.0003 of Section 3. Partial pivoting compares only the absolute sizes of the entries and misses this, whereas the comparison divided by the scale factors reveals it.

---

<br>

## 3. Example Demonstrating Pivoting Necessity

Consider the system:

```math
0.0003x_1 + 3.0000x_2 = 2.0001
```

```math
1.0000x_1 + 1.0000x_2 = 1.0000
```

The exact solution is $`x_1 = 1/3`$ and $`x_2 = 2/3`$.

### 3.1 Without Pivoting (4-Digit Arithmetic)

Using 4-digit arithmetic with the equations in the given order:

- The pivot is $`a_{11} = 0.0003`$, which is very small.
- The factor is $`1.0000 / 0.0003 = 3333.3`$, which is very large.
- After elimination: $`-9999x_2 = -6666`$.
- Back substitution gives $`x_2 = 0.6667`$, then $`x_1 = (2.0001 - 3.0000 \times 0.6667) / 0.0003 = (2.0001 - 2.0001) / 0.0003 = 0.0000`$.

$`x_2`$ is accurate, but with limited precision $`x_1`$ is **wrong**: it comes out as 0 instead of the true value $`1/3`$, a 100% relative error. The product $`3.0000 \times 0.6667`$ nearly equals 2.0001, so the subtraction cancels all significant digits (catastrophic cancellation), and the result is divided by the tiny pivot 0.0003. The large factor of 3333.3 magnifies the small perturbations introduced by 4-digit rounding to the point where significant digits are lost.

### 3.2 With Pivoting

Swap the two rows first (since $`|1.0000| > |0.0003|`$ in column 1):

```math
1.0000x_1 + 1.0000x_2 = 1.0000
```

```math
0.0003x_1 + 3.0000x_2 = 2.0001
```

Now:

- The pivot is $`a_{11} = 1.0000`$, which is well-scaled.
- The factor is $`0.0003 / 1.0000 = 0.0003`$, which is small.
- After elimination: $`2.9997x_2 = 1.9998`$.
- Back substitution gives $`x_2 = 0.6667`$, then $`x_1 = (1.0000 - 1.0000 \times 0.6667) / 1.0000 = 0.3333`$.

The result is much more accurate because the small factor (0.0003 instead of 3333.3) does not amplify roundoff errors.

> **Note:** This example dramatically shows why pivoting matters. A factor of 3333.3 amplifies every roundoff error by that amount. Pivoting keeps factors $`\leq 1`$, controlling error growth. In practice, partial pivoting is always used, and naive Gauss elimination is taught only for pedagogical purposes.

### 3.3 Results for Different Numbers of Significant Digits

The reason why $`x_2`$ was accurate in Section 3.1 while only $`x_1`$ was badly wrong appears in the back-substitution formula. Without pivoting, $`x_1`$ is computed as follows.

```math
x_1 = \frac{2.0001 - 3x_2}{0.0003}
```

If $`x_2`$ has an error $`\delta`$, the numerator changes by $`3\delta`$, and dividing it by 0.0003 makes the error in $`x_1`$ equal to $`3\delta/0.0003 = 10000\delta`$. After the row interchange, $`x_1 = (1 - x_2)/1`$, so an error $`\delta`$ in $`x_2`$ passes into $`x_1`$ only as $`\delta`$. The table below rounds $`x_2 = 2/3`$ to $`k`$ significant digits and then computes $`x_1`$ with each formula.

| Significant digits | $`x_2`$ | $`x_1`$ (no pivoting) | Relative error (%) | $`x_1`$ (pivoting) | Relative error (%) |
|:-----|:-----|:-----|-----:|:-----|-----:|
| 3 | 0.667 | -3 | 1000 | 0.333 | 0.1 |
| 4 | 0.6667 | 0 | 100 | 0.3333 | 0.01 |
| 5 | 0.66667 | 0.3 | 10 | 0.33333 | 0.001 |
| 6 | 0.666667 | 0.33 | 1 | 0.333333 | 0.0001 |
| 7 | 0.6666667 | 0.333 | 0.1 | 0.3333333 | 0.00001 |

*Table 2. Solution and relative error of the 0.0003 example for different numbers of significant digits*

The value of $`x_1`$ computed without pivoting improves tenfold with each extra digit, yet it is always 10000 times worse than the pivoted result with the same number of digits. A small pivot therefore throws away about four significant digits ($`\log_{10} 10000 = 4`$). In double precision (about 16 significant digits) the loss in this example is not noticeable, but with a pivot as small as $`10^{-12}`$ the same mechanism destroys most of the digits.

---

<br>

## 4. Gauss-Jordan Method

The **Gauss-Jordan method** extends standard Gauss elimination by eliminating entries **both above and below** the diagonal, producing a **diagonal matrix** rather than an upper triangular one.

The procedure is:

1. **Forward elimination**: Eliminate below the diagonal, as in standard Gauss elimination.
2. **Backward elimination**: Additionally eliminate entries **above** the diagonal.
3. **Normalization**: Divide each row by its diagonal element to produce the identity matrix on the left.

The result is that the augmented matrix $`[A|b]`$ is transformed into $`[I|x]`$, where $`x`$ is the solution. No back substitution is needed because the matrix is already diagonal (or identity after normalization).

**Key properties:**

- **No back substitution needed**: The solution is read directly from the transformed right-hand side.
- **Useful for computing the matrix inverse**: Applying Gauss-Jordan to $`[A|I]`$ yields $`[I|A^{-1}]`$.
- **More expensive**: The operation count is $`O(n^3)`$ with a larger constant than standard Gauss elimination, approximately $`n^3/2`$ multiplications versus $`n^3/3`$ for standard Gauss with back substitution.

> **Note:** **Matrix inversion and cost.** The Gauss-Jordan method is the standard textbook algorithm for computing matrix inverses. While it is about 50% more expensive than Gauss elimination with back substitution for solving a single system, it becomes competitive when solving for multiple right-hand sides simultaneously (as in computing the inverse).

**Example.** The following system is solved by the Gauss-Jordan method. In this example each step first divides the pivot row by the pivot to make the pivot 1 and then eliminates every other entry in that column.

```math
\left[\begin{array}{ccc|c} 2 & 1 & 1 & 3 \\ 4 & 3 & 3 & 7 \\ 8 & 7 & 9 & 19 \end{array}\right]
```

In column 1, row 1 is divided by 2 to give $`(1,\ 1/2,\ 1/2,\ 3/2)`$, and $`R_2 - 4R_1 = (0,\ 1,\ 1,\ 1)`$ and $`R_3 - 8R_1 = (0,\ 3,\ 5,\ 7)`$ are computed. In column 2 the pivot is already 1, so no division is needed, and row 2 is used to eliminate both row 1 above and row 3 below, giving $`R_1 - (1/2)R_2 = (1,\ 0,\ 0,\ 1)`$ and $`R_3 - 3R_2 = (0,\ 0,\ 2,\ 4)`$. In column 3, row 3 is divided by 2 to give $`(0,\ 0,\ 1,\ 2)`$, and $`R_2 - R_3 = (0,\ 1,\ 0,\ -1)`$ is computed. The column-3 entry of row 1 is already zero.

```math
\left[\begin{array}{ccc|c} 1 & 0 & 0 & 1 \\ 0 & 1 & 0 & -1 \\ 0 & 0 & 1 & 2 \end{array}\right]
```

The left side is now the identity matrix, so the solution $`x = (1,\ -1,\ 2)`$ is read directly from the right column.

**Operation count.** At step $`k`$, all $`n-1`$ rows other than the pivot row are eliminated. Each row needs one division for the multiplier and $`n-k+1`$ multiplications to update columns $`k+1`$ through $`n`$ and the right-hand side, so each row costs $`n-k+2`$ multiplications and divisions. Summing over all steps gives the following.

```math
\sum_{k=1}^{n}(n-1)(n-k+2) = (n-1)\,\frac{n^2+3n}{2} = \frac{n^3}{2} + O(n^2)
```

Gauss elimination eliminates only the $`n-k`$ rows below the diagonal, so it needs $`n^3/3 + O(n^2)`$ multiplications and divisions, whereas the Gauss-Jordan method also eliminates the rows above the diagonal and thus handles a constant $`n-1`$ rows per step. This difference produces the cost ratio of about 1.5 mentioned above.

---

<br>

## 5. Tridiagonal Systems: Thomas Algorithm

### 5.1 Tridiagonal Matrix Structure

A **tridiagonal matrix** has non-zero entries only on the main diagonal, the diagonal immediately above it (super-diagonal), and the diagonal immediately below it (sub-diagonal). All other entries are zero:

```math
\begin{bmatrix} f_1 & g_1 & & \\ e_2 & f_2 & g_2 & \\ & \ddots & \ddots & \ddots \\ & & e_n & f_n \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{bmatrix} = \begin{bmatrix} r_1 \\ r_2 \\ \vdots \\ r_n \end{bmatrix}
```

where $`e_i`$ is the sub-diagonal, $`f_i`$ is the main diagonal, and $`g_i`$ is the super-diagonal.

Tridiagonal systems are extremely common in scientific computing:

- **Discretized ODEs**: They arise from central finite differences for second-order boundary value problems.
- **Spline interpolation**: The system for cubic spline coefficients is tridiagonal.
- **Implicit time-stepping**: They arise in implicit methods for parabolic PDEs (e.g., the heat equation).

### 5.2 Thomas Algorithm

The **Thomas algorithm** is a specialized form of Gauss elimination that exploits the tridiagonal structure. It consists of two phases:

**Forward sweep (elimination):** For $`i = 2`$ to $`n`$:

```math
\text{factor} = \frac{e_i}{f_{i-1}}
```

```math
f_i = f_i - \text{factor} \cdot g_{i-1}
```

```math
r_i = r_i - \text{factor} \cdot r_{i-1}
```

**Back substitution:**

```math
x_n = \frac{r_n}{f_n}
```

For $`i = n-1`$ down to $`1`$:

```math
x_i = \frac{r_i - g_i \cdot x_{i+1}}{f_i}
```

Each step of the forward sweep only touches the current row and the row above it, and each step of back substitution only requires the next solution value. The algorithm only accesses the three diagonals and never fills in entries outside the band. Because it does not pivot, however, it is reliable when the $`f_i`$ cannot become nearly zero during elimination, as for diagonally dominant matrices.

### 5.3 Computational Cost

The Thomas algorithm requires $`3(n-1)`$ multiplications/divisions and $`2(n-1)`$ additions/subtractions in the forward sweep, plus $`2n-1`$ multiplications/divisions and $`n-1`$ subtractions in back substitution. This totals $`8n-7`$ operations, giving a cost of **$`O(n)`$**.

Counting the operations formula by formula gives the numbers above. In one row of the forward sweep, computing the multiplier $`e_i/f_{i-1}`$ is one division, updating $`f_i`$ is one multiplication and one subtraction, and updating $`r_i`$ is also one multiplication and one subtraction, so each row costs 5 operations and the $`n-1`$ rows cost $`5(n-1)`$. Back substitution needs one division for $`x_n`$ and one multiplication, one subtraction, and one division for each of the remaining $`n-1`$ values $`x_i`$, which is $`3(n-1) + 1 = 3n - 2`$ operations. The total is $`5(n-1) + 3n - 2 = 8n - 7`$.

> **Note:** **The count $`9n - 10`$.** If the multiplier $`e_i/f_{i-1}`$ is computed separately in the $`f_i`$ formula and in the $`r_i`$ formula, the forward sweep costs 6 operations per row, $`6(n-1)`$ in total. If back substitution then counts $`x_1`$ as a single division, $`x_1 = r_1/f_1`$, it costs $`3(n-2) + 2 = 3n - 4`$ operations and the total becomes $`9n - 10`$. However, $`x_1 = (r_1 - g_1x_2)/f_1`$, so $`x_1`$ also needs a multiplication and a subtraction, and with this correction the total is $`9n - 8`$. Computing the multiplier once and reusing it in both formulas reduces the count to $`8n - 7`$. Whichever way the operations are counted, the conclusion that the cost is $`O(n)`$ is the same.

This is a dramatic improvement over the $`O(n^3)`$ cost of general Gauss elimination. For a system of size $`n = 100{,}000`$, the ratio of operation counts is $`(2n^3/3)/(8n) = n^2/12 \approx 8 \times 10^8`$, so the Thomas algorithm is about $`10^9`$ times faster (comparing orders alone gives $`n^2 = 10^{10}`$).

> **[Data Structures]** Sparse matrix storage formats save memory by storing only the nonzero entries together with their positions and skip operations on the zero entries. The Thomas algorithm follows the same idea by storing and processing only the three diagonals as the vectors $`e`$, $`f`$, and $`g`$; since it handles about $`3n`$ values instead of $`n^2`$, both memory and work are $`O(n)`$.

### 5.4 Normalized Form and Stability Condition

After the forward sweep, dividing each row by its updated diagonal entry $`f_i'`$ gives a form with ones on the diagonal. The quotients are written as follows.

```math
h_i = \frac{g_i}{f_i'}, \qquad s_i = \frac{r_i'}{f_i'}
```

The system then becomes $`x_i + h_i x_{i+1} = s_i`$ ($`i < n`$) and $`x_n = s_n`$. In the update formulas $`f_i' = f_i - (e_i/f_{i-1}')g_{i-1}`$ and $`r_i' = r_i - (e_i/f_{i-1}')r_{i-1}'`$, we have $`g_{i-1}/f_{i-1}' = h_{i-1}`$ and $`r_{i-1}'/f_{i-1}' = s_{i-1}`$, so $`f_i' = f_i - e_ih_{i-1}`$ and $`r_i' = r_i - e_is_{i-1}`$. Hence $`h_i`$ and $`s_i`$ are computed in turn from the previous values alone.

```math
h_1 = \frac{g_1}{f_1}, \quad s_1 = \frac{r_1}{f_1}, \qquad h_i = \frac{g_i}{f_i - e_ih_{i-1}}, \quad s_i = \frac{r_i - e_is_{i-1}}{f_i - e_ih_{i-1}} \quad (i = 2, \ldots, n)
```

$`h_n`$ is not used. Back substitution finishes without any division: $`x_n = s_n`$ and $`x_i = s_i - h_ix_{i+1}`$ ($`i = n-1, \ldots, 1`$).

This form makes clear where the algorithm can fail. When the denominator $`d_i = f_i - e_ih_{i-1}`$ approaches zero, $`h_i`$ and $`s_i`$ become very large, and when it is zero the computation stops. The Thomas algorithm does not pivot, so the only safeguard is a property of the matrix. If the matrix is **strictly row diagonally dominant**, that is, $`|f_i| > |e_i| + |g_i|`$ in every row (with $`e_1 = g_n = 0`$), the following induction shows that no denominator vanishes.

1. $`|h_1| = |g_1|/|f_1| < 1`$.
2. Assuming $`|h_{i-1}| < 1`$, we get $`|d_i| \geq |f_i| - |e_i|\,|h_{i-1}| \geq |f_i| - |e_i| > |g_i| \geq 0`$, so $`d_i \neq 0`$.
3. Therefore $`|h_i| = |g_i|/|d_i| < 1`$, and by induction $`|h_i| < 1`$ and $`d_i \neq 0`$ for every $`i`$.

The bound $`|h_i| < 1`$ also matters in back substitution. If $`x_{i+1}`$ carries an error $`\delta`$, the error passed to $`x_i = s_i - h_ix_{i+1}`$ is $`|h_i|\delta < \delta`$, so errors do not grow as they move up the rows. Symmetric positive definite matrices also keep every pivot positive during elimination, so they too can be solved safely without pivoting.

### 5.5 Thomas Algorithm Example

The $`4 \times 4`$ tridiagonal system with $`f_i = 4`$ and $`e_i = g_i = -1`$ is solved for $`r = (2,\ 4,\ 6,\ 13)`$. Every row satisfies $`4 > 1 + 1`$, so the matrix is strictly diagonally dominant. The multipliers and updated values of the forward sweep and the normalized $`h_i`$ and $`s_i`$ are as follows.

| $`i`$ | Multiplier $`e_i/f_{i-1}'`$ | $`f_i'`$ | $`r_i'`$ | $`h_i`$ | $`s_i`$ |
|:-----|-----:|-----:|-----:|-----:|-----:|
| 1 | - | 4 | 2 | -0.25 | 0.5 |
| 2 | -0.25 | 3.75 | 4.5 | -0.266667 | 1.2 |
| 3 | -0.266667 | 3.733333 | 7.2 | -0.267857 | 1.928571 |
| 4 | -0.267857 | 3.732143 | 14.928571 | - | 4 |

*Table 3. Updated values of the forward sweep and normalized coefficients in the Thomas algorithm example*

For example, at $`i = 2`$ the multiplier is $`-1/4 = -0.25`$, so $`f_2' = 4 - (-0.25)(-1) = 3.75`$ and $`r_2' = 4 - (-0.25)(2) = 4.5`$. Normalizing gives $`h_2 = -1/3.75 = -4/15`$ and $`s_2 = 4.5/3.75 = 1.2`$. Back substitution proceeds as follows.

```math
x_4 = s_4 = 4, \quad x_3 = \frac{27}{14} - \left(-\frac{15}{56}\right)(4) = 3, \quad x_2 = 1.2 - \left(-\frac{4}{15}\right)(3) = 2, \quad x_1 = 0.5 - (-0.25)(2) = 1
```

The solution is $`(1,\ 2,\ 3,\ 4)`$. In the table every $`|h_i|`$ is below 1, and the values $`-0.25`$, $`-0.2667`$, $`-0.2679`$ quickly approach a single value. That value is the fixed point $`h = -2 + \sqrt{3} \approx -0.26795`$ of the recurrence $`h = -1/(4 + h)`$.

### 5.6 Tridiagonal Systems from a Heat Conduction Problem

How tridiagonal systems arise in practice can be seen in the steady-state temperature distribution of a thin rod of length $`L`$. If the temperatures at the two ends are fixed at $`T(0) = T_0`$ and $`T(L) = T_L`$ and heat escapes by convection from the lateral surface to the surrounding air at temperature $`T_a`$, the temperature $`T(x)`$ satisfies the following boundary-value problem.

```math
\frac{d^2T}{dx^2} + h'(T_a - T) = 0
```

Here $`h'`$ (in $`\text{m}^{-2}`$) is a coefficient that expresses how much heat is lost by convection. Dividing the rod into interior points $`x_1, \ldots, x_n`$ with spacing $`\Delta x = L/(n+1)`$ and replacing the second derivative with the centered difference $`(T_{i-1} - 2T_i + T_{i+1})/\Delta x^2`$ gives the following equation at each interior point.

```math
-T_{i-1} + (2 + h'\Delta x^2)\,T_i - T_{i+1} = h'\Delta x^2\,T_a, \qquad i = 1, \ldots, n
```

Each equation involves only the three values $`T_{i-1}`$, $`T_i`$, and $`T_{i+1}`$, so the coefficient matrix is tridiagonal with $`e_i = g_i = -1`$ and $`f_i = 2 + h'\Delta x^2`$. The values $`T_0`$ and $`T_{n+1} = T_L`$ in the first and last equations are already known, so moving them to the right-hand side gives $`r_1 = h'\Delta x^2T_a + T_0`$, $`r_n = h'\Delta x^2T_a + T_L`$, and $`r_i = h'\Delta x^2T_a`$ for the rest. When $`h' > 0`$, $`|f_i| = 2 + h'\Delta x^2 > 2 = |e_i| + |g_i|`$, so by the condition of Section 5.4 the Thomas algorithm works safely without pivoting.

With $`L = 10`$, $`n = 4`$ ($`\Delta x = 2`$), $`h' = 0.05`$, $`T_a = 20`$, $`T_0 = 40`$, and $`T_L = 200`$, we get $`f_i = 2 + 0.05 \times 4 = 2.2`$ and $`r = (44,\ 4,\ 4,\ 204)`$. The analytical solution of the differential equation is $`T(x) = T_a + C_1e^{\lambda x} + C_2e^{-\lambda x}`$ with $`\lambda = \sqrt{h'}`$, and the two boundary conditions $`C_1 + C_2 = T_0 - T_a`$ and $`C_1e^{\lambda L} + C_2e^{-\lambda L} = T_L - T_a`$ give $`C_1 \approx 19.2292`$ and $`C_2 \approx 0.7708`$. The two solutions compare as follows.

| $`x`$ | Numerical solution | Analytical solution | Error |
|-----:|-----:|-----:|-----:|
| 2 | 50.7866 | 50.5664 | 0.2203 |
| 4 | 67.7306 | 67.3485 | 0.3820 |
| 6 | 94.2206 | 93.7593 | 0.4613 |
| 8 | 135.5548 | 135.1695 | 0.3854 |

*Table 4. Finite-difference and analytical solutions for the heated rod*

The errors lie between 0.2 and 0.46 and arise because the truncation error of the centered difference is proportional to $`\Delta x^2`$. Improving the accuracy requires a smaller $`\Delta x`$ and correspondingly more unknowns $`n`$, but the cost of the Thomas algorithm is $`O(n)`$, so thousands or tens of thousands of grid points can be solved without difficulty. Applying the same discretization to general boundary-value problems is the subject of the finite-difference method in Chapter 24.

---

<br>

## 6. Python Implementation

Python's NumPy and SciPy libraries provide robust, production-quality implementations that use LAPACK routines with partial pivoting internally:

```python
import numpy as np

x = np.linalg.solve(A, b)  # Uses LAPACK with pivoting internally
```

For tridiagonal systems, SciPy provides a specialized solver:

```python
from scipy.linalg import solve_banded

# ab stores the diagonals in banded format
x = solve_banded((1, 1), ab, b)  # Exploits tridiagonal structure
```

> **Note:** In production code, always use `np.linalg.solve()` or `scipy.linalg.solve()` rather than implementing Gauss elimination from scratch. These functions use highly optimized LAPACK routines (DGESV) that include partial pivoting, cache-efficient blocking, and error handling.

### 6.1 Implementing Gauss Elimination with Partial Pivoting

Even when a library is available, implementing the algorithm once by hand makes the meaning of each step clear. The following function translates the pseudocode of Section 1.3 into NumPy.

```python
import numpy as np

def gauss_pivot(A, b):
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    n = len(b)
    Aug = np.hstack([A, b.reshape(-1, 1)])     # augmented matrix [A | b]
    for k in range(n - 1):
        p = k + np.argmax(np.abs(Aug[k:, k]))  # row of the largest |a_ik|, i >= k
        if Aug[p, k] == 0:
            raise ValueError("matrix is singular")
        if p != k:
            Aug[[k, p]] = Aug[[p, k]]          # swap rows k and p
        for i in range(k + 1, n):
            factor = Aug[i, k] / Aug[k, k]
            Aug[i, k:] -= factor * Aug[k, k:]
    x = np.zeros(n)
    x[n - 1] = Aug[n - 1, n] / Aug[n - 1, n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = (Aug[i, n] - Aug[i, i + 1:n] @ x[i + 1:]) / Aug[i, i]
    return x

A = np.array([[0, 2, 3], [4, 6, 7], [2, -3, 6]])
b = np.array([8, -3, 5])
x = gauss_pivot(A, b)
print(x)
print(np.allclose(x, np.linalg.solve(A, b)))
```

```text
[-5.42391304  0.02173913  2.65217391]
True
```

Each part of the code works as follows.

- `np.asarray(A, dtype=float)` converts the input to a floating-point array. If an integer array were used as is, an in-place update such as `Aug[i, k:] -= factor * ...`, which assigns a floating-point result into an integer array, would raise a `UFuncTypeError`, and other ways of writing it would truncate the fractional part.
- `np.hstack` appends $`b`$ as a column to the right of $`A`$ to form the augmented matrix. Row interchanges and row operations are then applied to the coefficients and the right-hand side together, which prevents the mistake of forgetting to swap the right-hand side.
- `np.argmax(np.abs(Aug[k:, k]))` returns the position within the subarray that starts at row $`k`$, so `k` must be added to obtain the row index in the full matrix.
- If the largest absolute value is zero, every remaining candidate in the column is zero and the matrix is singular. Production code often compares against a small tolerance instead of zero.
- `Aug[[k, p]] = Aug[[p, k]]` swaps the two rows at once with NumPy fancy indexing. Writing `Aug[k], Aug[p] = Aug[p], Aug[k]` instead can leave both rows equal because the views overlap.
- The elimination loop updates `Aug[i, k:]`, column $`k`$ through the end, in one operation. Columns to the left of $`k`$ are already zero and are not computed.
- In back substitution, `Aug[i, i + 1:n] @ x[i + 1:]` computes $`\sum_{j>i} a_{ij}x_j`$ as a dot product.

The result agrees with the hand computation of Section 1.4, and `np.allclose` confirms that it matches the library solution.

### 6.2 Implementing the Thomas Algorithm and the Banded Solver

The Thomas algorithm takes the three diagonals as length-$`n`$ vectors `e`, `f`, and `g` and translates the formulas of Section 5.2 directly. Indices start at 0, so `e[0]` and `g[n-1]` are not used.

```python
import numpy as np

def tridiag(e, f, g, r):
    n = len(f)
    f = np.array(f, dtype=float)               # copies, inputs stay intact
    r = np.array(r, dtype=float)
    for k in range(1, n):                      # forward sweep
        factor = e[k] / f[k - 1]
        f[k] = f[k] - factor * g[k - 1]
        r[k] = r[k] - factor * r[k - 1]
    x = np.zeros(n)
    x[n - 1] = r[n - 1] / f[n - 1]             # back substitution
    for k in range(n - 2, -1, -1):
        x[k] = (r[k] - g[k] * x[k + 1]) / f[k]
    return x

e = [0, -1, -1, -1]          # e[0] unused
f = [4, 4, 4, 4]
g = [-1, -1, -1, 0]          # g[n-1] unused
r = [2, 4, 6, 13]
print(tridiag(e, f, g, r))
```

```text
[1. 2. 3. 4.]
```

`f` and `r` are copied with `np.array` because the forward sweep overwrites both vectors. If the caller's arrays were modified in place, solving another right-hand side with the same coefficients would use values that have already changed. The multiplier `factor` is computed once and used in both update formulas.

SciPy's `scipy.linalg.solve_banded` solves general banded matrices, including tridiagonal ones. Its first argument `(1, 1)` gives the number of diagonals below and above the main diagonal, and the matrix is passed as a $`3 \times n`$ array `ab` stored by the rule `ab[u + i - j, j] = a[i, j]` (with $`u`$ the number of upper diagonals) [1]. For a tridiagonal matrix, row 0 of `ab` holds the super-diagonal shifted one place to the right, row 1 holds the main diagonal, and row 2 holds the sub-diagonal shifted one place to the left.

```python
import numpy as np
from scipy.linalg import solve_banded

n = 4
ab = np.zeros((3, n))
ab[0, 1:] = -1.0      # super-diagonal g, shifted right
ab[1, :] = 4.0        # main diagonal f
ab[2, :-1] = -1.0     # sub-diagonal e, shifted left
r = np.array([2.0, 4.0, 6.0, 13.0])
print(solve_banded((1, 1), ab, r))
```

```text
[1. 2. 3. 4.]
```

`ab[0, 0]` and `ab[2, n-1]` correspond to no matrix entry, so their values are not used. The solution is the same as that of the hand-written `tridiag`.

---

<br>

## 7. Determinant Evaluation with Gauss Elimination

### 7.1 Determinant of a Triangular Matrix and Row Operations

Gauss elimination is also an efficient way to compute determinants. The starting point is the fact that the determinant of an upper triangular matrix is the product of its diagonal entries. Expanding a $`3 \times 3`$ upper triangular matrix by cofactors along row 1, with $`a_{21} = a_{31} = a_{32} = 0`$, the minors of the second and third terms have a zero first column and vanish, leaving only the first term.

```math
\begin{vmatrix} a_{11} & a_{12} & a_{13} \\ 0 & a_{22} & a_{23} \\ 0 & 0 & a_{33} \end{vmatrix} = a_{11}\begin{vmatrix} a_{22} & a_{23} \\ 0 & a_{33} \end{vmatrix} - a_{12}\begin{vmatrix} 0 & a_{23} \\ 0 & a_{33} \end{vmatrix} + a_{13}\begin{vmatrix} 0 & a_{22} \\ 0 & 0 \end{vmatrix} = a_{11}a_{22}a_{33}
```

For a general $`n \times n`$ upper triangular matrix, expanding along column 1 leaves only $`a_{11}`$ times an upper triangular minor one size smaller, so by induction $`\det U = u_{11}u_{22}\cdots u_{nn}`$.

What remains is to know how elimination changes the determinant. Row operations act on the determinant as follows.

| Row operation | Effect on the determinant |
|:-----|:-----|
| Add a multiple of one row to another row | Unchanged |
| Interchange two rows | Sign changes |
| Multiply a row by a constant $`c`$ | Multiplied by $`c`$ |

*Table 5. Effect of row operations on the determinant*

The forward elimination of Gauss elimination uses only the first operation, and pivoting adds the second. Therefore the upper triangular matrix $`U`$ produced by forward elimination and the number of row interchanges $`p`$ give the determinant of the original matrix.

```math
\det A = (-1)^p\, u_{11}u_{22}\cdots u_{nn} = (-1)^p\, a_{11}a_{22}^{(1)}a_{33}^{(2)}\cdots a_{nn}^{(n-1)}
```

If a step divides rows by their pivots, as in the Gauss-Jordan method, the third operation enters, and the product of the divisors must be multiplied back.

### 7.2 Determinant Examples

**Without row interchanges.** For the example matrix of Section 4, $`A = \begin{bmatrix} 2 & 1 & 1 \\ 4 & 3 & 3 \\ 8 & 7 & 9 \end{bmatrix}`$, computing $`R_2 - 2R_1 = (0,\ 1,\ 1)`$, $`R_3 - 4R_1 = (0,\ 3,\ 5)`$, and then $`R_3 - 3R_2 = (0,\ 0,\ 2)`$ gives an upper triangular matrix with diagonal entries 2, 1, and 2. There is no interchange, so $`\det A = 2 \times 1 \times 2 = 4`$.

**With row interchanges.** In the example of Section 1.4, two row interchanges occurred and the diagonal entries were 4, $`-6`$, and $`23/6`$. Therefore:

```math
\det A = (-1)^2 \times 4 \times (-6) \times \frac{23}{6} = -92
```

Checking with a cofactor expansion along row 1 gives $`0 \cdot (36 + 21) - 2(24 - 14) + 3(-12 - 12) = -20 - 72 = -92`$, the same value. Omitting the sign of the interchanges gives an answer of the wrong sign whenever the number of interchanges is odd.

```python
import numpy as np

A = np.array([[0, 2, 3], [4, 6, 7], [2, -3, 6]], dtype=float)
print(np.linalg.det(A))
```

```text
-92.00000000000001
```

NumPy's `np.linalg.det` also computes the determinant through LU factorization (Chapter 10) [2]. The difference of about $`10^{-14}`$ in the last digits is roundoff error.

### 7.3 Cost Compared with Cofactor Expansion

Cofactor expansion turns an $`n \times n`$ determinant into $`n`$ determinants of size $`(n-1) \times (n-1)`$, so its number of multiplications $`M(n)`$ satisfies $`M(n) = n\,(M(n-1) + 1)`$ with $`M(1) = 0`$ and grows faster than $`n!`$. Gauss elimination, in contrast, needs only about $`2n^3/3`$ operations for forward elimination plus $`n-1`$ multiplications for the product of the diagonal entries.

| $`n`$ | Cofactor expansion multiplications $`M(n)`$ | Gauss elimination $`2n^3/3`$ |
|-----:|-----:|-----:|
| 3 | 9 | 18 |
| 5 | 205 | 83 |
| 10 | 6,235,300 | 667 |
| 20 | about $`4.2 \times 10^{18}`$ | 5,333 |

*Table 6. Cost of computing a determinant by cofactor expansion and by Gauss elimination*

For a small matrix such as $`n = 3`$, cofactor expansion is actually cheaper, but at $`n = 10`$ the gap is already about ten thousand, and at $`n = 20`$ cofactor expansion is practically impossible. This cost is also why Cramer's rule is not used for large systems.

> **Caution:** The magnitude of a determinant does not tell whether a system is ill-conditioned. The $`20 \times 20`$ matrix $`0.1I`$ has the very small determinant $`10^{-20}`$, yet it is perfectly well-conditioned. Conversely, a matrix with a large determinant can be ill-conditioned. Moreover, a product of many diagonal entries can overflow or underflow, so for large matrices `np.linalg.slogdet`, which returns the sign and the logarithm of the determinant separately, is used instead [3].

---

<br>

## 8. Deriving the Operation Count

This section derives step by step the operation count $`2n^3/3 + O(n^2)`$ of naive Gauss elimination summarized in Chapter 8. Each floating-point operation counts as one, with multiplications and divisions counted in one group and additions and subtractions in another.

> **Definition:** A **flop** is a unit that denotes one floating-point operation, that is, one addition, subtraction, multiplication, or division of floating-point numbers. Expressing the cost of an algorithm in flops allows algorithms to be compared independently of the speed of the computer.

### 8.1 Operation Count of Forward Elimination

Let the index of the outer loop be $`i`$ ($`i = 1, \ldots, n-1`$). Step $`i`$ uses pivot row $`i`$ to eliminate $`x_i`$ from the $`n-i`$ rows $`j = i+1, \ldots, n`$. The work in one row $`j`$ is as follows.

1. Compute the multiplier $`a_{ji}/a_{ii}`$: one division.
2. Update the $`n-i+1`$ entries $`a_{j,i+1}, \ldots, a_{jn}`$ and $`b_j`$ of row $`j`$ of the augmented matrix by $`a_{jk} \leftarrow a_{jk} - (\text{multiplier})\,a_{ik}`$: $`n-i+1`$ multiplications and $`n-i+1`$ subtractions.

$`a_{ji}`$ itself is certain to become zero, so it is not computed. Step $`i`$ therefore costs $`(n-i)(n+2-i)`$ multiplications and divisions and $`(n-i)(n+1-i)`$ additions and subtractions. Summing over all $`i`$ uses the following formulas.

```math
\sum_{i=1}^{m} 1 = m, \qquad \sum_{i=1}^{m} i = \frac{m(m+1)}{2}, \qquad \sum_{i=1}^{m} i^2 = \frac{m(m+1)(2m+1)}{6}
```

**Additions and subtractions.** Expanding the product gives $`(n-i)(n+1-i) = (n^2 + n) - (2n+1)i + i^2`$, and the formulas are applied with $`m = n-1`$.

```math
\sum_{i=1}^{n-1}(n-i)(n+1-i) = (n^2+n)(n-1) - (2n+1)\frac{(n-1)n}{2} + \frac{(n-1)n(2n-1)}{6} = \frac{n^3 - n}{3}
```

**Multiplications and divisions.** The substitution $`k = n - i`$ makes $`k`$ run from $`n-1`$ down to 1 and turns $`(n-i)(n+2-i)`$ into $`k(k+2)`$, which simplifies the computation.

```math
\sum_{k=1}^{n-1}k(k+2) = \frac{(n-1)n(2n-1)}{6} + (n-1)n = \frac{2n^3 + 3n^2 - 5n}{6} = \frac{n^3}{3} + \frac{n^2}{2} - \frac{5n}{6}
```

Adding the two results, forward elimination needs $`\frac{2}{3}n^3 + \frac{1}{2}n^2 - \frac{7}{6}n`$ operations in total.

### 8.2 Operation Count of Back Substitution

$`x_n = b_n/a_{nn}`$ is one division. For $`i = n-1, \ldots, 1`$, computing $`x_i = \bigl(b_i - \sum_{j=i+1}^{n}a_{ij}x_j\bigr)/a_{ii}`$ needs $`n-i`$ multiplications and $`n-i-1`$ additions for the sum, one subtraction from $`b_i`$, and finally one division. That is, each row $`i`$ costs $`n-i+1`$ multiplications and divisions and $`n-i`$ additions and subtractions (the same formula also holds for $`i = n`$).

```math
\sum_{i=1}^{n}(n-i+1) = \frac{n(n+1)}{2}, \qquad \sum_{i=1}^{n}(n-i) = \frac{n(n-1)}{2}
```

Adding the two sums, back substitution needs exactly $`n^2`$ operations.

### 8.3 Total Cost and the Extra Cost of Pivoting

Adding forward elimination and back substitution gives the following.

```math
\text{Total} = \frac{2}{3}n^3 + \frac{3}{2}n^2 - \frac{7}{6}n = \frac{2}{3}n^3 + O(n^2)
```

The following function follows the loops of the algorithm, counts the operations one by one, and compares the count with the formula.

```python
def count_gauss(n):
    md = asub = 0
    for i in range(1, n):              # pivot column i
        for j in range(i + 1, n + 1):  # rows below the pivot
            md += 1 + (n - i + 1)      # factor and n-i+1 products
            asub += n - i + 1          # n-i+1 subtractions
    for i in range(n, 0, -1):          # back substitution
        md += (n - i) + 1              # products and one division
        asub += n - i                  # additions and one subtraction
    return md, asub

for n in (3, 10, 100):
    md, asub = count_gauss(n)
    total = md + asub
    formula = 2 * n**3 / 3 + 3 * n**2 / 2 - 7 * n / 6
    print(n, md, asub, total, round(formula), round(2 * n**3 / 3))
```

```text
3 17 11 28 28 18
10 430 375 805 805 667
100 343300 338250 681550 681550 666667
```

The output columns are, in order, $`n`$, multiplications and divisions, additions and subtractions, the total, the exact formula, and the leading term $`2n^3/3`$. The counted total agrees exactly with the formula. When $`n`$ is small, the lower-order terms carry a large share and $`2n^3/3`$ underestimates the total noticeably, but as $`n`$ grows the gap becomes negligible.

| $`n`$ | Forward elimination | Back substitution | Total | $`2n^3/3`$ | Share of forward elimination (%) |
|-----:|-----:|-----:|-----:|-----:|-----:|
| 10 | 705 | 100 | 805 | 667 | 87.58 |
| 100 | 671,550 | 10,000 | 681,550 | 666,667 | 98.53 |
| 1000 | 667,165,500 | 1,000,000 | 668,165,500 | 666,666,667 | 99.85 |

*Table 7. Exact operation counts of Gauss elimination for several sizes*

Because the cost is proportional to $`n^3`$, doubling $`n`$ multiplies the computing time by about 8, and multiplying $`n`$ by 10 multiplies it by about 1000. The share of back substitution shrinks quickly as $`n`$ grows and is only 0.15% at $`n = 1000`$. This fact motivates the LU factorization of Chapter 10, which performs forward elimination only once when several right-hand sides share the same coefficient matrix.

**Extra cost of pivoting.** At step $`i`$, partial pivoting compares the absolute values of $`n-i+1`$ candidates, which takes $`n-i`$ comparisons, or $`\sum_{i=1}^{n-1}(n-i) = n(n-1)/2`$ in total. At most $`n-1`$ row interchanges occur, each moving $`n+1`$ entries, which is also $`O(n^2)`$. Both are negligible next to the $`O(n^3)`$ elimination, so, as stated in Section 1.3, there is no reason to omit pivoting.

---

<br>

## 9. Common Mistakes

The following mistakes are common when Gauss elimination and its variants are carried out by hand or implemented in code.

1. **Choosing the pivot without absolute values.** The pivot is the candidate with the largest absolute value, regardless of sign. If column 1 is $`(2,\ 4,\ -8)`$, the pivot is $`-8`$, not 4.
2. **Interchanging rows of the coefficient matrix only.** When rows are interchanged, the corresponding entries of the right-hand side $`b`$ must be interchanged as well. Working with the augmented matrix rules this mistake out structurally.
3. **Confusing a subarray index with a full index.** `np.argmax(np.abs(Aug[k:, k]))` returns a position counted from row $`k`$ as 0, so `k` must be added.
4. **Computing with integer arrays.** If the input is an integer array, in-place updates either raise an error or discard the fractional part. Convert to a floating-point array before computing.
5. **Omitting the sign of row interchanges in a determinant.** Forgetting $`p`$ in $`\det A = (-1)^p\prod u_{ii}`$ gives the wrong sign when the number of interchanges is odd.
6. **Actually dividing the coefficients by the scale factors.** In scaled partial pivoting, $`s_i`$ is used only in the comparison that selects the pivot; elimination uses the original coefficients.
7. **Applying the Thomas algorithm to a tridiagonal matrix that is not diagonally dominant.** $`\begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}`$ is nonsingular, yet $`f_1 = 0`$ causes a division by zero at the first multiplier. Such matrices are solved with a general solver with pivoting or with a banded solver.
8. **Believing that pivoting also cures ill-conditioned problems.** Pivoting only prevents error amplification by the algorithm. Problems in which small errors in the input data change the solution greatly are judged by the condition number of Chapter 11.

---

<br>

## Summary

| Concept | Key Summary |
|:-------|:-------|
| Naive Gauss | It does not pivot, costs $`O(n^3)`$, and is used only for teaching. |
| Gauss with Partial Pivoting (GEPP) | It picks the pivot with the largest absolute value at each step, costs $`O(n^3)`$, and is used for general systems. |
| Scaled partial pivoting | It picks the pivot by size relative to each row's scale factor, preventing wrong choices caused by rows of very different magnitudes. |
| Gauss-Jordan | It uses pivoting and costs $`O(n^3)`$ with a higher constant, and it is used for computing matrix inverses. |
| Thomas Algorithm | It solves tridiagonal systems at $`O(n)`$ cost without pivoting. |
| Pivoting rule | Always use pivoting for general systems, since the cost is negligible and the stability improvement is essential. |
| Exploiting structure | Exploit matrix structure whenever possible; the Thomas algorithm turns an $`O(n^3)`$ problem into $`O(n)`$ for tridiagonal systems. |
| Library use | Use library implementations in production, since they are optimized, tested, and handle edge cases. |
| Growth factor | The largest entry size during elimination divided by the original largest entry size; with partial pivoting it is at most $`2^{n-1}`$ and in practice it is mostly small. |
| Normalized Thomas algorithm | $`h_i = g_i/(f_i - e_ih_{i-1})`$, $`s_i = (r_i - e_is_{i-1})/(f_i - e_ih_{i-1})`$, $`x_i = s_i - h_ix_{i+1}`$; strict diagonal dominance guarantees $`\vert h_i \vert < 1`$. |
| Determinant evaluation | Multiply the diagonal entries of the upper triangular matrix from forward elimination and multiply by $`(-1)^p`$ for $`p`$ row interchanges. |
| Operation count | Forward elimination costs $`\frac{2}{3}n^3 + \frac{1}{2}n^2 - \frac{7}{6}n`$ and back substitution $`n^2`$; the $`n(n-1)/2`$ comparisons of pivoting are negligible. |
| Gauss-Jordan cost | Eliminating $`n-1`$ rows at every step needs about $`n^3/2`$ multiplications and divisions, about 1.5 times Gauss elimination. |

---

<br>

## Self-Check Questions

1. **Row selection in partial pivoting:** When partial pivoting is applied to $`A = \begin{bmatrix} 2 & 1 & 1 \\ 4 & -6 & 0 \\ -8 & 7 & 2 \end{bmatrix}`$, which row becomes the pivot row in the first step? Find the factors used in that step.

   <details>
   <summary>Answer</summary>

   The absolute values in column 1 are 2, 4, and 8, so row 3 becomes the pivot row and is swapped with row 1. The pivot is $`-8`$, and after the swap the factor for row 2 is $`4/(-8) = -0.5`$ and the factor for row 3 (formerly row 1) is $`2/(-8) = -0.25`$. Both factors have absolute value at most 1.

   </details>

2. **Scaled partial pivoting:** For $`2x_1 + 100000x_2 = 100000`$ and $`x_1 + x_2 = 2`$, which row does partial pivoting choose as the pivot row, and which does scaled partial pivoting choose?

   <details>
   <summary>Answer</summary>

   Partial pivoting chooses row 1 because $`|2| > |1|`$. Scaled partial pivoting uses $`s_1 = 100000`$ and $`s_2 = 1`$, giving ratios $`2/100000 = 2 \times 10^{-5}`$ and $`1/1 = 1`$, so it chooses row 2. In row 1 the coefficient of $`x_2`$ is so large that the coefficient 2 of $`x_1`$ is relatively small within that row.

   </details>

3. **Thomas algorithm:** Solve the tridiagonal system with $`f = (2, 2, 2)`$, $`e_2 = e_3 = -1`$, $`g_1 = g_2 = -1`$, and $`r = (1, 0, 1)`$ by the Thomas algorithm.

   <details>
   <summary>Answer</summary>

   In the forward sweep, for $`i = 2`$ the factor is $`-1/2 = -0.5`$, giving $`f_2 = 2 - (-0.5)(-1) = 1.5`$ and $`r_2 = 0 - (-0.5)(1) = 0.5`$. For $`i = 3`$ the factor is $`-1/1.5 = -2/3`$, giving $`f_3 = 2 - (-2/3)(-1) = 4/3`$ and $`r_3 = 1 - (-2/3)(0.5) = 4/3`$. Back substitution gives $`x_3 = (4/3)/(4/3) = 1`$, $`x_2 = (0.5 - (-1)(1))/1.5 = 1`$, and $`x_1 = (1 - (-1)(1))/2 = 1`$, so the solution is $`(1, 1, 1)`$.

   </details>

4. **Operation counts:** For $`n = 300`$, compare the multiplication counts of Gauss elimination with back substitution and of the Gauss-Jordan method, and find the total operation count of the Thomas algorithm for a tridiagonal system of the same size.

   <details>
   <summary>Answer</summary>

   Gauss elimination needs about $`n^3/3 = 9.0 \times 10^6`$ multiplications and Gauss-Jordan about $`n^3/2 = 1.35 \times 10^7`$, about 1.5 times as many. For a tridiagonal system, the Thomas algorithm needs only $`8n - 7 = 2393`$ operations.

   </details>

5. **Determinant evaluation:** Eliminate $`A = \begin{bmatrix} 1 & 1 & 1 \\ 2 & 2 & 5 \\ 4 & 6 & 8 \end{bmatrix}`$, interchanging rows only when a pivot is zero, and find $`\det A`$.

   <details>
   <summary>Answer</summary>

   $`R_2 - 2R_1 = (0,\ 0,\ 3)`$ and $`R_3 - 4R_1 = (0,\ 2,\ 4)`$. The pivot in column 2 is zero, so interchanging rows 2 and 3 gives an upper triangular matrix with diagonal entries 1, 2, and 3. There is one interchange, so $`\det A = (-1)^1 \times 1 \times 2 \times 3 = -6`$. This agrees with the cofactor expansion along row 1, $`1(16 - 30) - 1(16 - 20) + 1(12 - 8) = -6`$.

   </details>

6. **Operation count:** Find the exact operation count of naive Gauss elimination for $`n = 50`$, and state roughly how many times longer the computation takes when $`n`$ is increased to 100.

   <details>
   <summary>Answer</summary>

   Forward elimination costs $`\frac{2}{3}(125000) + \frac{1}{2}(2500) - \frac{7}{6}(50) = 83333.3 + 1250 - 58.3 = 84525`$ operations and back substitution $`50^2 = 2500`$, for a total of 87,025. The cost is proportional to $`n^3`$, so doubling $`n`$ makes it about $`2^3 = 8`$ times larger (the exact ratio is $`681550/87025 \approx 7.83`$).

   </details>

7. **Growth factor:** Find the growth factor of $`\begin{bmatrix} 0.0003 & 3 \\ 1 & 1 \end{bmatrix}`$ when it is eliminated without pivoting and when it is eliminated after a row interchange, and explain how this relates to the accuracy of the solution.

   <details>
   <summary>Answer</summary>

   Without pivoting the multiplier is $`1/0.0003 \approx 3333.3`$ and $`a_{22} = 1 - 3333.3 \times 3 \approx -9999`$, so $`\rho \approx 9999/3 = 3333`$. With the interchange the multiplier is 0.0003 and $`a_{22} = 3 - 0.0003 = 2.9997`$, so $`\rho = 1`$. Roundoff errors are proportional to the size of the numbers handled, so $`\rho = 3333`$ can cost about $`\log_{10}3333 \approx 3.5`$ additional significant digits; this is the same phenomenon as in Table 2, where the unpivoted result was 10000 times worse for the same number of digits.

   </details>

8. **Normalized Thomas algorithm:** Solve the $`3 \times 3`$ system with $`f_i = 3`$, $`e_i = g_i = -1`$, and $`r = (2,\ 1,\ 2)`$ with the $`h_i`$, $`s_i`$ recurrences of Section 5.4.

   <details>
   <summary>Answer</summary>

   $`h_1 = -1/3`$ and $`s_1 = 2/3`$. Since $`d_2 = 3 - (-1)(-1/3) = 8/3`$, $`h_2 = -1/(8/3) = -3/8`$ and $`s_2 = (1 + 2/3)/(8/3) = 5/8`$. Since $`d_3 = 3 - (-1)(-3/8) = 21/8`$, $`s_3 = (2 + 5/8)/(21/8) = 1`$. Back substitution gives $`x_3 = 1`$, $`x_2 = 5/8 + (3/8)(1) = 1`$, and $`x_1 = 2/3 + (1/3)(1) = 1`$, so the solution is $`(1,\ 1,\ 1)`$.

   </details>

---

<br>

## References

1. SciPy Developers, "scipy.linalg.solve_banded," *SciPy Documentation*. [Online]. Available: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_banded.html
2. NumPy Developers, "numpy.linalg.det," *NumPy Documentation*. [Online]. Available: https://numpy.org/doc/stable/reference/generated/numpy.linalg.det.html
3. NumPy Developers, "numpy.linalg.slogdet," *NumPy Documentation*. [Online]. Available: https://numpy.org/doc/stable/reference/generated/numpy.linalg.slogdet.html
