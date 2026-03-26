# Chapter 9 Lecture — Gauss Elimination: Improved Methods

> **Last Updated:** 2026-03-26

---

<br>

## Table of Contents

- [1. Partial Pivoting](#1-partial-pivoting)
  - [1.1 The Problem with Naive Gauss Elimination](#11-the-problem-with-naive-gauss-elimination)
  - [1.2 Partial Pivoting Strategy](#12-partial-pivoting-strategy)
  - [1.3 Modified Algorithm Pseudocode](#13-modified-algorithm-pseudocode)
- [2. Scaled Partial Pivoting](#2-scaled-partial-pivoting)
- [3. Example Demonstrating Pivoting Necessity](#3-example-demonstrating-pivoting-necessity)
  - [3.1 Without Pivoting (4-Digit Arithmetic)](#31-without-pivoting-4-digit-arithmetic)
  - [3.2 With Pivoting](#32-with-pivoting)
- [4. Gauss-Jordan Method](#4-gauss-jordan-method)
- [5. Tridiagonal Systems — Thomas Algorithm](#5-tridiagonal-systems--thomas-algorithm)
  - [5.1 Tridiagonal Matrix Structure](#51-tridiagonal-matrix-structure)
  - [5.2 Thomas Algorithm](#52-thomas-algorithm)
  - [5.3 Computational Cost](#53-computational-cost)
- [6. Python Implementation](#6-python-implementation)
- [Summary](#summary)

---

<br>

## 1. Partial Pivoting

### 1.1 The Problem with Naive Gauss Elimination

In naive Gauss elimination, the pivot element $a_{ii}$ is used directly as the divisor during each elimination step. This causes two critical problems:

1. **Division by zero**: If $a_{ii} = 0$, the algorithm fails outright
2. **Large round-off errors**: If $a_{ii}$ is very small (close to zero), the multiplier $a_{ki}/a_{ii}$ becomes very large, amplifying every round-off error in subsequent operations by that factor

Both problems stem from the same root cause: the pivot element is too small relative to other elements in the column.

### 1.2 Partial Pivoting Strategy

The solution is straightforward: before performing elimination on column $i$, search the rows from $i$ to $n$ to find the row with the **largest absolute value** in column $i$, and swap that row with row $i$. This ensures the pivot is always the largest available element in the column, which minimizes the magnitude of the multipliers.

Since the multiplier is computed as $\text{factor} = a_{ki}/a_{ii}$, making $|a_{ii}|$ as large as possible guarantees that all factors satisfy $|\text{factor}| \leq 1$. This bounds the error amplification at each step.

### 1.3 Modified Algorithm Pseudocode

The only modification to naive Gauss elimination is the addition of a pivot search and row swap before each elimination step:

```
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

> **Note:** Row swapping is an $O(n)$ operation per elimination step, so the total overhead of pivoting is $O(n^2)$ — negligible compared to the $O(n^3)$ cost of the elimination itself. There is essentially no reason to skip pivoting.

---

<br>

## 2. Scaled Partial Pivoting

Sometimes partial pivoting alone is insufficient. This happens when the rows of the matrix have **very different magnitudes**. A row with all large entries might dominate the pivot selection even when it is not the best choice relative to its own scale.

**Scaled partial pivoting** addresses this by normalizing each row before comparing pivot candidates:

1. Compute the **scale factor** for each row: $s_i = \max_j |a_{ij}|$
2. When selecting the pivot for column $i$, choose the row $k$ that maximizes $|a_{ki}| / s_k$ (the ratio of the candidate pivot to the row's scale factor)
3. Swap and proceed with standard elimination

The scale factors $s_i$ are computed once at the beginning and updated only when rows are swapped (the scale factors follow their rows). This normalization ensures that the pivot selection accounts for the overall magnitude of each row, preventing a row with uniformly large entries from always being chosen.

> **[Linear Algebra]** Scaled partial pivoting is equivalent to implicitly scaling the rows of the matrix before applying partial pivoting. The key insight is that the "best" pivot is not necessarily the largest in absolute value, but the largest relative to the other entries in its row.

---

<br>

## 3. Example Demonstrating Pivoting Necessity

Consider the system:

$$0.0003x_1 + 3.0000x_2 = 2.0001$$

$$1.0000x_1 + 1.0000x_2 = 1.0000$$

The exact solution is $x_1 = 1/3$ and $x_2 = 2/3$.

### 3.1 Without Pivoting (4-Digit Arithmetic)

Using 4-digit arithmetic with the equations in the given order:

- The pivot is $a_{11} = 0.0003$ (very small!)
- Factor = $1.0000 / 0.0003 = 3333.3$ (very large!)
- After elimination: $-9999x_2 = -6666$
- Back substitution: $x_2 = 0.6667$, then $x_1 = (2.0001 - 3.0000 \times 0.6667) / 0.0003 = 0.3333$

With limited precision, the answer is **wrong**. The enormous factor of 3333.3 amplifies every round-off error by that amount. The small perturbations introduced by 4-digit rounding are magnified to the point where significant digits are lost.

### 3.2 With Pivoting

Swap the two rows first (since $|1.0000| > |0.0003|$ in column 1):

$$1.0000x_1 + 1.0000x_2 = 1.0000$$

$$0.0003x_1 + 3.0000x_2 = 2.0001$$

Now:

- The pivot is $a_{11} = 1.0000$ (well-scaled)
- Factor = $0.0003 / 1.0000 = 0.0003$ (small factor!)
- After elimination: $2.9997x_2 = 1.9998$
- Back substitution: $x_2 = 0.6667$, then $x_1 = (1.0000 - 1.0000 \times 0.6667) / 1.0000 = 0.3333$

The result is much more accurate because the small factor (0.0003 instead of 3333.3) does not amplify round-off errors.

> **Note:** This example dramatically shows why pivoting matters. A factor of 3333.3 amplifies every round-off error by that amount. Pivoting keeps factors $\leq 1$, controlling error growth. In practice, partial pivoting is always used — naive Gauss elimination is only taught for pedagogical purposes.

---

<br>

## 4. Gauss-Jordan Method

The **Gauss-Jordan method** extends standard Gauss elimination by eliminating entries **both above and below** the diagonal, producing a **diagonal matrix** rather than an upper triangular one.

The procedure is:

1. **Forward elimination**: Same as standard Gauss elimination — eliminate below the diagonal
2. **Backward elimination**: Additionally eliminate entries **above** the diagonal
3. **Normalization**: Divide each row by its diagonal element to produce the identity matrix on the left

The result is that the augmented matrix $[A|b]$ is transformed into $[I|x]$, where $x$ is the solution. No back substitution is needed because the matrix is already diagonal (or identity after normalization).

**Key properties:**

- **No back substitution needed**: The solution is read directly from the transformed right-hand side
- **Useful for computing the matrix inverse**: Apply Gauss-Jordan to $[A|I]$ to obtain $[I|A^{-1}]$
- **More expensive**: The operation count is $O(n^3)$ with a larger constant than standard Gauss elimination — approximately $n^3/2$ multiplications vs. $n^3/3$ for standard Gauss with back substitution

> **[Linear Algebra]** The Gauss-Jordan method is the standard algorithm for computing matrix inverses. While it is about 50% more expensive than Gauss elimination with back substitution for solving a single system, it becomes competitive when solving for multiple right-hand sides simultaneously (as in computing the inverse).

---

<br>

## 5. Tridiagonal Systems — Thomas Algorithm

### 5.1 Tridiagonal Matrix Structure

A **tridiagonal matrix** has non-zero entries only on the main diagonal, the diagonal immediately above it (super-diagonal), and the diagonal immediately below it (sub-diagonal). All other entries are zero:

$$\begin{bmatrix} f_1 & g_1 & & \\ e_2 & f_2 & g_2 & \\ & \ddots & \ddots & \ddots \\ & & e_n & f_n \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{bmatrix} = \begin{bmatrix} r_1 \\ r_2 \\ \vdots \\ r_n \end{bmatrix}$$

where $e_i$ is the sub-diagonal, $f_i$ is the main diagonal, and $g_i$ is the super-diagonal.

Tridiagonal systems are extremely common in scientific computing:

- **Discretized ODEs**: Central finite differences for second-order boundary value problems
- **Spline interpolation**: Cubic spline coefficient systems
- **Implicit time-stepping**: Implicit methods for parabolic PDEs (e.g., heat equation)

### 5.2 Thomas Algorithm

The **Thomas algorithm** is a specialized form of Gauss elimination that exploits the tridiagonal structure. It consists of two phases:

**Forward sweep** (elimination): For $i = 2$ to $n$:

$$\text{factor} = \frac{e_i}{f_{i-1}}$$

$$f_i = f_i - \text{factor} \cdot g_{i-1}$$

$$r_i = r_i - \text{factor} \cdot r_{i-1}$$

**Back substitution**:

$$x_n = \frac{r_n}{f_n}$$

For $i = n-1$ down to $1$:

$$x_i = \frac{r_i - g_i \cdot x_{i+1}}{f_i}$$

Each step of the forward sweep only touches the current row and the row above it, and each step of back substitution only requires the next solution value. The algorithm only accesses the three diagonals and never fills in entries outside the band.

### 5.3 Computational Cost

The Thomas algorithm requires exactly $3(n-1)$ multiplications/divisions and $3(n-1)$ additions/subtractions in the forward sweep, plus $n$ operations in back substitution — giving a total cost of **$O(n)$**.

This is a dramatic improvement over the $O(n^3)$ cost of general Gauss elimination. For a system of size $n = 100{,}000$, the Thomas algorithm is approximately $10^{10}$ times faster.

> **[Data Structures]** The Thomas algorithm exploits the banded structure of the matrix. Only 3 diagonals are stored and processed, analogous to how sparse matrix formats save memory. This is a recurring theme: exploiting structure in the problem leads to dramatic efficiency gains.

---

<br>

## 6. Python Implementation

Python's NumPy/SciPy libraries provide robust, production-quality implementations that use LAPACK routines with partial pivoting internally:

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

---

<br>

## Summary

| Method | Pivot? | Cost | Best Use |
|:-------|:-------|:-----|:---------|
| Naive Gauss | No | $O(n^3)$ | Teaching only |
| Gauss with Partial Pivoting (GEPP) | Yes | $O(n^3)$ | General systems |
| Gauss-Jordan | Yes | $O(n^3)$ (higher constant) | Matrix inverse |
| Thomas Algorithm | No | $O(n)$ | Tridiagonal systems |

Key takeaways:

- **Always use pivoting** for general systems — the cost is negligible and the stability improvement is essential
- **Exploit structure** whenever possible — the Thomas algorithm turns an $O(n^3)$ problem into $O(n)$ for tridiagonal systems
- **Use library implementations** in production — they are optimized, tested, and handle edge cases

---
