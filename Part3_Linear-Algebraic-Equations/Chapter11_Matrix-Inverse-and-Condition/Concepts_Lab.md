# Chapter 11 Lab — Matrix Inverse and Condition Number

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 11

> **Prerequisites**: [Programming Language] MATLAB/Python. [Linear Algebra] LU factorization (Ch 8-10).
>
> **Learning Objectives**:
> 1. Compute matrix inverses numerically
> 2. Calculate and interpret condition numbers
> 3. Assess solution sensitivity to input perturbations

---

<br>

## Table of Contents

- [1. Matrix Inverse](#1-matrix-inverse)
  - [1.1 Computing the Inverse via LU Factorization](#11-computing-the-inverse-via-lu-factorization)
  - [1.2 Golden Rule: Solve, Don't Invert!](#12-golden-rule-solve-dont-invert)
- [2. Condition Number and Ill-Conditioning](#2-condition-number-and-ill-conditioning)
  - [2.1 Computing the Condition Number](#21-computing-the-condition-number)
  - [2.2 Geometric Interpretation](#22-geometric-interpretation)
  - [2.3 Hilbert Matrix — Exponential Ill-Conditioning](#23-hilbert-matrix--exponential-ill-conditioning)
  - [2.4 Error Amplification Demo](#24-error-amplification-demo)
  - [2.5 Diagnostic Exercise: Suspicious Stiffness Matrix](#25-diagnostic-exercise-suspicious-stiffness-matrix)
- [Summary](#summary)

---

<br>

## 1. Matrix Inverse

### 1.1 Computing the Inverse via LU Factorization

The inverse of a matrix $[A]$ can be computed column-by-column using LU factorization: decompose $[A]$ once, then solve $[A]\{x_k\} = \{e_k\}$ for each column $\{e_k\}$ of the identity matrix. Some useful identities for matrix inverses:

- $(AB)^{-1} = B^{-1}A^{-1}$
- $(A^T)^{-1} = (A^{-1})^T$
- For a $2 \times 2$ matrix: $\begin{bmatrix} a & b \\ c & d \end{bmatrix}^{-1} = \frac{1}{ad - bc} \begin{bmatrix} d & -b \\ -c & a \end{bmatrix}$

```python
import numpy as np
import scipy as sc

def lu_inverse(A):
    """Compute A^{-1} using LU factorization."""
    n = A.shape[0]
    A_inv = np.zeros((n, n))
    lu_obj = sc.linalg.lu_factor(A)
    I = np.eye(n)

    for j in range(n):
        A_inv[:, j] = sc.linalg.lu_solve(lu_obj, I[:, j])

    return A_inv

A = np.array([[3, -0.1, -0.2], [0.1, 7, -0.3], [0.3, -0.2, 10]])
A_inv = lu_inverse(A)
print(f'A^{{-1}} = \n{A_inv}')
print(f'||A A^{{-1}} - I|| = {np.linalg.norm(A @ A_inv - np.eye(3)):.2e}')
```

The LU-based approach decomposes $[A]$ into lower and upper triangular factors once, then reuses them for all $n$ right-hand sides. This is more efficient than computing $n$ independent Gauss eliminations.

---

<br>

### 1.2 Golden Rule: Solve, Don't Invert!

Computing $[A]^{-1}\{b\}$ explicitly is almost always the **wrong** approach to solving $[A]\{x\} = \{b\}$. Direct solving via LU factorization is both faster and more numerically accurate:

```python
import time

n = 500
A = np.random.rand(n, n) + n * np.eye(n)
b = np.random.rand(n)

# Method 1: Invert then multiply (SLOW, LESS ACCURATE)
start = time.time()
for _ in range(100):
    x_inv = np.linalg.inv(A) @ b
t_inv = (time.time() - start) / 100

# Method 2: Direct solve (FAST, MORE ACCURATE)
start = time.time()
for _ in range(100):
    x_solve = np.linalg.solve(A, b)
t_solve = (time.time() - start) / 100

print(f'Inverse + multiply: {t_inv:.6f}s')
print(f'Direct solve:       {t_solve:.6f}s')
print(f'Speedup: {t_inv/t_solve:.1f}x')
print(f'Error (inv):   {np.linalg.norm(A @ x_inv - b):.2e}')
print(f'Error (solve): {np.linalg.norm(A @ x_solve - b):.2e}')
```

> **Note:** Never compute $A^{-1}b$ when you can solve $Ax = b$ directly. Direct solving is approximately 3x faster and produces more accurate results. The only reason to compute the full inverse is when you need the matrix entries themselves (e.g., for physical interpretation or covariance analysis).

---

<br>

## 2. Condition Number and Ill-Conditioning

### 2.1 Computing the Condition Number

The condition number $\kappa(A) = \|A\| \cdot \|A^{-1}\|$ can be computed using different norms. NumPy provides a single function for all variants:

```python
# Computing condition number
A = np.array([[1, 2], [1.0001, 2]])
print(f'κ_2(A) = {np.linalg.cond(A, 2):.2e}')
print(f'κ_1(A) = {np.linalg.cond(A, 1):.2e}')
print(f'κ_∞(A) = {np.linalg.cond(A, np.inf):.2e}')

# 2-norm condition = ratio of singular values
U, S, Vt = np.linalg.svd(A)
print(f'σ_max/σ_min = {S[0]/S[-1]:.2e}')  # equals κ_2
```

The 2-norm condition number has a clean interpretation via the singular value decomposition (SVD): $\kappa_2(A) = \sigma_{max} / \sigma_{min}$, the ratio of the largest to smallest singular value. A matrix with one very small singular value is nearly singular, producing a large condition number.

---

<br>

### 2.2 Geometric Interpretation

Ill-conditioning has a clear geometric meaning for $2 \times 2$ systems. Each equation defines a line in the plane. A well-conditioned system has lines that cross at a steep angle, so the intersection point is well-defined. An ill-conditioned system has nearly parallel lines, so even a tiny perturbation shifts the intersection point dramatically:

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Well-conditioned
x_range = np.linspace(-2, 4, 100)
axes[0].plot(x_range, (5 - 2*x_range), label='2x + y = 5')
axes[0].plot(x_range, (x_range - 1), label='x - y = 1')
axes[0].set_title(f'Well-conditioned (κ ≈ {np.linalg.cond(np.array([[2,1],[1,-1]])):.1f})')
axes[0].legend()
axes[0].grid(True)

# Ill-conditioned
axes[1].plot(x_range, (5 - 2*x_range), label='2x + y = 5')
axes[1].plot(x_range, (5.001 - 2*x_range), label='2x + y = 5.001')
axes[1].set_title(f'Ill-conditioned (nearly parallel)')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()
```

In the well-conditioned case, the two lines have different slopes and intersect at a clear angle. The intersection point is stable under small perturbations. In the ill-conditioned case, the two lines are nearly identical (almost parallel), so a tiny change in the coefficients shifts the intersection from one end of the plot to the other.

---

<br>

### 2.3 Hilbert Matrix — Exponential Ill-Conditioning

The Hilbert matrix $H_{ij} = 1/(i + j - 1)$ is the standard textbook example of ill-conditioning. Its condition number grows **exponentially** with size, making it practically unsolvable beyond modest dimensions:

```python
from scipy.linalg import hilbert

print(f'{"n":>4} {"κ(H_n)":>15} {"digits lost":>12}')
print('-' * 35)
for n in range(2, 16):
    H = hilbert(n)
    kappa = np.linalg.cond(H)
    digits = int(np.log10(kappa)) if kappa > 1 else 0
    print(f'{n:4d} {kappa:15.2e} {digits:12d}')
```

For double precision ($\approx 16$ significant digits), the table shows that by $n = 12$, the Hilbert matrix has a condition number around $10^{16}$, meaning the computed solution has essentially **zero** reliable digits. Hilbert matrices arise naturally in polynomial least-squares fitting, which is one reason polynomial fitting beyond modest degrees is numerically dangerous.

---

<br>

### 2.4 Error Amplification Demo

The following experiment demonstrates how the condition number controls error amplification. We solve $[H]\{x\} = \{b\}$ where $\{x_{true}\} = \{1, 1, \ldots, 1\}$, then perturb $\{b\}$ by a tiny amount ($10^{-10}$) and re-solve:

```python
for n in [4, 8, 12]:
    H = hilbert(n)
    x_true = np.ones(n)
    b = H @ x_true

    # Solve
    x_computed = np.linalg.solve(H, b)

    # Perturb b
    b_pert = b + 1e-10 * np.random.randn(n)
    x_pert = np.linalg.solve(H, b_pert)

    rel_err = np.linalg.norm(x_computed - x_true) / np.linalg.norm(x_true)
    pert_err = np.linalg.norm(x_pert - x_computed) / np.linalg.norm(x_computed)

    print(f'n={n:2d}: κ={np.linalg.cond(H):.2e}, '
          f'solve error={rel_err:.2e}, perturbation effect={pert_err:.2e}')
```

The results show a clear pattern: as $n$ increases, the condition number grows exponentially, and both the solve error (even with exact arithmetic intent) and the perturbation sensitivity grow in proportion. For $n = 12$, a perturbation of $10^{-10}$ in $\{b\}$ can produce changes of order $10^6$ or more in $\{x\}$ — a magnification factor of $10^{16}$, matching the condition number.

---

<br>

### 2.5 Diagnostic Exercise: Suspicious Stiffness Matrix

In practice, ill-conditioning often appears in finite element method (FEM) stiffness matrices when boundary conditions are improperly applied or when the mesh has nearly degenerate elements. This example shows how to diagnose a suspicious matrix:

```python
# A suspicious stiffness matrix from a FEM simulation
K_susp = np.array([
    [1e6, -1e6, 0, 0],
    [-1e6, 2e6, -1e6, 0],
    [0, -1e6, 2e6, -1e6],
    [0, 0, -1e6, 1e6+1e-2]
])

f_susp = np.array([0, 100, 0, 0])

kappa = np.linalg.cond(K_susp)
digits_lost = int(np.log10(kappa))
x_susp = np.linalg.solve(K_susp, f_susp)

print(f'Condition number: {kappa:.2e}')
print(f'Digits of accuracy lost: ~{digits_lost}')
print(f'Solution: {x_susp}')
print(f'Residual: {np.linalg.norm(K_susp @ x_susp - f_susp):.2e}')
```

> **Note:** A large condition number doesn't necessarily mean the computed answer is wrong — it means you should be careful and verify. Check the residual $\|Ax - b\|$: a small residual with high $\kappa$ still indicates a reasonably good solution for the given data.

---

<br>

## Summary

| Concept | Function / Formula |
|:--------|:------------------|
| Compute inverse | `np.linalg.inv(A)` or `lu_inverse(A)` |
| Solve (preferred!) | `np.linalg.solve(A, b)` |
| Condition number | `np.linalg.cond(A)`, `np.linalg.cond(A, 2)` |
| Norms | `np.linalg.norm(A, ord)` with ord=1, 2, np.inf, 'fro' |
| Singular values | `np.linalg.svd(A)` — $\kappa_2 = \sigma_{max} / \sigma_{min}$ |

Key rules:

1. **Solve, don't invert** — `solve(A, b)` is faster and more accurate than `inv(A) @ b`
2. Always check $\kappa(A)$ before trusting results
3. Digits lost $\approx \log_{10}\kappa(A)$
4. Hilbert matrices: condition grows exponentially — classic ill-conditioning example
5. Small residual $\|Ax - b\|$ does not guarantee small error $\|x - x_{true}\|$ when $\kappa$ is large

---
