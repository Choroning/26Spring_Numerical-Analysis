# Chapter 9 Lab — Pivoting, Condition Number, and Tridiagonal Systems

> **Last Updated:** 2026-03-26

---

<br>

## Table of Contents

- [1. Partial Pivoting Implementation](#1-partial-pivoting-implementation)
  - [1.1 The Problem with Naive Elimination](#11-the-problem-with-naive-elimination)
  - [1.2 Gauss Elimination with Partial Pivoting](#12-gauss-elimination-with-partial-pivoting)
- [2. Condition Number](#2-condition-number)
  - [2.1 Well-Conditioned vs. Ill-Conditioned Systems](#21-well-conditioned-vs-ill-conditioned-systems)
  - [2.2 Perturbation Sensitivity of Ill-Conditioned Systems](#22-perturbation-sensitivity-of-ill-conditioned-systems)
- [3. Thomas Algorithm and Heat Conduction](#3-thomas-algorithm-and-heat-conduction)
  - [3.1 Thomas Algorithm Implementation](#31-thomas-algorithm-implementation)
  - [3.2 1D Steady-State Heat Conduction Application](#32-1d-steady-state-heat-conduction-application)
  - [3.3 Performance Comparison: Thomas vs. General Solver](#33-performance-comparison-thomas-vs-general-solver)
- [Summary](#summary)

---

<br>

## 1. Partial Pivoting Implementation

### 1.1 The Problem with Naive Elimination

Naive Gauss elimination fails when a pivot element is zero. Consider the following system where the $(1,1)$ entry is zero:

```python
import numpy as np

# Problem: naive Gauss fails when pivot is zero
A_bad = np.array([[0, 2, 1], [1, 1, 2], [2, 1, 1]], dtype=float)
b_bad = np.array([1, 1, 1], dtype=float)

# gauss_naive(A_bad, b_bad) would fail — division by zero!
```

The first pivot $a_{11} = 0$ causes immediate failure. Even if the pivot is not exactly zero but very small, the resulting large multipliers amplify round-off errors catastrophically.

### 1.2 Gauss Elimination with Partial Pivoting

The fix is to search for the largest element in the current column (from the pivot row downward) and swap rows before each elimination step:

```python
def gauss_pivot(A, b):
    """Solve Ax = b using Gaussian elimination with partial pivoting."""
    A = A.astype(float).copy()
    b = b.astype(float).copy()
    n = len(b)
    Aug = np.hstack([A, b.reshape(-1, 1)])
    nb = n + 1

    for k in range(n - 1):
        # Partial pivoting: find row with max |value| in column k
        imax = np.argmax(np.abs(Aug[k:, k])) + k
        if imax != k:
            Aug[[k, imax]] = Aug[[imax, k]]  # swap rows

        # Forward elimination (same as naive)
        for i in range(k + 1, n):
            factor = Aug[i, k] / Aug[k, k]
            Aug[i, k:nb] = Aug[i, k:nb] - factor * Aug[k, k:nb]

    # Back substitution
    x = np.zeros(n)
    x[n-1] = Aug[n-1, nb-1] / Aug[n-1, n-1]
    for i in range(n - 2, -1, -1):
        x[i] = (Aug[i, nb-1] - Aug[i, i+1:n] @ x[i+1:n]) / Aug[i, i]

    return x

# Test: now works!
x = gauss_pivot(A_bad, b_bad)
print(f'Solution: {x}')
print(f'Residual: {np.linalg.norm(A_bad @ x - b_bad):.2e}')
```

> **Note:** The key line is `imax = np.argmax(np.abs(Aug[k:, k])) + k` — this searches only the rows below (and including) the current pivot row. The `+ k` offset converts the local index back to the global row index. Without this offset, the wrong row would be swapped.

The row swap `Aug[[k, imax]] = Aug[[imax, k]]` uses NumPy's advanced indexing to exchange two rows in a single operation. The forward elimination and back substitution are identical to the naive version — the only addition is the pivot search and swap.

---

<br>

## 2. Condition Number

### 2.1 Well-Conditioned vs. Ill-Conditioned Systems

The **condition number** $\kappa(A)$ of a matrix quantifies how sensitive the solution of $Ax = b$ is to perturbations in the input data. A large condition number means small changes in $A$ or $b$ can cause large changes in $x$:

```python
# Well-conditioned system
A_good = np.array([[1, 0], [0, 1]], dtype=float)
print(f'Identity κ = {np.linalg.cond(A_good):.1f}')  # 1.0

# Ill-conditioned: Hilbert matrix
from scipy.linalg import hilbert
H = hilbert(5)
print(f'Hilbert(5) κ = {np.linalg.cond(H):.1f}')  # ~476607
```

The identity matrix has $\kappa = 1$ (the best possible), meaning perturbations are not amplified at all. The Hilbert matrix is a classic example of an ill-conditioned matrix — even at size $5 \times 5$, its condition number is nearly $500{,}000$.

> **[Linear Algebra]** The condition number $\kappa(A)$ is the ratio of the largest to smallest singular value. It bounds the worst-case error amplification: a relative perturbation of $\varepsilon$ in the input can cause a relative perturbation of up to $\kappa(A) \cdot \varepsilon$ in the output.

### 2.2 Perturbation Sensitivity of Ill-Conditioned Systems

The following experiment demonstrates the practical impact of ill-conditioning. A tiny perturbation in the right-hand side $b$ causes a disproportionately large change in the solution:

```python
# Perturbation sensitivity
n = 10
H10 = hilbert(n)
b_true = H10 @ np.ones(n)
x_exact = np.linalg.solve(H10, b_true)

# Perturb b slightly
b_perturbed = b_true + 1e-10 * np.random.randn(n)
x_perturbed = np.linalg.solve(H10, b_perturbed)

print(f'Hilbert(10) κ = {np.linalg.cond(H10):.2e}')
print(f'Input perturbation: {np.linalg.norm(b_perturbed - b_true):.2e}')
print(f'Output perturbation: {np.linalg.norm(x_perturbed - x_exact):.2e}')
# Output perturbation >> input perturbation!
```

The output perturbation is orders of magnitude larger than the input perturbation, with the amplification factor bounded by the condition number. For `hilbert(10)`, $\kappa \approx 10^{13}$, meaning a perturbation of $10^{-10}$ in $b$ can produce a change of up to $10^{3}$ in $x$. This makes the system essentially unsolvable in double precision.

> **Note:** When $\kappa(A) \cdot \varepsilon_{\text{mach}} \gtrsim 1$, the system is so ill-conditioned that no numerical method can produce a meaningful solution in standard double precision. The Hilbert matrix of size $n \geq 13$ reaches this threshold.

---

<br>

## 3. Thomas Algorithm and Heat Conduction

### 3.1 Thomas Algorithm Implementation

The **Thomas algorithm** solves tridiagonal systems in $O(n)$ time by exploiting the banded structure. Only the three diagonals need to be stored and processed:

```python
def thomas(e, f, g, r):
    """Solve tridiagonal system using Thomas algorithm.

    Args:
        e: sub-diagonal (n elements, e[0] unused)
        f: main diagonal (n elements)
        g: super-diagonal (n elements, g[n-1] unused)
        r: right-hand side (n elements)
    Returns:
        x: solution vector
    """
    e = e.astype(float).copy()
    f = f.astype(float).copy()
    g = g.astype(float).copy()
    r = r.astype(float).copy()
    n = len(f)

    # Forward sweep
    for k in range(1, n):
        factor = e[k] / f[k-1]
        f[k] = f[k] - factor * g[k-1]
        r[k] = r[k] - factor * r[k-1]

    # Back substitution
    x = np.zeros(n)
    x[n-1] = r[n-1] / f[n-1]
    for k in range(n - 2, -1, -1):
        x[k] = (r[k] - g[k] * x[k+1]) / f[k]

    return x
```

The `.copy()` calls are essential — without them, the function would modify the caller's arrays because NumPy arrays are passed by reference. The forward sweep eliminates the sub-diagonal entries one by one, and the back substitution recovers the solution from bottom to top.

> **[Python]** The `.astype(float).copy()` pattern ensures two things: (1) the computation uses floating-point arithmetic even if integer arrays are passed, and (2) the original arrays are not modified. This is a common defensive programming pattern in numerical Python.

### 3.2 1D Steady-State Heat Conduction Application

A classic application of tridiagonal systems is the 1D steady-state heat equation:

$$\frac{d^2T}{dx^2} = 0 \quad \text{with} \quad T(0) = T_L, \quad T(L) = T_R$$

Discretizing with central finite differences on a uniform grid with spacing $\Delta x$:

$$\frac{T_{i-1} - 2T_i + T_{i+1}}{(\Delta x)^2} = 0$$

This simplifies to the tridiagonal system $T_{i-1} - 2T_i + T_{i+1} = 0$, where the boundary conditions enter through the right-hand side:

```python
import matplotlib.pyplot as plt

def heat_conduction_1d(T_left, T_right, n_interior):
    """Solve 1D steady heat conduction with Dirichlet BCs."""
    n = n_interior

    e = np.ones(n)        # sub-diagonal
    f = -2 * np.ones(n)   # main diagonal
    g = np.ones(n)        # super-diagonal
    r = np.zeros(n)       # RHS

    # Boundary conditions
    r[0] = -T_left
    r[-1] = -T_right

    T_interior = thomas(e, f, g, r)

    # Full temperature profile
    x = np.linspace(0, 1, n + 2)
    T = np.concatenate([[T_left], T_interior, [T_right]])
    return x, T

# Solve for different mesh sizes
for n in [3, 6, 15]:
    x, T = heat_conduction_1d(100, 20, n)
    plt.plot(x, T, 'o-', label=f'n={n}')

plt.plot([0, 1], [100, 20], 'k--', label='Exact (linear)')
plt.xlabel('Position x')
plt.ylabel('Temperature T')
plt.title('1D Steady Heat Conduction')
plt.legend()
plt.grid(True)
plt.show()
```

For the Laplace equation ($d^2T/dx^2 = 0$), the exact solution is linear: $T(x) = T_L + (T_R - T_L) \cdot x / L$. The numerical solution matches exactly (up to round-off) for any mesh size because linear functions are represented exactly by the central difference approximation.

> **[Differential Equations]** The central difference approximation $T''(x_i) \approx (T_{i-1} - 2T_i + T_{i+1}) / (\Delta x)^2$ is second-order accurate. For the Laplace equation, it gives exact results because the truncation error involves $T^{(4)}(x)$, which is zero for linear (and quadratic) functions.

### 3.3 Performance Comparison: Thomas vs. General Solver

The $O(n)$ vs. $O(n^3)$ difference between the Thomas algorithm and general Gauss elimination is dramatic at large problem sizes:

```python
import time

sizes = [100, 1000, 10000, 100000]
t_thomas = []
t_numpy = []

for n in sizes:
    e = np.ones(n)
    f = -2 * np.ones(n)
    g = np.ones(n)
    r = np.random.rand(n)

    # Thomas
    start = time.time()
    for _ in range(100):
        thomas(e, f, g, r)
    t_thomas.append((time.time() - start) / 100)

    # numpy (build full matrix)
    A = np.diag(f) + np.diag(g[:-1], 1) + np.diag(e[1:], -1)
    start = time.time()
    for _ in range(10):
        np.linalg.solve(A, r)
    t_numpy.append((time.time() - start) / 10)

print(f'{"n":>8} {"Thomas (s)":>12} {"NumPy (s)":>12} {"Speedup":>10}')
for i, n in enumerate(sizes):
    speedup = t_numpy[i] / t_thomas[i] if t_thomas[i] > 0 else float('inf')
    print(f'{n:8d} {t_thomas[i]:12.6f} {t_numpy[i]:12.6f} {speedup:10.1f}x')
```

At $n = 100{,}000$, the Thomas algorithm completes in milliseconds while `np.linalg.solve` (which builds and factors a full $100{,}000 \times 100{,}000$ matrix) would require terabytes of memory and hours of computation. This demonstrates why exploiting matrix structure is essential in large-scale scientific computing.

> **Note:** In production, use `scipy.linalg.solve_banded()` for tridiagonal systems — it provides an optimized LAPACK implementation (DGTSV) that is faster than a pure Python Thomas algorithm while maintaining $O(n)$ complexity.

---

<br>

## Summary

| Method | Implementation | When to Use |
|:-------|:--------------|:------------|
| `gauss_pivot()` | Partial pivoting Gauss | General dense systems |
| `thomas()` | Thomas algorithm $O(n)$ | Tridiagonal systems |
| `np.linalg.solve()` | Optimized LAPACK | Production use (always) |
| `np.linalg.cond()` | Condition number | Diagnose ill-conditioning |

Key implementation takeaways:

- **Partial pivoting** adds minimal code (pivot search + row swap) but prevents division by zero and controls error growth
- **Condition number** should always be checked before trusting a solution — if $\kappa(A) \cdot \varepsilon_{\text{mach}} \approx 1$, the results are unreliable
- **Thomas algorithm** achieves $O(n)$ by exploiting tridiagonal structure — always use specialized solvers when the matrix has special structure

---
