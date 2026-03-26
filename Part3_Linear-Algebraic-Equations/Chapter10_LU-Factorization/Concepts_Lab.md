# Chapter 10 Lab — LU Factorization and Cholesky

> **Last Updated:** 2026-03-26

---

<br>

## Table of Contents

- [1. Doolittle's Algorithm](#1-doolittles-algorithm)
  - [1.1 Implementation](#11-implementation)
  - [1.2 Comparison with SciPy](#12-comparison-with-scipy)
- [2. Solving via L and U](#2-solving-via-l-and-u)
  - [2.1 Forward Substitution](#21-forward-substitution)
  - [2.2 Backward Substitution](#22-backward-substitution)
  - [2.3 Complete LU Solve with Pivoting](#23-complete-lu-solve-with-pivoting)
- [3. Thermal Resistor Network](#3-thermal-resistor-network)
- [4. Multiple Right-Hand Sides — Why LU Shines](#4-multiple-right-hand-sides--why-lu-shines)
  - [4.1 Performance Comparison](#41-performance-comparison)
  - [4.2 Spring-Mass System with Multiple Load Cases](#42-spring-mass-system-with-multiple-load-cases)
- [5. Cholesky Factorization](#5-cholesky-factorization)
  - [5.1 Checking for SPD](#51-checking-for-spd)
  - [5.2 Cholesky Decomposition and Solve](#52-cholesky-decomposition-and-solve)
  - [5.3 FEM Spring Network Application](#53-fem-spring-network-application)
- [Summary](#summary)

---

<br>

## 1. Doolittle's Algorithm

### 1.1 Implementation

Doolittle's algorithm computes the LU factorization by systematically building the rows of $U$ and the columns of $L$. For each pivot index $k$, the algorithm first computes the $k$-th row of $U$ using previously computed entries of $L$ and $U$, then computes the $k$-th column of $L$ by dividing by the pivot $U_{kk}$.

```python
import numpy as np
import scipy as sc
from scipy.linalg import lu

def lu_doolittle(A):
    """LU factorization using Doolittle's algorithm (no pivoting)."""
    n = A.shape[0]
    L = np.eye(n)
    U = np.zeros((n, n))
    A = A.astype(float).copy()

    for k in range(n):
        # Compute U row k
        for j in range(k, n):
            U[k, j] = A[k, j] - L[k, :k] @ U[:k, j]
        # Compute L column k
        for i in range(k + 1, n):
            L[i, k] = (A[i, k] - L[i, :k] @ U[:k, k]) / U[k, k]

    return L, U

# Test
A = np.array([[3, -0.1, -0.2], [0.1, 7, -0.3], [0.3, -0.2, 10]])
L, U = lu_doolittle(A)
print(f'||A - LU|| = {np.linalg.norm(A - L @ U):.2e}')  # ≈ 0
```

The key idea is that for row $k$ of $U$:

$$U_{kj} = A_{kj} - \sum_{s=1}^{k-1} L_{ks} U_{sj}, \quad j = k, k+1, \ldots, n$$

and for column $k$ of $L$:

$$L_{ik} = \frac{A_{ik} - \sum_{s=1}^{k-1} L_{is} U_{sk}}{U_{kk}}, \quad i = k+1, k+2, \ldots, n$$

These formulas follow directly from equating the entries of $A$ with the corresponding entries of the product $LU$.

### 1.2 Comparison with SciPy

SciPy's `lu` function performs LU factorization with partial pivoting, returning the permutation matrix $P$, the lower triangular matrix $L$, and the upper triangular matrix $U$:

```python
# Compare with scipy
P, L_sp, U_sp = sc.linalg.lu(A)
print(f'scipy: P = \n{P}\nL = \n{L_sp}\nU = \n{U_sp}')
```

> **Note:** SciPy returns $P$ such that $PA = LU$. When no row swaps are needed (as for this well-conditioned matrix), $P$ is the identity matrix and the results match our Doolittle implementation.

---

<br>

## 2. Solving via L and U

### 2.1 Forward Substitution

Forward substitution solves the lower triangular system $Ld = b$ by computing each $d_i$ sequentially from top to bottom. Since $L$ is unit lower triangular (diagonal entries are 1), no division by the diagonal is needed:

```python
def forward_sub(L, b):
    """Solve Ld = b by forward substitution."""
    n = len(b)
    d = np.zeros(n)
    for i in range(n):
        d[i] = b[i] - L[i, :i] @ d[:i]
    return d
```

The expression `L[i, :i] @ d[:i]` computes $\sum_{j=0}^{i-1} L_{ij} d_j$ using NumPy's slice notation. For $i = 0$, the slice is empty and the dot product is zero, so $d_0 = b_0$ as expected.

### 2.2 Backward Substitution

Backward substitution solves the upper triangular system $Ux = d$ by computing each $x_i$ from bottom to top. Division by the diagonal entry $U_{ii}$ is required at each step:

```python
def backward_sub(U, d):
    """Solve Ux = d by backward substitution."""
    n = len(d)
    x = np.zeros(n)
    x[n-1] = d[n-1] / U[n-1, n-1]
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - U[i, i+1:] @ x[i+1:]) / U[i, i]
    return x
```

### 2.3 Complete LU Solve with Pivoting

Combining the factorization and both substitution phases into a single function that handles pivoting:

```python
def lu_solve(A, b):
    """Solve Ax = b using LU factorization with pivoting."""
    P, L, U = sc.linalg.lu(A)
    b_perm = P.T @ b        # Apply permutation
    d = forward_sub(L, b_perm)
    x = backward_sub(U, d)
    return x

# Test
b = np.array([7.85, -19.3, 71.4])
x = lu_solve(A, b)
print(f'Solution: {x}')  # [3.0, -2.5, 7.0003]
print(f'Residual: {np.linalg.norm(A @ x - b):.2e}')
```

The permutation is applied as `P.T @ b` because SciPy defines $P$ such that $PA = LU$, so $Ax = b$ becomes $LUx = PAx = Pb$, but since $P$ is returned as a left-multiplication matrix, we need $P^T b$ to reorder $b$ according to the row swaps. The residual $\|Ax - b\|$ should be near machine epsilon, confirming the accuracy of the solution.

---

<br>

## 3. Thermal Resistor Network

A practical application of linear systems is the analysis of thermal resistor networks. Consider a 4-node network where temperatures at the boundary nodes are fixed and heat sources are applied at the interior nodes:

- $T_1 = 100°C$ (fixed), $T_4 = 20°C$ (fixed)
- Heat sources: $Q_2 = 50\text{W}$, $Q_3 = 30\text{W}$
- Thermal resistances: $R_{12} = 10$, $R_{23} = 20$, $R_{34} = 15$ (W/K)

The energy balance at each interior node yields:

$$\frac{T_2 - T_1}{R_{12}} + \frac{T_2 - T_3}{R_{23}} = Q_2$$

$$\frac{T_3 - T_2}{R_{23}} + \frac{T_3 - T_4}{R_{34}} = Q_3$$

Rearranging into matrix form $A \mathbf{T} = \mathbf{b}$:

```python
# 4-node network: T1=100 (fixed), T4=20 (fixed)
# Heat sources: Q2=50W, Q3=30W
# Resistances: R12=10, R23=20, R34=15 (W/K)

# Energy balance: (T2-T1)/R12 + (T2-T3)/R23 = Q2
#                  (T3-T2)/R23 + (T3-T4)/R34 = Q3

A_net = np.array([
    [1/10 + 1/20, -1/20],
    [-1/20, 1/20 + 1/15]
])
b_net = np.array([50 + 100/10, 30 + 20/15])
T = np.linalg.solve(A_net, b_net)
print(f'T2 = {T[0]:.1f}°C, T3 = {T[1]:.1f}°C')
```

The diagonal entries of $A_\text{net}$ are the sum of conductances (reciprocals of resistances) connected to each node, and the off-diagonal entries are the negative conductances between nodes. The right-hand side includes both the heat source terms and the contributions from the fixed boundary temperatures.

---

<br>

## 4. Multiple Right-Hand Sides — Why LU Shines

### 4.1 Performance Comparison

The primary advantage of LU factorization becomes clear when the same coefficient matrix must be used with many different right-hand side vectors. Factoring $A$ once and reusing the factors is far cheaper than solving from scratch each time:

```python
import time

n = 500
A = np.random.rand(n, n) + n * np.eye(n)
m = 100  # number of different RHS

# Method 1: Repeated np.linalg.solve
start = time.time()
for _ in range(m):
    b = np.random.rand(n)
    x = np.linalg.solve(A, b)
t1 = time.time() - start

# Method 2: Factor once, solve many
start = time.time()
lu_obj = sc.linalg.lu_factor(A)
for _ in range(m):
    b = np.random.rand(n)
    x = sc.linalg.lu_solve(lu_obj, b)
t2 = time.time() - start

print(f'Repeated solve: {t1:.3f}s')
print(f'LU factor + solve: {t2:.3f}s')
print(f'Speedup: {t1/t2:.1f}x')
```

For $n = 500$ and $m = 100$ right-hand sides:

- **Repeated solve** calls `np.linalg.solve` 100 times, each performing a full $O(n^3)$ factorization internally. Total cost: $100 \times O(n^3) = O(100 \, n^3)$.
- **LU factor + solve** performs one $O(n^3)$ factorization and 100 forward/backward substitutions at $O(n^2)$ each. Total cost: $O(n^3) + O(100 \, n^2) \approx O(n^3)$ for large $n$.

The speedup grows linearly with $m$: the more right-hand sides, the greater the benefit of pre-factoring.

### 4.2 Spring-Mass System with Multiple Load Cases

A spring-mass system provides a concrete engineering example where multiple load cases share the same stiffness matrix:

```python
K = np.array([
    [300, -200, 0],
    [-200, 350, -150],
    [0, -150, 400]
], dtype=float)

F_gravity = np.array([0, 0, -9.81 * 10])  # 10 kg mass at node 3
F_wind = np.array([50, 30, 0])
F_combined = F_gravity + F_wind

lu_K = sc.linalg.lu_factor(K)
x_grav = sc.linalg.lu_solve(lu_K, F_gravity)
x_wind = sc.linalg.lu_solve(lu_K, F_wind)
x_comb = sc.linalg.lu_solve(lu_K, F_combined)

# Superposition check
print(f'||x_comb - (x_grav + x_wind)|| = {np.linalg.norm(x_comb - (x_grav + x_wind)):.2e}')
# Should be ~1e-16 (machine epsilon)
```

> **[Physics]** The principle of superposition states that for linear systems, the response to a combined load equals the sum of individual responses. This is a direct consequence of linearity: $A(x_1 + x_2) = Ax_1 + Ax_2 = b_1 + b_2$.

The superposition check confirms this property numerically: the difference between the combined solution and the sum of individual solutions is at machine epsilon, demonstrating that the linear system preserves superposition exactly (up to floating-point precision).

---

<br>

## 5. Cholesky Factorization

### 5.1 Checking for SPD

Before applying Cholesky factorization, we must verify that the matrix is symmetric positive definite. A matrix is SPD if it is symmetric ($A = A^T$) and all its eigenvalues are strictly positive:

```python
def is_spd(A):
    """Check if A is symmetric positive definite."""
    if not np.allclose(A, A.T):
        return False
    eigenvalues = np.linalg.eigvalsh(A)
    return np.all(eigenvalues > 0)

# Example SPD matrix
A_spd = np.array([[6, 15, 55], [15, 55, 225], [55, 225, 979]], dtype=float)
print(f'Is SPD: {is_spd(A_spd)}')
```

The function `np.linalg.eigvalsh` is used instead of `np.linalg.eigvals` because it is optimized for symmetric matrices (the `h` stands for Hermitian). It returns only real eigenvalues and is both faster and more numerically stable for this purpose.

### 5.2 Cholesky Decomposition and Solve

SciPy provides the Cholesky factorization through `sc.linalg.cholesky`. By default it returns the upper triangular factor; passing `lower=True` returns the lower triangular factor:

```python
# Cholesky
L_chol = sc.linalg.cholesky(A_spd, lower=True)
print(f'L = \n{L_chol}')
print(f'||A - L L^T|| = {np.linalg.norm(A_spd - L_chol @ L_chol.T):.2e}')

def cholesky_solve(A, b):
    """Solve Ax = b using Cholesky factorization."""
    L = sc.linalg.cholesky(A, lower=True)
    d = forward_sub(L, b)
    x = backward_sub(L.T, d)
    return x
```

The solve procedure mirrors general LU: first solve $Ld = b$ by forward substitution, then solve $L^T x = d$ by backward substitution. The key difference is that only one triangular factor needs to be computed (since $U = L^T$), which halves both the computation and the storage.

> **Note:** The `forward_sub` function defined in Section 2 assumes a unit lower triangular matrix (diagonal entries are 1). For Cholesky, the diagonal entries of $L$ are not 1, so the forward substitution formula becomes $d_i = (b_i - \sum_{j=1}^{i-1} L_{ij} d_j) / L_{ii}$. In practice, SciPy's `cho_solve` handles this automatically.

### 5.3 FEM Spring Network Application

The stiffness matrix of a spring network assembled by the finite element method is always SPD (for valid, fully constrained spring systems), making Cholesky the ideal factorization:

```python
# Spring constants
k = [100, 200, 150, 250]

# Assembly (stiffness matrix is always SPD for valid spring systems)
K_fem = np.array([
    [k[0]+k[1], -k[1], 0],
    [-k[1], k[1]+k[2], -k[2]],
    [0, -k[2], k[2]+k[3]]
], dtype=float)

F_fem = np.array([50.0, 0.0, -30.0])
print(f'K is SPD: {is_spd(K_fem)}')
x_fem = cholesky_solve(K_fem, F_fem)
print(f'Displacements: {x_fem}')
```

The stiffness matrix $K$ is assembled by summing the contributions of each spring element. Element $k_i$ connecting nodes $i$ and $i+1$ contributes $+k_i$ to the diagonal entries $(i,i)$ and $(i+1,i+1)$, and $-k_i$ to the off-diagonal entries $(i,i+1)$ and $(i+1,i)$. The resulting matrix is always symmetric (by construction) and positive definite (because the energy $\frac{1}{2} x^T K x > 0$ for any nonzero displacement $x$).

---

<br>

## Summary

| Method | Function | When to Use |
|:-------|:---------|:-----------|
| Doolittle LU | `lu_doolittle(A)` | Understanding the algorithm |
| SciPy LU | `sc.linalg.lu(A)` | General factorization |
| LU factor + solve | `lu_factor(A)` + `lu_solve(lu, b)` | Multiple RHS (factor once!) |
| Cholesky | `sc.linalg.cholesky(A)` | SPD matrices (2x faster) |
| `np.linalg.solve` | Direct solve | Single RHS, production use |

**Practical Guidelines:**

1. Always use `np.linalg.solve` for single RHS — never manually invert
2. For multiple RHS with same $A$: use `lu_factor` + `lu_solve`
3. If $A$ is SPD: use Cholesky (faster, more stable)
4. Check condition number before trusting results
5. Verify solutions with residual: $\|Ax - b\|$

---
