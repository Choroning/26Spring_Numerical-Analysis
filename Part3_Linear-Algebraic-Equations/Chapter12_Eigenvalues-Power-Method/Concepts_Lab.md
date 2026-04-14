# Chapter 12 Lab — Power Method for Eigenvalue Computation

> **Last Updated:** 2026-04-14
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 13 (Power Method)

> **Prerequisites**: [Programming Language] Python, NumPy, SciPy. [Linear Algebra] Eigenvalues / eigenvectors, LU factorization (Ch 10).
>
> **Learning Objectives**:
> 1. Compute eigenpairs with `numpy.linalg.eig` and verify the trace / determinant identities
> 2. Implement the Power Method to find the dominant eigenvalue of a matrix
> 3. Apply the Inverse Power Method to find the smallest eigenvalue
> 4. Use the Shifted Power Method to find an eigenvalue near any target value
> 5. Analyze convergence rate as a function of the spectral gap $|\lambda_2 / \lambda_1|$

---

<br>

## Table of Contents

- [1. Eigenvalues and Eigenvectors](#1-eigenvalues-and-eigenvectors)
  - [1.1 Characteristic Polynomial and Key Identities](#11-characteristic-polynomial-and-key-identities)
  - [1.2 Geometric Interpretation](#12-geometric-interpretation)
  - [1.3 Using `numpy.linalg.eig`](#13-using-numpylinalgeig)
- [2. Power Method](#2-power-method)
  - [2.1 Algorithm](#21-algorithm)
  - [2.2 Step-by-Step Trace](#22-step-by-step-trace)
  - [2.3 Python Implementation](#23-python-implementation)
  - [2.4 Convergence Analysis](#24-convergence-analysis)
- [3. Inverse and Shifted Power Methods](#3-inverse-and-shifted-power-methods)
  - [3.1 Inverse Power Method — Smallest Eigenvalue](#31-inverse-power-method--smallest-eigenvalue)
  - [3.2 Shifted Power Method — Any Eigenvalue](#32-shifted-power-method--any-eigenvalue)
  - [3.3 Python Implementation](#33-python-implementation)
  - [3.4 Finding All Eigenvalues with Different Shifts](#34-finding-all-eigenvalues-with-different-shifts)
- [Summary](#summary)

---

<br>

## 1. Eigenvalues and Eigenvectors

### 1.1 Characteristic Polynomial and Key Identities

An **eigenvalue** $\lambda$ and **eigenvector** $\mathbf{v}$ of a square matrix $A$ satisfy $A\mathbf{v} = \lambda\mathbf{v}$ with $\mathbf{v} \neq \mathbf{0}$. Rearranging gives $(A - \lambda I)\mathbf{v} = \mathbf{0}$, which has a nontrivial solution only when the matrix is singular:

$$\det(A - \lambda I) = 0 \quad \leftarrow \text{characteristic equation}$$

For an $n \times n$ matrix this is a degree-$n$ polynomial in $\lambda$, yielding exactly $n$ eigenvalues (counting multiplicity, possibly complex). Several algebraic identities are useful for verification:

| Property | Formula |
|:---------|:--------|
| Trace identity | $\text{tr}(A) = \sum_i a_{ii} = \sum_i \lambda_i$ |
| Determinant identity | $\det(A) = \prod_i \lambda_i$ |
| Rank | equals the number of nonzero eigenvalues |

For a symmetric matrix, the eigenvectors returned by NumPy are orthonormal: $V^T V = I$.

---

<br>

### 1.2 Geometric Interpretation

For most vectors $\mathbf{u}$, the image $A\mathbf{u}$ points in a different direction from $\mathbf{u}$. An eigenvector is special: $A\mathbf{v}$ is exactly parallel to $\mathbf{v}$ and only scaled by $\lambda$. Geometrically, the map $\mathbf{x} \mapsto A\mathbf{x}$ transforms the unit circle into an ellipse whose axes are aligned with the eigenvectors:

```python
import numpy as np
import matplotlib.pyplot as plt

A_geo = np.array([[3.0, 1.0], [1.0, 3.0]])
eigvals_geo, eigvecs_geo = np.linalg.eig(A_geo)

theta   = np.linspace(0, 2 * np.pi, 300)
circle  = np.array([np.cos(theta), np.sin(theta)])
ellipse = A_geo @ circle

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
colors = ['crimson', 'steelblue']

ax = axes[0]
ax.plot(circle[0], circle[1], 'b-', lw=1.5)
for i, (lam, c) in enumerate(zip(eigvals_geo, colors)):
    v = eigvecs_geo[:, i]
    ax.annotate('', xy=v, xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=c, lw=2.5))
ax.set(xlim=[-1.5, 1.5], ylim=[-1.5, 1.5], aspect='equal',
       title='Before: unit circle + eigenvectors')

ax = axes[1]
ax.plot(ellipse[0], ellipse[1], 'b-', lw=1.5)
for i, (lam, c) in enumerate(zip(eigvals_geo, colors)):
    av = A_geo @ eigvecs_geo[:, i]
    ax.annotate('', xy=av, xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=c, lw=2.5))
ax.set(xlim=[-5, 5], ylim=[-5, 5], aspect='equal',
       title='After A·: ellipse + scaled eigenvectors')

plt.tight_layout()
plt.show()
print(f'Eigenvalues: {eigvals_geo}')
```

The eigenvectors $\mathbf{v}_1, \mathbf{v}_2$ are the only directions along which the transformation acts as a pure scaling. Every other direction gets rotated.

---

<br>

### 1.3 Using `numpy.linalg.eig`

`numpy.linalg.eig(A)` returns all $n$ eigenpairs in $O(n^3)$ time. The columns of the returned eigenvector matrix correspond (in order) to the entries of the eigenvalue array:

```python
# Lecture-notes example  (Ch. 13, Problem 13.1)
A = np.array([[20.0,  3.0,  2.0],
              [ 3.0,  9.0,  4.0],
              [ 2.0,  4.0, 12.0]])

eigvals, eigvecs = np.linalg.eig(A)

print('Eigenvalues:')
for i, lam in enumerate(eigvals):
    print(f'  lambda_{i+1} = {lam:.8f}')

print('\nEigenvectors (columns):')
print(np.round(eigvecs, 6))

# Verify  A v = lambda v
print('\nVerification  max|A v - lambda v|:')
for i in range(len(eigvals)):
    res = A @ eigvecs[:, i] - eigvals[i] * eigvecs[:, i]
    print(f'  lambda_{i+1}: {np.max(np.abs(res)):.2e}')
```

After computing eigenpairs you should always verify the trace and determinant identities as a sanity check:

```python
print('=== Key Properties ===')
print(f'trace(A)               = {np.trace(A):.6f}')
print(f'sum of eigenvalues     = {np.sum(eigvals):.6f}   <- must match')
print(f'det(A)                 = {np.linalg.det(A):.6f}')
print(f'product of eigenvalues = {np.prod(eigvals):.6f}   <- must match')

# For a symmetric matrix, eigenvectors are orthonormal
VTV = eigvecs.T @ eigvecs
print('\nV^T V (= I for symmetric A):')
print(np.round(VTV, 8))
```

If `trace ≠ sum(lambda)` or `det ≠ prod(lambda)`, something has gone wrong — either the matrix was entered incorrectly or a numerical issue has arisen.

---

<br>

## 2. Power Method

### 2.1 Algorithm

`np.linalg.eig` computes all $n$ eigenpairs — an $O(n^3)$ operation. For large matrices, only the **dominant** eigenvalue (the one with the largest $|\lambda|$) is often needed. The **Power Method** finds it in $O(n^2)$ per iteration using only matrix–vector products.

**Why does it work?** Expand the starting vector in the eigenbasis:

$$\mathbf{x}^{(0)} = c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_n\mathbf{v}_n$$

Then after $k$ multiplications by $A$:

$$A^k\mathbf{x}^{(0)} = \lambda_1^k\!\left[c_1\mathbf{v}_1 + c_2\!\left(\frac{\lambda_2}{\lambda_1}\right)^{\!k}\mathbf{v}_2 + \cdots\right]$$

If $|\lambda_1| > |\lambda_2| \ge \cdots \ge |\lambda_n|$, the ratio $(\lambda_2 / \lambda_1)^k \to 0$ and $A^k\mathbf{x}^{(0)}$ aligns with $\mathbf{v}_1$. The convergence rate per iteration is approximately $|\lambda_2 / \lambda_1|$.

**Algorithm** (normalize by the $\infty$-norm to prevent overflow). Given $\mathbf{x}^{(0)} = \mathbf{1}$, for $k = 1, 2, \ldots$:

1. $\mathbf{y}^{(k)} = A\,\mathbf{x}^{(k-1)}$
2. $\lambda^{(k)} = y^{(k)}_{i^*}$ where $i^* = \arg\max_i |y_i^{(k)}|$ — eigenvalue estimate
3. $\mathbf{x}^{(k)} = \mathbf{y}^{(k)} / \lambda^{(k)}$ — normalize so that $\|\mathbf{x}^{(k)}\|_\infty = 1$

---

<br>

### 2.2 Step-by-Step Trace

A hand trace makes the iteration concrete. Using the same $A$ as above with $\mathbf{x}^{(0)} = \mathbf{1}$:

```python
x_k = np.ones(3)
print(f'{"Iter":>5}  {"lambda est":>12}  {"x1":>10}  {"x2":>10}  {"x3":>10}')
print('-' * 57)
for k in range(1, 8):
    y     = A @ x_k
    lam_k = y[np.argmax(np.abs(y))]
    x_k   = y / lam_k
    print(f'{k:>5}  {lam_k:>12.6f}  {x_k[0]:>10.6f}  {x_k[1]:>10.6f}  {x_k[2]:>10.6f}')

print(f'\nTrue dominant lambda (numpy): {np.max(np.abs(eigvals)):.6f}')
```

Within a handful of iterations, `lambda est` settles near the true dominant eigenvalue and the normalized vector `x_k` approaches the corresponding eigenvector (with $\|\mathbf{x}\|_\infty = 1$).

---

<br>

### 2.3 Python Implementation

Wrap the iteration in a reusable function with a relative-change stopping criterion and a history for plotting:

```python
def power_method(A, x0=None, tol=1e-8, max_iter=500):
    """
    Find the dominant eigenvalue and eigenvector of A using the Power Method.

    Returns
    -------
    lam      : float        — dominant eigenvalue
    v        : (n,) ndarray — eigenvector (inf-norm normalized)
    lam_hist : list         — eigenvalue estimate history
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]
    x = np.ones(n) if x0 is None else np.array(x0, dtype=float)
    lam_hist = []
    lam_old  = np.inf

    for _ in range(max_iter):
        y   = A @ x
        lam = y[np.argmax(np.abs(y))]   # element with largest |value|
        x   = y / lam                    # normalize

        lam_hist.append(float(lam))
        if abs(lam - lam_old) / abs(lam) < tol:
            break
        lam_old = lam

    return float(lam), x, lam_hist
```

> **Note:** Normalizing by the element of largest absolute value (the $\infty$-norm) prevents overflow and keeps `lam` the eigenvalue estimate itself — no separate Rayleigh quotient is needed.

---

<br>

### 2.4 Convergence Analysis

The theoretical convergence rate of the Power Method is $|\lambda_2 / \lambda_1|$: the closer this ratio is to 1, the slower the convergence.

```python
lam_pm, v_pm, hist_pm = power_method(A, tol=1e-12)
lam_true = float(eigvals[np.argmax(np.abs(eigvals))])

print(f'Converged in {len(hist_pm)} iterations')
print(f'Dominant eigenvalue (Power Method): {lam_pm:.12f}')
print(f'True value (numpy)                : {lam_true:.12f}')
print(f'Absolute error                    : {abs(lam_pm - lam_true):.2e}')

sorted_lam = np.sort(np.abs(eigvals))[::-1]
rate = sorted_lam[1] / sorted_lam[0]
print(f'\nTheoretical convergence rate |λ₂/λ₁| = {rate:.6f}')
print(f'  -> error halves every {-np.log(2)/np.log(rate):.1f} iterations')

# Plot convergence
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
iters = range(1, len(hist_pm) + 1)

ax1.plot(iters, hist_pm, 'b.-', ms=6)
ax1.axhline(lam_true, color='r', ls='--', lw=2, label=f'True λ = {lam_true:.4f}')
ax1.set(xlabel='Iteration', ylabel='λ estimate', title='Eigenvalue Convergence')
ax1.legend(); ax1.grid(True, alpha=0.4)

errors = np.abs(np.array(hist_pm) - lam_true)
ax2.semilogy(iters, errors, 'r.-', ms=6, label='Actual error')
ref = errors[0] * rate ** np.array(range(len(errors)))
ax2.semilogy(iters, ref, 'k--', lw=1.5, label='|λ₂/λ₁|^k reference')
ax2.set(xlabel='Iteration', ylabel='|λ^(k) - λ_true|',
        title='Absolute Error (log scale)')
ax2.legend(); ax2.grid(True, which='both', ls=':', alpha=0.5)

plt.tight_layout()
plt.show()
```

The semi-log error plot should be a straight line with slope $\log|\lambda_2/\lambda_1|$. Well-separated eigenvalues give rapid linear convergence; a small spectral gap (say $|\lambda_2/\lambda_1| = 0.99$) requires thousands of iterations to converge.

---

<br>

## 3. Inverse and Shifted Power Methods

### 3.1 Inverse Power Method — Smallest Eigenvalue

If $A\mathbf{v} = \lambda\mathbf{v}$ and $A$ is invertible, then $A^{-1}\mathbf{v} = (1/\lambda)\mathbf{v}$. The smallest $|\lambda_i|$ of $A$ becomes the largest $|1/\lambda_i|$ of $A^{-1}$. Applying the Power Method to $A^{-1}$ therefore converges to $1/\lambda_{\min}$, and $\lambda_{\min} = 1 / \mu_{\max}$.

**Implementation:** never form $A^{-1}$ explicitly. Each iteration solves $A\mathbf{y} = \mathbf{x}$ using the LU factorization of $A$, computed once. Per-iteration cost is $O(n^2)$ after the $O(n^3)$ factorization.

---

<br>

### 3.2 Shifted Power Method — Any Eigenvalue

Apply the Inverse Power Method to $(A - \sigma I)$ for a user-chosen **shift** $\sigma$. The eigenvalues of $(A - \sigma I)^{-1}$ are $1/(\lambda_i - \sigma)$. The eigenvalue of $A$ closest to $\sigma$ has the largest $|1/(\lambda_i - \sigma)|$, so the iteration converges to it:

$$\lambda = \sigma + \frac{1}{\mu}, \qquad \mu = \text{dominant eigenvalue of } (A - \sigma I)^{-1}$$

Setting $\sigma = 0$ recovers the Inverse Power Method (finds smallest $|\lambda|$). Setting $\sigma$ near any target value picks out the interior eigenvalue closest to it.

---

<br>

### 3.3 Python Implementation

A single function handles both the Inverse Power Method ($\sigma = 0$) and the general Shifted variant:

```python
import scipy as sc

def shifted_power_method(A, sigma=0.0, x0=None, tol=1e-8, max_iter=500):
    """
    Find the eigenvalue of A nearest to `sigma` (Shifted Inverse Power Method).

    Algorithm:
      1. Build B = A - sigma * I
      2. Factor B = LU  (once)
      3. Power iteration on B^{-1}: solve B y = x, normalize
      4. True eigenvalue:  lambda = sigma + 1/mu   where mu = max(|y|)

    sigma = 0  ->  Inverse Power Method  (finds smallest |lambda|)
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]
    x = np.ones(n) if x0 is None else np.array(x0, dtype=float)

    B      = A - sigma * np.eye(n)
    lu_obj = sc.linalg.lu_factor(B)

    lam_hist = []
    lam_old  = np.inf

    for _ in range(max_iter):
        y   = sc.linalg.lu_solve(lu_obj, x)
        mu  = y[np.argmax(np.abs(y))]     # = 1 / (lambda - sigma)
        x   = y / mu
        lam = sigma + 1.0 / mu

        lam_hist.append(float(lam))
        if abs(lam - lam_old) / (abs(lam) + 1e-30) < tol:
            break
        lam_old = lam

    return float(lam), x, lam_hist
```

> **Note:** The LU factorization is computed once outside the loop, so each iteration only needs an $O(n^2)$ triangular solve. Never use `np.linalg.inv(B) @ x` inside the loop — it is both slower and less numerically accurate.

---

<br>

### 3.4 Finding All Eigenvalues with Different Shifts

By choosing shifts near each eigenvalue of interest, the Shifted Power Method can recover every eigenvalue of a small matrix without calling `eig`:

```python
lam_sorted = np.sort(np.real(eigvals))
print(f'All eigenvalues (numpy, ascending): {np.round(lam_sorted, 6)}\n')

tests = [
    ('Power Method   (dominant)', power_method,         dict(),            np.max(np.abs(eigvals))),
    ('sigma=0        (smallest)', shifted_power_method, dict(sigma=0.0),   lam_sorted[0]),
    ('sigma=10       (middle)',   shifted_power_method, dict(sigma=10.0),  lam_sorted[1]),
    ('sigma=24       (largest)',  shifted_power_method, dict(sigma=24.0),  lam_sorted[2]),
]

print(f'{"Method":<34}  {"Found λ":>14}  {"True λ":>12}  {"Error":>10}  {"Iters":>6}')
print('-' * 84)
for label, func, kwargs, lam_t in tests:
    lam_f, _, hist = func(A, tol=1e-12, **kwargs)
    print(f'{label:<34}  {lam_f:>14.8f}  {lam_t:>12.8f}  {abs(lam_f-lam_t):>10.2e}  {len(hist):>6}')
```

Shifts placed close to an eigenvalue lead to dramatically faster convergence: $|1/(\lambda_i - \sigma)|$ dominates the other reciprocals, so the effective spectral gap is very large. In practice, a rough estimate of $\sigma$ (e.g., from Gershgorin discs) is enough to converge in a handful of iterations.

---

<br>

## Summary

| Method | Finds | Cost per iteration | Convergence rate |
|:-------|:------|:------------------|:-----------------|
| `np.linalg.eig(A)` | All $n$ eigenpairs | $O(n^3)$ total | — |
| **Power Method** | Dominant $\lambda$ (largest $\|\lambda\|$) | $O(n^2)$ | $\|\lambda_2 / \lambda_1\|$ |
| **Inverse Power** | Smallest $\|\lambda\|$ | $O(n^2)$ after LU | $\|\lambda_{\min} / \lambda_{\text{next-smallest}}\|$ |
| **Shifted Power** | $\lambda$ nearest to $\sigma$ | $O(n^2)$ after LU | depends on shift quality |

Key formulas:

1. **Characteristic equation**: $\det(A - \lambda I) = 0$
2. **Trace identity**: $\text{tr}(A) = \sum_i \lambda_i$ (use as sanity check)
3. **Determinant identity**: $\det(A) = \prod_i \lambda_i$
4. **Power step**: $\lambda^{(k)} = \max\|A\mathbf{x}^{(k-1)}\|_\infty, \quad \mathbf{x}^{(k)} = A\mathbf{x}^{(k-1)} / \lambda^{(k)}$
5. **Shifted recovery**: $\lambda = \sigma + 1/\mu$ where $\mu$ is the dominant eigenvalue of $(A - \sigma I)^{-1}$

Practical tips:

- Always verify eigenpairs via the trace/determinant identities and the residual $\|A\mathbf{v} - \lambda\mathbf{v}\|$
- Factor $(A - \sigma I)$ once outside the iteration loop — never invert explicitly
- A small spectral gap means slow Power Method convergence; shift to amplify the effective gap
- Use the Shifted Power Method to probe interior eigenvalues without computing the whole spectrum

---
