# Chapter 13 Lab — Eigenvalues and Eigenvectors

[Repository Guide](../../README.md) | [English](Concepts_Lab.md) | [한국어](Concepts_Lab.ko.md)

> **Last Updated:** 2026-11-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 13

> **Learning Objectives:**
> 1. Compute eigenpairs with `numpy.linalg.eig` and verify them with the trace and determinant identities
> 2. Find the dominant eigenvalue of a matrix with the power method and analyze the convergence rate as a function of the eigenvalue ratio $`|\lambda_2 / \lambda_1|`$
> 3. Find the smallest eigenvalue and an eigenvalue near any target value with the inverse and shifted power methods
> 4. Exploit the special properties of symmetric matrices using `np.linalg.eigh` and build low-rank approximations from the spectral decomposition $`A = V\Lambda V^T`$
> 5. Solve the generalized eigenvalue problem $`K\mathbf{v} = \lambda M\mathbf{v}`$ with `scipy.linalg.eigh(K, M)` to compute the natural frequencies and mode shapes of a mass-spring system
> 6. Analyze the stability of a linear ODE system $`\dot{\mathbf{y}} = A\mathbf{y}`$ from its eigenvalues

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
  - [3.1 Inverse Power Method for the Smallest Eigenvalue](#31-inverse-power-method-for-the-smallest-eigenvalue)
  - [3.2 Shifted Power Method for Any Eigenvalue](#32-shifted-power-method-for-any-eigenvalue)
  - [3.3 Python Implementation](#33-python-implementation)
  - [3.4 Finding All Eigenvalues with Different Shifts](#34-finding-all-eigenvalues-with-different-shifts)
- [4. Method Comparison and Practical Guidelines](#4-method-comparison-and-practical-guidelines)
  - [4.1 Method Comparison](#41-method-comparison)
  - [4.2 Key Formulas](#42-key-formulas)
  - [4.3 Practical Tips](#43-practical-tips)
- [5. Symmetric Matrices and Spectral Decomposition](#5-symmetric-matrices-and-spectral-decomposition)
  - [5.1 Special Properties](#51-special-properties)
  - [5.2 `np.linalg.eigh` and Reconstruction](#52-nplinalgeigh-and-reconstruction)
  - [5.3 Low-Rank Approximation](#53-low-rank-approximation)
- [6. Application: Vibration Analysis](#6-application-vibration-analysis)
  - [6.1 Equations of Motion and the Generalized Eigenproblem](#61-equations-of-motion-and-the-generalized-eigenproblem)
  - [6.2 Natural Frequencies via `scipy.linalg.eigh(K, M)`](#62-natural-frequencies-via-scipylinalgeighk-m)
  - [6.3 Visualizing Mode Shapes](#63-visualizing-mode-shapes)
  - [6.4 Free-Vibration Time Response](#64-free-vibration-time-response)
- [7. Application: ODE Systems and Stability](#7-application-ode-systems-and-stability)
  - [7.1 General Solution via Eigenpairs](#71-general-solution-via-eigenpairs)
  - [7.2 Stability Classification](#72-stability-classification)
  - [7.3 Example: Problem 13.13](#73-example-problem-1313)
- [8. Method Comparison and Engineering Context](#8-method-comparison-and-engineering-context)
  - [8.1 Methods Overview](#81-methods-overview)
  - [8.2 Key Formulas](#82-key-formulas)
  - [8.3 Engineering Context](#83-engineering-context)
  - [8.4 Practical Tips](#84-practical-tips)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Eigenvalues and Eigenvectors

### 1.1 Characteristic Polynomial and Key Identities

An **eigenvalue** $`\lambda`$ and **eigenvector** $`\mathbf{v}`$ of a square matrix $`A`$ satisfy $`A\mathbf{v} = \lambda\mathbf{v}`$ with $`\mathbf{v} \neq \mathbf{0}`$. Rearranging gives $`(A - \lambda I)\mathbf{v} = \mathbf{0}`$, which has a nontrivial solution only when the matrix is singular:

```math
\det(A - \lambda I) = 0 \quad \leftarrow \text{characteristic equation}
```

For an $`n \times n`$ matrix this is a degree-$`n`$ polynomial in $`\lambda`$, yielding exactly $`n`$ eigenvalues (counting multiplicity, possibly complex). Several algebraic identities are useful for verification:

| Property | Formula |
|:---------|:--------|
| Trace identity | $`\text{tr}(A) = \sum_i a_{ii} = \sum_i \lambda_i`$ |
| Determinant identity | $`\det(A) = \prod_i \lambda_i`$ |
| Rank | For a diagonalizable matrix (e.g., a symmetric matrix), equals the number of nonzero eigenvalues (counting multiplicity) |

*Table 1. Algebraic identities for verifying eigenvalues*

For a symmetric matrix with distinct eigenvalues, the unit eigenvectors returned by `eig` are orthonormal: $`V^T V = I`$. To obtain orthonormal eigenvectors even when eigenvalues repeat, use `eigh` (Section 5.2).

---

<br>

### 1.2 Geometric Interpretation

For most vectors $`\mathbf{u}`$, the image $`A\mathbf{u}`$ points in a different direction from $`\mathbf{u}`$. An eigenvector is special: $`A\mathbf{v}`$ is exactly parallel to $`\mathbf{v}`$ and only scaled by $`\lambda`$. Geometrically, for a symmetric matrix $`A`$, the map $`\mathbf{x} \mapsto A\mathbf{x}`$ transforms the unit circle into an ellipse whose axes are aligned with the eigenvectors (for a general matrix, the axes follow the singular vectors):

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

In this example, the eigenvectors $`\mathbf{v}_1, \mathbf{v}_2`$ are the only directions along which the transformation acts as a pure scaling, and every other direction gets rotated. The code prints the eigenvalues 4 and 2.

---

<br>

### 1.3 Using `numpy.linalg.eig`

`numpy.linalg.eig(A)` returns all $`n`$ eigenpairs in $`O(n^3)`$ time. The columns of the returned eigenvector matrix correspond (in order) to the entries of the eigenvalue array:

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

Running it gives the eigenvalues $`21.72915618`$, $`13.18232626`$, and $`6.08851756`$. After computing eigenpairs, it is good practice to verify the result with the trace and determinant identities as a sanity check:

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

If `trace ≠ sum(lambda)` or `det ≠ prod(lambda)`, something has gone wrong: either the matrix was entered incorrectly or a numerical issue has arisen. For this matrix, the trace 41 and the determinant 1744 match the sum and the product of the eigenvalues.

---

<br>

## 2. Power Method

### 2.1 Algorithm

`np.linalg.eig` computes all $`n`$ eigenpairs, an $`O(n^3)`$ operation. For large matrices, only the **dominant** eigenvalue (the one with the largest $`|\lambda|`$) is often needed. The **power method** finds it in $`O(n^2)`$ per iteration using only matrix-vector products.

**Reason it works.** Assume that $`A`$ has $`n`$ linearly independent eigenvectors and expand the starting vector in the eigenbasis:

```math
\mathbf{x}^{(0)} = c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots + c_n\mathbf{v}_n
```

Then after $`k`$ multiplications by $`A`$:

```math
A^k\mathbf{x}^{(0)} = \lambda_1^k\!\left[c_1\mathbf{v}_1 + c_2\!\left(\frac{\lambda_2}{\lambda_1}\right)^{\!k}\mathbf{v}_2 + \cdots\right]
```

If $`|\lambda_1| > |\lambda_2| \ge \cdots \ge |\lambda_n|`$ and $`c_1 \neq 0`$, the ratio $`(\lambda_2 / \lambda_1)^k \to 0`$ and the direction of $`A^k\mathbf{x}^{(0)}`$ aligns with $`\mathbf{v}_1`$. The convergence rate per iteration is approximately $`|\lambda_2 / \lambda_1|`$.

**Algorithm** (normalize by the $`\infty`$-norm to prevent overflow). Given $`\mathbf{x}^{(0)} = \mathbf{1}`$, repeat the following for $`k = 1, 2, \ldots`$:

1. Compute $`\mathbf{y}^{(k)} = A\,\mathbf{x}^{(k-1)}`$.
2. Take the eigenvalue estimate $`\lambda^{(k)} = y^{(k)}_{i^*}`$, where $`i^* = \arg\max_i |y_i^{(k)}|`$.
3. Normalize with $`\mathbf{x}^{(k)} = \mathbf{y}^{(k)} / \lambda^{(k)}`$ so that $`\|\mathbf{x}^{(k)}\|_\infty = 1`$.

---

<br>

### 2.2 Step-by-Step Trace

A hand trace makes the iteration concrete. Using the same $`A`$ as above with $`\mathbf{x}^{(0)} = \mathbf{1}`$:

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

Within a handful of iterations, `lambda est` approaches the true dominant eigenvalue and the normalized vector `x_k` approaches the corresponding eigenvector (with $`\|\mathbf{x}\|_\infty = 1`$). The estimate at iteration 7 is 21.838011, about 0.11 away from the true value 21.729156.

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

> **Note:** Normalizing by the element of largest absolute value (the $`\infty`$-norm) prevents overflow and keeps `lam` the eigenvalue estimate itself, so no separate Rayleigh quotient is needed.

---

<br>

### 2.4 Convergence Analysis

The theoretical convergence rate of the power method is $`|\lambda_2 / \lambda_1|`$: the closer this ratio is to 1, the slower the convergence.

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

The semi-log error plot should approach a straight line with slope $`\log|\lambda_2/\lambda_1|`$. For this matrix, $`|\lambda_2/\lambda_1| = 0.606665`$, so the error halves about every 1.4 iterations and the method converges in 51 iterations with `tol=1e-12`. Well-separated eigenvalues give rapid linear convergence; a small spectral gap that brings the ratio close to 1 (say $`|\lambda_2/\lambda_1| = 0.99`$) requires thousands of iterations to converge.

---

<br>

## 3. Inverse and Shifted Power Methods

### 3.1 Inverse Power Method for the Smallest Eigenvalue

If $`A\mathbf{v} = \lambda\mathbf{v}`$ and $`A`$ is invertible, then $`A^{-1}\mathbf{v} = (1/\lambda)\mathbf{v}`$. The smallest $`|\lambda_i|`$ of $`A`$ becomes the largest $`|1/\lambda_i|`$ of $`A^{-1}`$. Applying the power method to $`A^{-1}`$ therefore converges to $`1/\lambda_{\min}`$, and $`\lambda_{\min} = 1 / \mu_{\max}`$.

**Implementation:** never form $`A^{-1}`$ explicitly. Each iteration solves $`A\mathbf{y} = \mathbf{x}`$ using the LU factorization of $`A`$, computed once. Per-iteration cost is $`O(n^2)`$ after the $`O(n^3)`$ factorization.

---

<br>

### 3.2 Shifted Power Method for Any Eigenvalue

Apply the inverse power method to $`(A - \sigma I)`$ for a user-chosen **shift** $`\sigma`$, so the method is precisely a shifted inverse power method. The eigenvalues of $`(A - \sigma I)^{-1}`$ are $`1/(\lambda_i - \sigma)`$. The eigenvalue of $`A`$ closest to $`\sigma`$ has the largest $`|1/(\lambda_i - \sigma)|`$, so the iteration converges to it:

```math
\lambda = \sigma + \frac{1}{\mu}, \qquad \mu = \text{dominant eigenvalue of } (A - \sigma I)^{-1}
```

Setting $`\sigma = 0`$ recovers the inverse power method (finds smallest $`|\lambda|`$). Setting $`\sigma`$ near any target value picks out the eigenvalue closest to it, so interior eigenvalues, which are neither the largest nor the smallest, can also be found.

---

<br>

### 3.3 Python Implementation

A single function handles both the inverse power method ($`\sigma = 0`$) and the general shifted variant:

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

> **Note:** The LU factorization is computed once outside the loop, so each iteration only needs an $`O(n^2)`$ triangular solve. Do not use `np.linalg.inv(B) @ x` inside the loop, because it is both slower and less numerically accurate.

---

<br>

### 3.4 Finding All Eigenvalues with Different Shifts

By choosing shifts near each eigenvalue of interest, the shifted power method can recover every eigenvalue of a small matrix without calling `eig`:

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

The closer a shift is to an eigenvalue, the more $`|1/(\lambda_i - \sigma)|`$ dominates the other reciprocals, so the effective eigenvalue ratio shrinks and convergence speeds up. The convergence ratio is $`|\lambda_j - \sigma| / |\lambda_k - \sigma|`$, where $`\lambda_j`$ is the eigenvalue closest to $`\sigma`$ and $`\lambda_k`$ the second closest. In the output above, $`\sigma = 24`$ gives a ratio of about 0.21 and converges in 18 iterations, whereas $`\sigma = 10`$ lies almost midway between the eigenvalues 6.09 and 13.18, gives a ratio of about 0.81, and needs 127 iterations. In practice, a rough $`\sigma`$ is chosen with tools such as Gershgorin discs, and convergence takes only a handful of iterations when $`\sigma`$ is much closer to the target eigenvalue than to any other.

---

<br>

## 4. Method Comparison and Practical Guidelines

### 4.1 Method Comparison

| Method | Finds | Cost per iteration | Convergence rate |
|:-------|:------|:------------------|:-----------------|
| `np.linalg.eig(A)` | All $`n`$ eigenpairs | $`O(n^3)`$ total | Not applicable (direct computation) |
| **Power Method** | Dominant $`\lambda`$ (largest $`\vert\lambda\vert`$) | $`O(n^2)`$ | $`\vert\lambda_2 / \lambda_1\vert`$ |
| **Inverse Power** | Smallest $`\vert\lambda\vert`$ | $`O(n^2)`$ after LU | $`\vert\lambda_{\min} / \lambda_{\text{next-smallest}}\vert`$ |
| **Shifted Power** | $`\lambda`$ nearest to $`\sigma`$ | $`O(n^2)`$ after LU | depends on shift quality |

*Table 2. Comparison of eigenvalue computation methods*

### 4.2 Key Formulas

1. **Characteristic equation:** $`\det(A - \lambda I) = 0`$
2. **Trace identity:** $`\text{tr}(A) = \sum_i \lambda_i`$ (use as sanity check)
3. **Determinant identity:** $`\det(A) = \prod_i \lambda_i`$
4. **Power step:** $`\lambda^{(k)} = y^{(k)}_{i^*}`$ (the entry of $`\mathbf{y}^{(k)} = A\mathbf{x}^{(k-1)}`$ with the largest absolute value, sign included), $`\mathbf{x}^{(k)} = A\mathbf{x}^{(k-1)} / \lambda^{(k)}`$
5. **Shifted recovery:** $`\lambda = \sigma + 1/\mu`$ where $`\mu`$ is the dominant eigenvalue of $`(A - \sigma I)^{-1}`$

### 4.3 Practical Tips

- Always verify eigenpairs via the trace and determinant identities and the residual $`\|A\mathbf{v} - \lambda\mathbf{v}\|`$.
- Factor $`(A - \sigma I)`$ once outside the iteration loop, and never invert it explicitly.
- A small spectral gap means slow power method convergence, so shift to amplify the effective gap.
- Use the shifted power method to probe interior eigenvalues without computing the whole spectrum.

---

<br>

## 5. Symmetric Matrices and Spectral Decomposition

### 5.1 Special Properties

When $`A = A^T`$ (symmetric), its eigenstructure has several guarantees that the general case lacks:

| Property | Statement |
|:---------|:----------|
| Real eigenvalues | All $`\lambda_i \in \mathbb{R}`$ |
| Orthogonal eigenvectors | $`\mathbf{v}_i^T \mathbf{v}_j = 0`$ for $`i \neq j`$ (can be chosen orthogonal for repeated eigenvalues) |
| Orthonormal eigenbasis | $`V^T V = V V^T = I`$ |
| Spectral decomposition | $`A = V \Lambda V^T = \displaystyle\sum_{i=1}^n \lambda_i\, \mathbf{v}_i\mathbf{v}_i^T`$ |

*Table 3. Properties of the eigenstructure of a real symmetric matrix*

These guarantees are why specialized symmetric-eigenvalue algorithms (the Jacobi eigenvalue method, QR on tridiagonal form) are much more efficient than general ones.

Use **`np.linalg.eigh`** (not `eig`) whenever you know the matrix is symmetric (Hermitian in the complex case). It:

- Returns only real eigenvalues (no spurious imaginary parts from roundoff).
- Sorts eigenvalues in ascending order.
- Runs about 2× faster than `eig`.
- Returns an orthonormal eigenvector matrix.

---

<br>

### 5.2 `np.linalg.eigh` and Reconstruction

The spectral decomposition $`A = V \Lambda V^T`$ lets you reconstruct a symmetric matrix exactly from its eigenpairs:

```python
import numpy as np

A = np.array([[20.0,  3.0,  2.0],
              [ 3.0,  9.0,  4.0],
              [ 2.0,  4.0, 12.0]])

# eigh: sorted, real, orthonormal
eigvals_h, V = np.linalg.eigh(A)   # ascending order
print('eigh eigenvalues (ascending):', np.round(eigvals_h, 8))

print('\nOrthonormality check  V^T V:')
print(np.round(V.T @ V, 10))

# Spectral decomposition  A = V Lambda V^T
Lambda = np.diag(eigvals_h)
A_reconstructed = V @ Lambda @ V.T
print(f'\nmax|A - V Lambda V^T| = {np.max(np.abs(A - A_reconstructed)):.2e}')
```

The reconstruction error should be close to the level of machine epsilon ($`10^{-14}`$ or smaller). This identity underlies many applications: once you have $`V`$ and $`\Lambda`$, expressions like $`A^k`$, $`e^A`$, and $`A^{-1}`$ reduce to $`V f(\Lambda) V^T`$, a diagonal operation sandwiched between cheap orthogonal multiplications.

---

<br>

### 5.3 Low-Rank Approximation

With the eigenvalues ordered as $`|\lambda_1| \ge |\lambda_2| \ge \cdots \ge |\lambda_n|`$, the spectral sum $`A = \sum_{i=1}^n \lambda_i \mathbf{v}_i \mathbf{v}_i^T`$ can be truncated to the $`k`$ terms largest in absolute value:

```math
A_k = \sum_{i=1}^{k} \lambda_i\, \mathbf{v}_i\mathbf{v}_i^T, \qquad
\|A - A_k\|_F = \sqrt{\sum_{i=k+1}^{n} \lambda_i^2}
```

This is the **best rank-$`k`$ approximation** in the Frobenius norm (Eckart-Young theorem for symmetric $`A`$). It underlies principal component analysis, image compression, and noise removal:

```python
print('Low-rank approximation error ||A - A_k||_F:')
for k in [1, 2, 3]:
    idx = np.argsort(np.abs(eigvals_h))[::-1][:k]
    A_k = sum(eigvals_h[i] * np.outer(V[:, i], V[:, i]) for i in idx)
    err_F = np.linalg.norm(A - A_k, 'fro')
    err_theory = np.sqrt(np.sum(eigvals_h[np.argsort(np.abs(eigvals_h))[:-k]] ** 2))
    print(f'  k={k}: ||A - A_k||_F = {err_F:.6f}  (theory: {err_theory:.6f})')
```

As $`k`$ grows, the approximation error drops because only the smallest-magnitude eigenvalues are discarded. For this matrix the error is 14.520460 for $`k = 1`$ and 6.088518 for $`k = 2`$, and for $`k = n`$ the reconstruction is exact.

---

<br>

## 6. Application: Vibration Analysis

### 6.1 Equations of Motion and the Generalized Eigenproblem

Consider three masses connected by four springs with fixed walls at both ends:

```text
Wall ─[k₁]─ m₁ ─[k₂]─ m₂ ─[k₃]─ m₃ ─[k₄]─ Wall
```

*Figure 1. Three-mass, four-spring system fixed to walls at both ends*

Newton's second law gives the undamped free-vibration equation:

```math
M\ddot{\mathbf{x}} + K\mathbf{x} = \mathbf{0}
```

Assuming harmonic motion $`\mathbf{x}(t) = \mathbf{v}\sin(\omega t)`$ leads to the **generalized eigenvalue problem**:

```math
K\mathbf{v} = \omega^2 M\mathbf{v}
```

where $`\lambda = \omega^2`$ is the squared **natural frequency** and $`\mathbf{v}`$ is the corresponding **mode shape**. The stiffness and mass matrices are:

```math
K = \begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix},
\qquad
M = \begin{bmatrix} m_1 & & \\ & m_2 & \\ & & m_3 \end{bmatrix}
```

Both $`K`$ and $`M`$ are symmetric positive definite (SPD), so `scipy.linalg.eigh(K, M)` is the right tool: it returns real eigenvalues in ascending order and $`M`$-normalized mode shapes ($`\phi_i^T M \phi_i = 1`$).

---

<br>

### 6.2 Natural Frequencies via `scipy.linalg.eigh(K, M)`

From $`\omega^2`$ you obtain the angular frequency $`\omega = \sqrt{\lambda}`$ (rad/s), the linear frequency $`f = \omega / (2\pi)`$ (Hz), and the period $`T = 1/f`$ (s):

```python
import scipy as sc

# System parameters
k1, k2, k3, k4 = 40.0, 30.0, 20.0, 10.0   # N/m
m1, m2, m3     =  2.0,  1.5,  1.0          # kg

K_vib = np.array([
    [k1+k2,   -k2,      0],
    [  -k2, k2+k3,    -k3],
    [    0,   -k3,  k3+k4]
])
M_vib = np.diag([m1, m2, m3])

# Solve  K v = omega^2 M v
omega_sq, modes = sc.linalg.eigh(K_vib, M_vib)   # eigenvalues ascending

omega_n = np.sqrt(omega_sq)
freq_hz = omega_n / (2 * np.pi)

print('=== Natural Frequencies ===')
print(f'{"Mode":>5}  {"omega^2":>14}  {"omega (rad/s)":>14}  {"f (Hz)":>10}  {"T (s)":>10}')
print('-' * 60)
for i in range(3):
    print(f'{i+1:>5}  {omega_sq[i]:>14.4f}  {omega_n[i]:>14.4f}  {freq_hz[i]:>10.4f}  {1/freq_hz[i]:>10.4f}')

# Verify generalized eigenvalue equation
print('\nVerification  max|K v - omega^2 M v|:')
for i in range(3):
    v   = modes[:, i]
    res = K_vib @ v - omega_sq[i] * M_vib @ v
    print(f'  Mode {i+1}: {np.max(np.abs(res)):.2e}')
```

Running it gives the natural frequencies $`\omega = 3.0078`$, $`5.6889`$, and $`7.5447`$ rad/s (0.4787, 0.9054, and 1.2008 Hz).

> **Note:** The generalized eigenproblem can be reduced to a standard one through a transformation such as $`\tilde{K} = M^{-1/2} K M^{-1/2}`$. `scipy.linalg.eigh(K, M)` performs this kind of reduction internally with a Cholesky factorization of $`M`$ and returns mode shapes that are $`M`$-orthonormal: $`\Phi^T M \Phi = I`$.

---

<br>

### 6.3 Visualizing Mode Shapes

Each column of `modes` is a vector describing how the three masses deflect in that vibration mode. Plotting the deflections against node position makes the mode structure visible:

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(13, 5))
x_pos    = [0, 1, 2, 3, 4]   # wall, m1, m2, m3, wall
x_masses = [1, 2, 3]

for i, ax in enumerate(axes):
    v      = modes[:, i]
    v_norm = v / np.max(np.abs(v))   # scale to max amplitude = 1
    disp   = np.array([0.0, v_norm[0], v_norm[1], v_norm[2], 0.0])

    ax.axhline(0, color='gray', ls='--', lw=1)
    ax.plot(x_pos, disp, 'b-o', ms=13, lw=2,
            markerfacecolor='steelblue', markeredgecolor='navy')
    ax.plot([0, 4], [0, 0], 'ks', ms=14)   # walls

    ax.set(xlim=[-0.5, 4.5], ylim=[-1.7, 1.7],
           title=f'Mode {i+1}\nω = {omega_n[i]:.2f} rad/s  ({freq_hz[i]:.3f} Hz)',
           xlabel='Node', ylabel='Normalized displacement')
    ax.set_xticks(x_pos, ['Wall', 'm₁', 'm₂', 'm₃', 'Wall'])
    ax.grid(True, alpha=0.3)

plt.suptitle('Three-Mass Spring System — Mode Shapes', fontsize=13)
plt.tight_layout()
plt.show()
```

The lowest-frequency mode typically has all masses moving in the same direction, and higher modes introduce sign changes (nodes) along the chain. The number of sign changes grows with the mode number, a pattern familiar from vibrating-string harmonics.

---

<br>

### 6.4 Free-Vibration Time Response

Given initial conditions $`\mathbf{x}(0) = \mathbf{x}_0`$ and $`\dot{\mathbf{x}}(0) = \mathbf{v}_0`$, the time response is a superposition of the modes:

```math
\mathbf{x}(t) = \sum_{i=1}^{n}\! \left[c_i \cos(\omega_i t) + s_i \sin(\omega_i t)\right] \boldsymbol\phi_i
```

The modal amplitudes come from $`M`$-projections of the initial conditions onto the mode shapes: $`c_i = \boldsymbol\phi_i^T M \mathbf{x}_0`$ and $`s_i = \boldsymbol\phi_i^T M \mathbf{v}_0 / \omega_i`$.

```python
# Initial condition: displace m1 by 50 mm, release from rest
x0_vib  = np.array([0.05, 0.0, 0.0])
dx0_vib = np.zeros(3)

c_amp = modes.T @ M_vib @ x0_vib
s_amp = modes.T @ M_vib @ dx0_vib / (omega_n + 1e-30)

t = np.linspace(0, 3.0, 1000)
x_t = np.zeros((3, len(t)))
for i in range(3):
    x_t += np.outer(modes[:, i],
                    c_amp[i] * np.cos(omega_n[i] * t)
                    + s_amp[i] * np.sin(omega_n[i] * t))

fig, ax = plt.subplots(figsize=(10, 4))
for j in range(3):
    ax.plot(t, x_t[j] * 1000, lw=1.8, label=f'm_{j+1}')
ax.set(xlabel='Time (s)', ylabel='Displacement (mm)',
       title='Free Vibration: Initial Displacement of m₁ = 50 mm')
ax.legend(); ax.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()
```

The response is a sum of three pure sinusoids at the three natural frequencies. With no damping, the motion continues indefinitely; in a real system, energy would slowly leak out through damping, and the higher modes typically decay faster.

---

<br>

## 7. Application: ODE Systems and Stability

### 7.1 General Solution via Eigenpairs

For a linear first-order system $`\dot{\mathbf{y}} = A\mathbf{y}`$ with $`\mathbf{y}(0) = \mathbf{y}_0`$, if $`A`$ has $`n`$ linearly independent eigenvectors (is diagonalizable), the general solution is a sum of eigen-modes:

```math
\mathbf{y}(t) = \sum_{i=1}^{n} c_i\, \mathbf{v}_i\, e^{\lambda_i t}, \qquad
\mathbf{c} = V^{-1}\mathbf{y}_0
```

Each mode evolves independently along its eigenvector direction $`\mathbf{v}_i`$, scaled exponentially by $`e^{\lambda_i t}`$. The coefficient vector $`\mathbf{c}`$ is obtained by solving $`V\mathbf{c} = \mathbf{y}_0`$.

---

<br>

### 7.2 Stability Classification

The long-term behavior of every mode is determined by the real part of its eigenvalue:

| Condition | Behavior |
|:----------|:---------|
| All $`\text{Re}(\lambda_i) < 0`$ | Solution **decays** → **asymptotically stable** |
| Any $`\text{Re}(\lambda_i) > 0`$ | Solution **grows** → **unstable** |
| All $`\text{Re}(\lambda_i) \leq 0`$, some $`= 0`$ | **Marginally stable** (bounded but non-decaying) |

*Table 4. Stability classification by the real parts of the eigenvalues*

The third row assumes that $`A`$ is diagonalizable, as in Section 7.1. If an eigenvalue with zero real part lacks enough linearly independent eigenvectors, the solution can grow in proportion to $`t`$.

Imaginary parts of $`\lambda_i`$ produce oscillations at angular frequency $`\text{Im}(\lambda_i)`$, modulated by the exponential envelope $`e^{\text{Re}(\lambda_i) t}`$.

The **stiffness ratio** $`|\lambda_{\max} / \lambda_{\min}|`$ (of the real parts) measures how disparate the decay timescales are. Large stiffness ratios make a system numerically **stiff**. Explicit ODE integrators need very small steps to follow the fastest mode, even though the slow mode dominates the answer.

---

<br>

### 7.3 Example: Problem 13.13

Consider the system:

```math
\begin{cases} \dot{y}_1 = -5 y_1 + 3 y_2 \\ \dot{y}_2 = 100 y_1 - 301 y_2 \end{cases},
\qquad \mathbf{y}(0) = [50,\; 100]^T
```

```python
A_ode = np.array([[ -5.0,    3.0],
                  [100.0, -301.0]])
y0    = np.array([50.0, 100.0])

lam_ode, V_ode = np.linalg.eig(A_ode)
print('Eigenvalues:')
for i, lam in enumerate(lam_ode):
    status = 'decays' if np.real(lam) < 0 else 'GROWS'
    print(f'  lambda_{i+1} = {lam:.6f}   -> mode {status}')
stable = np.all(np.real(lam_ode) < 0)
print(f'System is {"STABLE" if stable else "UNSTABLE"}')

# Modal constants  c = V^{-1} y0
c_ode = np.linalg.solve(V_ode, y0)
print(f'\nModal constants: c1 = {c_ode[0]:.4f},  c2 = {c_ode[1]:.4f}')

# Time solution
t = np.linspace(0, 1.5, 800)
y_t = np.zeros((2, len(t)))
for i in range(2):
    y_t += np.real(c_ode[i] * np.outer(V_ode[:, i], np.exp(lam_ode[i] * t)))

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(t, y_t[0], 'b-', lw=2, label='y₁(t)')
ax.plot(t, y_t[1], 'r-', lw=2, label='y₂(t)')
ax.set(xlabel='Time t', ylabel='y', title='ODE Solution  ẏ = A y')
ax.legend(); ax.grid(True, alpha=0.4)

ax = axes[1]
ax.plot(y_t[0], y_t[1], 'k-', lw=2)
ax.plot(y0[0], y0[1], 'go', ms=12, label='Start y(0)', zorder=5)
ax.plot(y_t[0, -1], y_t[1, -1], 'rs', ms=12, label='End', zorder=5)
ax.set(xlabel='y₁', ylabel='y₂', title='Phase Portrait')
ax.legend(); ax.grid(True, alpha=0.4)

plt.tight_layout()
plt.show()

print(f'\nTime constants:  tau_1 = 1/|lambda_1| = {1/abs(lam_ode[0]):.4f} s  (slow mode)')
print(f'                 tau_2 = 1/|lambda_2| = {1/abs(lam_ode[1]):.6f} s  (fast mode)')
print(f'Stiffness ratio  |lambda_2 / lambda_1| = {abs(lam_ode[1]/lam_ode[0]):.1f}')
```

Both eigenvalues of this matrix are real, negative, and widely separated ($`\lambda_1 \approx -4`$, $`\lambda_2 \approx -302`$), so the system is stable but stiff. The printed time constants are $`\tau_1 \approx 0.2506`$ s and $`\tau_2 \approx 0.0033`$ s, and the stiffness ratio is about 75.7. The fast mode decays within a few hundredths of a second, and the slow mode dominates the visible transient. This is a classic warm-up for the stiff-ODE methods covered in Chapter 23.

---

<br>

## 8. Method Comparison and Engineering Context

### 8.1 Methods Overview

| Method | Finds | Cost | Notes |
|:-------|:------|:-----|:------|
| `np.linalg.eig(A)` | All $`n`$ eigenpairs | $`O(n^3)`$ | General matrices; may give complex values |
| `np.linalg.eigh(A)` | All $`n`$ eigenpairs | $`O(n^3)`$ | Symmetric only; real, sorted ascending |
| **Power Method** (Section 2) | Dominant $`\lambda`$ | $`O(n^2)`$/iter | Rate $`\approx \vert\lambda_2/\lambda_1\vert`$ |
| **Inverse and Shifted Power** (Section 3) | Smallest or nearest-to-$`\sigma`$ $`\lambda`$ | $`O(n^2)`$/iter (after LU) | Use for interior eigenvalues |
| `sc.linalg.eigh(K, M)` | Generalized $`K\mathbf{v} = \lambda M\mathbf{v}`$ | $`O(n^3)`$ | SPD pair; natural frequencies |

*Table 5. Overview of eigenvalue computation methods*

### 8.2 Key Formulas

| Concept | Formula |
|:--------|:--------|
| Spectral decomposition | $`A = V \Lambda V^T`$ (symmetric $`A`$, $`V^T V = I`$) |
| Low-rank approximation | $`A_k = \sum_{i=1}^k \lambda_i\, \mathbf{v}_i \mathbf{v}_i^T`$ |
| Frobenius error | $`\Vert A - A_k\Vert_F = \sqrt{\sum_{i=k+1}^n \lambda_i^2}`$ (eigenvalues sorted by decreasing absolute value) |
| Natural frequency | $`\omega_i = \sqrt{\lambda_i}`$ from $`K\mathbf{v} = \lambda M\mathbf{v}`$ |
| ODE solution | $`\mathbf{y}(t) = \sum_i c_i\, \mathbf{v}_i\, e^{\lambda_i t}, \quad \mathbf{c} = V^{-1}\mathbf{y}_0`$ |
| Stability criterion | Asymptotically stable $`\Leftrightarrow`$ all $`\text{Re}(\lambda_i) < 0`$ |

*Table 6. Key formulas for symmetric eigenproblems and applications*

### 8.3 Engineering Context

| Application | Eigenvalue meaning | Eigenvector meaning |
|:------------|:-------------------|:--------------------|
| Structural vibration | $`\lambda = \omega^2`$ (squared natural freq.) | Mode shape |
| ODE system | Decay / growth rate $`e^{\lambda t}`$ | Mode direction in state space |
| Iterative solvers (Chapter 12) | Spectral radius $`\rho(M)`$ of the iteration matrix $`M`$ | Convergence direction |
| Condition number (Chapter 11) | $`\kappa_2 = \vert\lambda\vert_{\max} / \vert\lambda\vert_{\min}`$ (symmetric matrix, 2-norm) | Not applicable |

*Table 7. Meaning of eigenvalues and eigenvectors in engineering applications*

### 8.4 Practical Tips

1. Always prefer `eigh` over `eig` for symmetric matrices. It is faster, real-valued, and sorted.
2. For large matrices where only a few eigenvalues are needed, use the power method (Section 2) or the shifted power method (Section 3) instead of the full `eig`.
3. Verify generalized eigenproblems by checking $`K\mathbf{v} - \lambda M\mathbf{v} \approx \mathbf{0}`$.
4. A large stiffness ratio $`|\lambda_{\max} / \lambda_{\min}|`$ signals a stiff ODE, for which implicit integrators are needed.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Characteristic equation | The eigenvalues are the roots of $`\det(A - \lambda I) = 0`$, and an $`n \times n`$ matrix has $`n`$ eigenvalues counting multiplicity. |
| Trace and determinant identities | The identities $`\text{tr}(A) = \sum_i \lambda_i`$ and $`\det(A) = \prod_i \lambda_i`$ verify computed eigenvalues. |
| `np.linalg.eig` | It computes all eigenpairs in $`O(n^3)`$, and the columns of the eigenvector matrix follow the order of the eigenvalue array. |
| Power method | Repeated multiplication by $`A`$ with normalization by the largest-magnitude entry gives the dominant eigenvalue and eigenvector at $`O(n^2)`$ per iteration. |
| Convergence ratio | The power method error shrinks by about $`\vert\lambda_2/\lambda_1\vert`$ per iteration, and a ratio close to 1 means slow convergence. |
| Inverse power method | The power method is applied to $`A^{-1}`$ by solving $`A\mathbf{y} = \mathbf{x}`$ with an LU factorization, which gives the smallest $`\vert\lambda\vert`$. |
| Shifted power method | The power method is applied to $`(A - \sigma I)^{-1}`$ and $`\lambda = \sigma + 1/\mu`$ recovers the eigenvalue closest to $`\sigma`$. |
| Properties of symmetric matrices | A real symmetric matrix has real eigenvalues and orthonormal eigenvectors and decomposes as $`A = V \Lambda V^T`$. |
| `np.linalg.eigh` | This function for symmetric matrices returns real eigenvalues sorted in ascending order and orthonormal eigenvectors. |
| Low-rank approximation | Keeping the $`k`$ terms largest in absolute value, $`A_k = \sum_{i=1}^k \lambda_i \mathbf{v}_i \mathbf{v}_i^T`$ is the best rank-$`k`$ approximation in the Frobenius norm, with error equal to the square root of the sum of the squared discarded eigenvalues. |
| Generalized eigenvalue problem | The vibration problem $`K\mathbf{v} = \omega^2 M\mathbf{v}`$ is solved with `scipy.linalg.eigh(K, M)`, and the natural frequencies are $`\omega_i = \sqrt{\lambda_i}`$. |
| Modal superposition | The free-vibration response is a sum of modal sinusoids whose amplitudes come from $`M`$-projections of the initial conditions. |
| ODE solution | For a diagonalizable $`A`$, the solution of $`\dot{\mathbf{y}} = A\mathbf{y}`$ is $`\mathbf{y}(t) = \sum_i c_i \mathbf{v}_i e^{\lambda_i t}`$ with $`\mathbf{c} = V^{-1}\mathbf{y}_0`$. |
| Stability | The system is asymptotically stable if every eigenvalue has a negative real part and unstable if any has a positive real part. |
| Stiffness ratio | A large $`\vert\lambda_{\max}/\lambda_{\min}\vert`$ of the real parts indicates a stiff system, for which explicit integrators need very small steps. |

---

<br>

## Self-Check Questions

1. **Power method by hand:** Apply three iterations of the power method to $`A = \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}`$ with $`\mathbf{x}^{(0)} = [1, 0]^T`$ and give the eigenvalue estimates. Where do the estimates converge, and what is the convergence ratio per iteration?

   <details>
   <summary>Answer</summary>

   Since $`A\mathbf{x}^{(0)} = [2, 1]^T`$, $`\lambda^{(1)} = 2`$ and $`\mathbf{x}^{(1)} = [1, 0.5]^T`$. Since $`A\mathbf{x}^{(1)} = [2.5, 2]^T`$, $`\lambda^{(2)} = 2.5`$ and $`\mathbf{x}^{(2)} = [1, 0.8]^T`$. Since $`A\mathbf{x}^{(2)} = [2.8, 2.6]^T`$, $`\lambda^{(3)} = 2.8`$. The eigenvalues of $`A`$ are 3 and 1, so the estimates converge to 3 with convergence ratio $`\vert\lambda_2/\lambda_1\vert = 1/3`$ per iteration.

   </details>

2. **Trace and determinant identities:** For the matrix $`A`$ of Section 1.3, given the two eigenvalues 21.729156 and 13.182326, find the remaining eigenvalue with the trace identity and check it with the determinant identity.

   <details>
   <summary>Answer</summary>

   Since $`\text{tr}(A) = 20 + 9 + 12 = 41`$, the remaining eigenvalue is $`41 - 21.729156 - 13.182326 = 6.088518`$. The product of the three eigenvalues is $`21.729156 \times 13.182326 \times 6.088518 \approx 1744.0`$, which matches $`\det(A) = 1744`$.

   </details>

3. **Choosing a shift:** If the shifted power method is applied to the same matrix with $`\sigma = 15`$, which eigenvalue does it converge to? Find $`\mu`$ at convergence and the convergence ratio per iteration.

   <details>
   <summary>Answer</summary>

   The distances from 15 are $`\vert 21.729 - 15\vert = 6.729`$, $`\vert 13.182 - 15\vert = 1.818`$, and $`\vert 6.089 - 15\vert = 8.911`$, so the method converges to the closest eigenvalue, 13.182326. Then $`\mu = 1/(13.182326 - 15) \approx -0.5502`$, and the convergence ratio is $`1.818/6.729 \approx 0.270`$.

   </details>

4. **Condition on the starting vector:** What happens in theory if the starting vector of the power method has no component along the dominant eigenvector $`\mathbf{v}_1`$ ($`c_1 = 0`$)? What happens in actual floating-point computation?

   <details>
   <summary>Answer</summary>

   In the expansion of Section 2.1, $`c_1 = 0`$ removes the $`\lambda_1^k`$ term, so in exact arithmetic the iteration converges to the eigenvalue of largest magnitude among the remaining components (e.g., $`\lambda_2`$ if $`c_2 \neq 0`$ and $`\vert\lambda_2\vert > \vert\lambda_3\vert`$). In floating-point computation, roundoff error gradually introduces a $`\mathbf{v}_1`$ component, so the iteration usually converges to $`\lambda_1`$ in the end, although it may need many more iterations.

   </details>

5. **Cost of the inverse power method:** For $`m`$ iterations of the inverse power method on a dense $`n \times n`$ matrix, compare the total cost of factoring the matrix once by LU with the cost of factoring it anew in every iteration.

   <details>
   <summary>Answer</summary>

   Factoring once costs one $`O(n^3)`$ factorization plus an $`O(n^2)`$ triangular solve per iteration, for a total of $`O(n^3 + mn^2)`$. Factoring in every iteration costs $`O(mn^3)`$ in total, and the gap widens as $`m`$ grows.

   </details>

6. **Low-rank approximation error:** The matrix of Section 5.2 has eigenvalues 21.729156, 13.182326, and 6.088518. Find the Frobenius errors of its rank-1 and rank-2 approximations.

   <details>
   <summary>Answer</summary>

   $`\Vert A - A_1\Vert_F = \sqrt{13.182326^2 + 6.088518^2} = \sqrt{173.774 + 37.070} = \sqrt{210.844} \approx 14.520`$ and $`\Vert A - A_2\Vert_F = \sqrt{6.088518^2} = 6.088518`$. These match the code output of Section 5.3, 14.520460 and 6.088518.

   </details>

7. **Generalized eigenvalue problem:** Two masses $`m_1 = m_2 = 2`$ kg are connected between two walls by three springs of stiffness $`k = 10`$ N/m. Write $`K`$ and $`M`$ and find the natural frequencies and mode shapes.

   <details>
   <summary>Answer</summary>

   $`K = \begin{bmatrix} 20 & -10 \\ -10 & 20 \end{bmatrix}`$ and $`M = \begin{bmatrix} 2 & 0 \\ 0 & 2 \end{bmatrix}`$. From $`\det(K - \omega^2 M) = (20 - 2\omega^2)^2 - 100 = 0`$, $`20 - 2\omega^2 = \pm 10`$, so $`\omega^2 = 5`$ or $`15`$. Hence $`\omega_1 = \sqrt{5} \approx 2.236`$ rad/s ($`f_1 \approx 0.356`$ Hz) and $`\omega_2 = \sqrt{15} \approx 3.873`$ rad/s ($`f_2 \approx 0.616`$ Hz). The mode shapes are $`[1, 1]^T`$ (in phase) for $`\omega_1`$ and $`[1, -1]^T`$ (out of phase) for $`\omega_2`$, the same result as Example 13.4 in Chapter 13.

   </details>

8. **Stability classification:** Determine the stability of $`\dot{\mathbf{y}} = A\mathbf{y}`$ for each matrix. (a) $`\begin{bmatrix} -2 & 1 \\ 0 & -5 \end{bmatrix}`$ (b) $`\begin{bmatrix} 0 & 1 \\ -4 & 0 \end{bmatrix}`$ (c) $`\begin{bmatrix} 1 & 2 \\ 2 & 1 \end{bmatrix}`$

   <details>
   <summary>Answer</summary>

   (a) is triangular, so its eigenvalues are the diagonal entries $`-2`$ and $`-5`$; both are negative, so the system is asymptotically stable (stiffness ratio 2.5). (b) has $`\lambda^2 + 4 = 0`$, so $`\lambda = \pm 2i`$; the system is marginally stable and oscillates without decay at angular frequency 2 (period $`\pi`$). (c) has $`(1 - \lambda)^2 - 4 = 0`$, so its eigenvalues are 3 and $`-1`$, and the positive eigenvalue 3 makes the system unstable.

   </details>

9. **Stiffness ratio and time constants:** For the system of Section 7.3, find the time constants $`\tau_i = 1/|\lambda_i|`$ and the stiffness ratio, and explain the difficulty they cause for explicit integrators.

   <details>
   <summary>Answer</summary>

   The trace is $`-306`$ and the determinant is $`1505 - 300 = 1205`$, so $`\lambda = \frac{-306 \pm \sqrt{306^2 - 4 \cdot 1205}}{2}`$ gives $`\lambda_1 \approx -3.990`$ and $`\lambda_2 \approx -302.010`$. The time constants are $`\tau_1 \approx 0.2506`$ s and $`\tau_2 \approx 0.00331`$ s, and the stiffness ratio is $`302.010/3.990 \approx 75.7`$. The fast mode vanishes within $`5\tau_2 \approx 0.017`$ s, but an explicit integrator must keep using small steps on the scale of $`\tau_2`$ for stability. Integrating up to $`5\tau_1 \approx 1.25`$ s, when the slow mode dies out, therefore takes a very large number of steps.

   </details>
