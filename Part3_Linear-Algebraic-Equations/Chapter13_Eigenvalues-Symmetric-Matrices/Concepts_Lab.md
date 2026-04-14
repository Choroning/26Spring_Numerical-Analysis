# Chapter 13 Lab — Symmetric Eigenproblems and Applications

> **Last Updated:** 2026-04-14
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 13

> **Prerequisites**: [Programming Language] Python, NumPy, SciPy, Matplotlib. [Linear Algebra] Power Method (Ch 12 Lab), eigen-decomposition basics.
>
> **Learning Objectives**:
> 1. Exploit the special properties of symmetric matrices using `np.linalg.eigh`
> 2. Use the spectral decomposition $A = V\Lambda V^T$ to build low-rank approximations
> 3. Solve the generalized eigenvalue problem $K\mathbf{v} = \lambda M\mathbf{v}$ with `scipy.linalg.eigh(K, M)`
> 4. Compute natural frequencies and mode shapes of a mass–spring system
> 5. Analyze the stability of a linear ODE system $\dot{\mathbf{y}} = A\mathbf{y}$ from its eigenvalues

---

<br>

## Table of Contents

- [1. Symmetric Matrices and Spectral Decomposition](#1-symmetric-matrices-and-spectral-decomposition)
  - [1.1 Special Properties](#11-special-properties)
  - [1.2 `np.linalg.eigh` and Reconstruction](#12-nplinalgeigh-and-reconstruction)
  - [1.3 Low-Rank Approximation](#13-low-rank-approximation)
- [2. Application: Vibration Analysis](#2-application-vibration-analysis)
  - [2.1 Equations of Motion and the Generalized Eigenproblem](#21-equations-of-motion-and-the-generalized-eigenproblem)
  - [2.2 Natural Frequencies via `scipy.linalg.eigh(K, M)`](#22-natural-frequencies-via-scipylinalgeighk-m)
  - [2.3 Visualizing Mode Shapes](#23-visualizing-mode-shapes)
  - [2.4 Free-Vibration Time Response](#24-free-vibration-time-response)
- [3. Application: ODE Systems and Stability](#3-application-ode-systems-and-stability)
  - [3.1 General Solution via Eigenpairs](#31-general-solution-via-eigenpairs)
  - [3.2 Stability Classification](#32-stability-classification)
  - [3.3 Example: Lecture Problem 13.13](#33-example-lecture-problem-1313)
- [Summary](#summary)

---

<br>

## 1. Symmetric Matrices and Spectral Decomposition

### 1.1 Special Properties

When $A = A^T$ (symmetric), its eigenstructure has several guarantees that the general case lacks:

| Property | Statement |
|:---------|:----------|
| Real eigenvalues | All $\lambda_i \in \mathbb{R}$ |
| Orthogonal eigenvectors | $\mathbf{v}_i^T \mathbf{v}_j = 0$ for $i \neq j$ |
| Orthonormal eigenbasis | $V^T V = V V^T = I$ |
| Spectral decomposition | $A = V \Lambda V^T = \displaystyle\sum_{i=1}^n \lambda_i\, \mathbf{v}_i\mathbf{v}_i^T$ |

These guarantees are why specialized symmetric-eigenvalue algorithms (Jacobi, QR on tridiagonal form) are much more efficient than general ones.

Use **`np.linalg.eigh`** (not `eig`) whenever you know the matrix is symmetric (Hermitian in the complex case). It:

- Returns only real eigenvalues (no spurious imaginary parts from roundoff)
- Sorts eigenvalues in ascending order
- Runs about 2× faster than `eig`
- Returns an orthonormal eigenvector matrix

---

<br>

### 1.2 `np.linalg.eigh` and Reconstruction

The spectral decomposition $A = V \Lambda V^T$ lets you reconstruct a symmetric matrix exactly from its eigenpairs:

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

The reconstruction error should be on the order of machine precision ($10^{-14}$ or smaller). This identity underlies many applications: once you have $V$ and $\Lambda$, expressions like $A^k$, $e^A$, and $A^{-1}$ reduce to $V f(\Lambda) V^T$ — a diagonal operation sandwiched between cheap orthogonal multiplications.

---

<br>

### 1.3 Low-Rank Approximation

The spectral sum $A = \sum_{i=1}^n \lambda_i \mathbf{v}_i \mathbf{v}_i^T$ can be truncated to the $k$ dominant terms:

$$A_k = \sum_{i=1}^{k} \lambda_i\, \mathbf{v}_i\mathbf{v}_i^T, \qquad
\|A - A_k\|_F = \sqrt{\sum_{i=k+1}^{n} \lambda_i^2}$$

This is the **best rank-$k$ approximation** in the Frobenius norm (Eckart–Young theorem for symmetric $A$). It underlies principal component analysis, image compression, and noise removal:

```python
print('Low-rank approximation error ||A - A_k||_F:')
for k in [1, 2, 3]:
    idx = np.argsort(np.abs(eigvals_h))[::-1][:k]
    A_k = sum(eigvals_h[i] * np.outer(V[:, i], V[:, i]) for i in idx)
    err_F = np.linalg.norm(A - A_k, 'fro')
    err_theory = np.sqrt(np.sum(eigvals_h[np.argsort(np.abs(eigvals_h))[:-k]] ** 2))
    print(f'  k={k}: ||A - A_k||_F = {err_F:.6f}  (theory: {err_theory:.6f})')
```

As $k$ grows, the approximation error drops by discarding the smallest-magnitude eigenvalues. For $k = n$, the reconstruction is exact.

---

<br>

## 2. Application: Vibration Analysis

### 2.1 Equations of Motion and the Generalized Eigenproblem

Consider three masses connected by four springs with fixed walls at both ends:

```
Wall ─[k₁]─ m₁ ─[k₂]─ m₂ ─[k₃]─ m₃ ─[k₄]─ Wall
```

Newton's second law gives the undamped free-vibration equation:

$$M\ddot{\mathbf{x}} + K\mathbf{x} = \mathbf{0}$$

Assuming harmonic motion $\mathbf{x}(t) = \mathbf{v}\sin(\omega t)$ leads to the **generalized eigenvalue problem**:

$$K\mathbf{v} = \omega^2 M\mathbf{v}$$

where $\lambda = \omega^2$ is the squared **natural frequency** and $\mathbf{v}$ is the corresponding **mode shape**. The stiffness and mass matrices are

$$K = \begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix},
\qquad
M = \begin{bmatrix} m_1 & & \\ & m_2 & \\ & & m_3 \end{bmatrix}$$

Both $K$ and $M$ are symmetric positive definite (SPD), so `scipy.linalg.eigh(K, M)` is the right tool: it returns real eigenvalues in ascending order and $M$-normalized mode shapes ($\phi_i^T M \phi_i = 1$).

---

<br>

### 2.2 Natural Frequencies via `scipy.linalg.eigh(K, M)`

From $\omega^2$ you obtain the angular frequency $\omega = \sqrt{\lambda}$ (rad/s), the linear frequency $f = \omega / (2\pi)$ (Hz), and the period $T = 1/f$ (s):

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

> **Note:** The generalized eigenproblem reduces to a standard one via $\tilde{K} = M^{-1/2} K M^{-1/2}$, but `scipy.linalg.eigh(K, M)` does this internally and returns mode shapes that are $M$-orthonormal: $\Phi^T M \Phi = I$.

---

<br>

### 2.3 Visualizing Mode Shapes

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

The lowest-frequency mode typically has all masses moving in the same direction; higher modes introduce sign changes (nodes) along the chain. The number of sign changes grows with the mode number — a pattern familiar from vibrating-string harmonics.

---

<br>

### 2.4 Free-Vibration Time Response

Given initial conditions $\mathbf{x}(0) = \mathbf{x}_0$ and $\dot{\mathbf{x}}(0) = \mathbf{v}_0$, the time response is a superposition of the modes:

$$\mathbf{x}(t) = \sum_{i=1}^{n}\! \left[c_i \cos(\omega_i t) + s_i \sin(\omega_i t)\right] \boldsymbol\phi_i$$

The modal amplitudes come from $M$-projections of the initial conditions onto the mode shapes: $c_i = \boldsymbol\phi_i^T M \mathbf{x}_0$ and $s_i = \boldsymbol\phi_i^T M \mathbf{v}_0 / \omega_i$.

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

The response is a sum of three pure sinusoids at the three natural frequencies. With no damping, the motion continues indefinitely; in a real system, energy would slowly leak out through damping and the higher modes decay fastest.

---

<br>

## 3. Application: ODE Systems and Stability

### 3.1 General Solution via Eigenpairs

For a linear first-order system $\dot{\mathbf{y}} = A\mathbf{y}$ with $\mathbf{y}(0) = \mathbf{y}_0$, the general solution is a sum of eigen-modes:

$$\mathbf{y}(t) = \sum_{i=1}^{n} c_i\, \mathbf{v}_i\, e^{\lambda_i t}, \qquad
\mathbf{c} = V^{-1}\mathbf{y}_0$$

Each mode evolves independently along its eigenvector direction $\mathbf{v}_i$, scaled exponentially by $e^{\lambda_i t}$. The coefficient vector $\mathbf{c}$ is obtained by solving $V\mathbf{c} = \mathbf{y}_0$.

---

<br>

### 3.2 Stability Classification

The long-term behavior of every mode is determined by the real part of its eigenvalue:

| Condition | Behavior |
|:----------|:---------|
| All $\text{Re}(\lambda_i) < 0$ | Solution **decays** → **asymptotically stable** |
| Any $\text{Re}(\lambda_i) > 0$ | Solution **grows** → **unstable** |
| All $\text{Re}(\lambda_i) \leq 0$, some $= 0$ | **Marginally stable** (bounded but non-decaying) |

Imaginary parts of $\lambda_i$ produce oscillations at frequency $\text{Im}(\lambda_i)$, modulated by the exponential envelope $e^{\text{Re}(\lambda_i) t}$.

The **stiffness ratio** $|\lambda_{\max} / \lambda_{\min}|$ (of the real parts) measures how disparate the decay timescales are. Large stiffness ratios make a system numerically "stiff": explicit ODE integrators need very small steps to follow the fastest mode, even though the slow mode dominates the answer.

---

<br>

### 3.3 Example: Lecture Problem 13.13

Consider the system

$$\begin{cases} \dot{y}_1 = -5 y_1 + 3 y_2 \\ \dot{y}_2 = 100 y_1 - 301 y_2 \end{cases},
\qquad \mathbf{y}(0) = [50,\; 100]^T$$

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

Both eigenvalues of this matrix are real, negative, and widely separated ($\lambda_1 \approx -4$, $\lambda_2 \approx -302$), so the system is stable but stiff. The fast mode decays within a few hundredths of a second; the slow mode dominates the visible transient. This is a classic warm-up for the stiff-ODE methods covered later in the course (Part 6, Ch 23).

---

<br>

## Summary

### Methods Overview

| Method | Finds | Cost | Notes |
|:-------|:------|:-----|:------|
| `np.linalg.eig(A)` | All $n$ eigenpairs | $O(n^3)$ | General matrices; may give complex values |
| `np.linalg.eigh(A)` | All $n$ eigenpairs | $O(n^3)$ | Symmetric only; real, sorted ascending |
| **Power Method** (Ch 12 Lab) | Dominant $\lambda$ | $O(n^2)$/iter | Rate $\approx \|\lambda_2/\lambda_1\|$ |
| **Inverse / Shifted Power** (Ch 12 Lab) | Smallest or nearest-to-$\sigma$ $\lambda$ | $O(n^2)$/iter (after LU) | Use for interior eigenvalues |
| `sc.linalg.eigh(K, M)` | Generalized $K\mathbf{v} = \lambda M\mathbf{v}$ | $O(n^3)$ | SPD pair; natural frequencies |

### Key Formulas

| Concept | Formula |
|:--------|:--------|
| Spectral decomposition | $A = V \Lambda V^T$ (symmetric $A$, $V^T V = I$) |
| Low-rank approximation | $A_k = \sum_{i=1}^k \lambda_i\, \mathbf{v}_i \mathbf{v}_i^T$ |
| Frobenius error | $\|A - A_k\|_F = \sqrt{\sum_{i=k+1}^n \lambda_i^2}$ |
| Natural frequency | $\omega_i = \sqrt{\lambda_i}$ from $K\mathbf{v} = \lambda M\mathbf{v}$ |
| ODE solution | $\mathbf{y}(t) = \sum_i c_i\, \mathbf{v}_i\, e^{\lambda_i t}, \quad \mathbf{c} = V^{-1}\mathbf{y}_0$ |
| Stability criterion | Stable $\Leftrightarrow$ all $\text{Re}(\lambda_i) < 0$ |

### Engineering Context

| Application | Eigenvalue meaning | Eigenvector meaning |
|:------------|:-------------------|:--------------------|
| Structural vibration | $\lambda = \omega^2$ (squared natural freq.) | Mode shape |
| ODE system | Decay / growth rate $e^{\lambda t}$ | Mode direction in state space |
| Iterative solvers (Ch 12) | Spectral radius $\rho(T)$ | Convergence direction |
| Condition number (Ch 11) | $\kappa = \lambda_{\max} / \lambda_{\min}$ | — |

Practical tips:

1. Always prefer `eigh` over `eig` for symmetric matrices — faster, real-valued, sorted
2. For large matrices use Power / Shifted Power (Ch 12 Lab) instead of the full `eig`
3. Verify generalized eigenproblems by checking $K\mathbf{v} - \lambda M\mathbf{v} \approx \mathbf{0}$
4. A large stiffness ratio $|\lambda_{\max} / \lambda_{\min}|$ signals a stiff ODE — implicit integrators needed

---
