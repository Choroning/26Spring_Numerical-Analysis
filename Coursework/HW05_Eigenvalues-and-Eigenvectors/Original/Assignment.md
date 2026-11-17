# HW5 Assignment Handout — Eigenvalues and Eigenvectors

[Repository Guide](../../../README.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **Last Updated:** 2026-11-17

This is the Week 6 homework of Numerical Analysis. The task is to fill in the blanks of the lab notebook `Lab06-Eigen.ipynb` distributed by the instructor and run it. The notebook has six parts: eigenvalues and eigenvectors, the power method, the inverse and shifted power methods, symmetric matrices and spectral decomposition, an application to vibration analysis, and an application to ODE systems and stability. There are five blanks in total, each between two comment lines marked `COMPLETE`. The first code cell of the notebook imports `numpy` as `np`, `scipy` as `sc`, and `matplotlib.pyplot` as `plt`.

Each code block below reproduces a cell of the notebook, and the places where code is missing are the blanks to fill. The distributed blank notebook is [Lab06-Eigen.ipynb](Lab06-Eigen.ipynb).

- **Deadline:** one week after distribution
- **Submission:** Submit the notebook online after filling in the blanks and running it.

---

<br>

## Problem 1. Computing Eigenvalues and Eigenvectors (Notebook Section 1.3)

An eigenvalue $`\lambda`$ and an eigenvector $`\mathbf{v}`$ of a square matrix $`A`$ satisfy the following.

$$
A\mathbf{v} = \lambda\mathbf{v}, \qquad \mathbf{v} \neq \mathbf{0}
$$

Geometrically, $`A`$ only scales the eigenvector and does not rotate it. Rearranging gives $`(A - \lambda I)\mathbf{v} = \mathbf{0}`$, and a nontrivial solution requires $`A - \lambda I`$ to be singular, which yields the characteristic equation $`\det(A - \lambda I) = 0`$. For an $`n \times n`$ matrix this is a polynomial of degree $`n`$ in $`\lambda`$, yielding exactly $`n`$ eigenvalues counting multiplicity, possibly complex.

Section 1.3 of the notebook uses `numpy.linalg.eig` to find the eigenvalues and eigenvectors of the lecture-notes example matrix (Chapter 13, Problem 13.1). Fill in the blank of the next cell.

```python
# Lecture-notes example  (Ch. 13, Problem 13.1)
A = np.array([[20.0,  3.0,  2.0],
              [ 3.0,  9.0,  4.0],
              [ 2.0,  4.0, 12.0]])

# --- COMPLETE ---
eigvals, eigvecs = 
# --- COMPLETE ---

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

The distributed notebook shows the following output for this cell.

```text
Eigenvalues:
  lambda_1 = 21.72915618
  lambda_2 = 13.18232626
  lambda_3 = 6.08851756

Eigenvectors (columns):
[[ 0.89819   0.427012  0.104481]
 [ 0.30972  -0.44601  -0.839732]
 [ 0.311976 -0.786598  0.532855]]

Verification  max|A v - lambda v|:
  lambda_1: 3.55e-15
  lambda_2: 7.11e-15
  lambda_3: 1.33e-15
```

## Problem 2. Power Method (Notebook Section 2.1)

`np.linalg.eig` computes all $`n`$ eigenpairs, an $`O(n^3)`$ operation. For large matrices often only the dominant eigenvalue (largest $`|\lambda|`$) is needed, and the **power method** computes it in $`O(n^2)`$ per iteration.

The reason it works is as follows. Expanding $`\mathbf{x}^{(0)}`$ in the eigenbasis as $`\mathbf{x}^{(0)} = c_1\mathbf{v}_1 + c_2\mathbf{v}_2 + \cdots`$ gives the following.

```math
A^k\mathbf{x}^{(0)} = \lambda_1^k\!\left[c_1\mathbf{v}_1 + c_2\!\left(\frac{\lambda_2}{\lambda_1}\right)^{\!k}\mathbf{v}_2 + \cdots\right]
```

If $`|\lambda_1| > |\lambda_2|`$, then $`(\lambda_2/\lambda_1)^k \to 0`$, so $`A^k\mathbf{x}^{(0)}`$ aligns with $`\mathbf{v}_1`$. The convergence rate is about $`|\lambda_2/\lambda_1|`$ per iteration.

The algorithm normalizes by the infinity norm to prevent overflow. Starting from $`\mathbf{x}^{(0)} = \mathbf{1}`$, it repeats the following for $`k = 1, 2, \ldots`$.

1. $`\mathbf{y}^{(k)} = A\,\mathbf{x}^{(k-1)}`$
2. $`\lambda^{(k)} = y^{(k)}_{i^*}`$, where $`i^* = \arg\max_i |y_i^{(k)}|`$ (eigenvalue estimate)
3. $`\mathbf{x}^{(k)} = \mathbf{y}^{(k)} / \lambda^{(k)}`$ (normalize so that $`\|\mathbf{x}^{(k)}\|_\infty = 1`$)

Fill in the blanks of the following function to complete the power method.

```python
def power_method(A, x0=None, tol=1e-8, max_iter=500):
    """
    Find the dominant eigenvalue and eigenvector of A using the Power Method.

    Parameters
    ----------
    A        : (n, n) array
    x0       : (n,)   initial guess  (default: ones)
    tol      : float  relative change in lambda as stopping criterion
    max_iter : int

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
        # --- COMPLETE ---
        lam =                                         
        x   =                                         
        # --- COMPLETE ---

        lam_hist.append(float(lam))
        if abs(lam - lam_old) / abs(lam) < tol:
            break
        lam_old = lam

    return float(lam), x, lam_hist
```

The next cell applies the completed function to the matrix `A` of Problem 1 and analyzes its convergence with plots. The plotting part is omitted here.

```python
# Test power_method and analyze convergence
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
```

The distributed notebook shows the following output for this cell.

```text
Converged in 51 iterations
Dominant eigenvalue (Power Method): 21.729156182199
True value (numpy)                : 21.729156182169
Absolute error                    : 3.02e-11

Theoretical convergence rate |λ₂/λ₁| = 0.606665
  -> error halves every 1.4 iterations
```

## Problem 3. Inverse and Shifted Power Methods (Notebook Sections 3.1 and 3.2)

If $`A\mathbf{v} = \lambda\mathbf{v}`$, then $`A^{-1}\mathbf{v} = (1/\lambda)\mathbf{v}`$. The smallest $`|\lambda_i|`$ of $`A`$ therefore becomes the largest $`|1/\lambda_i|`$ of $`A^{-1}`$, and applying the power method to $`A^{-1}`$ converges to $`1/\lambda_{\min}`$. The implementation never forms $`A^{-1}`$ explicitly. Each iteration solves $`A\mathbf{y} = \mathbf{x}`$, which costs only $`O(n^2)`$ once the LU factorization has been computed.

Applying the inverse power method to $`(A - \sigma I)`$ for a user-supplied **shift** $`\sigma`$ gives the following.

$$
\text{eigenvalues of }(A-\sigma I)^{-1}: \quad \frac{1}{\lambda_i - \sigma}
$$

The eigenvalue closest to $`\sigma`$ has the largest $`|1/(\lambda_i - \sigma)|`$, so the method converges to it. Setting $`\sigma = 0`$ recovers the inverse power method. Fill in the blanks of the following function.

```python
def shifted_power_method(A, sigma=0.0, x0=None, tol=1e-8, max_iter=500):
    """
    Find the eigenvalue of A nearest to `sigma` (Shifted Inverse Power Method).

    Algorithm:
      1. Build B = A - sigma*I
      2. Factor B = LU  (once)
      3. Power iteration on B^{-1}: solve B y = x, normalize
      4. True eigenvalue:  lambda = sigma + 1/mu   where mu = max(|y|)

    sigma = 0  ->  Inverse Power Method  (finds smallest |lambda|)
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]
    x = np.ones(n) if x0 is None else np.array(x0, dtype=float)

    # --- COMPLETE ---
    B      = A - sigma * np.eye(n)                   # shifted matrix
    lu_obj = 
    # --- COMPLETE ---

    lam_hist = []
    lam_old  = np.inf

    for _ in range(max_iter):
        # --- COMPLETE ---
        y   = 
        mu  = y[np.argmax(np.abs(y))]                # largest element of y  (= 1 / (lambda - sigma))
        x   = y / mu
        lam = 
        # --- COMPLETE ---

        lam_hist.append(float(lam))
        if abs(lam - lam_old) / (abs(lam) + 1e-30) < tol:
            break
        lam_old = lam

    return float(lam), x, lam_hist
```

The next cell finds all three eigenvalues with different shifts.

```python
# Find all three eigenvalues using different shifts
lam_sorted = np.sort(np.real(eigvals))
print(f'All eigenvalues (numpy, ascending): {np.round(lam_sorted, 6)}')
print()

tests = [
    ('Power Method   (dominant)', power_method,           dict(),           np.max(np.abs(eigvals))),
    ('sigma=0        (smallest)', shifted_power_method,   dict(sigma=0.0),  lam_sorted[0]),
    ('sigma=10       (middle)',   shifted_power_method,   dict(sigma=10.0), lam_sorted[1]),
    ('sigma=24       (largest)', shifted_power_method,   dict(sigma=24.0), lam_sorted[2]),
]

print(f'{"Method":<34}  {"Found lambda":>14}  {"True lambda":>12}  {"Error":>10}  {"Iters":>6}')
print('-' * 84)
for label, func, kwargs, lam_t in tests:
    lam_f, _, hist = func(A, tol=1e-12, **kwargs)
    print(f'{label:<34}  {lam_f:>14.8f}  {lam_t:>12.8f}  {abs(lam_f-lam_t):>10.2e}  {len(hist):>6}')
```

The distributed notebook shows the following output for this cell.

```text
All eigenvalues (numpy, ascending): [ 6.088518 13.182326 21.729156]

Method                                Found lambda   True lambda       Error   Iters
------------------------------------------------------------------------------------
Power Method   (dominant)              21.72915618   21.72915618    3.02e-11      51
sigma=0        (smallest)               6.08851756    6.08851756    2.68e-12      38
sigma=10       (middle)                13.18232626   13.18232626    5.06e-12     127
sigma=24       (largest)               21.72915618   21.72915618    1.34e-12      18
```

## Problem 4. Natural Frequencies of a Vibrating System (Notebook Section 5.1)

Three masses are connected by four springs, with fixed walls at both ends.

```text
Wall ─[k₁]─ m₁ ─[k₂]─ m₂ ─[k₃]─ m₃ ─[k₄]─ Wall
```

Newton's second law gives the free-vibration equation $`M\ddot{\mathbf{x}} + K\mathbf{x} = \mathbf{0}`$. Assuming harmonic motion $`\mathbf{x}(t) = \mathbf{v}\sin(\omega t)`$ leads to the following **generalized eigenvalue problem**.

$$
K\mathbf{v} = \omega^2 M\mathbf{v}
$$

Here $`\lambda = \omega^2`$ is the squared natural frequency and $`\mathbf{v}`$ is the mode shape. The stiffness and mass matrices are as follows.

```math
K = \begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix},
\qquad
M = \begin{bmatrix} m_1 & & \\ & m_2 & \\ & & m_3 \end{bmatrix}
```

`scipy.linalg.eigh(K, M)` solves the generalized eigenvalue problem for a symmetric positive definite pair and returns the eigenvalues in ascending order. Fill in the angular frequency `omega_n` (rad/s) and the frequency `freq_hz` (Hz) in the blanks of the next cell.

```python
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

# --- COMPLETE ---
omega_n =                                        #  — rad/s
freq_hz =                                        #  — Hz
# --- COMPLETE ---

print('=== Natural Frequencies ===')
print(f'{"Mode":>5}  {"omega^2":>14}  {"omega (rad/s)":>14}  {"f (Hz)":>10}  {"T (s)":>10}')
print('-' * 60)
for i in range(3):
    print(f'{i+1:>5}  {omega_sq[i]:>14.4f}  {omega_n[i]:>14.4f}  {freq_hz[i]:>10.4f}  {1/freq_hz[i]:>10.4f}')

print('\nMode shapes (columns = eigenvectors, M-normalized):')
print(np.round(modes, 6))

# Verify generalized eigenvalue equation
print('\nVerification  max|K v - omega^2 M v|:')
for i in range(3):
    v   = modes[:, i]
    res = K_vib @ v - omega_sq[i] * M_vib @ v
    print(f'  Mode {i+1}: {np.max(np.abs(res)):.2e}')
```

The distributed notebook shows the following output for this cell.

```text
=== Natural Frequencies ===
 Mode         omega^2   omega (rad/s)      f (Hz)       T (s)
------------------------------------------------------------
    1          9.0470          3.0078      0.4787      2.0889
    2         32.3636          5.6889      0.9054      1.1045
    3         56.9227          7.5447      1.2008      0.8328

Mode shapes (columns = eigenvectors, M-normalized):
[[-0.329371  0.484607 -0.395817]
 [-0.569877  0.085174  0.578492]
 [-0.543958 -0.720715 -0.429743]]

Verification  max|K v - omega^2 M v|:
  Mode 1: 4.44e-15
  Mode 2: 7.11e-15
  Mode 3: 1.42e-14
```

## Problem 5. Modal Constants of an ODE System (Notebook Sections 6.1 and 6.2)

The solution of $`\dot{\mathbf{y}} = A\mathbf{y}`$ with $`\mathbf{y}(0) = \mathbf{y}_0`$ is as follows.

```math
\mathbf{y}(t) = \sum_{i=1}^{n} c_i\, \mathbf{v}_i\, e^{\lambda_i t}, \qquad \mathbf{c} = V^{-1}\mathbf{y}_0
```

Stability is determined by the real parts of the eigenvalues.

| Condition | Behavior |
|:-----|:-----|
| All $`\text{Re}(\lambda_i) < 0`$ | The solution decays and is asymptotically stable |
| Any $`\text{Re}(\lambda_i) > 0`$ | The solution grows and is unstable |
| All $`\text{Re}(\lambda_i) \leq 0`$ with some $`= 0`$ | The system is marginally stable |

*Table 1. Stability of a linear ODE system by the real parts of the eigenvalues*

The example (Lecture Problem 13.13) is the following system.

```math
\begin{cases} \dot{y}_1 = -5y_1 + 3y_2 \\ \dot{y}_2 = 100y_1 - 301y_2 \end{cases},
\qquad \mathbf{y}(0) = [50,\; 100]^T
```

Fill in the blank of the next cell with code that finds the modal constants $`\mathbf{c} = V^{-1}\mathbf{y}_0`$.

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

# --- COMPLETE: find modal constants c = V^{-1} y0 ---
c_ode = 
# --- COMPLETE ---

print(f'\nModal constants: c1 = {c_ode[0]:.4f},  c2 = {c_ode[1]:.4f}')

# Time solution
t = np.linspace(0, 1.5, 800)
y_t = np.zeros((2, len(t)))
for i in range(2):
    y_t += np.real(c_ode[i] * np.outer(V_ode[:, i], np.exp(lam_ode[i] * t)))

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

ax = axes[0]
ax.plot(t, y_t[0], 'b-', lw=2, label='$y_1(t)$')
ax.plot(t, y_t[1], 'r-', lw=2, label='$y_2(t)$')
ax.set(xlabel='Time $t$', ylabel='$y$', title='ODE Solution $\\dot{\\mathbf{y}} = A\\mathbf{y}$')
ax.legend(); ax.grid(True, alpha=0.4)

ax = axes[1]
ax.plot(y_t[0], y_t[1], 'k-', lw=2)
ax.plot(y0[0], y0[1], 'go', ms=12, label='Start $\\mathbf{y}(0)$', zorder=5)
ax.plot(y_t[0,-1], y_t[1,-1], 'rs', ms=12, label='End', zorder=5)
ax.set(xlabel='$y_1$', ylabel='$y_2$', title='Phase Portrait')
ax.legend(); ax.grid(True, alpha=0.4)

plt.tight_layout()
plt.show()

print(f'\nTime constants:  tau_1 = 1/|lambda_1| = {1/abs(lam_ode[0]):.4f} s  (slow mode)')
print(f'                 tau_2 = 1/|lambda_2| = {1/abs(lam_ode[1]):.6f} s  (fast mode)')
print(f'Stiffness ratio  |lambda_2 / lambda_1| = {abs(lam_ode[1]/lam_ode[0]):.1f}')
```

The distributed notebook shows the following text output for this cell. Between the two printed parts there are plots of the solution over time and of the phase portrait.

```text
Eigenvalues:
  lambda_1 = -3.989933   -> mode decays
  lambda_2 = -302.010067   -> mode decays
System is STABLE

Modal constants: c1 = 53.6413,  c2 = 82.8879

Time constants:  tau_1 = 1/|lambda_1| = 0.2506 s  (slow mode)
                 tau_2 = 1/|lambda_2| = 0.003311 s  (fast mode)
Stiffness ratio  |lambda_2 / lambda_1| = 75.7
```

---

<br>

After solving the problems, see the [assignment explanation](../Assignment-Explanation.md) for the answers and solutions and the [completed notebook](../Solution/Lab06-Eigen_Completed.ipynb) for the full code.
