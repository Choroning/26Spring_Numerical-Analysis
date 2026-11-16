# HW2 Assignment Handout — Linear Systems and Gaussian Elimination

[Repository Guide](../../../README.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **Last Updated:** 2026-11-16

This is the Week 3 homework of Numerical Analysis. The task is to fill in the blanks of the lab notebook `Lab03-LinearSystem.ipynb` distributed by the instructor and run it. The notebook has four parts: setting up and solving $`A\mathbf{x} = \mathbf{b}`$, naive Gaussian elimination, partial pivoting, and tridiagonal systems with the Thomas algorithm. There are four blanks in total, one in each part, each between two comment lines marked `COMPLETE`. The first code cell of the notebook imports `numpy` as `np`, `matplotlib.pyplot` as `plt`, and `time`.

Each code block below reproduces a cell of the notebook, and the places where code is missing are the blanks to fill. The comments left next to the blanks are descriptions from the distributed notebook. The distributed blank notebook is [Lab03-LinearSystem.ipynb](Lab03-LinearSystem.ipynb).

- **Deadline:** one week after distribution
- **Submission:** Submit the notebook online after filling in the blanks and running it.

---

<br>

## Problem 1. Right-Hand Side of the Spring-Mass System (Notebook Section 1.3)

Three masses are connected by four springs between two fixed walls.

```text
Wall ─[k₁]─ m₁ ─[k₂]─ m₂ ─[k₃]─ m₃ ─[k₄]─ Wall
              ↓F₁         ↓F₂         ↓F₃
```

*Figure 1. Three masses connected by four springs between two walls and the force acting on each mass*

Applying force equilibrium at each mass, with positive $`x`$ to the right, gives the following equations.

```math
\begin{aligned}
m_1: &\quad (k_1 + k_2)\,x_1 - k_2\,x_2 = F_1 \\
m_2: &\quad -k_2\,x_1 + (k_2 + k_3)\,x_2 - k_3\,x_3 = F_2 \\
m_3: &\quad -k_3\,x_2 + (k_3 + k_4)\,x_3 = F_3
\end{aligned}
```

In matrix form $`A\mathbf{x} = \mathbf{b}`$, they read as follows.

```math
\begin{bmatrix} k_1+k_2 & -k_2 & 0 \\ -k_2 & k_2+k_3 & -k_3 \\ 0 & -k_3 & k_3+k_4 \end{bmatrix}
\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} =
\begin{bmatrix} F_1 \\ F_2 \\ F_3 \end{bmatrix}
```

The matrix $`A`$ is symmetric and tridiagonal, a structure that becomes important in Section 4 of the notebook. Fill in the blank of the next cell to assemble the right-hand side vector $`\mathbf{b} = [F_1, F_2, F_3]`$.

```python
# Spring constants (N/m)
k1, k2, k3, k4 = 40.0, 30.0, 20.0, 10.0

# Applied forces (N)
F1, F2, F3 = 20.0, 40.0, 30.0

# Coefficient matrix (derived from equilibrium equations above)
A = np.array([
    [k1+k2,   -k2,      0],
    [  -k2, k2+k3,    -k3],
    [    0,   -k3,  k3+k4]
])

# --- COMPLETE: assemble the right-hand side vector b = [F1, F2, F3] ---
b = np.array([ ,  ,  ])
# --- COMPLETE ---

print('Coefficient matrix A:')
print(A)
print('\nRight-hand side b:', b)
```

The cell that follows solves the system with NumPy and verifies the solution through the residual.

```python
# Solve with NumPy
x = np.linalg.solve(A, b)

print('Displacements:')
for i, xi in enumerate(x):
    print(f'  x{i+1} = {xi:.6f} m')

# Verify: compute the residual r = A*x - b (should be ~0)
residual = A @ x - b
print(f'\nResidual A*x - b = {residual}')
print(f'Max |residual|   = {np.max(np.abs(residual)):.2e}  (machine precision: ~2e-16)')
```

## Problem 2. Naive Gaussian Elimination (Notebook Section 2.4)

Gaussian elimination transforms $`A\mathbf{x} = \mathbf{b}`$ into an equivalent upper triangular system through elementary row operations and then solves it by back substitution. First, $`A`$ and $`\mathbf{b}`$ are combined into the augmented matrix $`[A \mid \mathbf{b}]`$.

```math
\left[\begin{array}{ccc|c}
a_{11} & a_{12} & a_{13} & b_1 \\
a_{21} & a_{22} & a_{23} & b_2 \\
a_{31} & a_{32} & a_{33} & b_3
\end{array}\right]
```

**Forward elimination:** For each pivot column $`k = 0, 1, \ldots, n-2`$ and each row below it, $`i = k+1, \ldots, n-1`$, do the following.

1. Compute the elimination factor $`f_{ik} = a_{ik}^{(k)} / a_{kk}^{(k)}`$.
2. Update row $`i`$ as $`\text{row}_i \leftarrow \text{row}_i - f_{ik} \cdot \text{row}_k`$.

This zeros out entry $`(i, k)`$ and leaves an upper triangular matrix.

**Back substitution:** Start from the last equation and work upward.

```math
x_{n-1} = \frac{b_{n-1}^{(n)}}{a_{n-1,n-1}^{(n)}}, \qquad
x_i = \frac{b_i^{(i+1)} - \displaystyle\sum_{j=i+1}^{n-1} a_{ij}^{(i+1)}\, x_j}{a_{ii}^{(i+1)}}, \quad i = n-2, \ldots, 0
```

The notebook traces the process step by step on the following system, whose true solution is $`(x_1, x_2, x_3) = (2, 3, -1)`$.

```math
\begin{bmatrix} 2 & 1 & -1 \\ -3 & -1 & 2 \\ -2 & 1 & 2 \end{bmatrix}
\begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix} =
\begin{bmatrix} 8 \\ -11 \\ -3 \end{bmatrix}
```

- **Step 1 (pivot column $`k = 0`$):** Eliminating below row 0 with $`f_{10} = -3/2`$ and $`f_{20} = -2/2 = -1`$ turns rows 1 and 2 of the augmented matrix into $`[0, 1/2, 1/2 \mid 1]`$ and $`[0, 2, 1 \mid 5]`$.
- **Step 2 (pivot column $`k = 1`$):** Eliminating below row 1 with $`f_{21} = 2/(1/2) = 4`$ turns row 2 into $`[0, 0, -1 \mid 1]`$.
- **Back substitution:** $`x_3 = 1/(-1) = -1`$, $`x_2 = (1 - (1/2)(-1))/(1/2) = 3`$, and $`x_1 = (8 - (1)(3) - (-1)(-1))/2 = 2`$.

Fill in the blanks of the function `gauss_naive` to complete the forward elimination and the back substitution.

```python
def gauss_naive(A, b):
    """
    Solve Ax = b using Naive Gaussian Elimination (no pivoting).

    Parameters
    ----------
    A : (n, n) array-like  — coefficient matrix
    b : (n,) array-like    — right-hand side vector

    Returns
    -------
    x   : (n,) ndarray — solution vector
    Aug : (n, n+1) ndarray — augmented matrix after elimination
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float).reshape(-1, 1)
    n  = A.shape[0]
    nb = n + 1                         # number of columns in augmented matrix

    # Build augmented matrix [A | b]
    Aug = np.hstack((A, b))

    # ── 1. Forward Elimination ────────────────────────────────────────────────
    for k in range(n - 1):             # pivot column k = 0, 1, ..., n-2
        for i in range(k + 1, n):     # rows to eliminate below pivot
            # --- COMPLETE: compute elimination factor and update row i ---
            factor =
            Aug[i, k:nb] =                                  # row i -= factor * row k 
            # --- COMPLETE ---

    # ── 2. Back Substitution ──────────────────────────────────────────────────
    x = np.zeros(n)
    # --- COMPLETE: solve for x[n-1] first, then work backwards ---
    x[n-1] =
    for i in range(n - 2, -1, -1):
        x[i] =
    # --- COMPLETE ---

    return x, Aug
```

The two cells that follow test the function on the $`3 \times 3`$ example above and on the spring-mass system of Problem 1.

```python
# Test gauss_naive on the 3x3 example traced above
A_test = np.array([[ 2,  1, -1],
                   [-3, -1,  2],
                   [-2,  1,  2]], dtype=float)
b_test = np.array([8.0, -11.0, -3.0])

x_sol, Aug_result = gauss_naive(A_test, b_test)

print('Upper triangular augmented matrix after elimination:')
print(np.round(Aug_result, 6))
print()
print(f'Solution x  = {x_sol}')
print(f'Expected    = [2, 3, -1]')
print(f'Residual    = {A_test @ x_sol - b_test}')
```

```python
# Apply gauss_naive to the spring-mass system
x_gauss, _ = gauss_naive(A, b)
x_numpy     = np.linalg.solve(A, b)

print('Spring-mass displacements:')
print(f'  Gauss Elimination : {x_gauss}')
print(f'  np.linalg.solve   : {x_numpy}')
print(f'  Max difference    : {np.max(np.abs(x_gauss - x_numpy)):.2e}')
```

## Problem 3. Partial Pivoting (Notebook Section 3.2)

Naive Gaussian elimination divides by the pivot element $`a_{kk}^{(k)}`$ at each step. If the pivot is zero, the algorithm fails, and if it is very small, roundoff errors are amplified in the same way as in catastrophic cancellation. For example, the following system has the zero leading pivot $`a_{11} = 0`$, so naive elimination divides by zero in its first step.

```math
\begin{bmatrix} 0 & -3 & 7 \\ 1 & 2 & -1 \\ 5 & -2 & 0 \end{bmatrix}\mathbf{x} =
\begin{bmatrix} 4 \\ 0 \\ 3 \end{bmatrix}
```

```python
A_bad = np.array([[0.0, -3.0,  7.0],
                  [1.0,  2.0, -1.0],
                  [5.0, -2.0,  0.0]])
b_bad = np.array([4.0, 0.0, 3.0])

print('True solution (via numpy):')
print('  x =', np.linalg.solve(A_bad, b_bad))

print('\nNaive Gaussian Elimination result:')
x_fail, Aug_fail = gauss_naive(A_bad, b_bad)
print('  x =', x_fail)
print('Augmented matrix after elimination:')
print(Aug_fail)
```

**Partial pivoting** reorders the rows before each elimination step. At pivot step $`k`$, it does the following.

1. Find the row $`r \geq k`$ with the largest absolute value in column $`k`$.
2. Swap rows $`k`$ and $`r`$ in the augmented matrix.
3. Proceed with forward elimination.

This guarantees that the pivot is never zero as long as the system has a unique solution, and the elimination factor satisfies $`|f_{ik}| = |a_{ik}| / |a_{kk}| \leq 1`$, which prevents error amplification. The resulting algorithm is called Gaussian elimination with partial pivoting (GEPP). Fill in the blank of the function `gauss_pivot` so that it swaps rows `k` and `imax`.

```python
def gauss_pivot(A, b):
    """
    Solve Ax = b using Gaussian Elimination with Partial Pivoting.

    Parameters
    ----------
    A : (n, n) array-like  — coefficient matrix
    b : (n,) array-like    — right-hand side vector

    Returns
    -------
    x   : (n,) ndarray — solution vector
    Aug : (n, n+1) ndarray — augmented matrix after elimination
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float).reshape(-1, 1)
    n  = A.shape[0]
    nb = n + 1
    Aug = np.hstack((A, b))

    for k in range(n - 1):
        # ── Partial pivoting ────────────────────────────────────────────────
        # Find the row index (>= k) with the largest |value| in column k
        imax = np.argmax(np.abs(Aug[k:, k])) + k
        if imax != k:
            # --- COMPLETE: swap row k and row imax in the augmented matrix ---
            Aug[[k, imax]] =                        # swap rows
            # --- COMPLETE ---

        # ── Forward elimination (same as naive) ─────────────────────────────
        for i in range(k + 1, n):
            factor = Aug[i, k] / Aug[k, k]
            Aug[i, k:nb] = Aug[i, k:nb] - factor * Aug[k, k:nb]

    # ── Back substitution (same as naive) ───────────────────────────────────
    x = np.zeros(n)
    x[n-1] = Aug[n-1, nb-1] / Aug[n-1, n-1]
    for i in range(n - 2, -1, -1):
        x[i] = (Aug[i, nb-1] - Aug[i, i+1:n] @ x[i+1:n]) / Aug[i, i]

    return x, Aug
```

The cell that follows tests the function on the example that failed above.

```python
# Test gauss_pivot on the previously-failing zero-pivot example
x_pivot, Aug_pivot = gauss_pivot(A_bad, b_bad)
x_ref              = np.linalg.solve(A_bad, b_bad)

print('Gauss with partial pivoting:')
print(f'  x       = {x_pivot}')
print(f'  x (ref) = {x_ref}')
print(f'  Max error: {np.max(np.abs(x_pivot - x_ref)):.2e}')
print()
print('Upper triangular augmented matrix after pivoting:')
print(np.round(Aug_pivot, 6))
```

## Problem 4. The Thomas Algorithm (Notebook Section 4.2)

Many engineering problems, such as 1D heat conduction, beam deflection, and finite-difference PDEs, produce tridiagonal systems in which only the main diagonal and the diagonals directly below and above it are nonzero.

```math
\begin{bmatrix}
f_0 & g_0 &        &        &     \\
e_1 & f_1 & g_1    &        &     \\
    & e_2 & f_2    & g_2    &     \\
    &     & \ddots & \ddots & \ddots \\
    &     &        & e_{n-1} & f_{n-1}
\end{bmatrix}\mathbf{x} = \mathbf{r}
```

Applying Gaussian elimination to this structure gives an algorithm with $`O(n)`$ operations instead of $`O(n^3)`$. The forward sweep of the Thomas algorithm eliminates the subdiagonal (the $`e`$ entries) for $`k = 1, 2, \ldots, n-1`$.

```math
\text{factor} = \frac{e_k}{f_{k-1}}, \qquad
f_k \leftarrow f_k - \text{factor} \cdot g_{k-1}, \qquad
r_k \leftarrow r_k - \text{factor} \cdot r_{k-1}
```

The back substitution for $`k = n-2, n-3, \ldots, 0`$ is as follows.

```math
x_{n-1} = \frac{r_{n-1}}{f_{n-1}}, \qquad x_k = \frac{r_k - g_k\, x_{k+1}}{f_k}
```

Fill in the blanks of the function `thomas` to complete the forward sweep.

```python
def thomas(e, f, g, r):
    """
    Solve a tridiagonal system using the Thomas Algorithm.

    Parameters (all 1-D arrays of length n)
    ----------
    e : subdiagonal    e[0] unused, e[1]..e[n-1] active
    f : main diagonal  f[0]..f[n-1]  (modified internally)
    g : superdiagonal  g[0]..g[n-2] active, g[n-1] unused
    r : right-hand side r[0]..r[n-1]  (modified internally)

    Returns
    -------
    x : (n,) ndarray — solution vector
    """
    f = np.array(f, dtype=float)   # work on copies
    r = np.array(r, dtype=float)
    n = len(f)

    # ── Forward sweep: eliminate subdiagonal ─────────────────────────────────
    for k in range(1, n):
        # --- COMPLETE: compute factor and update f[k] and r[k] ---
        factor =
        f[k] =
        r[k] =
        # --- COMPLETE ---

    # ── Back substitution ────────────────────────────────────────────────────
    x = np.zeros(n)
    x[-1] = r[-1] / f[-1]
    for k in range(n - 2, -1, -1):
        x[k] = (r[k] - g[k] * x[k+1]) / f[k]

    return x
```

Section 4.3 of the notebook applies this function to steady-state 1D heat conduction. A rod of length $`L = 1\,\text{m}`$ has end temperatures $`T(0) = T_L`$ and $`T(L) = T_R`$. Dividing the rod into $`n+1`$ equal segments of width $`\Delta x = L/(n+1)`$ and applying $`d^2T/dx^2 \approx (T_{i-1} - 2T_i + T_{i+1})/\Delta x^2 = 0`$ gives $`T_{i-1} - 2T_i + T_{i+1} = 0`$ $`(i = 1, \ldots, n)`$. This is a tridiagonal system with $`e_i = g_i = 1`$ and $`f_i = -2`$, and the boundary conditions enter through $`r_1 = -T_L`$ and $`r_n = -T_R`$ (all other $`r_i = 0`$).

```python
def heat_conduction_1d(n_interior, T_left=50.0, T_right=150.0):
    """Solve steady-state 1D heat conduction on n_interior interior nodes."""
    n = n_interior
    f =  -2.0 * np.ones(n)   # main diagonal
    e =   1.0 * np.ones(n)   # subdiagonal  (e[0] unused)
    g =   1.0 * np.ones(n)   # superdiagonal (g[n-1] unused)
    r =   np.zeros(n)

    # Boundary conditions modify the first and last RHS entries
    r[0]  -= T_left
    r[-1] -= T_right

    T_interior = thomas(e, f, g, r)

    x_all = np.linspace(0, 1, n + 2)
    T_all = np.concatenate([[T_left], T_interior, [T_right]])
    return x_all, T_all


# Exact solution: linear profile (for zero internal heat generation)
T_left, T_right = 50.0, 150.0
x_exact = np.linspace(0, 1, 500)
T_exact = T_left + (T_right - T_left) * x_exact

plt.figure(figsize=(8, 4))
for n in [3, 6, 15]:
    x, T = heat_conduction_1d(n)
    plt.plot(x, T, 'o-', ms=6, label=f'Thomas, {n} interior nodes')
plt.plot(x_exact, T_exact, 'k--', lw=2, label='Exact (linear)')
plt.xlabel('Position $x$ (m)')
plt.ylabel('Temperature $T$ (°C)')
plt.title('Steady-State 1D Heat Conduction — Thomas Algorithm')
plt.legend()
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()
```

---

<br>

After solving the problems, see the [assignment explanation](../Assignment-Explanation.md) for the answers and solutions and the [completed notebook](../Solution/Lab03-LinearSystem_Completed.ipynb) for the full code.
