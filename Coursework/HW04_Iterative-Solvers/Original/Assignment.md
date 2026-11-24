# HW4 Assignment Handout — Iterative Solvers

[Repository Guide](../../../README.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **Last Updated:** 2026-11-24

This is the Week 5 homework of Numerical Analysis. The task is to fill in the blanks of the lab notebook `Lab05-IterativeSolver.ipynb` distributed by the instructor and run it. The notebook has six parts: matrix splitting and the fixed-point framework, Jacobi iteration, the Gauss-Seidel method, convergence analysis, successive over-relaxation (SOR), and an application to 2D steady-state heat conduction. There are six blanks in total, each between two comment lines marked `COMPLETE`. The first code cell of the notebook imports `numpy` as `np`, `matplotlib.pyplot` as `plt`, and `time`.

Each code block below reproduces a cell of the notebook, and the places where code is missing are the blanks to fill. The comments left next to the blanks are descriptions from the distributed notebook. The distributed blank notebook is [Lab05-IterativeSolver.ipynb](Lab05-IterativeSolver.ipynb).

---

<br>

## Problem 1. Matrix Splitting (Notebook Section 1.3)

Every square matrix can be split into three parts.

$$
A = D + L + U
$$

Here $`D`$ is the diagonal part of $`A`$, $`L`$ the strictly lower triangular part with zeros on and above the diagonal, and $`U`$ the strictly upper triangular part with zeros on and below the diagonal.

```math
A = \begin{bmatrix} a_{11} & a_{12} & a_{13} \\ a_{21} & a_{22} & a_{23} \\ a_{31} & a_{32} & a_{33} \end{bmatrix}
= \underbrace{\begin{bmatrix} a_{11} & & \\ & a_{22} & \\ & & a_{33} \end{bmatrix}}_{D}
+ \underbrace{\begin{bmatrix} 0 & & \\ a_{21} & 0 & \\ a_{31} & a_{32} & 0 \end{bmatrix}}_{L}
+ \underbrace{\begin{bmatrix} 0 & a_{12} & a_{13} \\ & 0 & a_{23} \\ & & 0 \end{bmatrix}}_{U}
```

All classical iterative methods rewrite $`A\mathbf{x} = \mathbf{b}`$ as the fixed-point iteration $`\mathbf{x}^{(k+1)} = T\mathbf{x}^{(k)} + \mathbf{c}`$. $`T`$ is the iteration matrix and $`\mathbf{c}`$ a constant vector, and the iteration converges if and only if the spectral radius satisfies $`\rho(T) = \max_i |\lambda_i(T)| < 1`$. Fill in the blanks of the function `split_DLU` to compute $`L`$ and $`U`$. The latter part of the same cell tests the function on the example matrix of the notebook.

```python
def split_DLU(A):
    """
    Split matrix A into diagonal D, strictly lower triangular L,
    and strictly upper triangular U such that A = D + L + U.

    Returns
    -------
    D, L, U : each (n, n) ndarray
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]

    # --- COMPLETE ---
    D = np.diag(np.diag(A))                          # diagonal matrix (use np.diag twice)
    L =                                              # (strictly lower)
    U =                                              # (strictly upper)
    # --- COMPLETE ---

    return D, L, U


# Test with the example from the lecture notes
A_ex = np.array([[ 3.0, -0.1, -0.2],
                 [ 0.1,  7.0, -0.3],
                 [ 0.3, -0.2, 10.0]])
b_ex = np.array([7.85, -19.3, 71.4])   # true solution: x = [3, -2.5, 7]

D, L, U = split_DLU(A_ex)
print('D =\n', D)
print('L =\n', L)
print('U =\n', U)
print('\nSanity check  D + L + U - A =\n', D + L + U - A_ex)
```

The distributed notebook shows the following output for this cell.

```text
D =
 [[ 3.  0.  0.]
 [ 0.  7.  0.]
 [ 0.  0. 10.]]
L =
 [[ 0.   0.   0. ]
 [ 0.1  0.   0. ]
 [ 0.3 -0.2  0. ]]
U =
 [[ 0.  -0.1 -0.2]
 [ 0.   0.  -0.3]
 [ 0.   0.   0. ]]

Sanity check  D + L + U - A =
 [[0. 0. 0.]
 [0. 0. 0.]
 [0. 0. 0.]]
```

## Problem 2. Jacobi Iteration (Notebook Section 2.2)

Isolating $`x_i`$ from the $`i`$-th equation $`\sum_j a_{ij} x_j = b_i`$ gives the following.

```math
x_i^{(k+1)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j \neq i} a_{ij}\, x_j^{(k)}\right), \quad i = 0, 1, \ldots, n-1
```

All values on the right use the previous iterate $`\mathbf{x}^{(k)}`$, so the entire new vector can be computed simultaneously (easy to parallelize). In matrix form, splitting $`A = D + (L + U)`$ and writing $`D\mathbf{x} = \mathbf{b} - (L+U)\mathbf{x}`$ gives the following.

$$
\mathbf{x}^{(k+1)} = \underbrace{-D^{-1}(L+U)}_{T_J}\mathbf{x}^{(k)} + D^{-1}\mathbf{b}
$$

A common convergence test uses the relative approximate error, and the iteration stops when $`\varepsilon_a`$ falls below a user-specified tolerance $`\varepsilon_s`$.

```math
\varepsilon_a = \frac{\|\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}\|_2}{\|\mathbf{x}^{(k+1)}\|_2} \times 100\%
```

Fill in the blank of the function `jacobi` to compute `x_new` with one Jacobi iteration.

```python
def jacobi(A, b, x0=None, tol=1e-6, max_iter=200):
    """
    Solve Ax = b using the Jacobi iteration.

    Parameters
    ----------
    A        : (n, n) array — coefficient matrix (must be diagonally dominant)
    b        : (n,)   array — right-hand side
    x0       : (n,)   array — initial guess (default: zeros)
    tol      : float        — relative error tolerance (stopping criterion)
    max_iter : int          — maximum number of iterations

    Returns
    -------
    x        : (n,) ndarray — approximate solution
    errors   : list of float — relative approximate error at each iteration
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)

    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)
    errors = []

    d = np.diag(A)                  # diagonal entries a_{ii}
    R = A - np.diag(d)              # off-diagonal part  (L + U)

    for _ in range(max_iter):
        # --- COMPLETE ---
        x_new =                                       
        # --- COMPLETE ---

        # Relative approximate error
        if np.linalg.norm(x_new) > 0:
            ea = np.linalg.norm(x_new - x) / np.linalg.norm(x_new) * 100
        else:
            ea = np.inf
        errors.append(ea)

        x = x_new
        if ea < tol:
            break

    return x, errors
```

The cell that follows tests the function on the notebook's example and prints an iteration table.

```python
# Test Jacobi on the lecture-notes example  (true solution: x = [3, -2.5, 7])
x_true = np.array([3.0, -2.5, 7.0])

x_jac, err_jac = jacobi(A_ex, b_ex, tol=1e-6)
x_ref          = np.linalg.solve(A_ex, b_ex)

print(f'Jacobi converged in {len(err_jac)} iterations')
print(f'Solution : {np.round(x_jac, 8)}')
print(f'Reference: {np.round(x_ref, 8)}')
print(f'Max error: {np.max(np.abs(x_jac - x_ref)):.2e}')

# Print iteration table
print(f'\n{"Iter":>5}  {"x1":>12}  {"x2":>12}  {"x3":>12}  {"ea (%)": >12}')
print('-' * 60)
x_k = np.zeros(3)
d_A = np.diag(A_ex)
R_A = A_ex - np.diag(d_A)
for k in range(min(8, len(err_jac))):
    x_k = (b_ex - R_A @ x_k) / d_A
    print(f'{k+1:>5}  {x_k[0]:>12.6f}  {x_k[1]:>12.6f}  {x_k[2]:>12.6f}  {err_jac[k]:>12.6f}')
```

## Problem 3. The Gauss-Seidel Method (Notebook Section 3.1)

Jacobi iteration computes all components of $`\mathbf{x}^{(k+1)}`$ from old values only. The Gauss-Seidel method instead uses each newly computed component immediately.

```math
x_i^{(k+1)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j=0}^{i-1} a_{ij}\, x_j^{(k+1)} - \sum_{j=i+1}^{n-1} a_{ij}\, x_j^{(k)}\right)
```

The first sum uses values already updated in this iteration, and the second sum uses old values. In matrix form, splitting $`A = (D+L) + U`$ gives the following.

$$
\mathbf{x}^{(k+1)} = \underbrace{-(D+L)^{-1}U}_{T_{GS}}\mathbf{x}^{(k)} + (D+L)^{-1}\mathbf{b}
$$

| Feature | Jacobi iteration | Gauss-Seidel method |
|:-----|:-----|:-----|
| Values used | All old $`\mathbf{x}^{(k)}`$ | New for $`j < i`$, old for $`j > i`$ |
| Convergence rate | Slower | Typically about 2 times faster |
| Parallelism | Easy (components are independent) | Sequential (each depends on the previous ones) |
| Memory | Needs two copies of $`\mathbf{x}`$ | Can update in place |

*Table 1. Comparison of Jacobi iteration and the Gauss-Seidel method*

The notebook shows the first iteration of the example from $`\mathbf{x}^{(0)} = \mathbf{0}`$ as follows. The values in bold are those updated in this iteration.

```math
\begin{aligned}
x_1^{(1)} &= \frac{7.85 - (-0.1)(0) - (-0.2)(0)}{3} = 2.6167 \\
x_2^{(1)} &= \frac{-19.3 - (0.1)(\mathbf{2.6167}) - (-0.3)(0)}{7} = -2.7945 \\
x_3^{(1)} &= \frac{71.4 - (0.3)(\mathbf{2.6167}) - (-0.2)(\mathbf{-2.7945})}{10} = 7.0056
\end{aligned}
```

Fill in the blanks of the function `gauss_seidel` so that each component is updated with the most recently available values.

```python
def gauss_seidel(A, b, x0=None, tol=1e-6, max_iter=200):
    """
    Solve Ax = b using the Gauss-Seidel iteration.

    Parameters
    ----------
    A, b, x0, tol, max_iter  — same as jacobi()

    Returns
    -------
    x      : (n,) ndarray — approximate solution
    errors : list of float — relative approximate error per iteration
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)
    errors = []

    for _ in range(max_iter):
        x_old = x.copy()

        # Update each component using the most recently available values
        for i in range(n):
            # --- COMPLETE ---
            sigma =                                  
            x[i]  =                                  
            # --- COMPLETE ---

        # Relative approximate error
        if np.linalg.norm(x) > 0:
            ea = np.linalg.norm(x - x_old) / np.linalg.norm(x) * 100
        else:
            ea = np.inf
        errors.append(ea)

        if ea < tol:
            break

    return x, errors
```

The cell that follows compares the iteration counts of the two methods on the same example and plots their convergence histories.

```python
# Compare Jacobi vs Gauss-Seidel on the lecture-notes example
x_gs, err_gs = gauss_seidel(A_ex, b_ex, tol=1e-6)

print(f'Jacobi       converged in {len(err_jac):3d} iterations')
print(f'Gauss-Seidel converged in {len(err_gs):3d} iterations')
print(f'GS speedup (iterations): {len(err_jac)/len(err_gs):.1f}x')
print(f'\nGS solution: {np.round(x_gs, 8)}')
print(f'True value : {x_true}')

plt.figure(figsize=(8, 4))
plt.semilogy(range(1, len(err_jac)+1), err_jac, 'b.-', ms=8, label='Jacobi')
plt.semilogy(range(1, len(err_gs)+1),  err_gs,  'r.-', ms=8, label='Gauss-Seidel')
plt.axhline(1e-6, color='gray', ls='--', label='Tolerance $10^{-6}$')
plt.xlabel('Iteration $k$')
plt.ylabel('Relative error $\\varepsilon_a$ (%)')
plt.title('Jacobi vs Gauss-Seidel Convergence')
plt.legend()
plt.grid(True, which='both', ls=':', alpha=0.5)
plt.tight_layout()
plt.show()
```

## Problem 4. Diagonal Dominance Check (Notebook Section 4.1)

A matrix $`A`$ is called strictly diagonally dominant if the following holds for every row $`i`$.

$$
|a_{ii}| > \sum_{j \neq i} |a_{ij}|
$$

If $`A`$ is strictly diagonally dominant, both Jacobi iteration and the Gauss-Seidel method converge for any initial guess $`\mathbf{x}^{(0)}`$. This is only a sufficient condition, so iterative methods may still converge for matrices that are not strictly diagonally dominant. The notebook checks the example matrix of Problem 1 row by row as follows.

- **Row 0:** $`|3| > |-0.1| + |-0.2| = 0.3`$
- **Row 1:** $`|7| > |0.1| + |-0.3| = 0.4`$
- **Row 2:** $`|10| > |0.3| + |-0.2| = 0.5`$

Fill in the blanks of the function `is_diag_dominant` to compute, for each row, the sum of the absolute values of the off-diagonal entries `offdiag_sum` and the margin `margins`. A row with a positive margin is dominant. The latter part of the same cell tests the function on three matrices.

```python
def is_diag_dominant(A, strict=True):
    """
    Check whether A is (strictly) diagonally dominant.

    Returns
    -------
    dominant : bool
    margins  : (n,) array — |a_ii| - sum_{j!=i}|a_ij| for each row
                             Positive margin means the row is dominant.
    """
    A = np.array(A, dtype=float)

    # --- COMPLETE ---
    diag_abs = np.abs(np.diag(A))                    # |a_ii|
    offdiag_sum =                                    # row sums of |A| minus the diagonal
    margins =
    # --- COMPLETE ---

    if strict:
        dominant = bool(np.all(margins > 0))
    else:
        dominant = bool(np.all(margins >= 0))
    return dominant, margins


# Test on well-conditioned and poorly-conditioned matrices
A_good = np.array([[ 4.0, -1.0,  0.5],
                   [ 0.2,  5.0, -1.0],
                   [-0.5,  0.3,  6.0]])

A_bad  = np.array([[ 1.0,  2.0, -3.0],
                   [ 1.5,  2.0,  0.5],
                   [-1.0,  1.0,  1.0]])

dom_good, mar_good = is_diag_dominant(A_good)
dom_bad,  mar_bad  = is_diag_dominant(A_bad)

print('A_good — strictly diagonally dominant:', dom_good)
print('  Row margins:', np.round(mar_good, 4))

print('\nA_bad  — strictly diagonally dominant:', dom_bad)
print('  Row margins:', np.round(mar_bad, 4))

# Check our lecture-notes example
dom_ex, mar_ex = is_diag_dominant(A_ex)
print(f'\nA_ex   — strictly diagonally dominant: {dom_ex},  margins: {np.round(mar_ex, 4)}')
```

The distributed notebook shows the following output for this cell.

```text
A_good — strictly diagonally dominant: True
  Row margins: [2.5 3.8 5.2]

A_bad  — strictly diagonally dominant: False
  Row margins: [-4.  0. -1.]

A_ex   — strictly diagonally dominant: True,  margins: [2.7 6.6 9.5]
```

## Problem 5. Successive Over-Relaxation (Notebook Section 5.1)

The Gauss-Seidel method takes a step from $`x_i^{(k)}`$ to a new value $`\tilde{x}_i^{(k+1)}`$. Successive over-relaxation (SOR) overshoots this update by a weighting factor $`\omega > 1`$ to accelerate convergence.

```math
x_i^{(k+1)} = (1 - \omega)\,x_i^{(k)}
+ \frac{\omega}{a_{ii}}\left(b_i - \sum_{j < i} a_{ij}\, x_j^{(k+1)} - \sum_{j > i} a_{ij}\, x_j^{(k)}\right)
```

| $`\omega`$ | Name | Effect |
|:-----|:-----|:-----|
| $`\omega = 1`$ | Gauss-Seidel method | No acceleration |
| $`1 < \omega < 2`$ | Over-relaxation (SOR) | Faster convergence |
| $`0 < \omega < 1`$ | Under-relaxation | Slower; useful to stabilize divergent cases |
| $`\omega \geq 2`$ | - | Diverges |

*Table 2. Kinds of relaxation by the range of the weighting factor*

The optimal $`\omega`$ minimizes $`\rho(T_{SOR}(\omega))`$. For the tridiagonal model problem (the 5-point Laplace stencil on an $`n \times n`$ grid), it is as follows.

$$
\omega_{\text{opt}} = \frac{2}{1 + \sqrt{1 - \rho(T_J)^2}}
$$

Fill in the blank of the function `sor` to complete the SOR update that blends the Gauss-Seidel value `x_gs_i` with the old value.

```python
def sor(A, b, omega=1.0, x0=None, tol=1e-6, max_iter=500):
    """
    Solve Ax = b using Successive Over-Relaxation (SOR).

    Parameters
    ----------
    A, b, x0, tol, max_iter — same as gauss_seidel()
    omega  : float — relaxation parameter (0 < omega < 2)

    Returns
    -------
    x      : (n,) ndarray
    errors : list of float
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)
    errors = []

    for _ in range(max_iter):
        x_old = x.copy()

        for i in range(n):
            # Gauss-Seidel update (x_gs_i)
            sigma  = A[i, :i] @ x[:i] + A[i, i+1:] @ x_old[i+1:]
            x_gs_i = (b[i] - sigma) / A[i, i]

            # --- COMPLETE: apply SOR blending ---
            x[i] =
            # --- COMPLETE ---

        if np.linalg.norm(x) > 0:
            ea = np.linalg.norm(x - x_old) / np.linalg.norm(x) * 100
        else:
            ea = np.inf
        errors.append(ea)

        if ea < tol:
            break

    return x, errors
```

The cell that follows compares Jacobi iteration, the Gauss-Seidel method, and SOR with several values of $`\omega`$ on the heat network of Section 3.2 of the notebook. Section 3.2 sets up the system for the three interior nodes, $`K = \begin{bmatrix} 4 & -1 & 0 \\ -1 & 4 & -1 \\ 0 & -1 & 4 \end{bmatrix}`$ and $`\mathbf{f} = [200, 0, 50]^T`$, as `K_heat` and `f_heat`, and stores the error history of the Gauss-Seidel solve to the tolerance $`10^{-8}`$ in `err_heat_gs`.

```python
# Compare Jacobi, GS, and SOR (ω = 1.1, 1.3, 1.5) on the heat network
results = {}
x_jac_h, err_jac_h = jacobi(K_heat, f_heat, tol=1e-8)
results['Jacobi'] = err_jac_h
results['GS (ω=1.0)'] = err_heat_gs

for w in [1.1, 1.3, 1.5, 1.7]:
    _, err = sor(K_heat, f_heat, omega=w, tol=1e-8)
    results[f'SOR ω={w}'] = err

plt.figure(figsize=(9, 5))
styles = ['b--', 'r-', 'g-', 'm-', 'c-', 'orange']
for (label, errs), style in zip(results.items(), styles):
    plt.semilogy(range(1, len(errs)+1), errs, style, lw=1.8, label=f'{label} ({len(errs)} iters)')
plt.axhline(1e-8, color='gray', ls=':', label='Tolerance')
plt.xlabel('Iteration $k$')
plt.ylabel('Relative error $\\varepsilon_a$ (%)')
plt.title('Jacobi / Gauss-Seidel / SOR Comparison (Heat Network)')
plt.legend(fontsize=8)
plt.grid(True, which='both', ls=':', alpha=0.5)
plt.tight_layout()
plt.show()

print(f'\n{"Method":<14}  {"Iterations":>12}')
print('-' * 28)
for label, errs in results.items():
    print(f'{label:<14}  {len(errs):>12}')
```

## Problem 6. Grid SOR for 2D Heat Conduction (Notebook Section 6.2)

The boundary conditions of the square plate $`[0,1]\times[0,1]`$ are $`T = 0`$°C at the bottom ($`y=0`$), $`T = 100`$°C at the top ($`y=1`$), and $`T = 0`$°C on the left ($`x=0`$) and right ($`x=1`$). The steady-state temperature satisfies the 2D Laplace equation.

$$
\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0
$$

Discretizing with finite differences on an $`(n+2)\times(n+2)`$ grid with spacing $`h = 1/(n+1)`$ gives the following equation.

$$
T_{i-1,j} + T_{i+1,j} + T_{i,j-1} + T_{i,j+1} - 4T_{i,j} = 0
$$

This 5-point stencil produces a large, sparse, symmetric system with $`n^2`$ unknowns and only $`5n^2`$ nonzeros, which suits iterative solvers. Instead of assembling the $`n^2 \times n^2`$ matrix, the notebook sweeps directly over the grid with the stencil. Fill in the blank of the function `solve_laplace_2d` to complete the update for `method='sor'`.

```python
def solve_laplace_2d(n, T_top=100.0, T_bot=0.0, T_left=0.0, T_right=0.0,
                     method='gs', omega=1.0, tol=1e-6, max_iter=5000):
    """
    Solve the 2D Laplace equation on an (n+2) x (n+2) grid using
    a point iterative method (Gauss-Seidel or SOR).

    Interior nodes: indices 1..n, 1..n
    Returns T : (n+2, n+2) temperature array, iters : iteration count
    """
    T = np.zeros((n + 2, n + 2))

    # Boundary conditions
    T[0,  :]  = T_bot
    T[-1, :]  = T_top
    T[:, 0]   = T_left
    T[:, -1]  = T_right

    iters = 0
    for it in range(max_iter):
        T_old = T.copy()

        # Sweep interior nodes
        for i in range(1, n + 1):          # row index
            for j in range(1, n + 1):      # column index
                # Gauss-Seidel update for the 5-point stencil
                T_gs = 0.25 * (T[i-1, j] + T[i+1, j] + T[i, j-1] + T[i, j+1])

                if method == 'sor':
                    # --- COMPLETE: apply SOR blending ---
                    T[i, j] =
                    # --- COMPLETE ---
                else:
                    T[i, j] = T_gs

        iters += 1

        # Convergence check on interior nodes
        diff = np.max(np.abs(T[1:-1, 1:-1] - T_old[1:-1, 1:-1]))
        if diff < tol:
            break

    return T, iters
```

The cell that follows compares the Gauss-Seidel method with SOR at $`\omega = 1.7`$ for $`n = 30`$ and plots the temperature distributions.

```python
# Solve and visualize for n = 30
n_grid = 30

T_gs_2d,  it_gs  = solve_laplace_2d(n_grid, method='gs',  tol=1e-6)
T_sor_2d, it_sor = solve_laplace_2d(n_grid, method='sor', omega=1.7, tol=1e-6)

print(f'2D Laplace (n={n_grid}, {n_grid**2} unknowns):')
print(f'  Gauss-Seidel : {it_gs:4d} iterations')
print(f'  SOR (ω=1.7)  : {it_sor:4d} iterations')
print(f'  SOR speedup  : {it_gs / it_sor:.1f}x')
print(f'  Max |GS - SOR|: {np.max(np.abs(T_gs_2d - T_sor_2d)):.2e}')

# Plot
x = y = np.linspace(0, 1, n_grid + 2)
X, Y = np.meshgrid(x, y)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (T_plot, title) in zip(axes,
    [(T_gs_2d,  f'Gauss-Seidel ({it_gs} iters)'),
     (T_sor_2d, f'SOR ω=1.7 ({it_sor} iters)')]):
    cf = ax.contourf(X, Y, T_plot, 20, cmap='hot')
    fig.colorbar(cf, ax=ax, label='T (°C)')
    ax.contour(X, Y, T_plot, 10, colors='k', linewidths=0.5, alpha=0.4)
    ax.set(title=title, xlabel='x', ylabel='y')

plt.suptitle('2D Steady-State Heat Conduction: $\\nabla^2 T = 0$', fontsize=13)
plt.tight_layout()
plt.show()
```

---

<br>

After solving the problems, see the [assignment explanation](../Assignment-Explanation.md) for the answers and solutions and the [completed notebook](../Solution/Lab05-IterativeSolver_Completed.ipynb) for the full code.
