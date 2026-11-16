# HW7 Assignment Handout — Optimization Methods

[Repository Guide](../../../README.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **Last Updated:** 2026-11-16

This is the Week 9 homework of Numerical Analysis. The task is to fill in the blanks of the lab notebook `Lab08-Optimization.ipynb` distributed by the instructor and run it. The notebook has six parts: local and global minima, golden-section search, parabolic interpolation, steepest descent, `scipy.optimize` with engineering applications, and constrained optimization (optional). There are three blanks in total, each between two comment lines marked `COMPLETE`. The Exercise cell of Section 5.4 that extends the cylinder design problem was distributed with its code filled in, and the two Exercise cells distributed with empty code are reproduced in [Additional Exercises in the Notebook](#additional-exercises-in-the-notebook) at the end of this handout. The first code cell of the notebook imports `numpy` as `np`, `matplotlib.pyplot` as `plt`, and `minimize_scalar` and `minimize` from `scipy.optimize`.

Each code block below reproduces a cell of the notebook, and the places where code is missing are the blanks to fill. The distributed blank notebook is [Lab08-Optimization.ipynb](Lab08-Optimization.ipynb).

- **Deadline:** one week after distribution
- **Submission:** Submit the notebook online after filling in the blanks and running it.

---

<br>

## Example Function of the Notebook (Notebook Section 1)

**Optimization** seeks the values of decision variables that minimize (or maximize) an objective function. For smooth $`f`$, a **necessary condition** for a local minimum at $`x^*`$ is $`f'(x^*) = 0`$ in one dimension and $`\nabla f(x^*) = \mathbf{0}`$ in several dimensions. A **sufficient condition** for a local **minimum** in one dimension is $`f'(x^*) = 0`$ together with $`f''(x^*) > 0`$. Most numerical methods find **local** minima, and finding the **global** minimum is much harder.

The notebook uses the following function throughout.

$$
f(x) = \frac{x^2}{10} - 2\sin(x)
$$

```python
# Running example: f(x) = x²/10 - 2sin(x)
f1d = lambda x: x**2 / 10.0 - 2.0 * np.sin(x)
```

The notebook marks the global minimum in $`[0, 4]`$ and the local minimum in $`[5, 8]`$ on a plot of this function and explains that numerical methods started in $`[0, 4]`$ find the global minimum, while methods started in $`[5, 8]`$ find only a local minimum.

## Problem 1. Golden-Section Search (Notebook Sections 2.1-2.3)

Analogous to bisection for root-finding, golden-section search shrinks a bracket $`[x_l, x_u]`$ known to contain a minimum of a **unimodal** function. At each step, two interior points $`x_1 > x_2`$ are placed symmetrically.

$$
d = (\varphi - 1)(x_u - x_l), \qquad x_1 = x_l + d, \quad x_2 = x_u - d
$$

Here $`\varphi = \frac{1+\sqrt{5}}{2} \approx 1.618`$ is the **golden ratio**. For minimization, $`f(x_1)`$ and $`f(x_2)`$ are compared and the bracket is updated as follows.

| Condition | Keep interval | Rationale |
|:-----|:-----|:-----|
| $`f(x_1) < f(x_2)`$ | $`[x_2, x_u]`$ | The minimum is to the right of $`x_2`$ |
| $`f(x_1) \geq f(x_2)`$ | $`[x_l, x_1]`$ | The minimum is to the left of $`x_1`$ |

*Table 1. Bracket update rule of golden-section search*

After the update, one of the old interior points becomes a new interior point, so only **one new function evaluation** is needed per iteration. Each iteration reduces the interval by the factor $`(\varphi - 1) \approx 0.618`$, which is linear convergence, and the approximate relative error is $`\varepsilon_a = (2 - \varphi) \left|\dfrac{\Delta x}{x_{\rm opt}}\right|`$ with $`\Delta x = x_u - x_l`$. Fill in the two blanks of the following function.

```python
def goldmin(f, xl, xu, es=1e-7, maxit=100):
    """
    Find the minimum of a unimodal f on [xl, xu] using Golden-Section Search.

    Parameters
    ----------
    f       : callable — unimodal objective function
    xl, xu  : float   — initial bracket
    es      : float   — stopping criterion (approx. relative error)
    maxit   : int     — maximum iterations

    Returns
    -------
    xopt    : float  — minimizer estimate
    fopt    : float  — f(xopt)
    ea      : float  — final approximate relative error
    n_iter  : int    — iterations performed
    history : list of (iter, xopt, f(xopt), ea, interval_width)
    """
    phi = (1.0 + np.sqrt(5)) / 2.0

    d  = (phi - 1) * (xu - xl)
    x1 = xl + d;  f1 = f(x1)
    x2 = xu - d;  f2 = f(x2)

    history = []
    xopt = x1 if f1 < f2 else x2

    for i in range(maxit):
        xint = xu - xl   # current interval width

        if f1 < f2:
            # minimum is to the right of x2: new bracket = [x2, xu]
            xopt = x1
            # --- COMPLETE: update bracket and interior points ---
            xl = x2          # shift left bound to old x2
            x2 = x1          # old x1 becomes new left interior point
            f2 = f1
            x1 =              #   — new right interior point
            f1 = f(x1)
            # --- COMPLETE ---
        else:
            # minimum is to the left of x1: new bracket = [xl, x1]
            xopt = x2
            # --- COMPLETE: update bracket and interior points ---
            xu = x1          # shift right bound to old x1
            x1 = x2          # old x2 becomes new right interior point
            f1 = f2
            x2 =              #  — new left interior point
            f2 = f(x2)
            # --- COMPLETE ---

        ea = (2.0 - phi) * abs(xint / xopt) if xopt != 0 else (2.0 - phi) * abs(xint)
        history.append((i + 1, xopt, f(xopt), ea, xu - xl))

        if ea < es:
            break

    return xopt, f(xopt), ea, i + 1, history
```

The next cell applies the completed function to the example function.

```python
# Test goldmin on the running example
x_true = 1.4275517797083348   # known global minimum (from scipy)

xopt, fopt, ea, n_iter, history = goldmin(f1d, xl=0.0, xu=4.0, es=1e-9, maxit=100)

print(f'Minimum found : x* = {xopt:.10f}')
print(f'f(x*)         : {fopt:.10f}')
print(f'True minimum  : {x_true:.10f}')
print(f'|error|       : {abs(xopt - x_true):.3e}')
print(f'Iterations    : {n_iter}')

# Iteration table
print(f'\n{"Iter":>5}  {"xopt":>14}  {"f(xopt)":>14}  {"ea":>12}  {"interval":>10}')
print('-' * 62)
for row in history[::3]:   # print every 3rd row to keep table short
    print(f'{row[0]:>5}  {row[1]:>14.8f}  {row[2]:>14.8f}  {row[3]:>12.6e}  {row[4]:>10.6f}')
```

The distributed notebook shows the following output for this cell.

```text
Minimum found : x* = 1.4275517674
f(x*)         : -1.7757256531
True minimum  : 1.4275517797
|error|       : 1.229e-08
Iterations    : 45

 Iter            xopt         f(xopt)            ea    interval
--------------------------------------------------------------
    1      1.52786405     -1.76472025  1.000000e+00    2.472136
    4      1.52786405     -1.76472025  2.360680e-01    0.583592
    7      1.44271910     -1.77547480  5.901699e-02    0.137767
   10      1.42261911     -1.77569914  1.412887e-02    0.032522
   13      1.42736407     -1.77572561  3.324285e-03    0.007678
   16      1.42736407     -1.77572561  7.847573e-04    0.001812
   19      1.42752750     -1.77572565  1.852349e-04    0.000428
   22      1.42756608     -1.77572565  4.372684e-05    0.000101
   25      1.42755134     -1.77572565  1.032261e-05    0.000024
   28      1.42755134     -1.77572565  2.436838e-06    0.000006
   31      1.42755216     -1.77572565  5.752591e-07    0.000001
   34      1.42755185     -1.77572565  1.358003e-07    0.000000
   37      1.42755177     -1.77572565  3.205810e-08    0.000000
   40      1.42755177     -1.77572565  7.567891e-09    0.000000
   43      1.42755177     -1.77572565  1.786537e-09    0.000000
```

## Problem 2. Parabolic Interpolation (Notebook Sections 3.1-3.3)

Instead of using only the bracket endpoints, parabolic interpolation fits a **quadratic polynomial** through three points $`(x_1, f_1)`$, $`(x_2, f_2)`$, $`(x_3, f_3)`$ and jumps directly to the vertex of the parabola.

```math
x_4 = x_2 - \frac{1}{2}\cdot
\frac{(x_2 - x_1)^2[f_2 - f_3] - (x_2 - x_3)^2[f_2 - f_1]}
{(x_2 - x_1)[f_2 - f_3] - (x_2 - x_3)[f_2 - f_1]}
```

If $`f`$ is well approximated by a parabola near the minimum, this can be much faster than golden-section search. After computing $`x_4`$, the worst of the three current points is replaced and the process repeats, terminating when $`|x_4^{\rm new} - x_4^{\rm old}|`$ is small enough. Near a smooth minimum the convergence is superlinear with order about 1.32, but the method can fail if the three points are nearly collinear (denominator $`\approx 0`$). Fill in the blank of the following function.

```python
def parabolic(f, x1, x2, x3, es=1e-7, maxit=50):
    """
    Minimize f using Parabolic Interpolation.

    Requires three initial points x1 < x2 < x3 with f(x2) < f(x1) and f(x2) < f(x3).

    Parameters
    ----------
    f       : callable
    x1, x2, x3 : float — initial three points  (x1 < x2 < x3)
    es      : float   — stopping criterion
    maxit   : int

    Returns
    -------
    x4, f(x4), ea, n_iter, history
    """
    history = []
    x4_old = x2

    for i in range(maxit):
        y1, y2, y3 = f(x1), f(x2), f(x3)

        num = (x2 - x1)**2 * (y2 - y3) - (x2 - x3)**2 * (y2 - y1)
        den = (x2 - x1) * (y2 - y3)    - (x2 - x3) * (y2 - y1)

        if abs(den) < 1e-30:
            break   # denominator near zero — cannot proceed

        # --- COMPLETE: parabolic interpolation formula ---
        x4 = 
        # --- COMPLETE ---

        ea = abs((x4 - x4_old) / x4) if x4 != 0 else abs(x4 - x4_old)
        history.append((i + 1, x4, f(x4), ea))

        if ea < es:
            break

        # Update: replace the worst of the three points
        pts = sorted([(f(p), p) for p in [x1, x2, x3, x4]], key=lambda t: t[0])
        x1, x2, x3 = sorted([pts[0][1], pts[1][1], pts[2][1]])

        x4_old = x4

    return x4, f(x4), ea, i + 1, history
```

The next cell compares golden-section search and parabolic interpolation. The plotting part is omitted here.

```python
# Compare Golden-Section vs Parabolic Interpolation
xopt_gs, fopt_gs, _, n_gs, hist_gs = goldmin(f1d, 0.0, 4.0, es=1e-10, maxit=100)
xopt_pi, fopt_pi, _, n_pi, hist_pi = parabolic(f1d, 0.0, 1.0, 4.0, es=1e-10, maxit=50)

print(f'Golden-Section:         {n_gs:3d} iterations,  final |error| = {abs(xopt_gs - x_true):.3e}')
print(f'Parabolic Interpolation:{n_pi:3d} iterations,  final |error| = {abs(xopt_pi - x_true):.3e}')
print('\nNote: Parabolic converges faster when f is smooth near the minimum.')
print('\nUse: curvature information of a function')
```

The distributed notebook shows the following output for this cell.

```text
Golden-Section:          49 iterations,  final |error| = 1.229e-08
Parabolic Interpolation:  9 iterations,  final |error| = 2.222e-09

Note: Parabolic converges faster when f is smooth near the minimum.

Use: curvature information of a function
```

## Problem 3. Steepest Descent (Notebook Sections 4.1-4.4)

For a function $`f(\mathbf{x})`$ with $`\mathbf{x} \in \mathbb{R}^n`$, the **gradient** $`\nabla f`$ points in the direction of steepest ascent. Moving in the **negative gradient direction** therefore decreases $`f`$ most rapidly. Starting from $`\mathbf{x}^{(0)}`$, the algorithm repeats the following for $`k = 0, 1, 2, \ldots`$.

1. Compute the gradient $`\mathbf{g}^{(k)} = \nabla f(\mathbf{x}^{(k)})`$.
2. Choose the step size $`\alpha^{(k)}`$ (line search).
3. Update with $`\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} - \alpha^{(k)} \mathbf{g}^{(k)}`$.
4. Stop if $`\|\nabla f(\mathbf{x}^{(k)})\| < \varepsilon`$ or $`\|\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}\| < \delta`$.

A step size that is too large makes $`f`$ increase (overshooting the minimum), and one that is too small makes convergence slow. The backtracking (Armijo) rule starts with $`\alpha = 1`$ and sets $`\alpha \leftarrow c\,\alpha`$ $`(c < 1)`$ while $`f(\mathbf{x} - \alpha\mathbf{g}) \geq f(\mathbf{x})`$. Steepest descent converges linearly, and for strongly convex $`f`$, $`\|\mathbf{x}^{(k)} - \mathbf{x}^*\| \leq \left(\frac{\kappa - 1}{\kappa + 1}\right)^k \|\mathbf{x}^{(0)} - \mathbf{x}^*\|`$. Here $`\kappa = \lambda_{\max}/\lambda_{\min}`$ is the **condition number** of the Hessian, and a high condition number leads to slow, zigzagging convergence.

The notebook uses the two-dimensional example $`f(x, y) = 2 + x - y + 2x^2 + 2xy + y^2`$, whose analytical minimum is $`(-1, 1.5)`$. The next cell defines the function and its gradient. The plotting part is omitted here.

```python
# Running 2D example: f(x, y) = 2 + x - y + 2x² + 2xy + y²
# Analytical minimum: ∂f/∂x = 1+4x+2y=0, ∂f/∂y = -1+2x+2y=0 → x=-1, y=3/2

def f2d(x):
    return 2.0 + x[0] - x[1] + 2*x[0]**2 + 2*x[0]*x[1] + x[1]**2

def grad_f2d(x):
    return np.array([1.0 + 4*x[0] + 2*x[1],
                    -1.0 + 2*x[0] + 2*x[1]])

x2_true = np.array([-1.0, 1.5])   # exact minimum
```

Fill in the update formula of steepest descent in the blank of the next cell.

```python
def backtrack_linesearch(f, x, grad, c=0.5, alpha0=1.0, max_ls=50):
    """
    Backtracking line search: find alpha so that f(x - alpha*grad) < f(x).

    Starts from alpha0 and multiplies by c until descent condition is satisfied.
    """
    alpha = alpha0
    fx = f(x)
    for _ in range(max_ls):
        if f(x - alpha * grad) < fx:
            break
        alpha *= c
    return alpha


def steepest_descent(f, grad_f, x0, tol=1e-6, maxit=1000, record=True):
    """
    Minimize f using steepest descent with backtracking line search.

    Parameters
    ----------
    f      : callable — objective function, accepts 1-D array
    grad_f : callable — gradient function, returns 1-D array
    x0     : array-like — initial guess
    tol    : float    — stopping criterion: ||grad f|| < tol
    maxit  : int      — maximum iterations
    record : bool     — record iteration history

    Returns
    -------
    x       : array  — minimizer estimate
    f(x)    : float
    n_iter  : int
    history : list of (iter, x.copy(), f(x), ||grad||, alpha)
    """
    x = np.array(x0, dtype=float)
    history = []

    for i in range(maxit):
        g = grad_f(x)
        grad_norm = np.linalg.norm(g)

        if record:
            history.append((i, x.copy(), f(x), grad_norm, None))

        if grad_norm < tol:
            break

        # --- COMPLETE: gradient descent update ---
        alpha = backtrack_linesearch(f, x, g)   # backtracking line search
        x =                                    # — steepest descent step
        # --- COMPLETE ---

        if record and history:
            history[-1] = (history[-1][0], history[-1][1],
                           history[-1][2], history[-1][3], alpha)

    return x, f(x), i + 1, history
```

The next cell applies the completed function to the two-dimensional example. The plotting part is omitted here.

```python
# Test steepest descent on the 2D example
x0 = np.array([-2.0, 3.0])
xopt2d, fopt2d, n_sd, hist_sd = steepest_descent(f2d, grad_f2d, x0, tol=1e-8, maxit=200)

print(f'Starting point    : {x0}')
print(f'Minimizer found   : {np.round(xopt2d, 8)}')
print(f'True minimum      : {x2_true}')
print(f'f(minimizer)      : {fopt2d:.10f}')
print(f'f(true min)       : {f2d(x2_true):.10f}')
print(f'|error|           : {np.linalg.norm(xopt2d - x2_true):.3e}')
print(f'Iterations        : {n_sd}')
```

The distributed notebook shows the following output for this cell. After the output come plots of the descent path and the gradient norm.

```text
Starting point    : [-2.  3.]
Minimizer found   : [-1.00000001  1.50000001]
True minimum      : [-1.   1.5]
f(minimizer)      : 0.7500000000
f(true min)       : 0.7500000000
|error|           : 2.107e-08
Iterations        : 200
```

## Additional Exercises in the Notebook

Besides the blanks, the notebook contains two Exercise cells whose code cells are empty.

**Exercise A (end of notebook Section 2).** The following exercise comes after golden-section search is completed. Apply `goldmin` to find the **local** minimum of $`f(x) = x^2/10 - 2\sin(x)`$ in the interval $`[5, 8]`$. How does its value compare to the global minimum found above?

```python
# --- YOUR CODE HERE ---
# Apply goldmin to find the local minimum in [5, 8]
```

**Exercise B (end of notebook Section 4).** The following exercise comes after steepest descent is tested. Try steepest descent with a **fixed** step size (for example, $`\alpha = 0.01, 0.1, 0.5`$) instead of backtracking, and observe how the choice of $`\alpha`$ affects convergence. Hint: replace `backtrack_linesearch(...)` with a constant `alpha`.

```python
# --- YOUR CODE HERE ---
# Try fixed step sizes alpha = 0.01, 0.1, 0.5 on f2d
# Plot convergence curves for each
```

---

<br>

After solving the problems, see the [assignment explanation](../Assignment-Explanation.md) for the answers and solutions and the [completed notebook](../Solution/Lab08-Optimization_Completed.ipynb) for the full code.
