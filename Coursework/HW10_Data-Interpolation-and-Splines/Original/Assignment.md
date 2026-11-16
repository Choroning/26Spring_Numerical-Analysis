# HW10 Assignment Handout — Data Interpolation and Splines

[Repository Guide](../../../README.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **Last Updated:** 2026-11-16

This is the Week 12 homework of Numerical Analysis. The task is to fill in the blank of the lab notebook `Lab11-Spline.ipynb` distributed by the instructor and run it. The notebook has seven parts: why interpolation is needed, Newton's divided differences, Lagrange interpolation, Runge's phenomenon, quadratic splines, cubic splines, and `scipy.interpolate.CubicSpline`. The other cells, including the cells of Exercises 2.1 and 3.1, were distributed with their code filled in, and the only blank is the one between the two comment lines marked `COMPLETE` in the function `lagrange_basis` of Section 3. The first code cell of the notebook imports `numpy` as `np` and `matplotlib.pyplot` as `plt`, and imports `CubicSpline` and `interp1d` from `scipy.interpolate`.

The code blocks below reproduce cells of the notebook, and the place where code is missing is the blank to fill. The distributed blank notebook is [Lab11-Spline.ipynb](Lab11-Spline.ipynb).

- **Deadline:** one week after distribution
- **Submission:** Submit the notebook online after filling in the blank and running it.

---

<br>

## Problem. Lagrange Basis Functions (Notebook Section 3.1)

The cell in Section 1.2 of the notebook defines Dataset A as follows.

```python
# Dataset A: Air density vs. temperature (from textbook Table 17.1)
T_air  = np.array([300., 400., 500.])        # temperature [K]
rho_air = np.array([0.616, 0.525, 0.457])   # density [kg/m^3]
```

Section 2 of the notebook builds Newton's divided-difference table `b_air = newtdd_table(T_air, rho_air)` from this data and evaluates Newton's interpolating polynomial with `newtdd_eval`. The Lagrange form of Section 3 writes the interpolating polynomial directly in terms of the basis functions $`L_i(x)`$.

$$
p_n(x) = \sum_{i=0}^{n} y_i L_i(x)
$$

Each basis function satisfies $`L_i(x_j) = \delta_{ij}`$, that is, it equals 1 at its own node and 0 at all other nodes, and is given by the following.

```math
L_i(x) = \prod_{\substack{j=0 \\ j \neq i}}^{n} \frac{x - x_j}{x_i - x_j}
```

The notebook compares the two methods as follows and adds that both produce the unique polynomial of degree at most $`n`$ passing through $`n+1`$ points.

| Feature | Newton | Lagrange |
|:-----|:-----|:-----|
| Adding a data point | Easy (one new column) | Must recompute all |
| Central concept | Divided differences | Basis functions |
| Result | Identical polynomial | Identical polynomial |

*Table 1. Comparison of Newton and Lagrange interpolation*

Fill in the blank of the function `lagrange_basis` to evaluate the $`i`$th Lagrange basis function. The function `lagrange_interp` in the same cell uses it to compute interpolated values.

```python
def lagrange_basis(x, i, xx):
    """
    Evaluate the i-th Lagrange basis polynomial at xx.

    Parameters
    ----------
    x  : array of interpolation nodes
    i  : int — which basis polynomial
    xx : scalar or array — evaluation point(s)

    Returns
    -------
    Li : value(s) of L_i(xx)
    """
    x  = np.array(x, dtype=float)
    xx = np.atleast_1d(np.array(xx, dtype=float))
    n  = len(x)
    Li = np.ones(xx.shape)

    for j in range(n):
        if j != i:
            # - - - - - - - COMPLETE - - - - - - - 
            Li =                # Lagrange basis
            # - - - - - - - COMPLETE - - - - - - - 
    return Li


def lagrange_interp(x, y, xx):
    """
    Evaluate the Lagrange interpolating polynomial at xx.

    Parameters
    ----------
    x, y : array — interpolation nodes and values
    xx   : scalar or array — evaluation point(s)

    Returns
    -------
    yy : interpolated value(s)
    """
    x  = np.array(x, dtype=float)
    y  = np.array(y, dtype=float)
    xx = np.atleast_1d(np.array(xx, dtype=float))
    n  = len(x)
    yy = np.zeros(xx.shape)

    for i in range(n):
        yy = yy + y[i] * lagrange_basis(x, i, xx)
    return yy if yy.size > 1 else yy[0]
```

The next cell plots the three basis functions of Dataset A, compares the Lagrange and Newton interpolation curves, and prints both values at 350 K. The last part of the cell is as follows.

```python
rho_L = lagrange_interp(T_air, rho_air, 350.)
rho_N = newtdd_eval(T_air, b_air, 350.)
print(f'Lagrange at 350 K = {rho_L:.8f}')
print(f'Newton   at 350 K = {rho_N:.8f}')
print(f'Difference        = {abs(rho_L - rho_N):.2e}  (should be ~machine precision)')
```

The distributed notebook shows the following output for this cell.

```text
Lagrange at 350 K = 0.56762500
Newton   at 350 K = 0.56762500
Difference        = 0.00e+00  (should be ~machine precision)
```

---

<br>

After solving the problem, see the [assignment explanation](../Assignment-Explanation.md) for the answer and solution and the [completed notebook](../Solution/Lab11-Spline_Completed.ipynb) for the full code.
