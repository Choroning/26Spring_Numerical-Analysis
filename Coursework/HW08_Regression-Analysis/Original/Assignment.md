# HW8 Assignment Handout — Regression Analysis

[Repository Guide](../../../README.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **Last Updated:** 2026-11-24

This is the Week 10 homework of Numerical Analysis. The task is to fill in the blanks of the lab notebook `Lab09-Regression.ipynb` distributed by the instructor and run it. The notebook has six parts: why regression is needed, simple linear regression, linearizable nonlinear regression, general linear least squares, nonlinear regression with `scipy.optimize.curve_fit`, and model comparison and selection. There are five blanks in total, each between two comment lines marked `COMPLETE`. The first code cell of the notebook imports `numpy` as `np`, `matplotlib.pyplot` as `plt`, and `scipy.stats` as `stats`, and imports `curve_fit` from `scipy.optimize`.

Each code block below reproduces a cell of the notebook, and the places where code is missing are the blanks to fill. The distributed blank notebook is [Lab09-Regression.ipynb](Lab09-Regression.ipynb).

The cell in Section 1.2 of the notebook defines the datasets used throughout the lab as follows and draws their scatter plots.

```python
# Dataset 1: Beam deflection (clearly linear)
load       = np.array([0, 10, 20, 30, 40, 50, 60, 70], dtype=float)   # kN
deflection = np.array([0.0, 2.1, 4.3, 6.2, 8.5, 10.4, 12.6, 14.8])   # mm

# Dataset 2: Stress-response data (power law behavior)
stress   = np.linspace(1, 8, 8) * 10.0                                # MPa
response = np.array([25, 70, 380, 550, 610, 1220, 830, 1450], dtype=float)

# Dataset 3: Borosilicate glass refractive index (Cauchy equation)
wavelength = np.array([0.6563, 0.6439, 0.5890, 0.5338,
                       0.5086, 0.4861, 0.4340, 0.3988])  # micrometers
ridx       = np.array([1.50883, 1.50917, 1.51124, 1.51386,
                       1.51534, 1.51690, 1.52136, 1.52546])
```

---

<br>

## Problem 1. Coefficients of Simple Linear Regression (Notebook Section 2.1)

Fit $`\hat{y} = a_0 + a_1 x`$ by minimizing $`S_r = \sum(y_i - a_0 - a_1 x_i)^2`$. Setting the partial derivatives to zero gives the following.

$$
a_1 = \frac{n\sum x_i y_i - \sum x_i \sum y_i}{n\sum x_i^2 - \left(\sum x_i\right)^2} = \frac{\sum x_i y_i - n\bar{x}\bar{y}}{\sum x_i^2 - n\bar{x}^2}, \qquad a_0 = \bar{y} - a_1 \bar{x}
$$

Fill in the blanks of the function `linreg` to compute the slope `a1` and the intercept `a0`.

```python
def linreg(xi, yi):
    """
    Simple linear regression: fit y = a0 + a1*x via least squares.

    Returns
    -------
    a0 : intercept
    a1 : slope
    """
    xi = np.array(xi, dtype=float)
    yi = np.array(yi, dtype=float)
    n  = len(xi)

    # --- COMPLETE ---
    a1 =                  
    a0 =                  
    # --- COMPLETE ---

    return a0, a1
```

The next cell applies this function to the beam-deflection data and plots the fitted line and the residuals. Its computation and printing part is as follows.

```python
a0, a1 = linreg(load, deflection)
print(f'Slope     a1 = {a1:.6f} mm/kN')
print(f'Intercept a0 = {a0:.6f} mm')
print(f'Stiffness  k = 1/a1 = {1/a1:.4f} kN/mm')
```

The distributed notebook shows the following output for this cell.

```text
Slope     a1 = 0.210357 mm/kN
Intercept a0 = 0.000000 mm
Stiffness  k = 1/a1 = 4.7538 kN/mm
```

## Problem 2. Goodness-of-Fit Statistics (Notebook Section 2.2)

The notebook summarizes the goodness-of-fit statistics as follows.

| Quantity | Formula | Meaning |
|:-----|:-----|:-----|
| $`S_T`$ | $`\sum(y_i - \bar{y})^2`$ | Total sum of squares |
| $`S_r`$ (SSE) | $`\sum(y_i - \hat{y}_i)^2`$ | Residual (unexplained) sum of squares |
| SSR | $`S_T - S_r`$ | Explained sum of squares |
| $`R^2`$ | $`1 - S_r/S_T`$ | Fraction of variance explained |
| $`s_e`$ | $`\sqrt{S_r/(n-2)}`$ | Standard error of the estimate |

*Table 1. Goodness-of-fit statistics*

Fill in the blanks of the function `goodness_of_fit` to compute `SST`, `SSE`, `SSR`, and `R2`.

```python
def goodness_of_fit(yi, yi_hat, n_params=2):
    """
    Compute goodness-of-fit statistics.

    Returns dict: SST, SSE, SSR, R2, R2adj, se
    """
    n     = len(yi)
    y_bar = np.mean(yi)

    # --- COMPLETE ---
    SST =             # total sum of squares
    SSE =             # residual sum of squares
    SSR =             # regression sum of squares
    R2  =             # coefficient of determination
    # --- COMPLETE ---

    se    = np.sqrt(SSE / (n - n_params))
    R2adj = 1 - (SSE / (n - n_params)) / (SST / (n - 1))
    return dict(SST=SST, SSE=SSE, SSR=SSR, R2=R2, R2adj=R2adj, se=se)
```

The next cell applies this function to the line of Problem 1.

```python
yhat_lin = a0 + a1 * load
gof_lin  = goodness_of_fit(deflection, yhat_lin, n_params=2)

print('=== Goodness of Fit: Beam Deflection (Linear) ===')
print(f'SST    = {gof_lin["SST"]:.6f}   (total variation)')
print(f'SSE    = {gof_lin["SSE"]:.6f}   (residual/unexplained)')
print(f'SSR    = {gof_lin["SSR"]:.6f}   (explained by model)')
print(f'R2     = {gof_lin["R2"]:.8f}')
print(f'Adj R2 = {gof_lin["R2adj"]:.8f}')
print(f'se     = {gof_lin["se"]:.6f} mm')
print(f'SST = SSR + SSE? {np.isclose(gof_lin["SST"], gof_lin["SSR"] + gof_lin["SSE"])}')
```

The distributed notebook shows the following output for this cell.

```text
=== Goodness of Fit: Beam Deflection (Linear) ===
SST    = 185.898750   (total variation)
SSE    = 0.048214   (residual/unexplained)
SSR    = 185.850536   (explained by model)
R2     = 0.99974064
Adj R2 = 0.99969742
se     = 0.089642 mm
SST = SSR + SSE? True
```

## Problem 3. Linearizing the Power Model (Notebook Section 3.1)

Taking $`\log_{10}`$ of both sides of the power model $`y = \alpha x^\beta`$ gives a linear equation.

$$
\underbrace{\log_{10} y}_{Y} = \underbrace{\log_{10} \alpha}_{A_0} + \beta \underbrace{\log_{10} x}_{X} \quad\Rightarrow\quad Y = A_0 + \beta X
$$

The parameters are recovered as $`\beta = A_1`$ and $`\alpha = 10^{A_0}`$. The notebook cautions that log-linearization minimizes the SSE in the transformed (log) space rather than in the original $`y`$ space, so for large noise `curve_fit` of Section 5 should be used for exact nonlinear least squares. Fill in the blanks of the function `power_fit` to transform the data to log space and recover the power-model parameters.

```python
def power_fit(xi, yi):
    """
    Fit y = alpha * x^beta using log10 linearization.

    Returns
    -------
    alpha, beta : power law parameters
    A0, A1      : intercept and slope in log10 space
    """
    xi = np.array(xi, dtype=float)
    yi = np.array(yi, dtype=float)

    # --- COMPLETE: transform to log10 space ---
    X =              
    Y =              
    # --- COMPLETE ---

    A0, A1 = linreg(X, Y)      # fit in log space

    # --- COMPLETE: recover power law parameters ---
    alpha =          
    beta  =          
    # --- COMPLETE ---

    return alpha, beta, A0, A1
```

The next cell applies this function to the stress-response data, compares it in plots with the straight-line fit of the same data from the preceding cell (`gof_bad`), and prints $`R^2`$ and $`s_e`$ of both fits. The beginning of the cell is as follows.

```python
alpha_pw, beta_pw, A0, A1 = power_fit(stress, response)
print(f'Power law fit:  y = {alpha_pw:.4f} * x^{beta_pw:.4f}')
print(f'Exponent beta ~ {beta_pw:.2f}  (response ~ stress squared)')

yhat_pw = alpha_pw * stress**beta_pw
gof_pw  = goodness_of_fit(response, yhat_pw, n_params=2)
```

The distributed notebook shows the following output for this cell.

```text
Power law fit:  y = 0.2741 * x^1.9842
Exponent beta ~ 1.98  (response ~ stress squared)

Linear:    R2=0.88049,  se=189.79
Power law: R2=0.80882,  se=240.04
```

## Problem 4. Normal Equations for Polynomial Regression (Notebook Section 4.1)

Any model linear in the parameters $`\mathbf{c}`$ can be written as $`\hat{y} = \mathbf{X}\mathbf{c}`$, where the rows of the design matrix $`\mathbf{X}`$ ($`n\times(m+1)`$) are $`[z_0(x_i), z_1(x_i), \ldots, z_m(x_i)]`$. The least-squares solution satisfies the following normal equations.

$$
\mathbf{X}^T\mathbf{X}\mathbf{c} = \mathbf{X}^T\mathbf{y} \quad\Rightarrow\quad \mathbf{c} = (\mathbf{X}^T\mathbf{X})^{-1}\mathbf{X}^T\mathbf{y}
$$

| Basis functions $`z_j(x)`$ | Model |
|:-----|:-----|
| $`1, x`$ | Simple linear |
| $`1, x, x^2, \ldots, x^m`$ | Polynomial of degree $`m`$ |
| $`1, 1/\lambda^2, 1/\lambda^4`$ | Cauchy's equation |
| $`1, x_1, x_2`$ | Multiple linear regression |

*Table 2. Basis functions and models*

When $`n > m+1`$ the system is overdetermined, so use the normal equations or `np.linalg.lstsq`, which is more numerically stable. Fill in the blanks of the function `poly_lstsq` to build the design matrix and solve the normal equations.

```python
def poly_lstsq(xi, yi, m):
    """
    Polynomial regression: fit y = c0 + c1*x + ... + cm*x^m.

    Parameters
    ----------
    xi, yi : array-like
    m      : int  — polynomial degree

    Returns
    -------
    c : array  — coefficients [c0, c1, ..., cm]
    X : array  — design matrix  (n x (m+1))
    """
    xi = np.array(xi, dtype=float)
    yi = np.array(yi, dtype=float)
    n  = len(xi)

    # --- COMPLETE: build design matrix ---
    X = np.zeros((n, m + 1))
    for j in range(m + 1):
        X[:, j] =            # (j-th column = x to the power j)
    # --- COMPLETE ---

    # --- COMPLETE: solve normal equations X^T X c = X^T y ---
    c =                      
    # --- COMPLETE ---

    return c, X
```

The next cell fits polynomials of degrees 1, 2, and 3 to velocity data and plots them.

```python
# Polynomial regression on velocity data (quadratic model expected)
x_vel = np.array([0., 1., 2., 3., 4., 5.], dtype=float)  # time [s]
y_vel = np.array([2.1, 7.7, 13.6, 27.2, 40.9, 61.1])     # velocity [m/s]

fig, ax = plt.subplots(figsize=(8, 4))
ax.scatter(x_vel, y_vel, c='k', marker='s', s=80, zorder=5, label='Data')

colors = ['blue', 'red', 'green']
for m, col in zip([1, 2, 3], colors):
    c, X = poly_lstsq(x_vel, y_vel, m)
    x_plt = np.linspace(-0.2, 5.5, 300)
    y_plt = sum(c[j] * x_plt**j for j in range(m + 1))
    yhat  = X @ c
    gof_p = goodness_of_fit(y_vel, yhat, n_params=m + 1)
    ax.plot(x_plt, y_plt, color=col, lw=2,
            label=f'Deg {m}: $R^2={gof_p["R2"]:.5f}$, $s_e={gof_p["se"]:.3f}$')

ax.set(xlabel='Time [s]', ylabel='Velocity [m/s]',
       title='Polynomial Regression: Degrees 1, 2, 3')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
plt.tight_layout(); plt.show()

c2, _ = poly_lstsq(x_vel, y_vel, 2)
print(f'Quadratic: y = {c2[0]:.4f} + {c2[1]:.4f}t + {c2[2]:.4f}t²')
```

The distributed notebook shows the following output for this cell.

```text
Quadratic: y = 2.4786 + 2.3593t + 1.8607t²
```

## Problem 5. Nonlinear Regression of the Antoine Equation (Notebook Section 5.1)

Some models cannot be linearized, such as the Michaelis-Menten equation $`y = V_{\max} x / (K_m + x)`$ and the Antoine equation $`\log_{10} P = A - B/(C + T)`$. `curve_fit` minimizes $`\sum(y_i - f(x_i;\mathbf{c}))^2`$ directly with the Levenberg-Marquardt algorithm, a damped Newton method. The notebook shows its usage as follows.

```python
popt, pcov = curve_fit(f, xdata, ydata, p0=initial_guess)
# popt : optimal parameters
# pcov : covariance matrix  ->  perr = np.sqrt(np.diag(pcov))  gives std errors
```

Fill in the blank of the next cell to fit the Antoine equation to the vapor-pressure data of water.

```python
# Antoine equation: log10(P) = A - B/(C+T)  -> P = 10^(A - B/(C+T))
T_ant = np.array([40, 50, 60, 70, 80, 90, 100], dtype=float)    # degC
P_ant = np.array([55.3, 92.5, 149.4, 233.7, 355.1, 525.8, 760.0])  # mmHg

def antoine_model(T, A, B, C):
    return 10.0**(A - B / (C + T))

# --- COMPLETE: fit the Antoine equation ---
p0_ant = [8.0, 1500.0, 200.0]
popt_ant, pcov_ant =     
# --- COMPLETE ---

A_opt, B_opt, C_opt = popt_ant
perr_ant = np.sqrt(np.diag(pcov_ant))
gof_ant  = goodness_of_fit(P_ant, antoine_model(T_ant, *popt_ant), n_params=3)

print('=== Antoine Equation Fit (water vapor pressure) ===')
print(f'A = {A_opt:.6f} +/- {perr_ant[0]:.6f}')
print(f'B = {B_opt:.4f} +/- {perr_ant[1]:.4f}')
print(f'C = {C_opt:.4f} +/- {perr_ant[2]:.4f}')
print(f'R2 = {gof_ant["R2"]:.8f}')
print(f'Literature: A=8.07131, B=1730.63, C=233.426')

T_fit = np.linspace(35, 105, 300)
fig, ax = plt.subplots(figsize=(8, 4))
ax.scatter(T_ant, P_ant, c='k', marker='s', s=80, zorder=5, label='Data')
ax.plot(T_fit, antoine_model(T_fit, *popt_ant), 'r-', lw=2,
        label=f'Antoine fit $R^2={gof_ant["R2"]:.6f}$')
ax.set(xlabel='Temperature [°C]', ylabel='Vapor pressure [mmHg]',
       title='Antoine Equation: Vapor Pressure of Water')
ax.legend(); ax.grid(True, alpha=0.3)
plt.tight_layout(); plt.show()
```

The distributed notebook shows the following output for this cell.

```text
=== Antoine Equation Fit (water vapor pressure) ===
A = 7.997574 +/- 0.007401
B = 1687.0509 +/- 4.6041
C = 229.7113 +/- 0.4237
R2 = 0.99999999
Literature: A=8.07131, B=1730.63, C=233.426
```

---

<br>

After solving the problems, see the [assignment explanation](../Assignment-Explanation.md) for the answers and solutions and the [completed notebook](../Solution/Lab09-Regression_Completed.ipynb) for the full code.
