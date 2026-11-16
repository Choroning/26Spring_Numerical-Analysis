# HW1 Assignment Handout — Sources of Numerical Errors

[Repository Guide](../../../README.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **Last Updated:** 2026-11-16

This is the Week 2 homework of Numerical Analysis. The task is to fill in the blanks of the lab notebook `Lab02-SourceOfErrors.ipynb` distributed by the instructor and run it. The notebook has four parts: error definitions and metrics, truncation error and Taylor series, roundoff error and floating-point arithmetic, and finite differences with combined error analysis. There are five blanks in total, each between two comment lines marked `COMPLETE` or `Complete`. The first code cell of the notebook imports `numpy` as `np`, `matplotlib.pyplot` as `plt`, and `math`.

Each code block below reproduces a cell of the notebook, and the places where code is missing are the blanks to fill. The distributed blank notebook is [Lab02-SourceOfErrors.ipynb](Lab02-SourceOfErrors.ipynb).

- **Deadline:** one week after distribution
- **Submission:** Submit the notebook online after filling in the blanks and running it.

---

<br>

## Problem 1. Maclaurin Partial Sums of the Exponential (Notebook Section 2.2)

The exponential function has the especially clean Maclaurin expansion below.

$$
e^x = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \cdots = \sum_{k=0}^{\infty} \frac{x^k}{k!}
$$

The notebook approximates $`e^{0.5}`$ by adding terms one at a time and watches the error decrease. The preceding cell sets the true value as follows.

```python
# True value
x = 0.5
f_true = np.exp(x)
print(f'True value e^{x} = {f_true:.15f}')
```

Fill in the blanks of the next cell to compute the remaining partial sums `f_approx[3]`, `f_approx[4]`, and `f_approx[5]`.

```python
# Build approximations term by term  (k=0,1,2,...,5 → 6 terms)
# f_approx[k] = sum_{i=0}^{k} x^i / i!

# --- COMPLETE: fill in the missing partial sums ---
f_approx = np.zeros(6)
f_approx[0] = 1                               # k=0: just the constant term
f_approx[1] = f_approx[0] + x                # k=1: + x
f_approx[2] = f_approx[1] + x**2 / math.factorial(2)   # k=2: + x^2/2!
f_approx[3] =                                # k=3: + x^3/3!
f_approx[4] =                                # k=4: + x^4/4!
f_approx[5] =                                # k=5: + x^5/5!
# --- COMPLETE ---

print(f'{"Terms":>6}  {"Approximation":>20}  {"True Rel Err (%)":>18}')
print('-' * 50)
for k in range(6):
    et = abs(f_true - f_approx[k]) / f_true * 100
    print(f'{k+1:>6}  {f_approx[k]:>20.15f}  {et:>18.10f}')
```

## Problem 2. Series Evaluation with a Stopping Criterion (Notebook Section 2.3)

Rather than fixing the number of terms in advance, a practical algorithm keeps adding terms until $`\varepsilon_a < \varepsilon_s`$. Section 1 of the notebook defines the approximate relative error as follows.

```math
\varepsilon_a = \frac{| x_{\text{new}} - x_{\text{prev}} |}{|x_{\text{new}}|} \times 100\%
```

In practice the true value is unknown, and the Scarborough criterion tells how small $`\varepsilon_a`$ must be to guarantee $`n`$ correct significant figures in the result.

```math
\varepsilon_a < \varepsilon_s = (0.5 \times 10^{2-n})\%
```

The notebook defines this criterion as `scarborough = lambda n: 0.5 * 10**(2 - n)`. Fill in the blank of the function `exp_series` so that each iteration adds the next term to the partial sum.

```python
def exp_series(x, es=1e-4, maxit=50):
    """
    Approximate e^x using the Maclaurin series, stopping when
    approximate relative error < es (%) or after maxit iterations.

    Parameters
    ----------
    x     : evaluation point
    es    : stopping criterion in % (default 1e-4 %)
    maxit : maximum number of terms (default 50)

    Returns
    -------
    fx   : approximated value
    ea   : final approximate relative error (%)
    k    : number of terms used
    """
    sol = 1.0       # k=0 term
    ea  = 100.0

    print(f'{"k":>4}  {"Approximation":>20}  {"ea (%)":>15}')
    print('-' * 45)

    for k in range(1, maxit + 1):
        sol_old = sol
        # - - - Complete - - -
        sol = sol + 
        # - - - Complete - - -
        if sol != 0:
            ea = abs(sol - sol_old) / abs(sol) * 100
        print(f'{k:>4}  {sol:>20.15f}  {ea:>15.10f}')
        if ea < es:
            break

    return sol, ea, k
```

The next cell runs the function with four significant figures as the target.

```python
# Target: 4 significant figures → es = 0.005 %
es = scarborough(4)
print(f'Stopping criterion: {es} %\n')
fx, ea, k = exp_series(0.5, es=es)
print(f'\nResult : {fx:.10f}')
print(f'True   : {np.exp(0.5):.10f}')
print(f'Terms  : {k}')
```

The notebook gives the expected output of this cell as follows.

```text
Stopping criterion: 0.005 %

   k         Approximation           ea (%)
---------------------------------------------
   1     1.500000000000000    33.3333333333
   2     1.625000000000000     7.6923076923
   3     1.645833333333333     1.2658227848
   4     1.648437500000000     0.1579778831
   5     1.648697916666667     0.0157952930
   6     1.648719618055555     0.0013162571

Result : 1.6487196181
True   : 1.6487212707
Terms  : 6
```

## Problem 3. Quadratic Roots and Catastrophic Cancellation (Notebook Section 3.5)

When two nearly equal numbers are subtracted, the significant digits cancel and the result is dominated by roundoff. This is called **catastrophic (or subtractive) cancellation**. For example, in the quadratic formula $`x = \dfrac{-b \pm \sqrt{b^2 - 4ac}}{2a}`$ with $`b^2 \gg 4ac`$, one root involves subtracting two nearly equal large numbers. The notebook first compares the standard formula with a stable form on the equation $`x^2 - 100000.0001x + 1 = 0`$.

```python
# Equation: x^2 - 100000.0001 x + 1 = 0
# Roots: x1 ≈ 100000.0001, x2 ≈ 0.000010000
a, b, c = 1, -100000.0001, 1

discriminant = math.sqrt(b**2 - 4*a*c)

# Standard quadratic formula
x1_std = (-b + discriminant) / (2*a)
x2_std = (-b - discriminant) / (2*a)    # cancellation occurs here!

# Numerically stable form: use x2 = c / (a * x1)
x2_stable = c / (a * x1_std)

# True roots (high precision)
x1_true = 100000.00000999990
x2_true = 0.000010000000100000

print(f'Standard   x1 = {x1_std:.10f}   error = {abs(x1_std - x1_true):.2e}')
print(f'Standard   x2 = {x2_std:.10f}   error = {abs(x2_std - x2_true):.2e}')
print(f'Stable     x2 = {x2_stable:.10f}   error = {abs(x2_stable - x2_true):.2e}')
```

```text
Standard   x1 = 100000.0000900000   error = 8.00e-05
Standard   x2 = 0.0000100000   error = 3.99e-12
Stable     x2 = 0.0000100000   error = 1.09e-13
```

The exercise cell that follows investigates how the error of $`x_2`$ from the standard formula changes as $`|b|`$ grows. Fill in the blanks `x1_s` and `x2_stab` so that the cell computes both roots and the relative error in $`x_2`$ for each value of $`b`$ below.

```python
# --- Exercise ---
# Investigate how the error in x2_std changes as |b| grows.
# For each value of b below, compute both roots and the relative error in x2.

print(f'{"b":>12}  {"x2 standard":>16}  {"x2 stable":>16}  {"rel err std (%)":>17}  {"rel err stable (%)":>20}')
print('-' * 90)

for b_val in [1e2, 1e4, 1e6, 1e8, 1e10]:
    a_, b_, c_ = 1, -b_val, 1
    disc = math.sqrt(b_**2 - 4*a_*c_)
    # --- COMPLETE: compute x2_std and x2_stable ---
    x1_s = 
    x2_s = (-b_ - disc) / (2*a_)  # standard formula
    x2_stab =   # stable formula: c / (a * x1)
    # --- COMPLETE ---
    x2_ref = 1.0 / x1_s             # reference (stable)
    err_std  = abs(x2_s    - x2_ref) / abs(x2_ref) * 100
    err_stab = abs(x2_stab - x2_ref) / abs(x2_ref) * 100
    print(f'{b_val:>12.0e}  {x2_s:>16.10f}  {x2_stab:>16.10f}  {err_std:>17.6e}  {err_stab:>20.6e}')
```

The distributed notebook shows the following output for this cell.

```text
           b       x2 standard         x2 stable    rel err std (%)    rel err stable (%)
------------------------------------------------------------------------------------------
       1e+02      0.0100010002      0.0100010002       1.221123e-11          0.000000e+00
       1e+04      0.0001000000      0.0001000000       1.117663e-07          0.000000e+00
       1e+06      0.0000010000      0.0000010000       7.614492e-04          0.000000e+00
       1e+08      0.0000000075      0.0000000100       2.549419e+01          0.000000e+00
       1e+10      0.0000000000      0.0000000001       1.000000e+02          0.000000e+00
```

## Problem 4. Three Finite-Difference Formulas (Notebook Section 4.3)

Isolating $`f'(x)`$ from the Taylor expansion $`f(x + h) = f(x) + f'(x)h + \frac{f''(x)}{2!}h^2 + \cdots`$ gives the following approximations.

| Approximation | Formula | Error order |
|:-----|:-----|:-----|
| Forward difference | $`f'(x) \approx \dfrac{f(x+h)-f(x)}{h}`$ | $`O(h)`$ |
| Backward difference | $`f'(x) \approx \dfrac{f(x)-f(x-h)}{h}`$ | $`O(h)`$ |
| Central difference | $`f'(x) \approx \dfrac{f(x+h)-f(x-h)}{2h}`$ | $`O(h^2)`$ |

*Table 1. Three finite-difference formulas from the Taylor expansion*

The example function is $`f(x) = -0.1x^4 - 0.15x^3 - 0.5x^2 - 0.25x + 1.2`$, and its analytical derivative is $`f'(x) = -0.4x^3 - 0.45x^2 - x - 0.25`$. The preceding cell evaluates the needed function values at $`x = 0.5`$ with $`h = 0.5`$.

```python
f     = lambda x: -0.1*x**4 - 0.15*x**3 - 0.5*x**2 - 0.25*x + 1.2
dfdx  = lambda x: -0.4*x**3 - 0.45*x**2 -       x - 0.25

x  = 0.5
h  = 0.5

fm = f(x - h)   # f(x-h)
fi = f(x)       # f(x)
fp = f(x + h)   # f(x+h)

print(f'f(x-h) = {fm:8.4f}')
print(f'f(x)   = {fi:8.4f}')
print(f'f(x+h) = {fp:8.4f}')
print(f"True f'(x) = {dfdx(x):8.4f}")
```

```text
f(x-h) =   1.2000
f(x)   =   0.9250
f(x+h) =   0.2000
True f'(x) =  -0.9125
```

Fill in the three finite-difference formulas in the blanks of the next cell.

```python
# --- COMPLETE: fill in the three finite difference formulas ---
dfdx_forward  =                 
dfdx_backward =                 
dfdx_central  =                 
# --- COMPLETE ---

true_val = dfdx(x)
et_fwd = abs(true_val - dfdx_forward)  / abs(true_val) * 100
et_bwd = abs(true_val - dfdx_backward) / abs(true_val) * 100
et_cen = abs(true_val - dfdx_central)  / abs(true_val) * 100

print(f'{"Method":>12}  {"Estimate":>10}  {"True Rel Err (%)":>18}')
print('-' * 45)
print(f'{"Forward":>12}  {dfdx_forward:>10.5f}  {et_fwd:>18.4f}')
print(f'{"Backward":>12}  {dfdx_backward:>10.5f}  {et_bwd:>18.4f}')
print(f'{"Central":>12}  {dfdx_central:>10.5f}  {et_cen:>18.4f}')
print(f'{"True":>12}  {true_val:>10.5f}')
```

## Problem 5. Forward and Central Difference Errors versus Step Size (Notebook Section 4.5)

Section 4.4 of the notebook explains that as the step size $`h`$ decreases, the truncation error of the central difference shrinks in proportion to $`h^2`$ while the roundoff error grows in proportion to $`\varepsilon_{\text{mach}} / h`$, so an optimal step size $`h^* \approx \varepsilon_{\text{mach}}^{1/3}`$ minimizes the total error. It then performs a step-size analysis of the central difference with the following cell.

```python
def step_size_analysis(func, dfunc, x, n_steps=15):
    """
    Compute central-difference derivative for step sizes h = 10^0, 10^-1, ..., 10^{-(n-1)}
    and record the absolute error against the analytical derivative.
    """
    true_val = dfunc(x)
    H = np.zeros(n_steps)
    E = np.zeros(n_steps)

    h = 1.0
    for i in range(n_steps):
        H[i] = h
        approx = (func(x + h) - func(x - h)) / (2 * h)
        E[i]   = abs(true_val - approx)
        h /= 10

    return H, E

H, E = step_size_analysis(f, dfdx, x=0.5)

print(f'{"Step size":>14}  {"Approx deriv":>18}  {"Abs error":>14}')
print('-' * 52)
for i in range(len(H)):
    approx = (f(0.5 + H[i]) - f(0.5 - H[i])) / (2 * H[i])
    print(f'{H[i]:>14.1e}  {approx:>18.14f}  {E[i]:>14.4e}')
```

The exercise reads as follows. Repeat the step-size analysis for the forward and central differences and plot both on the same figure to visually confirm that central differences converge as $`O(h^2)`$ while forward differences converge as $`O(h)`$. Fill in the blanks `E_fwd[i]` and `E_cen[i]` of the next cell.

```python
n_steps = 15
H  = np.zeros(n_steps)
E_fwd = np.zeros(n_steps)
E_cen = np.zeros(n_steps)

true_val = dfdx(0.5)
h = 1.0
for i in range(n_steps):
    H[i] = h
    # - - - - Complete - - - 
    E_fwd[i] =  
    E_cen[i] =  
    # - - - - Complete - - - 
    h /= 10

plt.figure(figsize=(8, 5))
plt.loglog(H, E_fwd, 'rs-', markersize=6, label='Forward difference')
plt.loglog(H, E_cen, 'bo-', markersize=6, label='Central difference')

h_ref = np.logspace(0, -8, 80)
plt.loglog(h_ref, 0.5  * h_ref,    'r--', alpha=0.5, label=r'$O(h)$')
plt.loglog(h_ref, 0.05 * h_ref**2, 'b--', alpha=0.5, label=r'$O(h^2)$')

plt.xlabel('Step size $h$')
plt.ylabel('Absolute error')
plt.title('Forward vs. central difference: convergence rate')
plt.legend()
plt.grid(True, which='both', ls=':', alpha=0.5)
plt.tight_layout()
```

---

<br>

After solving the problems, see the [assignment explanation](../Assignment-Explanation.md) for the answers and solutions and the [completed notebook](../Solution/Lab02-SourceOfErrors_Completed.ipynb) for the full code.
