# HW6 Assignment Handout — Root-Finding Methods

[Repository Guide](../../../README.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **Last Updated:** 2026-11-17

This is the Week 7 homework of Numerical Analysis. The task is to fill in the blanks of the lab notebook `Lab07-RootSolving.ipynb` distributed by the instructor and run it. The notebook has seven parts: motivation and the graphical approach, incremental search, bisection, false position, Newton-Raphson, the secant method, and `scipy.optimize` with a method comparison. There are five blanks in total, each between two comment lines marked `COMPLETE`. The two Exercise cells that apply bisection and Newton-Raphson to the savings problem were distributed with their code filled in. The first code cell of the notebook imports `numpy` as `np`, `matplotlib.pyplot` as `plt`, `brentq` and `root_scalar` from `scipy.optimize`, and `partial` from `functools`.

Each code block below reproduces a cell of the notebook, and the places where code is missing are the blanks to fill. The distributed blank notebook is [Lab07-RootSolving.ipynb](Lab07-RootSolving.ipynb).

- **Deadline:** one week after distribution
- **Submission:** Submit the notebook online after filling in the blanks and running it.

---

<br>

## Example Problems of the Notebook (Notebook Section 1)

A **root** of a function $`f`$ is a value $`x^*`$ such that $`f(x^*) = 0`$. The notebook uses two engineering problems as examples.

- **Savings account:** Shin deposits $`p`$ dollars per month at annual interest rate $`r`$, and the total after $`n`$ payments is $`A = \frac{12p}{r}\left[\left(1+\frac{r}{12}\right)^n - 1\right]`$. The goal is to find $`r`$ such that $`A = 1500`$ dollars after $`n = 24`$ monthly deposits of 60 dollars, with $`f(r) = \frac{12p}{r}\left[(1+r/12)^n - 1\right] - A = 0`$.
- **Bungee jumper:** The velocity of a falling bungee jumper of mass $`m`$ at time $`t`$ is $`v(m) = \sqrt{mg/c_d} \tanh\left(\sqrt{c_d g/m}\, t\right)`$. The goal is to find the mass $`m`$ such that $`v = 36`$ m/s at $`t = 4`$ s with the drag coefficient $`c_d = 0.25`$ kg/m, with $`f(m) = \sqrt{mg/c_d} \tanh\left(\sqrt{c_d g/m}\, t\right) - v = 0`$.

The next cell defines the two functions. The plotting part is omitted here.

```python
# Problem parameters
cd, g_acc, v_target, t_fall = 0.25, 9.81, 36.0, 4.0
deposit, duration, target_amount = 60.0, 24, 1500.0

# Define the two root-finding functions
def f_bungee(m):
    """f(m) = v(m) - v_target  (bungee jumper mass problem)"""
    return np.sqrt(m * g_acc / cd) * np.tanh(np.sqrt(cd * g_acc / m) * t_fall) - v_target

def f_savings(r):
    """f(r) = A(r) - target  (savings account interest rate)"""
    ratio = 1.0 + r / 12.0
    return 12.0 * deposit / r * (ratio**duration - 1.0) - target_amount
```

The graph shows that the bungee jumper root lies in $`[50, 250]`$, with $`f(50) < 0`$ and $`f(250) > 0`$. If $`f`$ is continuous and $`f(x_l) \cdot f(x_u) < 0`$, then at least one root exists in $`[x_l, x_u]`$ (Intermediate Value Theorem).

## Problem 1. Incremental Search (Notebook Section 2)

Before applying any root-finding method, brackets that enclose the roots must be located. The incremental search divides $`[x_{\min}, x_{\max}]`$ into $`n_s - 1`$ equal subintervals and checks each for a sign change.

1. Create $`n_s`$ equally spaced points $`x_0, x_1, \ldots, x_{n_s-1}`$ in $`[x_{\min}, x_{\max}]`$.
2. For each pair $`(x_k, x_{k+1})`$, record the bracket if $`f(x_k) \cdot f(x_{k+1}) < 0`$.
3. Return all brackets found.

If $`n_s`$ is too small, a bracket may contain an even number of roots, which are then missed. Fill in the blank of the following function with the condition that detects a sign change.

```python
def incsearch(func, xmin, xmax, ns=50):
    """
    Find all sign-change brackets of func in [xmin, xmax].

    Parameters
    ----------
    func : callable
    xmin, xmax : float  — search interval
    ns   : int          — number of points (ns-1 subintervals)

    Returns
    -------
    nb   : int          — number of brackets found
    xb   : list of (xl, xu) tuples
    """
    x = np.linspace(xmin, xmax, ns)
    nb = 0
    xb = []

    for k in range(ns - 1):
        # --- COMPLETE: detect sign change between x[k] and x[k+1] ---
        if :
            nb += 1
            xb.append((x[k], x[k+1]))
        # --- COMPLETE ---

    return nb, xb
```

The next cell tests the completed function on a function with multiple roots. The plotting part is omitted here.

```python
# Test incsearch on a function with multiple roots
f_multi = lambda x: np.sin(10 * x) + np.cos(3 * x)
xmin, xmax = 3.0, 6.0

for ns in [20, 50, 100]:
    nb, xb = incsearch(f_multi, xmin, xmax, ns=ns)
    print(f'ns={ns:3d}  ->  {nb} bracket(s) found')

# Best result
nb, xb = incsearch(f_multi, xmin, xmax, ns=100)
print(f'\nBrackets (ns=100):')
for xl_b, xu_b in xb:
    print(f'  [{xl_b:.4f}, {xu_b:.4f}]')
```

The distributed notebook shows the following output for this cell.

```text
ns= 20  ->  7 bracket(s) found
ns= 50  ->  5 bracket(s) found
ns=100  ->  9 bracket(s) found

Brackets (ns=100):
  [3.2424, 3.2727]
  [3.3636, 3.3939]
  [3.7273, 3.7576]
  [4.2121, 4.2424]
  [4.2424, 4.2727]
  [4.6970, 4.7273]
  [5.1515, 5.1818]
  [5.1818, 5.2121]
  [5.6667, 5.6970]
```

## Problem 2. Bisection Method (Notebook Section 3.1)

Given a bracket $`[x_l, x_u]`$ with $`f(x_l) \cdot f(x_u) < 0`$, bisection repeats the following steps.

1. Compute the midpoint $`x_r = \dfrac{x_l + x_u}{2}`$.
2. Compute the approximate relative error $`\varepsilon_a = \left|\dfrac{x_r^{\rm new} - x_r^{\rm old}}{x_r^{\rm new}}\right| \times 100\%`$.
3. Update the bracket. If $`f(x_l) \cdot f(x_r) < 0`$, the root is in $`[x_l, x_r]`$, so set $`x_u = x_r`$; otherwise the root is in $`[x_r, x_u]`$, so set $`x_l = x_r`$.
4. Repeat until $`\varepsilon_a < \varepsilon_s`$ or the maximum number of iterations is reached.

The bracket width is halved at each step, so the error after $`n`$ iterations is at most $`(x_u - x_l)/2^n`$. This is linear convergence with rate 0.5, and the minimum number of iterations to achieve the absolute tolerance $`\varepsilon_{\rm abs}`$ is $`n \geq \left\lceil \log_2\!\left(\frac{x_u - x_l}{\varepsilon_{\rm abs}}\right) \right\rceil`$. The method always converges if the initial bracket is valid, but the convergence is slow. Fill in the blanks of the following function.

```python
def bisect(func, xl, xu, es=1e-7, maxit=50):
    """
    Find a root of func in [xl, xu] using the bisection method.

    Parameters
    ----------
    func : callable — continuous function
    xl, xu : float  — initial bracket  (must satisfy f(xl)*f(xu) < 0)
    es     : float  — stopping criterion (approx. relative error)
    maxit  : int    — maximum iterations

    Returns
    -------
    xr      : float — root estimate
    fxr     : float — f(xr)
    ea      : float — final approximate relative error
    n_iter  : int   — iterations performed
    history : list of (iter, xr, f(xr), ea)
    """
    if func(xl) * func(xu) > 0:
        raise ValueError('f(xl) and f(xu) must have opposite signs (bracket check failed)')

    xr_old = xl       # previous midpoint (initialized to xl)
    history = []

    for i in range(maxit):
        # --- COMPLETE ---
        xr =                                # — midpoint
        ea = abs((xr - xr_old) / xr) if xr != 0 else abs(xr - xr_old)  # approx. relative error

        if func(xl) * func(xr) < 0:     #  — sign test: root in left half?
            xu = xr
        else:
            xl =                             #      — root in right half
        # --- COMPLETE ---

        history.append((i + 1, xr, func(xr), ea))
        xr_old = xr

        if ea < es:
            break

    return xr, func(xr), ea, i + 1, history
```

The next cell applies the completed function to the bungee jumper problem.

```python
# Apply bisection to the bungee jumper problem
m_true = 142.7376  # known true root (for error analysis)

xr, fxr, ea, n_iter, history = bisect(f_bungee, 50, 250, es=1e-9, maxit=60)

print(f'Root found : m = {xr:.8f} kg')
print(f'f(root)    : {fxr:.3e}')
print(f'Rel. error : {ea:.3e}')
print(f'Iterations : {n_iter}')

# Print iteration table
print(f'\n{"Iter":>5}  {"xr":>14}  {"f(xr)":>14}  {"ea":>12}')
print('-' * 52)
for row in history:
    print(f'{row[0]:>5}  {row[1]:>14.8f}  {row[2]:>14.6e}  {row[3]:>12.6e}')
```

The distributed notebook shows the following output for this cell.

```text
Root found : m = 142.73763308 kg
f(root)    : -5.203e-10
Rel. error : 6.525e-10
Iterations : 31

 Iter              xr           f(xr)            ea
----------------------------------------------------
    1    150.00000000    1.420430e-01  6.666667e-01
    2    100.00000000   -1.197380e+00  5.000000e-01
    3    125.00000000   -4.086015e-01  2.000000e-01
    4    137.50000000   -1.108041e-01  9.090909e-02
    5    143.75000000    2.057713e-02  4.347826e-02
    6    140.62500000   -4.379968e-02  2.222222e-02
    7    142.18750000   -1.129258e-02  1.098901e-02
    8    142.96875000    4.720768e-03  5.464481e-03
    9    142.57812500   -3.266139e-03  2.739726e-03
   10    142.77343750    7.322383e-04  1.367989e-03
   11    142.67578125   -1.265717e-03  6.844627e-04
   12    142.72460938   -2.664315e-04  3.421143e-04
   13    142.74902344    2.329804e-04  1.710279e-04
   14    142.73681641   -1.670630e-05  8.552125e-05
   15    142.74291992    1.081418e-04  4.275880e-05
   16    142.73986816    4.571898e-05  2.137986e-05
   17    142.73834229    1.450664e-05  1.069004e-05
   18    142.73757935   -1.099755e-06  5.345050e-06
   19    142.73796082    6.703461e-06  2.672518e-06
   20    142.73777008    2.801858e-06  1.336261e-06
   21    142.73767471    8.510526e-07  6.681308e-07
   22    142.73762703   -1.243508e-07  3.340655e-07
   23    142.73765087    3.633510e-07  1.670327e-07
   24    142.73763895    1.195001e-07  8.351637e-08
   25    142.73763299   -2.425352e-09  4.175818e-08
   26    142.73763597    5.853737e-08  2.087909e-08
   27    142.73763448    2.805601e-08  1.043955e-08
   28    142.73763373    1.281533e-08  5.219773e-09
   29    142.73763336    5.194991e-09  2.609887e-09
   30    142.73763318    1.384819e-09  1.304943e-09
   31    142.73763308   -5.202665e-10  6.524716e-10
```

## Problem 3. False Position (Notebook Section 4.1)

Instead of bisecting at the **midpoint**, false position (regula falsi) finds the x-intercept of the **line** connecting $`(x_l, f(x_l))`$ and $`(x_u, f(x_u))`$.

```math
x_r = x_u - f(x_u)\,\frac{x_l - x_u}{f(x_l) - f(x_u)}
```

This is equivalent to linear interpolation, so if $`f`$ is approximately linear on $`[x_l, x_u]`$, it gives a better estimate than the midpoint. The bracket update rule (sign test) is the same as in bisection. The geometry is sketched below.

```text
  f(xl) *
          \                 <- True f(x)
           \    ___---
  0 --------\--/--------  x-axis
              \/ <- xr (line intercept)
          *  f(xu)
```

If $`f`$ is strongly curved, one bracket endpoint often stays fixed for many iterations, which can make false position slower than bisection in some cases. Fill in the blank of the following function.

```python
def false_position(func, xl, xu, es=1e-7, maxit=50):
    """
    Find a root of func in [xl, xu] using the False Position (Regula Falsi) method.

    Parameters
    ----------
    func : callable
    xl, xu : float  — initial bracket
    es     : float  — stopping criterion (approx. relative error)
    maxit  : int    — maximum iterations

    Returns
    -------
    xr, f(xr), ea, n_iter, history
    """
    if func(xl) * func(xu) > 0:
        raise ValueError('f(xl) and f(xu) must have opposite signs')

    xr_old = xl
    history = []

    for i in range(maxit):
        fl, fu = func(xl), func(xu)

        # --- COMPLETE: linear interpolation (secant line x-intercept) ---
        xr = 
        # --- COMPLETE ---

        fr = func(xr)
        ea = abs((xr - xr_old) / xr) if xr != 0 else abs(xr - xr_old)

        history.append((i + 1, xr, fr, ea))
        xr_old = xr

        if ea < es:
            break

        if fl * fr < 0:
            xu = xr
        else:
            xl = xr

    return xr, func(xr), ea, i + 1, history
```

The next cell compares bisection and false position on the bungee jumper problem. The plotting part is omitted here.

```python
# Compare bisection vs. false position on the bungee jumper problem
xr_bs, _, _, _, hist_bs = bisect(f_bungee, 50, 250, es=1e-10, maxit=60)
xr_fp, _, _, _, hist_fp = false_position(f_bungee, 50, 250, es=1e-10, maxit=60)

err_bs = [abs(row[1] - m_true) for row in hist_bs]
err_fp = [abs(row[1] - m_true) for row in hist_fp]

print(f'Bisection:      {len(hist_bs):2d} iterations,  final |error| = {err_bs[-1]:.3e}')
print(f'False Position: {len(hist_fp):2d} iterations,  final |error| = {err_fp[-1]:.3e}')
```

The distributed notebook shows the following output for this cell.

```text
Bisection:      34 iterations,  final |error| = 3.312e-05
False Position: 42 iterations,  final |error| = 3.313e-05
```

## Problem 4. Newton-Raphson Method (Notebook Section 5.1)

The Newton-Raphson method uses the derivative $`f'(x)`$ to construct the tangent line at the current guess and takes its x-intercept as the next guess.

$$
x_{i+1} = x_i - \frac{f(x_i)}{f'(x_i)}
$$

Drawing the tangent line $`y = f(x_i) + f'(x_i)(x - x_i)`$ at $`x_i`$, setting $`y = 0`$, and solving for $`x`$ gives $`x_{i+1}`$. Near the root $`x^*`$, the error $`e_i = x_i - x^*`$ satisfies $`e_{i+1} \approx -\frac{f''(x^*)}{2f'(x^*)} e_i^2`$, which is quadratic convergence: the number of correct decimal digits roughly doubles each iteration. The method needs $`f'(x)`$ analytically, however, and may diverge if $`f'(x_i) \approx 0`$ or the starting guess is poor.

| Failure mode | Cause | Remedy |
|:-----|:-----|:-----|
| Divergence | Initial guess too far from the root | Use a bracketing method first |
| Oscillation | Local extremum near the root | Use a smaller step or switch to bisection |
| Zero derivative | $`f'(x_i) \approx 0`$ | Add a small perturbation or switch methods |

*Table 1. Cases in which the Newton-Raphson method fails*

Fill in the blank of the following function.

```python
def newtraph(f, fp, x0, es=1e-7, maxit=30):
    """
    Find a root of f using the Newton-Raphson method.

    Parameters
    ----------
    f    : callable — function
    fp   : callable — derivative of f
    x0   : float   — initial guess
    es   : float   — stopping criterion (approx. relative error)
    maxit: int     — maximum iterations

    Returns
    -------
    x1, f(x1), ea, n_iter, history
    """
    history = []
    xi = x0

    for i in range(maxit):
        fxi  = f(xi)
        fpxi = fp(xi)

        if abs(fpxi) < 1e-30:
            raise ZeroDivisionError('Derivative near zero — Newton-Raphson fails')

        # --- COMPLETE ---
        x1 =                                # — Newton-Raphson update
        ea = abs((x1 - xi) / x1) if x1 != 0 else abs(x1 - xi)
        # --- COMPLETE ---

        history.append((i + 1, x1, f(x1), ea))

        if ea < es:
            break

        xi = x1

    return x1, f(x1), ea, i + 1, history
```

A preceding cell defines the test function as `f_nr = lambda x: x - np.exp(-x)` (root about 0.5671), `fp_nr = lambda x: 1.0 + np.exp(-x)`, and `x_true = 0.5671432904097838`. The next cell applies the completed function to it.

```python
# Test newtraph on f(x) = x - e^{-x}, root = 0.5671...
xr_nr, fxr_nr, ea_nr, n_nr, hist_nr = newtraph(f_nr, fp_nr, x0=0.0, es=1e-12, maxit=20)

print(f'Root found : {xr_nr:.15f}')
print(f'True root  : {x_true:.15f}')
print(f'|error|    : {abs(xr_nr - x_true):.3e}')
print(f'Iterations : {n_nr}')

print(f'\n{"Iter":>5}  {"xi":>18}  {"f(xi)":>14}  {"ea":>14}')
print('-' * 58)
for row in hist_nr:
    print(f'{row[0]:>5}  {row[1]:>18.12f}  {row[2]:>14.6e}  {row[3]:>14.6e}')
```

The distributed notebook shows the following output for this cell.

```text
Root found : 0.567143290409784
True root  : 0.567143290409784
|error|    : 1.110e-16
Iterations : 5

 Iter                  xi           f(xi)              ea
----------------------------------------------------------
    1      0.500000000000   -1.065307e-01    1.000000e+00
    2      0.566311003197   -1.304510e-03    1.170929e-01
    3      0.567143165035   -1.964805e-07    1.467287e-03
    4      0.567143290410   -4.440892e-15    2.210639e-07
    5      0.567143290410    1.110223e-16    5.089684e-15
```

## Problem 5. Secant Method (Notebook Section 6.1)

Computing $`f'(x)`$ analytically is sometimes difficult or expensive. The secant method approximates the derivative using the two previous iterates.

$$
f'(x_i) \approx \frac{f(x_i) - f(x_{i-1})}{x_i - x_{i-1}}
$$

Substituting this into the Newton-Raphson formula gives the following.

```math
x_{i+1} = x_i - f(x_i)\,\frac{x_{i-1} - x_i}{f(x_{i-1}) - f(x_i)}
```

The secant method converges with order equal to the golden ratio $`\phi = \frac{1+\sqrt{5}}{2} \approx 1.618`$. This is superlinear but slower than the quadratic convergence of Newton-Raphson.

| Method | Order | Evaluations per step | Effective rate |
|:-----|:-----|:-----|:-----|
| Newton-Raphson | 2 (quadratic) | $`f`$ and $`f'`$ | High |
| Secant | 1.618 | $`f`$ only | Medium |
| Bisection | 1 (linear) | $`f`$ only | Slow |

*Table 2. Convergence of Newton-Raphson, secant, and bisection compared*

Fill in the blank of the following function.

```python
def secant(f, x0, x1, es=1e-7, maxit=30):
    """
    Find a root of f using the Secant method.

    Parameters
    ----------
    f     : callable — function (no derivative needed)
    x0, x1: float   — two initial guesses (do NOT need to bracket the root)
    es    : float   — stopping criterion
    maxit : int     — maximum iterations

    Returns
    -------
    x2, f(x2), ea, n_iter, history
    """
    history = []

    for i in range(maxit):
        f0, f1 = f(x0), f(x1)

        if abs(f1 - f0) < 1e-30:
            raise ZeroDivisionError('f(x1) ≈ f(x0) — secant method fails (zero denominator)')

        # --- COMPLETE: secant update formula ---
        x2 = 
        # --- COMPLETE ---

        ea = abs((x2 - x1) / x2) if x2 != 0 else abs(x2 - x1)
        history.append((i + 1, x2, f(x2), ea))

        if ea < es:
            break

        x0, x1 = x1, x2

    return x2, f(x2), ea, i + 1, history
```

The next cell compares Newton-Raphson and the secant method on the same test function. The plotting part is omitted here.

```python
# Compare Newton-Raphson vs Secant on f(x) = x - e^{-x}
xr_sc, fxr_sc, ea_sc, n_sc, hist_sc = secant(f_nr, x0=0.0, x1=1.0, es=1e-12, maxit=30)

print('=== f(x) = x - e^{-x},  true root = 0.5671432904097838 ===')
print(f'Newton-Raphson:  {n_nr} iters,  |error| = {abs(xr_nr - x_true):.3e}')
print(f'Secant:          {n_sc} iters,  |error| = {abs(xr_sc - x_true):.3e}')
```

The distributed notebook shows the following output for this cell.

```text
=== f(x) = x - e^{-x},  true root = 0.5671432904097838 ===
Newton-Raphson:  5 iters,  |error| = 1.110e-16
Secant:          6 iters,  |error| = 0.000e+00
```

---

<br>

After solving the problems, see the [assignment explanation](../Assignment-Explanation.md) for the answers and solutions and the [completed notebook](../Solution/Lab07-RootSolving_Completed.ipynb) for the full code.
