# HW12 Assignment Handout — Numerical Differentiation

[Repository Guide](../../../README.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **Last Updated:** 2026-11-24

This is the Week 14 homework of Numerical Analysis. The task is to fill in the blank of the lab notebook `Lab13-diff.ipynb` distributed by the instructor and run it. The notebook has eight parts: the benchmark function and a review of Taylor expansion, the forward, backward, and central differences, higher-order stencils and Richardson extrapolation, second derivatives and mixed partial derivatives, step-size selection and the complex-step method, multivariate differentiation (Jacobians and Hessians), FEM derivative recovery, and engineering applications. All other code cells, including the cell of Section 7 marked `COMPLETE` and the Exercise cells, were distributed with their code filled in, and the only blank lies between the two comment lines marked `COMPLETE` in Section 2.1. The first code cell of the notebook imports `numpy` as `np` and `matplotlib.pyplot` as `plt`, imports `solve` from `scipy.linalg`, and turns off warning output.

Each code block below reproduces a cell of the notebook, and the place where code is missing is the blank to fill. The distributed blank notebook is [Lab13-diff.ipynb](Lab13-diff.ipynb).

---

<br>

## Problem. Central Difference Function (Notebook Section 2.1)

Section 1 of the notebook uses the same polynomial as the textbook (Chapra and Canale) as the benchmark function.

$$
f(x) = 0.1x^5 - 0.2x^4 + 3.0x^3 + 2.0x^2 - x + 0.5
$$

Its analytic derivatives are as follows.

$$
f'(x) = 0.5x^4 - 0.8x^3 + 9.0x^2 + 4.0x - 1
$$

$$
f''(x) = 2x^3 - 2.4x^2 + 18x + 4
$$

$$
f'''(x) = 6x^2 - 4.8x + 18
$$

Every finite difference formula comes from Taylor's theorem.

$$
f(x+h) = f(x) + h f'(x) + \frac{h^2}{2!} f''(x) + \frac{h^3}{3!} f'''(x) + \cdots
$$

$$
f(x-h) = f(x) - h f'(x) + \frac{h^2}{2!} f''(x) - \frac{h^3}{3!} f'''(x) + \cdots
$$

By rearranging these expansions, one can isolate $`f'(x)`$, $`f''(x)`$, and so on with controlled truncation error. The next cell sets the benchmark function and its derivatives, prints their values at $`x_0 = 0.5`$, and draws the plots.

```python
# Benchmark function and its analytic derivatives
f    = lambda x: 0.1*x**5 - 0.2*x**4 + 3.0*x**3 + 2.0*x**2 - x + 0.5
df   = lambda x: 0.5*x**4 - 0.8*x**3 + 9.0*x**2 + 4.0*x - 1.0
d2f  = lambda x: 2.0*x**3 - 2.4*x**2 + 18.0*x + 4.0
d3f  = lambda x: 6.0*x**2 - 4.8*x + 18.0

# Evaluation point and step size
x0 = 0.5
print(f'Benchmark at x = {x0}:')
print(f'  f({x0})    = {f(x0):.8f}')
print(f'  f\'({x0})   = {df(x0):.8f}')
print(f'  f\'\'({x0})  = {d2f(x0):.8f}')
print(f'  f\'\'\'({x0}) = {d3f(x0):.8f}')

# Visualize f and its derivatives
x_plot = np.linspace(-1, 2, 400)
fig, axes = plt.subplots(1, 2, figsize=(13, 4))

ax = axes[0]
ax.plot(x_plot, f(x_plot),   'b-',  lw=2, label='$f(x)$')
ax.plot(x_plot, df(x_plot),  'r-',  lw=2, label="$f'(x)$")
ax.plot(x_plot, d2f(x_plot), 'g-',  lw=2, label="$f''(x)$")
ax.axvline(x0, color='k', ls=':', lw=1.5, label=f'$x_0={x0}$')
ax.set(xlabel='x', ylabel='value', title='Benchmark Function and Derivatives')
ax.legend(); ax.grid(True, alpha=0.3)

# Taylor expansion visualization at x0
ax = axes[1]
h_demo = 0.4
x_expand = np.linspace(x0 - h_demo - 0.1, x0 + h_demo + 0.1, 200)
T1 = f(x0) + df(x0)*(x_expand - x0)
T2 = T1 + d2f(x0)/2*(x_expand - x0)**2
T3 = T2 + d3f(x0)/6*(x_expand - x0)**3
ax.plot(x_expand, f(x_expand),  'b-',  lw=2.5, label='$f(x)$')
ax.plot(x_expand, T1,           'r--', lw=1.5, label='Taylor $O(h)$')
ax.plot(x_expand, T2,           'g--', lw=1.5, label='Taylor $O(h^2)$')
ax.plot(x_expand, T3,           'm--', lw=1.5, label='Taylor $O(h^3)$')
ax.scatter([x0], [f(x0)], c='k', s=80, zorder=6)
ax.set(xlabel='x', ylabel='value',
       title=f'Taylor Approximation around $x_0={x0}$')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

plt.tight_layout(); plt.show()
```

```text
Benchmark at x = 0.5:
  f(0.5)    = 0.86562500
  f'(0.5)   = 3.18125000
  f''(0.5)  = 12.65000000
  f'''(0.5) = 17.10000000
```

Section 2.1 of the notebook solves the Taylor expansion $`f(x+h) = f(x) + hf'(x) + \frac{h^2}{2}f''(\xi)`$ for $`f'(x)`$ to obtain the three formulas below.

| Name | Formula | Truncation error | Uses |
|:-----|:-----|:-----|:-----|
| Forward | $`f'(x) \approx \dfrac{f(x+h) - f(x)}{h}`$ | $`O(h)`$, first order | $`f(x)`$, $`f(x+h)`$ |
| Backward | $`f'(x) \approx \dfrac{f(x) - f(x-h)}{h}`$ | $`O(h)`$, first order | $`f(x-h)`$, $`f(x)`$ |
| Central | $`f'(x) \approx \dfrac{f(x+h) - f(x-h)}{2h}`$ | $`O(h^2)`$, second order | $`f(x-h)`$, $`f(x+h)`$ |

*Table 1. Three basic finite difference formulas for the first derivative*

The central difference is second order for the following reason. Subtracting the two Taylor expansions cancels the $`h^2`$ term.

$$
f(x+h) - f(x-h) = 2h f'(x) + \frac{2h^3}{6} f'''(x) + \cdots
$$

$$
\Rightarrow f'(x) = \frac{f(x+h)-f(x-h)}{2h} - \frac{h^2}{6}f'''(x) + O(h^4)
$$

The leading error is $`O(h^2)`$, compared to $`O(h)`$ for the forward and backward differences. Fill in the blank of the next cell to complete the function `central_diff`. The functions `forward_diff` and `backward_diff` in the same cell were distributed filled in.

```python
def forward_diff(func, x, h):
    """
    First-order forward difference approximation of f'(x).
    Truncation error: O(h).

    Parameters
    ----------
    func : callable — function f
    x    : float   — point of differentiation
    h    : float   — step size

    Returns
    -------
    float : approximation of f'(x)
    """
    fdiff = (func(x + h) - func(x)) / h
    return fdiff 


def backward_diff(func, x, h):
    """
    First-order backward difference approximation of f'(x).
    Truncation error: O(h).
    """
    bdiff = (func(x) - func(x - h)) / h
    return bdiff


def central_diff(func, x, h):
    """
    Second-order central difference approximation of f'(x).
    Truncation error: O(h^2).
    """
    # - - - - - COMPLETE - - - - - 
    cdiff = 
    # - - - - - COMPLETE - - - - - 
    return cdiff
    


# Compare at x0 = 0.5
exact = df(x0)
print(f'Analytic f\'({x0}) = {exact:.10f}\n')
print(f'{"h":>8}  {"Forward":>14}  {"err_F":>10}  {"Backward":>14}  {"err_B":>10}  {"Central":>14}  {"err_C":>10}')
print('-' * 95)
for h in [0.5, 0.25, 0.1, 0.05, 0.01, 0.001]:
    fwd  = forward_diff(f, x0, h)
    bwd  = backward_diff(f, x0, h)
    cnt  = central_diff(f, x0, h)
    print(f'{h:>8.4f}  {fwd:>14.8f}  {abs(fwd-exact):>10.2e}  '
          f'{bwd:>14.8f}  {abs(bwd-exact):>10.2e}  '
          f'{cnt:>14.8f}  {abs(cnt-exact):>10.2e}')
```

The distributed notebook shows the following output for this cell.

```text
Analytic f'(0.5) = 3.1812500000

       h         Forward       err_F        Backward       err_B         Central       err_C
-----------------------------------------------------------------------------------------------
  0.5000      7.06875000    3.89e+00      0.73125000    2.45e+00      3.90000000    7.19e-01
  0.2500      4.94179687    1.76e+00      1.77773437    1.40e+00      3.35976563    1.79e-01
  0.1000      3.84231000    6.61e-01      2.57721000    6.04e-01      3.20976000    2.85e-02
  0.0500      3.50463188    3.23e-01      2.87211937    3.09e-01      3.18837563    7.13e-03
  0.0100      3.24478505    6.35e-02      3.11828495    6.30e-02      3.18153500    2.85e-04
  0.0010      3.18757785    6.33e-03      3.17492785    6.32e-03      3.18125285    2.85e-06
```

---

<br>

After solving the problem, see the [assignment explanation](../Assignment-Explanation.md) for the answer and solution and the [completed notebook](../Solution/Lab13-diff_Completed.ipynb) for the full code.
