# HW11 Assignment Handout — Numerical Integration

[Repository Guide](../../../README.md) | [English](Assignment.md) | [한국어](Assignment.ko.md)

> **Last Updated:** 2026-11-24

This is the Week 13 homework of Numerical Analysis. The task is to fill in the blank of the lab notebook `Lab12-Int.ipynb` distributed by the instructor and run it. The notebook has eight parts: the benchmark function, the trapezoidal rule, Simpson's rules, Richardson extrapolation and Romberg integration, Gauss-Legendre quadrature, Gauss-Lobatto quadrature, `scipy.integrate`, and engineering applications. All other code cells, including the Exercise cells, were distributed with their code filled in, and the only blank lies between the two comment lines marked `COMPLETE`. The first code cell of the notebook imports `numpy` as `np` and `matplotlib.pyplot` as `plt`, as well as `integrate` from `scipy` and `leggauss` from `numpy.polynomial.legendre`.

Each code block below reproduces a cell of the notebook, and the place where code is missing is the blank to fill. The distributed blank notebook is [Lab12-Int.ipynb](Lab12-Int.ipynb).

---

<br>

## Problem. Composite Trapezoidal Rule (Notebook Section 2.2)

Section 1 of the notebook uses a single **benchmark function** whose exact integral is known to compare all the methods.

```math
f(x) = 0.2 + 25x - 200x^2 + 675x^3 - 900x^4 + 400x^5, \quad x \in [0,\, 0.8]
```

```math
I_{\text{exact}} = \int_0^{0.8} f(x)\,dx = 1.640533\ldots
```

The percent true relative error is defined as follows.

```math
\varepsilon_t = \left|\frac{I_{\text{exact}} - I_{\text{approx}}}{I_{\text{exact}}}\right| \times 100\,\%
```

The next cell sets the benchmark function, the exact integral, and the error function, and plots the function.

```python
# Benchmark function and exact integral
f   = lambda x: 0.2 + 25*x - 200*x**2 + 675*x**3 - 900*x**4 + 400*x**5
a, b  = 0.0, 0.8
# Exact integral, computed in closed form:
#   F(x) = 0.2x + 12.5 x^2 - (200/3) x^3 + (675/4) x^4 - 180 x^5 + (200/3) x^6
#   F(0.8) - F(0) = 49216/30000 = 3076/1875 = 1.6405333333...
I_exact = 49216 / 30000     # = 1.6405333333... (full double precision)

def pct_err(I_approx):
    """Percent true relative error."""
    return abs(I_exact - I_approx) / abs(I_exact) * 100

# Visualise
x_plot = np.linspace(a, b, 300)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(x_plot, f(x_plot), 'b-', lw=2)
ax.fill_between(x_plot, f(x_plot), alpha=0.15, color='blue')
ax.set(xlabel='x', ylabel='f(x)', title='Benchmark Function  (shaded area = exact integral)')
ax.text(0.5, 1.5, f'$I_{{\\rm exact}} = {I_exact:.10f}$', fontsize=13)
ax.grid(True, alpha=0.3)
plt.tight_layout(); plt.show()
```

The single-segment trapezoidal rule of Section 2.1 approximates $`f`$ by a linear polynomial through $`(a, f(a))`$ and $`(b, f(b))`$.

```math
I \approx (b-a)\,\frac{f(a)+f(b)}{2}
```

The formula is exact for polynomials of degree at most 1, and its truncation error is as follows.

$$
E_t = -\frac{1}{12}f''(\xi)(b-a)^3, \quad \xi \in [a,b]
$$

The composite trapezoidal rule of Section 2.2 divides $`[a,b]`$ into $`n`$ equal sub-intervals of width $`h = (b-a)/n`$.

$$
I \approx \frac{h}{2}\left[f(x_0) + 2\sum_{i=1}^{n-1}f(x_i) + f(x_n)\right]
$$

The global error is proportional to $`h^2 \propto n^{-2}`$, so doubling $`n`$ reduces the error by a factor of about 4. Fill in the blank of the next cell to complete the function `composite_trap`.

```python
# ---- Single-segment trapezoidal rule ----
I_trap1 = (b - a) * (f(a) + f(b)) / 2.0
print(f'Single-segment Trapezoidal: I = {I_trap1:.6f},  error = {pct_err(I_trap1):.2f}%')

# ---- Composite trapezoidal rule ----
def composite_trap(func, a, b, n):
    """
    Composite trapezoidal rule.

    Parameters
    ----------
    func : callable
    a, b : integration limits
    n    : number of sub-intervals

    Returns
    -------
    I : float
    """
    x = np.linspace(a, b, n + 1)
    y = func(x)
    h = (b - a) / n
    # - - - - - - COMPLETE - - - - - - - 
    I = h/2 * 
    # - - - - - - COMPLETE - - - - - - - 
    return I

# Test with increasing n
print('\nComposite Trapezoidal Rule:')
print(f'{"n":>6}  {"I":>12}  {"error %":>10}')
for n in [1, 2, 4, 8, 16, 32, 64]:
    I = composite_trap(f, a, b, n)
    print(f'{n:>6}  {I:>12.8f}  {pct_err(I):>10.4f}')
```

The distributed notebook shows the following output for this cell.

```text
Single-segment Trapezoidal: I = 0.172800,  error = 89.47%

Composite Trapezoidal Rule:
     n             I     error %
     1    0.17280000     89.4668
     2    1.06880000     34.8505
     4    1.48480000      9.4928
     8    1.60080000      2.4220
    16    1.63055000      0.6085
    32    1.63803438      0.1523
    64    1.63990840      0.0381
```

---

<br>

After solving the problem, see the [assignment explanation](../Assignment-Explanation.md) for the answer and solution and the [completed notebook](../Solution/Lab12-Int_Completed.ipynb) for the full code.
