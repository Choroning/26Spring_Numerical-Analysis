# Chapter 3 Lab — Sources of Numerical Errors

> **Last Updated:** 2026-03-26

---

<br>

## Table of Contents

- [1. Error Definitions and Metrics](#1-error-definitions-and-metrics)
  - [1.1 Scarborough Criterion Table](#11-scarborough-criterion-table)
  - [1.2 sqrt(2) Approximation — Babylonian/Heron's Method](#12-sqrt2-approximation--babylonianherons-method)
- [2. Truncation Error and Taylor Series](#2-truncation-error-and-taylor-series)
  - [2.1 Maclaurin Series for e^x — Term-by-Term](#21-maclaurin-series-for-ex--term-by-term)
  - [2.2 Iterative Computation with Stopping Criterion](#22-iterative-computation-with-stopping-criterion)
  - [2.3 Visualization of Truncation Error Convergence](#23-visualization-of-truncation-error-convergence)
- [3. Round-Off Error and Floating-Point](#3-round-off-error-and-floating-point)
  - [3.1 Machine Epsilon](#31-machine-epsilon)
  - [3.2 Binary Representation of Decimals](#32-binary-representation-of-decimals)
  - [3.3 Accumulation of Round-Off Error](#33-accumulation-of-round-off-error)
  - [3.4 Kahan Summation Algorithm](#34-kahan-summation-algorithm)
  - [3.5 Catastrophic Cancellation — Quadratic Formula](#35-catastrophic-cancellation--quadratic-formula)
- [Summary](#summary)

---

<br>

## 1. Error Definitions and Metrics

### 1.1 Scarborough Criterion Table

The **Scarborough criterion** provides a direct mapping from the desired number of significant figures to the required approximate percent relative error. The criterion states that if $|\varepsilon_a| < (0.5 \times 10^{2-n})\%$, the result is correct to at least $n$ significant figures.

We can compute the threshold for each value of $n$ using a simple lambda function:

```python
import numpy as np
import matplotlib.pyplot as plt
import math

# Scarborough criterion: minimum epsilon_a for n significant figures
scarborough = lambda n: 0.5 * 10**(2 - n)

# Table of criteria
for n in range(1, 8):
    print(f'{n} sig. figs → |ε_a| < {scarborough(n):.4e}%')
# 1 → 5.0000e+00%, 2 → 5.0000e-01%, ..., 7 → 5.0000e-05%
```

This table is a practical reference: before starting an iterative computation, decide how many significant figures you need, look up the corresponding threshold, and use it as your stopping criterion $\varepsilon_s$.

### 1.2 sqrt(2) Approximation — Babylonian/Heron's Method

The **Babylonian method** (also known as **Heron's method**) is one of the oldest iterative algorithms, dating back to ancient Mesopotamia. It computes $\sqrt{S}$ using the recurrence:

$$x_{n+1} = \frac{1}{2}\left(x_n + \frac{S}{x_n}\right)$$

The intuition is elegant: if $x_n$ is an overestimate of $\sqrt{S}$, then $S/x_n$ is an underestimate (and vice versa). Averaging the two gives a better estimate that is closer to $\sqrt{S}$ than either value alone.

> **[Linear Algebra]** The Babylonian method is actually a special case of Newton's method applied to $f(x) = x^2 - S$. Newton's method: $x_{n+1} = x_n - f(x_n)/f'(x_n) = x_n - (x_n^2 - S)/(2x_n) = (x_n + S/x_n)/2$. This explains its quadratic convergence (errors square each iteration).

The following code demonstrates how quickly the method converges, tracking both the true percent relative error $\varepsilon_t$ and the approximate percent relative error $\varepsilon_a$:

```python
target = 2
x = 1.0  # initial guess
true_val = math.sqrt(target)

print(f'{"Iter":>4} {"x":>20} {"ε_t (%)":>12} {"ε_a (%)":>12}')
print('-' * 52)

for i in range(6):
    x_old = x
    x = 0.5 * (x + target / x)
    et = abs(true_val - x) / true_val * 100
    ea = abs(x - x_old) / abs(x) * 100 if i > 0 else float('inf')
    print(f'{i+1:4d} {x:20.15f} {et:12.6e} {ea:12.6e}')
```

Observe the **quadratic convergence**: the number of correct digits roughly doubles with each iteration. By iteration 5 or 6, the true error is at or near machine epsilon — the algorithm has converged to the best possible double-precision representation of $\sqrt{2}$.

This demonstrates a key principle: for well-designed iterative methods, $\varepsilon_a$ closely tracks $\varepsilon_t$, so we can use $\varepsilon_a$ as a reliable proxy even when the true value is unknown.

---

<br>

## 2. Truncation Error and Taylor Series

### 2.1 Maclaurin Series for e^x — Term-by-Term

The **Maclaurin series** (Taylor series expanded about $x = 0$) for $e^x$ is:

$$e^x = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \cdots = \sum_{k=0}^{\infty} \frac{x^k}{k!}$$

This is a convergent series for all real $x$, meaning we can approximate $e^x$ to any desired accuracy by including enough terms. Each additional term reduces the truncation error.

The following code builds the approximation term by term for $e^{0.5}$ and shows how the error decreases with each additional term:

```python
# Term-by-term approximation of e^0.5
x = 0.5
true_val = math.exp(x)

f_approx = np.zeros(6)
f_approx[0] = 1.0                           # 0th order: 1
f_approx[1] = f_approx[0] + x              # 1st: + x
f_approx[2] = f_approx[1] + x**2/2         # 2nd: + x²/2!
f_approx[3] = f_approx[2] + x**3/6         # 3rd: + x³/3!
f_approx[4] = f_approx[3] + x**4/24        # 4th: + x⁴/4!
f_approx[5] = f_approx[4] + x**5/120       # 5th: + x⁵/5!

for k in range(6):
    et = abs(true_val - f_approx[k]) / true_val * 100
    print(f'Order {k}: {f_approx[k]:.10f}, ε_t = {et:.6e}%')
```

The true value of $e^{0.5} \approx 1.6487212707$. As we include more terms, the approximation converges rapidly:

- 0th order (just $1$): error $\approx 39.3\%$
- 1st order ($1 + x = 1.5$): error $\approx 9.02\%$
- 5th order: error is negligibly small

The factorial in the denominator grows much faster than the power of $x$ in the numerator, which is why the series converges so quickly for moderate values of $x$.

### 2.2 Iterative Computation with Stopping Criterion

In practice, we do not hard-code each term. Instead, we write a loop that adds terms iteratively and stops when the approximate error drops below a threshold:

```python
def exp_series(x, es=1e-4, maxit=50):
    """Compute e^x via Taylor series with stopping criterion.

    Args:
        x: input value
        es: stopping criterion (% approximate error)
        maxit: maximum iterations
    Returns:
        (result, final_error_%, iterations)
    """
    sol = 1.0
    ea = 100.0

    for k in range(1, maxit + 1):
        sol_old = sol
        sol = sol + x**k / math.factorial(k)
        ea = abs(sol - sol_old) / abs(sol) * 100
        if ea < es:
            break

    return sol, ea, k

result, error, iters = exp_series(0.5)
print(f'e^0.5 ≈ {result:.15f} (after {iters} iterations, ε_a = {error:.2e}%)')
print(f'True:    {math.exp(0.5):.15f}')
```

Key observations about this implementation:

- The stopping criterion `es=1e-4` corresponds to $\varepsilon_s = 0.0001\%$. By the Scarborough criterion, this guarantees at least 6 significant figures ($0.5 \times 10^{2-6} = 0.00005\%$, and $0.0001\% < 0.0005\%$ which gives 5 sig figs).
- The `maxit` parameter provides a safety net: if the series fails to converge (e.g., for very large $|x|$), the loop terminates after a maximum number of iterations rather than running forever.
- The function returns all three pieces of information the caller needs: the result, the final error estimate, and the number of iterations used.

### 2.3 Visualization of Truncation Error Convergence

A semilog plot reveals the exponential rate at which the truncation error decreases as more terms are added:

```python
orders = range(1, 15)
errors = [abs(math.exp(0.5) - sum(0.5**k/math.factorial(k) for k in range(n+1))) for n in orders]

plt.semilogy(list(orders), errors, 'bo-')
plt.xlabel('Number of Terms')
plt.ylabel('|Truncation Error|')
plt.title('Taylor Series Convergence for e^0.5')
plt.grid(True)
plt.show()
```

The **straight line on the semilog plot** confirms exponential convergence: each additional term reduces the error by a roughly constant factor. This is a characteristic of Taylor series for well-behaved (analytic) functions.

For $e^{0.5}$, the error drops below machine epsilon ($\approx 2.22 \times 10^{-16}$) after approximately 10-12 terms. Adding more terms beyond this point provides no further improvement because round-off error dominates — this is the point where truncation error and round-off error meet.

---

<br>

## 3. Round-Off Error and Floating-Point

### 3.1 Machine Epsilon

**Machine epsilon** is the smallest floating-point number $\varepsilon$ such that $1 + \varepsilon \neq 1$ in the computer's arithmetic. It represents the fundamental precision limit of the floating-point format.

```python
# Machine epsilon
eps = np.finfo(float).eps
print(f'Machine epsilon: {eps}')  # 2.220446049250313e-16

# Demonstration: 1 + eps/2 == 1 is True!
print(1 + eps/2 == 1)   # True (eps/2 is below representable precision)
print(1 + eps == 1)      # False (eps is exactly at the boundary)
```

The first test (`1 + eps/2 == 1`) returns `True` because `eps/2` is too small to affect the representation of 1 in 64-bit floating point. The number $1 + \varepsilon/2$ cannot be distinguished from $1$ — there is no 64-bit floating-point number between $1$ and $1 + \varepsilon$.

The second test (`1 + eps == 1`) returns `False` because machine epsilon is precisely the smallest increment that produces a distinguishable result. The number $1 + \varepsilon$ is the next representable floating-point number after $1$.

> **Note:** Machine epsilon is a **relative** measure. Near $1$, the absolute precision is $\varepsilon \approx 2.22 \times 10^{-16}$. Near $10^{10}$, the absolute precision is $10^{10} \times \varepsilon \approx 2.22 \times 10^{-6}$. The relative precision is always the same.

### 3.2 Binary Representation of Decimals

Many decimal fractions that appear "simple" to humans (like $0.1$ or $0.3$) have **non-terminating representations in binary**, just as $1/3 = 0.333...$ is non-terminating in decimal. This is a fundamental source of round-off error in floating-point arithmetic.

```python
def decimal_to_binary(num, bits=52):
    """Show binary representation of a decimal fraction."""
    result = ''
    for _ in range(bits):
        num *= 2
        if num >= 1:
            result += '1'
            num -= 1
        else:
            result += '0'
    return '0.' + result

print(f'0.1 in binary: {decimal_to_binary(0.1, 20)}')  # non-terminating!
print(f'0.3 in binary: {decimal_to_binary(0.3, 20)}')  # also non-terminating!
# This is why 0.1 + 0.2 != 0.3 in floating-point!
```

The algorithm works by repeatedly doubling the fraction. If the result exceeds 1, the next binary digit is 1 (and we subtract 1); otherwise it is 0. For $0.1$:

$$0.1 \times 2 = 0.2 \to 0, \quad 0.2 \times 2 = 0.4 \to 0, \quad 0.4 \times 2 = 0.8 \to 0, \quad 0.8 \times 2 = 1.6 \to 1, \quad \ldots$$

The pattern $0.0\overline{0011}$ repeats infinitely. Since the computer can only store 52 bits of mantissa, the representation is truncated, introducing a small error. This is why the famous `0.1 + 0.2 != 0.3` surprise occurs in virtually every programming language that uses IEEE 754 floating point.

### 3.3 Accumulation of Round-Off Error

When many small round-off errors accumulate over thousands of operations, the total error can become significant. This is demonstrated by summing $0.0001$ ten thousand times — the mathematical result is exactly $1.0$, but floating-point arithmetic introduces a discrepancy:

```python
# Summing 0.0001 ten thousand times
total_naive = sum(0.0001 for _ in range(10000))
print(f'Naive sum: {total_naive}')
print(f'Error: {abs(1.0 - total_naive):.2e}')  # ~9.38e-14
```

Each addition of $0.0001$ introduces a tiny round-off error (because $0.0001$ cannot be represented exactly in binary). Over 10,000 additions, these errors accumulate to approximately $10^{-13}$ — small in absolute terms, but potentially problematic in applications that require high accuracy or that use this sum in further calculations (where the error may be amplified).

### 3.4 Kahan Summation Algorithm

The **Kahan summation algorithm** (also called **compensated summation**) is an elegant technique that dramatically reduces the accumulation of round-off errors in large summations. The key idea is to maintain a separate "compensation" variable that tracks the running round-off error and corrects for it in the next addition:

```python
def kahan_sum(values):
    """Kahan summation for reduced round-off error."""
    total = 0.0
    compensation = 0.0
    for val in values:
        y = val - compensation          # compensate for lost low-order bits
        temp = total + y
        compensation = (temp - total) - y  # recover what was lost
        total = temp
    return total

total_kahan = kahan_sum(0.0001 for _ in range(10000))
print(f'Kahan sum: {total_kahan}')
print(f'Error: {abs(1.0 - total_kahan):.2e}')  # 0.00e+00 — perfect!
```

Here is how the algorithm works step by step:

1. `y = val - compensation`: Adjust the next value to add by subtracting the accumulated compensation (the error from previous steps).
2. `temp = total + y`: Perform the addition (this introduces a new round-off error).
3. `compensation = (temp - total) - y`: Compute the round-off error introduced in step 2. Since `temp = total + y` in exact arithmetic, `(temp - total) - y` should be zero — but in floating point, it captures exactly the round-off error that was lost.
4. `total = temp`: Update the running total.

> **Note:** Kahan summation adds a "compensation" variable that tracks the accumulated round-off error and corrects for it in the next addition. This is one of the most important numerical techniques for maintaining accuracy in large summations.

The result is remarkable: while naive summation produces an error of $\sim 10^{-13}$, Kahan summation achieves **zero error** for this test case. The compensation mechanism effectively gives us double the precision of a single floating-point variable.

### 3.5 Catastrophic Cancellation — Quadratic Formula

The standard quadratic formula for $ax^2 + bx + c = 0$:

$$x = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}$$

suffers from **catastrophic cancellation** when $b^2 \gg 4ac$. In this case, $\sqrt{b^2 - 4ac} \approx |b|$, so one of the two roots involves subtracting two nearly equal numbers.

Specifically, when $b > 0$, the root $x_1 = \frac{-b + \sqrt{b^2 - 4ac}}{2a}$ involves subtracting $b$ from $\sqrt{b^2 - 4ac}$ — two numbers that are nearly identical when $b$ is large.

The stable alternative uses **Vieta's formulas**: for $ax^2 + bx + c = 0$, the product of roots is $x_1 \cdot x_2 = c/a$. So if we can compute one root accurately, we get the other as $x_2 = c/(a \cdot x_1)$ without any subtraction.

```python
# Standard quadratic formula for ax² + bx + c = 0:
# x = (-b ± sqrt(b²-4ac)) / (2a)
# When b >> sqrt(b²-4ac), one root suffers catastrophic cancellation

a, c = 1, 1

print(f'{"b":>12} {"x2_standard":>20} {"x2_stable":>20} {"rel_error":>15}')
print('-' * 70)

for exp in range(2, 11):
    b = 10**exp
    disc = math.sqrt(b**2 - 4*a*c)

    # Standard formula (cancellation in x2)
    x1_std = (-b + disc) / (2*a)
    x2_std = (-b - disc) / (2*a)

    # Stable alternative: x2 = c / (a * x1)
    x2_stable = c / (a * x1_std)

    rel_err = abs(x2_std - x2_stable) / abs(x2_stable) * 100
    print(f'{b:12.0e} {x2_std:20.10f} {x2_stable:20.10f} {rel_err:15.6e}%')
```

> **[Calculus]** The stable formula uses Vieta's formulas: for $ax^2 + bx + c = 0$, the product of roots is $x_1 \cdot x_2 = c/a$. So if we can compute $x_1$ accurately, we get $x_2 = c/(a \cdot x_1)$ without subtraction.

As $b$ increases, the standard formula loses more and more significant digits in $x_1$ (the root computed with the $-b + \sqrt{b^2 - 4ac}$ subtraction), while the stable alternative maintains full precision. For $b = 10^{10}$, the standard formula may have lost nearly all significant digits, while Vieta's approach remains accurate.

The lesson is general: whenever you encounter a subtraction of nearly equal quantities in a formula, look for an algebraically equivalent form that avoids the subtraction. Common techniques include:

- **Rationalizing**: multiply by the conjugate (e.g., $\frac{(\sqrt{a} - \sqrt{b})(\sqrt{a} + \sqrt{b})}{\sqrt{a} + \sqrt{b}} = \frac{a - b}{\sqrt{a} + \sqrt{b}}$)
- **Vieta's formulas**: use the product-of-roots relationship
- **Taylor expansion**: when terms nearly cancel, expand to find the leading non-cancelling term

---

<br>

## Summary

| Topic | Key Concept |
|:------|:-----------|
| Scarborough Criterion | $|\varepsilon_a| < 0.5 \times 10^{2-n}\%$ → n significant figures |
| Babylonian Method | $x_{n+1} = (x_n + S/x_n)/2$ — quadratic convergence |
| Taylor Series | $e^x = \sum x^k/k!$ — iterative approximation with stopping criterion |
| Machine Epsilon | `np.finfo(float).eps` $\approx 2.22 \times 10^{-16}$ |
| Kahan Summation | Compensation technique for accurate large sums |
| Catastrophic Cancellation | Use algebraic alternatives to avoid subtracting nearly-equal numbers |

---
