# Chapter 3 Lab — Sources of Numerical Errors

[Repository Guide](../../README.md) | [English](Concepts_Lab.md) | [한국어](Concepts_Lab.ko.md)

> **Last Updated:** 2026-10-27
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 4

> **Learning Objectives:**
> 1. Choose a stopping tolerance with the Scarborough criterion and apply it to iterative computations
> 2. Verify by computation how the truncation error of a Taylor series shrinks with the number of terms
> 3. Analyze roundoff errors in floating-point arithmetic using machine epsilon and the binary representation of decimals
> 4. Reduce accumulated roundoff error and catastrophic cancellation with Kahan summation and algebraic reformulation

---

## Table of Contents

- [1. Error Definitions and Metrics](#1-error-definitions-and-metrics)
  - [1.1 Scarborough Criterion Table](#11-scarborough-criterion-table)
  - [1.2 sqrt(2) Approximation: Babylonian/Heron's Method](#12-sqrt2-approximation-babylonianherons-method)
- [2. Truncation Error and Taylor Series](#2-truncation-error-and-taylor-series)
  - [2.1 Maclaurin Series for e^x: Term-by-Term](#21-maclaurin-series-for-ex-term-by-term)
  - [2.2 Iterative Computation with Stopping Criterion](#22-iterative-computation-with-stopping-criterion)
  - [2.3 Visualization of Truncation Error Convergence](#23-visualization-of-truncation-error-convergence)
- [3. Roundoff Error and Floating-Point](#3-roundoff-error-and-floating-point)
  - [3.1 Machine Epsilon](#31-machine-epsilon)
  - [3.2 Binary Representation of Decimals](#32-binary-representation-of-decimals)
  - [3.3 Accumulation of Roundoff Error](#33-accumulation-of-roundoff-error)
  - [3.4 Kahan Summation Algorithm](#34-kahan-summation-algorithm)
  - [3.5 Catastrophic Cancellation: Quadratic Formula](#35-catastrophic-cancellation-quadratic-formula)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Error Definitions and Metrics

### 1.1 Scarborough Criterion Table

The **Scarborough criterion** provides a direct mapping from the desired number of significant figures to the required approximate percent relative error. The criterion states that if $`|\varepsilon_a| < (0.5 \times 10^{2-n})\%`$, the result is correct to at least $`n`$ significant figures.

We can compute the threshold for each value of $`n`$ using a simple lambda function.

```python
import numpy as np
import matplotlib.pyplot as plt
import math

# Scarborough criterion: maximum epsilon_a for n significant figures
scarborough = lambda n: 0.5 * 10**(2 - n)

# Table of criteria
for n in range(1, 8):
    print(f'{n} sig. figs → |ε_a| < {scarborough(n):.4e}%')
# 1 → 5.0000e+00%, 2 → 5.0000e-01%, ..., 7 → 5.0000e-06%
```

This table is a practical reference: before starting an iterative computation, decide how many significant figures you need, look up the corresponding threshold, and use it as your stopping criterion $`\varepsilon_s`$.

### 1.2 sqrt(2) Approximation: Babylonian/Heron's Method

The **Babylonian method** is one of the oldest iterative algorithms, dating back to ancient Mesopotamia, and is also known as **Heron's method**. It computes $`\sqrt{S}`$ using the following recurrence.

```math
x_{n+1} = \frac{1}{2}\left(x_n + \frac{S}{x_n}\right)
```

The intuition is simple: if $`x_n`$ is an overestimate of $`\sqrt{S}`$, then $`S/x_n`$ is an underestimate, and vice versa. Averaging the two gives an estimate that is closer to $`\sqrt{S}`$ than either value alone.

> **Note:** **Relation to the Newton-Raphson method.** The Babylonian method is a special case of the Newton-Raphson method (Chapter 6) applied to $`f(x) = x^2 - S`$. The Newton-Raphson iteration gives $`x_{n+1} = x_n - f(x_n)/f'(x_n) = x_n - (x_n^2 - S)/(2x_n) = (x_n + S/x_n)/2`$. This relation explains its quadratic convergence (the error at each iteration is roughly proportional to the square of the previous error).

The following code demonstrates how quickly the method converges, tracking both the true percent relative error $`\varepsilon_t`$ and the approximate percent relative error $`\varepsilon_a`$.

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

The output shows **quadratic convergence**: the number of correct digits roughly doubles with each iteration. By iteration 5, the true error is $`\varepsilon_t \approx 1.57 \times 10^{-14}\%`$ (a relative error of about $`1.6 \times 10^{-16}`$), which is at the level of machine epsilon, and the algorithm has converged to the best double-precision representation of $`\sqrt{2}`$.

This demonstrates a key principle: for well-designed iterative methods, $`\varepsilon_a`$ closely tracks $`\varepsilon_t`$ (in this example, the $`\varepsilon_a`$ of each iteration is nearly equal to the $`\varepsilon_t`$ of the previous one), so we can use $`\varepsilon_a`$ as a reliable proxy even when the true value is unknown.

---

<br>

## 2. Truncation Error and Taylor Series

### 2.1 Maclaurin Series for e^x: Term-by-Term

The **Maclaurin series** is the Taylor series expanded about $`x = 0`$, and the Maclaurin series for $`e^x`$ is:

```math
e^x = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \cdots = \sum_{k=0}^{\infty} \frac{x^k}{k!}
```

This series converges for all real $`x`$, meaning we can approximate $`e^x`$ to any desired accuracy by including enough terms. Each additional term reduces the truncation error.

The following code builds the approximation term by term for $`e^{0.5}`$ and shows how the error decreases with each additional term.

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

The true value is $`e^{0.5} \approx 1.6487212707`$. As we include more terms, the approximation converges rapidly.

- The 0th order (just $`1`$) has an error of about $`39.3\%`$.
- The 1st order ($`1 + x = 1.5`$) has an error of about $`9.02\%`$.
- The 5th order has an error of about $`0.0014\%`$, which is negligibly small.

The factorial in the denominator grows much faster than the power of $`x`$ in the numerator, which is why the series converges so quickly for moderate values of $`x`$.

### 2.2 Iterative Computation with Stopping Criterion

In practice, we do not hard-code each term. Instead, we write a loop that adds terms iteratively and stops when the approximate error drops below a threshold.

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

Key observations about this implementation are as follows.

- The stopping criterion `es=1e-4` corresponds to $`\varepsilon_s = 0.0001\%`$. This value is below the 5-figure threshold $`0.5 \times 10^{2-5} = 0.0005\%`$ but above the 6-figure threshold $`0.5 \times 10^{2-6} = 0.00005\%`$, so the Scarborough criterion guarantees at least 5 significant figures. In fact, `exp_series(0.5)` stops after 7 iterations at $`\varepsilon_a \approx 9.40 \times 10^{-5}\%`$, with a true relative error of about $`6.2 \times 10^{-6}\%`$.
- The `maxit` parameter provides a safety net: if the criterion is not met within `maxit` terms, for example when $`|x|`$ is very large and many terms are needed, the loop terminates after a maximum number of iterations rather than running forever.
- The function returns all three pieces of information the caller needs: the result, the final error estimate, and the number of iterations used.

### 2.3 Visualization of Truncation Error Convergence

A semilog plot reveals how quickly the truncation error decreases as more terms are added.

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

**On the semilog plot, the error follows a curve that falls more steeply than a straight line.** The factor by which each additional term reduces the error shrinks from about 0.16 to 0.04, so the decrease is even faster than geometric decay at a constant ratio. The reason is that the next term $`x^{n+1}/(n+1)!`$ has a factorial in its denominator, which is characteristic of Taylor series of analytic functions such as $`e^x`$ that converge for all $`x`$.

The horizontal axis of the code holds the order $`n`$, and the number of terms is $`n + 1`$. For $`e^{0.5}`$, the absolute error reaches about $`6.7 \times 10^{-16}`$ at order 13 (14 terms), the level of machine epsilon ($`\approx 2.22 \times 10^{-16}`$), and at order 14 the computed value equals `math.exp(0.5)`, so the error is 0 (a zero cannot be drawn on a semilog plot, so this point does not appear). Adding more terms beyond this point provides no further improvement because roundoff error dominates. This is the point where truncation error and roundoff error meet.

---

<br>

## 3. Roundoff Error and Floating-Point

### 3.1 Machine Epsilon

**Machine epsilon** is the gap $`\varepsilon`$ between 1 and the next larger floating-point number, which is also the smallest power of 2 such that $`1 + \varepsilon \neq 1`$ in the computer's arithmetic. It represents the fundamental precision limit of the floating-point format.

```python
# Machine epsilon
eps = np.finfo(float).eps
print(f'Machine epsilon: {eps}')  # 2.220446049250313e-16

# Demonstration: 1 + eps/2 == 1 is True!
print(1 + eps/2 == 1)   # True (eps/2 is below representable precision)
print(1 + eps == 1)      # False (eps is exactly at the boundary)
```

The first test (`1 + eps/2 == 1`) returns `True`. The number $`1 + \varepsilon/2`$ lies exactly halfway between $`1`$ and $`1 + \varepsilon`$, and since there is no 64-bit floating-point number between them, round-half-to-even rounds it to 1. In other words, `eps/2` is too small to affect the representation of 1 in 64-bit floating point.

The second test (`1 + eps == 1`) returns `False`. The number $`1 + \varepsilon`$ is the next representable floating-point number after $`1`$, so it is represented exactly without rounding and is distinguishable from 1. Machine epsilon is the smallest power of 2 that produces a distinguishable result when added to 1, while the smallest number that does so is slightly larger than $`\varepsilon/2`$.

> **Note:** Machine epsilon is a **relative** measure. Near $`1`$, the gap between adjacent numbers (the absolute precision) is $`\varepsilon \approx 2.22 \times 10^{-16}`$. Near $`10^{10}`$, it is about $`10^{10} \times \varepsilon \approx 2.22 \times 10^{-6}`$ (the actual gap is about $`1.9 \times 10^{-6}`$). The absolute precision grows in proportion to the magnitude, while the relative precision stays nearly the same regardless of magnitude.

### 3.2 Binary Representation of Decimals

Many decimal fractions that appear simple to humans (like $`0.1`$ or $`0.3`$) have **non-terminating representations in binary**, just as $`1/3 = 0.333...`$ is non-terminating in decimal. This is a fundamental source of roundoff error in floating-point arithmetic.

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

The algorithm works by repeatedly doubling the fraction. If the result is at least 1, the next binary digit is 1 (and we subtract 1); otherwise it is 0. For $`0.1`$:

```math
0.1 \times 2 = 0.2 \to 0, \quad 0.2 \times 2 = 0.4 \to 0, \quad 0.4 \times 2 = 0.8 \to 0, \quad 0.8 \times 2 = 1.6 \to 1, \quad \ldots
```

The pattern $`0.0\overline{0011}`$ repeats infinitely. Since double precision stores only 52 bits of mantissa (53 bits including the implicit leading 1), the representation is rounded to a finite number of bits, introducing a small error. For example, the value stored for 0.1 is exactly 0.1000000000000000055511151231257827.... This is why `0.1 + 0.2 != 0.3` occurs in virtually every programming language that uses IEEE 754 floating point.

> **[Computer Architecture]** The IEEE 754 double-precision format stores a number as a binary bit pattern of 1 sign bit, 11 exponent bits, and 52 mantissa bits, so a fraction can be represented exactly with finitely many bits only if its denominator is a power of 2 (0.5, 0.25, 0.375, and so on). A decimal fraction such as 0.1 = 1/10, whose denominator has a prime factor other than 2, has a non-terminating binary expansion as shown above and is rounded at the last mantissa bit when stored; this rounding at storage time is the source of the error described in this section.

### 3.3 Accumulation of Roundoff Error

When many small roundoff errors accumulate over thousands of operations, the total error can become significant. This is demonstrated by summing $`0.0001`$ ten thousand times. The mathematical result is exactly $`1.0`$, but floating-point arithmetic introduces a discrepancy.

```python
# Summing 0.0001 ten thousand times
total_naive = sum(0.0001 for _ in range(10000))
print(f'Naive sum: {total_naive}')
print(f'Error: {abs(1.0 - total_naive):.2e}')  # ~9.38e-14
```

> **Note:** **Differences between Python versions.** The comment `~9.38e-14` in the code above is the result in Python 3.11 and earlier. Since Python 3.12, the built-in `sum()` uses compensated summation for floating-point sums, so this code prints `1.0` and an error of `0.00e+00`. To see the error of naive summation, start from `total = 0.0` and repeat `total += 0.0001` ten thousand times in a `for` loop; the result is then `0.9999999999999062`, with an error of about `9.38e-14`.

Each addition of $`0.0001`$ introduces a tiny roundoff error (because $`0.0001`$ cannot be represented exactly in binary). Over 10,000 additions, these errors accumulate to approximately $`10^{-13}`$. This is small in absolute terms, but potentially problematic in applications that require high accuracy or that use this sum in further calculations, where the error may be amplified.

### 3.4 Kahan Summation Algorithm

The **Kahan summation algorithm** is a technique that dramatically reduces the accumulation of roundoff errors in large summations, and is also called **compensated summation**. The key idea is to maintain a separate compensation variable that tracks the running roundoff error and corrects for it in the next addition.

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

The algorithm works step by step as follows.

1. `y = val - compensation`: Adjust the next value to add by subtracting the accumulated compensation (the error from previous steps).
2. `temp = total + y`: Perform the addition (this introduces a new roundoff error).
3. `compensation = (temp - total) - y`: Compute the roundoff error introduced in step 2. Since `temp = total + y` in exact arithmetic, `(temp - total) - y` should be zero, but in floating point it captures exactly the roundoff error that was lost.
4. `total = temp`: Update the running total.

> **Note:** Kahan summation adds a compensation variable that tracks the accumulated roundoff error and corrects for it in the next addition. This is one of the most important numerical techniques for maintaining accuracy in large summations.

The result is remarkable: in this test, naive summation produces an error of $`\sim 10^{-13}`$ (with the `for` loop of the note above), while Kahan summation achieves **zero error**. The compensation mechanism effectively gives close to double the precision of a single floating-point variable.

### 3.5 Catastrophic Cancellation: Quadratic Formula

The standard quadratic formula for $`ax^2 + bx + c = 0`$ is:

```math
x = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}
```

It suffers from **catastrophic cancellation** (the subtractive cancellation of Section 6.1 of the lecture note) when $`b^2 \gg 4ac`$. In this case, $`\sqrt{b^2 - 4ac} \approx |b|`$, so one of the two roots involves subtracting two nearly equal numbers.

Specifically, when $`b > 0`$, the root $`x_1 = \frac{-b + \sqrt{b^2 - 4ac}}{2a}`$ involves subtracting $`b`$ from $`\sqrt{b^2 - 4ac}`$, two numbers that are nearly identical when $`b`$ is large.

The stable alternative uses **Vieta's formulas**: for $`ax^2 + bx + c = 0`$, the product of roots is $`x_1 \cdot x_2 = c/a`$. So if we can compute one root accurately, we get the other as $`x_2 = c/(a \cdot x_1)`$ without any subtraction.

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

> **Note:** **Wrong choice of root in the code above.** When $`b > 0`$, the root that suffers cancellation is `x1_std = (-b + disc) / (2*a)`, while `x2_std = (-b - disc) / (2*a)` is accurate without cancellation. Contrary to its comment, however, the code above computes `x2_stable = c / (a * x1_std)` from the inaccurate `x1_std`, so the `x2_stable` column actually becomes less accurate as $`b`$ grows (`-134217728` at $`b = 10^8`$). At $`b = 10^9`$, `math.sqrt(b**2 - 4*a*c)` equals $`b`$, so `x1_std` becomes 0 and a `ZeroDivisionError` is raised. Vieta's formula should be used in the other direction, obtaining the small root as $`x_1 = c/(a \cdot x_2)`$ from the accurate large root $`x_2`$. The corrected code follows.

```python
# Corrected: compute the large root without cancellation, then the small root by Vieta
a, c = 1, 1

print(f'{"b":>12} {"x1_standard":>22} {"x1_stable":>22} {"rel_error":>15}')
print('-' * 76)

for exp in range(2, 11):
    b = 10**exp
    disc = math.sqrt(b**2 - 4*a*c)

    x1_std = (-b + disc) / (2*a)    # small root: cancellation
    x2 = (-b - disc) / (2*a)        # large root: no cancellation

    x1_stable = c / (a * x2)        # Vieta: x1 * x2 = c / a

    rel_err = abs(x1_std - x1_stable) / abs(x1_stable) * 100
    print(f'{b:12.0e} {x1_std:22.15e} {x1_stable:22.15e} {rel_err:15.6e}%')
```

> **Note:** **Direction of Vieta's formula.** The product of roots $`x_1 \cdot x_2 = c/a`$ gives either root from the other without subtraction, but the root used as the starting point must be one computed accurately without cancellation. For $`b > 0`$, the corrected code starts from the accurate large root $`x_2`$ and obtains the small root as $`x_1 = c/(a \cdot x_2)`$.

As $`b`$ increases, the standard formula loses more and more significant digits in $`x_1`$ (the root computed with the $`-b + \sqrt{b^2 - 4ac}`$ subtraction), while the stable alternative maintains full precision. In the corrected code, the standard $`x_1`$ has a relative error of about 25% at $`b = 10^8`$ and becomes 0 for $`b \geq 10^9`$, losing all significant digits, while the Vieta $`x_1`$ stays accurate at about $`-1/b`$ for every $`b`$.

The lesson is general: whenever we encounter a subtraction of nearly equal quantities in a formula, we look for an algebraically equivalent form that avoids the subtraction. Common techniques include the following.

- **Rationalizing:** Multiply by the conjugate, for example $`\frac{(\sqrt{a} - \sqrt{b})(\sqrt{a} + \sqrt{b})}{\sqrt{a} + \sqrt{b}} = \frac{a - b}{\sqrt{a} + \sqrt{b}}`$.
- **Vieta's formulas:** Use the product-of-roots relationship.
- **Taylor expansion:** When terms nearly cancel, expand to find the leading non-cancelling term.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Scarborough Criterion | If $`\vert\varepsilon_a\vert < (0.5 \times 10^{2-n})\%`$, the result is correct to at least $`n`$ significant figures. |
| Babylonian Method | It computes $`\sqrt{S}`$ by $`x_{n+1} = (x_n + S/x_n)/2`$ and converges quadratically. |
| Taylor Series | The terms of $`e^x = \sum x^k/k!`$ are added iteratively until the stopping criterion is met. |
| Machine Epsilon | `np.finfo(float).eps` $`\approx 2.22 \times 10^{-16}`$ is the gap between 1 and the next floating-point number. |
| Kahan Summation | It tracks the lost roundoff error in a compensation variable to make large sums accurate. |
| Catastrophic Cancellation | It is avoided with algebraic alternatives (rationalizing, Vieta's formulas, Taylor expansion) that remove the subtraction of nearly equal numbers. |

---

<br>

## Self-Check Questions

1. **Scarborough and stopping criteria:** If `exp_series` is called with `es=1e-6`, that is, $`\varepsilon_s = 10^{-6}\%`$, how many significant figures does the Scarborough criterion guarantee at least?

   <details>
   <summary>Answer</summary>

   The 7-figure threshold is $`0.5 \times 10^{2-7}\% = 5 \times 10^{-6}\%`$, and the 8-figure threshold is $`5 \times 10^{-7}\%`$. Since $`10^{-6}\%`$ is below the 7-figure threshold and above the 8-figure threshold, at least 7 figures are guaranteed.

   </details>

2. **Babylonian method:** Starting from $`S = 2`$ and $`x_0 = 1`$, compute $`x_1`$ and $`x_2`$ by hand and find the true percent relative error of $`x_2`$.

   <details>
   <summary>Answer</summary>

   $`x_1 = (1 + 2/1)/2 = 1.5`$ and $`x_2 = (1.5 + 2/1.5)/2 = (1.5 + 1.333333)/2 \approx 1.416667`$. Since $`\sqrt{2} \approx 1.414214`$, $`\varepsilon_t = |1.414214 - 1.416667| / 1.414214 \times 100\% \approx 0.173\%`$.

   </details>

3. **Machine epsilon:** In double precision, explain why `1 + eps/2 == 1` is `True` but `1 + 0.6*eps == 1` is `False`.

   <details>
   <summary>Answer</summary>

   There is no floating-point number between $`1`$ and $`1 + \varepsilon`$, so a value between them is rounded to the nearer one. $`1 + 0.6\varepsilon`$ is closer to $`1 + \varepsilon`$ than the midpoint, so it becomes $`1 + \varepsilon`$ and differs from 1. $`1 + \varepsilon/2`$ is exactly the midpoint, so round-half-to-even gives 1.

   </details>

4. **Avoiding catastrophic cancellation:** Explain the problem with computing $`\sqrt{x + 1} - \sqrt{x}`$ directly at $`x = 10^8`$, and give the rationalized form.

   <details>
   <summary>Answer</summary>

   $`\sqrt{10^8 + 1} \approx 10000.00005`$ and $`\sqrt{10^8} = 10000`$ are nearly equal, so the subtraction cancels the leading digits. Computing it directly in double precision gives $`5.000000055588 \times 10^{-5}`$, whose relative error with respect to the true value $`4.999999998750 \times 10^{-5}`$ is about $`1.4 \times 10^{-8}`$, a loss of about 8 of the 16 digits. Multiplying by the conjugate gives $`\frac{1}{\sqrt{x + 1} + \sqrt{x}}`$, which has no subtraction and keeps almost all digits accurate.

   </details>

5. **Stable computation of a quadratic:** Give a procedure to compute both roots of $`x^2 + 10^8 x + 1 = 0`$ accurately, together with their approximate values.

   <details>
   <summary>Answer</summary>

   First compute the large root without cancellation, $`x_2 = (-b - \sqrt{b^2 - 4ac})/(2a) \approx -10^8`$. Then obtain the small root by Vieta's formula, $`x_1 = c/(a \cdot x_2) \approx -1.0 \times 10^{-8}`$. Computing $`x_1`$ directly with the standard formula gives about $`-7.45 \times 10^{-9}`$ because of cancellation, a relative error of about 25%.

   </details>
