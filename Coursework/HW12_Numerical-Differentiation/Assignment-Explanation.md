# HW12 Assignment Explanation — Numerical Differentiation

[Repository Guide](../../README.md) | [English](Assignment-Explanation.md) | [한국어](Assignment-Explanation.ko.md)

> **Last Updated:** 2026-11-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 4, 21

> **Learning Objectives:**
> 1. Derive the central difference formula for the first derivative by subtracting two Taylor expansions and implement it as a function
> 2. Verify the results of the implemented function by hand and through the average of the forward and backward differences
> 3. Confirm from the error ratio under a shrinking step size that the forward and backward differences are $`O(h)`$ and the central difference is $`O(h^2)`$
> 4. Explain the error of the central difference for the benchmark function exactly by the terms of the Taylor expansion
> 5. Explain why applying Richardson extrapolation to the central difference gives the same result as the five-point central difference

---

## Table of Contents

- [1. Assignment Overview](#1-assignment-overview)
- [2. Requirements](#2-requirements)
- [3. Problem: Central Difference Function](#3-problem-central-difference-function)
  - [3.1 Approach](#31-approach)
  - [3.2 Code Written](#32-code-written)
  - [3.3 Results and Verification](#33-results-and-verification)
  - [3.4 Results Used by Later Cells](#34-results-used-by-later-cells)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Assignment Overview

> **Try It First:** This document contains the answers and solutions. Solve the problems in the [assignment handout](Original/Assignment.md) yourself before reading on.

The assignment is to fill in the single blank of the lab notebook `Lab13-diff.ipynb` distributed by the instructor and run it. The notebook has eight parts: the benchmark function and a review of Taylor expansion, the forward, backward, and central differences, higher-order stencils and Richardson extrapolation, second derivatives and mixed partial derivatives, step-size selection and the complex-step method, multivariate differentiation, FEM derivative recovery, and engineering applications, and the blank lies inside the function `central_diff` of Section 2.1. All other code cells, including the Exercise cells, were distributed with their code filled in. The problem statement is in the [assignment handout](Original/Assignment.md), the distributed blank notebook is [Original/Lab13-diff.ipynb](Original/Lab13-diff.ipynb), and the notebook with the blank filled in and run is [Solution/Lab13-diff_Completed.ipynb](Solution/Lab13-diff_Completed.ipynb).

| Problem | Notebook section | Blank to fill |
|:-----|:-----|:-----|
| Problem | 2.1 Three basic formulas | The computation of `cdiff` in the function `central_diff` |

*Table 1. The blank that makes up the assignment*

The derivations and error orders of the forward, backward, and central differences are covered in Section 10 of the [Chapter 4 lecture note](../../Lecture-Notes/Ch04_Roundoff-and-Truncation-Errors/Concepts_Lecture.md) and Section 2 of the [Chapter 21 lecture note](../../Lecture-Notes/Ch21_Numerical-Differentiation/Concepts_Lecture.md), and Richardson extrapolation for derivatives in Section 5 of the Chapter 21 lecture note. This document repeats only what the solution needs.

---

<br>

## 2. Requirements

The completed notebook must satisfy the following requirements. The rest of the cell stays as distributed, and only the blank is filled in.

1. **Formula:** `cdiff` must be the central difference approximation of the first derivative at the point `x` with step size `h`. Like `forward_diff` and `backward_diff` in the same cell, it uses only the arguments `func`, `x`, and `h`.
2. **Output:** The comparison table that the cell prints for six step sizes at $`x_0 = 0.5`$ must match the output in the distributed notebook.
3. **Reuse:** Later cells of the notebook reuse this function for the convergence plot, Exercise 2.1, Richardson extrapolation, the step-size analysis, and the comparison with the complex-step method, so it must work for any function, point, and step size.

---

<br>

## 3. Problem: Central Difference Function

### 3.1 Approach

Writing the Taylor expansions at the two points a distance $`h`$ to either side of $`x`$ gives the following.

$$
f(x \pm h) = f(x) \pm h f'(x) + \frac{h^2}{2} f''(x) \pm \frac{h^3}{6} f'''(x) + \cdots
$$

Subtracting the two expansions cancels the even-order terms and leaves only the odd-order terms, doubled.

$$
f(x+h) - f(x-h) = 2h f'(x) + \frac{h^3}{3} f'''(x) + \cdots
$$

Dividing both sides by $`2h`$ gives $`f'(x) \approx \frac{f(x+h) - f(x-h)}{2h}`$, and the leading term of the remaining error is $`\frac{h^2}{6} f'''(x)`$. The blank is therefore the difference of `func(x + h)` and `func(x - h)` divided by `2 * h`. Note that the denominator is $`2h`$, not $`h`$, because the two points are a distance $`2h`$ apart. The derivation of the central difference is covered in Section 2.3 of the [Chapter 21 lecture note](../../Lecture-Notes/Ch21_Numerical-Differentiation/Concepts_Lecture.md).

### 3.2 Code Written

```python
    cdiff = (func(x + h) - func(x - h)) / (2 * h)
```

### 3.3 Results and Verification

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

The output matches the output in the distributed notebook. Three checks confirm the result.

**Hand calculation.** For $`h = 0.5`$, $`f(1.0) = 0.1 - 0.2 + 3 + 2 - 1 + 0.5 = 4.4`$ and $`f(0) = 0.5`$, so the central difference is $`(4.4 - 0.5)/1 = 3.9`$. With $`f(0.5) = 0.865625`$, the forward difference is $`(4.4 - 0.865625)/0.5 = 7.06875`$ and the backward difference is $`(0.865625 - 0.5)/0.5 = 0.73125`$, and their average $`(7.06875 + 0.73125)/2 = 3.9`$ equals the central difference. This relation holds in every row of the table.

**Exact structure of the error.** The benchmark function is a quintic polynomial, so $`f^{(5)} = 12`$ is constant and all higher derivatives are 0. Let $`D_{\text{F}}(h)`$ denote the forward difference and, following the notation of Section 3.2 of the notebook, $`D(h)`$ the central difference. Writing the Taylor expansions out in full, the differences between the approximations and the true value at $`x_0 = 0.5`$ are as follows, where $`f''(0.5) = 12.65`$, $`f'''(0.5) = 17.1`$, and $`f^{(4)}(0.5) = 1.2`$.

$$
D_{\text{F}}(h) - f'(x_0) = \frac{f''}{2}h + \frac{f'''}{6}h^2 + \frac{f^{(4)}}{24}h^3 + \frac{f^{(5)}}{120}h^4 = 6.325h + 2.85h^2 + 0.05h^3 + 0.1h^4
$$

$$
D(h) - f'(x_0) = \frac{f'''}{6}h^2 + \frac{f^{(5)}}{120}h^4 = 2.85h^2 + 0.1h^4
$$

The backward difference differs from the forward difference only in the sign of the odd-order terms: $`-6.325h + 2.85h^2 - 0.05h^3 + 0.1h^4`$. The signed errors are as follows.

| $`h`$ | Forward difference error | Backward difference error | Central difference error | $`2.85h^2 + 0.1h^4`$ |
|:-----|:-----|:-----|:-----|:-----|
| 0.5 | 3.887500e+00 | -2.450000e+00 | 7.187500e-01 | 7.187500e-01 |
| 0.25 | 1.760547e+00 | -1.403516e+00 | 1.785156e-01 | 1.785156e-01 |
| 0.1 | 6.610600e-01 | -6.040400e-01 | 2.851000e-02 | 2.851000e-02 |
| 0.05 | 3.233819e-01 | -3.091306e-01 | 7.125625e-03 | 7.125625e-03 |
| 0.01 | 6.353505e-02 | -6.296505e-02 | 2.850010e-04 | 2.850010e-04 |
| 0.001 | 6.327850e-03 | -6.322150e-03 | 2.850000e-06 | 2.850000e-06 |

*Table 2. Signed errors of the three differences versus step size and the error formula of the central difference*

The central-difference error agrees with $`2.85h^2 + 0.1h^4`$ in every row. The central estimate also exceeds the true value because $`f'''(0.5) > 0`$.

**Convergence ratio of the error.** Halving $`h`$ from 0.5 to 0.25 shrinks the central-difference error by a factor of 4.03 and the forward-difference error by 2.21, and going from 0.1 to 0.05 shrinks them by 4.001 and 2.04. Reducing $`h`$ tenfold from 0.01 to 0.001 shrinks the central-difference error by 100 and the forward-difference error by 10.04. As $`h`$ shrinks, the ratios approach 4 and 100 for $`O(h^2)`$ and 2 and 10 for $`O(h)`$, because the higher-order terms carry less weight. The convergence-plot cell of the notebook also measures orders of 1.00 for the forward difference and 2.00 for the central difference from the errors with $`10^{-4} < h \le 10^{-1}`$.

### 3.4 Results Used by Later Cells

The function `richardson_extrap` of Section 3.2 of the notebook computes $`D(h)`$ and $`D(h/2)`$ with `central_diff` and combines them as $`\frac{4D(h/2) - D(h)}{3}`$. The output of the completed notebook is as follows.

```text
Derivative estimates at x0=0.5  (exact = 3.1812500000)

       h    Central O(h^2)         err       5-pt O(h^4)         err   Richardson O(h^4)         err
---------------------------------------------------------------------------------------------------------
  0.5000      3.9000000000    7.19e-01      3.1562500000    2.50e-02        3.1796875000    1.56e-03
  0.2500      3.3597656250    1.79e-01      3.1796875000    1.56e-03        3.1811523438    9.77e-05
  0.1000      3.2097600000    2.85e-02      3.1812100000    4.00e-05        3.1812475000    2.50e-06
  0.0500      3.1883756250    7.13e-03      3.1812475000    2.50e-06        3.1812498438    1.56e-07
  0.0250      3.1830312891    1.78e-03      3.1812498438    1.56e-07        3.1812499902    9.77e-09
  0.0100      3.1815350010    2.85e-04      3.1812499960    4.00e-09        3.1812499997    2.50e-10
```

- **Error of the extrapolation:** Substituting the error formula $`D(h) - f' = 2.85h^2 + 0.1h^4`$ of Section 3.3 gives $`\frac{4D(h/2) - D(h)}{3} - f' = -0.025h^4`$. This is $`-0.025 \times 0.0625 = -1.5625 \times 10^{-3}`$ at $`h = 0.5`$ and $`-2.5 \times 10^{-6}`$ at $`h = 0.1`$, matching the errors in the table. The $`h^2`$ term is removed exactly, and the error becomes $`O(h^4)`$.
- **Agreement with the five-point central difference:** The Richardson value at step size $`h`$ equals the five-point central difference at step size $`h/2`$ (for example, the extrapolated value 3.1796875000 at $`h = 0.5`$ and the five-point value at $`h = 0.25`$). Let $`k = h/2`$; then the two expressions coincide as follows.

$$
\frac{4D(k) - D(2k)}{3} = \frac{8\left[f(x+k) - f(x-k)\right] - \left[f(x+2k) - f(x-2k)\right]}{12k}
$$

The function is also used in the convergence plot of Section 2.1 of the notebook, the estimate of $`f'(1.5)`$ in Exercise 2.1, the comparison of higher-order formulas in Section 3, and the step-size analysis and comparison with the complex-step method in Section 5. The cells of Section 5 print the optimal step size of the central difference as $`2.24 \times 10^{-6}`$ and show a central-difference error of $`2.78 \times 10^{-10}`$ at $`h = 10^{-5}`$.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Central difference formula | Subtracting the Taylor expansions of $`f(x+h)`$ and $`f(x-h)`$ cancels the even-order terms and gives $`f'(x) \approx \frac{f(x+h) - f(x-h)}{2h}`$, whose denominator is the distance $`2h`$ between the two points. |
| Relation to forward and backward differences | The central difference is the average of the forward and backward differences, and averaging cancels the odd-order error terms of opposite sign. |
| Order of convergence | Halving the step size shrinks the central-difference error by about 4 and the forward and backward errors by about 2, and the measured orders were 2.00 and 1.00. |
| Error for the benchmark function | For the quintic benchmark function, the central-difference error at $`x_0 = 0.5`$ is exactly $`2.85h^2 + 0.1h^4`$. |
| Richardson extrapolation | $`\frac{4D(h/2) - D(h)}{3}`$ removes the $`h^2`$ term so that the error becomes $`O(h^4)`$, and it equals the five-point central difference with step size $`h/2`$. |

---

<br>

## Self-Check Questions

1. **Size of the denominator:** Why is the denominator of the central difference formula $`2h`$ rather than $`h`$, and what happens if one divides by $`h`$?

   <details>
   <summary>Answer</summary>

   The two points $`x - h`$ and $`x + h`$ are a distance $`2h`$ apart. Dividing by $`h`$ gives about $`2f'(x)`$, twice the slope. For the benchmark function with $`h = 0.5`$, the result would be 7.8 instead of 3.9.

   </details>

2. **Average and central difference:** Show with a formula that the average of the forward and backward differences equals the central difference, and explain how this relates to the order of the error.

   <details>
   <summary>Answer</summary>

   $`\frac{1}{2}\left(\frac{f(x+h) - f(x)}{h} + \frac{f(x) - f(x-h)}{h}\right) = \frac{f(x+h) - f(x-h)}{2h}`$. The error of the forward difference is $`\frac{f''}{2}h + \cdots`$ and that of the backward difference is $`-\frac{f''}{2}h + \cdots`$, so their first-order terms have opposite signs, and averaging removes the $`O(h)`$ term and leaves the $`O(h^2)`$ term.

   </details>

3. **Error ratio:** In Table 2, why is the error ratio of the central difference 4.03 rather than exactly 4 when $`h`$ goes from 0.5 to 0.25?

   <details>
   <summary>Answer</summary>

   The error is $`2.85h^2 + 0.1h^4`$, which includes an $`h^4`$ term. Halving $`h`$ shrinks the $`h^4`$ term by a factor of 16, faster than the $`h^2`$ term, so the ratio is slightly above 4. As $`h`$ shrinks, this term carries less weight and the ratio approaches 4.

   </details>

4. **Sign of the error:** In Table 2, explain why the central estimate exceeds the true value for every $`h`$.

   <details>
   <summary>Answer</summary>

   The central-difference error is $`\frac{f'''(x_0)}{6}h^2 + \frac{f^{(5)}}{120}h^4`$, and both $`f'''(0.5) = 17.1`$ and $`f^{(5)} = 12`$ are positive, so the error is always positive.

   </details>

5. **Error of Richardson extrapolation:** If the central-difference error is $`c_2 h^2 + c_4 h^4`$, find the error of $`\frac{4D(h/2) - D(h)}{3}`$.

   <details>
   <summary>Answer</summary>

   $`\frac{4(c_2 h^2/4 + c_4 h^4/16) - (c_2 h^2 + c_4 h^4)}{3} = \frac{c_4 h^4/4 - c_4 h^4}{3} = -\frac{c_4}{4}h^4`$. For the benchmark function $`c_4 = 0.1`$, so the error is $`-0.025h^4`$.

   </details>
