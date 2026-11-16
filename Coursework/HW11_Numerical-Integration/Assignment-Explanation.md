# HW11 Assignment Explanation — Numerical Integration

[Repository Guide](../../README.md) | [English](Assignment-Explanation.md) | [한국어](Assignment-Explanation.ko.md)

> **Last Updated:** 2026-11-16
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 19, 20

> **Learning Objectives:**
> 1. Implement the composite trapezoidal rule with NumPy array indexing and sums
> 2. Verify the results of the implemented function by hand, including the case of a single sub-interval
> 3. Confirm from the error ratio under doubling of the number of sub-intervals that the error is proportional to $`h^2`$
> 4. Compare the error estimate based on the average second derivative with the actual error
> 5. Explain that the first column of Romberg integration consists of composite trapezoidal values and interpret the extrapolated results

---

## Table of Contents

- [1. Assignment Overview](#1-assignment-overview)
- [2. Requirements](#2-requirements)
- [3. Problem: Composite Trapezoidal Rule](#3-problem-composite-trapezoidal-rule)
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

The assignment is to fill in the single blank of the lab notebook `Lab12-Int.ipynb` distributed by the instructor and run it. The notebook has eight parts: the benchmark function, the trapezoidal rule, Simpson's rules, Richardson extrapolation and Romberg integration, Gauss-Legendre quadrature, Gauss-Lobatto quadrature, `scipy.integrate`, and engineering applications, and the blank lies inside the function `composite_trap` of Section 2. All other code cells, including the Exercise cells, were distributed with their code filled in. The problem statement is in the [assignment handout](Original/Assignment.md), the distributed blank notebook is [Original/Lab12-Int.ipynb](Original/Lab12-Int.ipynb), and the notebook with the blank filled in and run is [Solution/Lab12-Int_Completed.ipynb](Solution/Lab12-Int_Completed.ipynb).

| Problem | Notebook section | Blank to fill |
|:-----|:-----|:-----|
| Problem | 2.2 Composite trapezoidal rule | The computation of `I` in the function `composite_trap` |

*Table 1. The blank that makes up the assignment*

The derivation and error analysis of the trapezoidal and composite trapezoidal rules are covered in Sections 2 and 5 of the [Chapter 19 lecture note](../../Lecture-Notes/Ch19_Numerical-Integration-Formulas/Concepts_Lecture.md), and Richardson extrapolation and Romberg integration in Section 2 of the [Chapter 20 lecture note](../../Lecture-Notes/Ch20_Numerical-Integration-of-Functions/Concepts_Lecture.md). This document repeats only what the solution needs.

---

<br>

## 2. Requirements

The completed notebook must satisfy the following requirements. The rest of the cell stays as distributed, and only the blank is filled in.

1. **Formula:** `I` must be the value of the composite trapezoidal rule over $`n`$ sub-intervals. The cell has already computed the nodes `x`, the function values `y`, and the sub-interval width `h`, and the blank is the parenthesized part after `I = h/2 * `.
2. **Single sub-interval:** For $`n = 1`$ the value must equal `I_trap1`, the single-segment trapezoidal rule in the same cell.
3. **Array operations:** `func` takes a NumPy array and returns an array, so the computation must work with array indexing and sums, without loops. Later cells of the notebook reuse this function for the convergence study, the comparison with Simpson's rules, Romberg integration, and the engineering applications.
4. **Output:** The output of the cell must match the output in the distributed notebook.

---

<br>

## 3. Problem: Composite Trapezoidal Rule

### 3.1 Approach

Let the nodes be $`x_i = a + ih`$ $`(i = 0, 1, \ldots, n)`$ and the function values $`f_i = f(x_i)`$. Applying the trapezoidal rule to each sub-interval $`[x_{i-1}, x_i]`$ and adding, each interior node is used once by each of its two neighboring sub-intervals and is added twice, while the two end nodes are added once.

$$
I \approx \frac{h}{2}\left(f_0 + 2\sum_{i=1}^{n-1} f_i + f_n\right)
$$

Since `y = func(x)` is the array $`(f_0, f_1, \ldots, f_n)`$ of length $`n + 1`$, each part of the formula corresponds to the following array operation.

| Term of the formula | Array operation |
|:-----|:-----|
| $`f_0`$ | `y[0]` |
| $`f_n`$ | `y[-1]` |
| $`\sum_{i=1}^{n-1} f_i`$ | `np.sum(y[1:-1])` |

*Table 2. Correspondence between the terms of the composite trapezoidal rule and array operations*

For $`n = 1`$, `y[1:-1]` is an empty array and `np.sum` returns 0, so the expression becomes $`\frac{h}{2}(f_0 + f_1)`$, which is the single-segment trapezoidal rule. The derivation of the composite trapezoidal rule is covered in Section 5.2 of the [Chapter 19 lecture note](../../Lecture-Notes/Ch19_Numerical-Integration-Formulas/Concepts_Lecture.md).

### 3.2 Code Written

```python
    I = h/2 * (y[0] + 2*np.sum(y[1:-1]) + y[-1])
```

### 3.3 Results and Verification

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

The output matches the output in the distributed notebook. Three checks confirm the result.

**Hand calculation.** $`f(0) = 0.2`$, $`f(0.4) = 2.456`$, and $`f(0.8) = 0.232`$. For $`n = 1`$, $`h = 0.8`$, so $`I = 0.4 \times (0.2 + 0.232) = 0.1728`$, which equals the single-segment trapezoidal value. For $`n = 2`$, $`h = 0.4`$, so $`I = 0.2 \times (0.2 + 2 \times 2.456 + 0.232) = 0.2 \times 5.344 = 1.0688`$. Both values match the output.

**Convergence ratio of the error.** Rerunning the extracted code to compute the error $`E = I_{\text{exact}} - I`$ and the error ratio when $`n`$ doubles gives the following.

| $`n`$ | $`h`$ | $`I`$ | Error $`E`$ | Ratio | $`4h^2`$ |
|:-----|:-----|:-----|:-----|:-----|:-----|
| 1 | 0.8 | 0.1728000000 | 1.467733e+00 | | 2.560000e+00 |
| 2 | 0.4 | 1.0688000000 | 5.717333e-01 | 2.5672 | 6.400000e-01 |
| 4 | 0.2 | 1.4848000000 | 1.557333e-01 | 3.6712 | 1.600000e-01 |
| 8 | 0.1 | 1.6008000000 | 3.973333e-02 | 3.9195 | 4.000000e-02 |
| 16 | 0.05 | 1.6305500000 | 9.983333e-03 | 3.9800 | 1.000000e-02 |
| 32 | 0.025 | 1.6380343750 | 2.498958e-03 | 3.9950 | 2.500000e-03 |
| 64 | 0.0125 | 1.6399083984 | 6.249349e-04 | 3.9987 | 6.250000e-04 |

*Table 3. Error and error ratio of the composite trapezoidal rule versus the number of sub-intervals*

The ratio approaches 4, in agreement with the theory that the error is proportional to $`h^2`$. The log-log slope that the notebook's convergence-study cell computes from the errors for $`n = 2`$ through $`256`$ is also $`-1.98`$, close to $`-2`$.

**Comparison with the error estimate.** The error estimate of the Chapter 19 lecture note is $`E_a \approx -\frac{b-a}{12}\,\overline{f''}\,h^2`$, where the average second derivative is $`\overline{f''} = \frac{f'(b) - f'(a)}{b - a}`$. From $`f'(x) = 25 - 400x + 2025x^2 - 3600x^3 + 2000x^4`$, $`f'(0) = 25`$ and $`f'(0.8) = -23`$, so $`\overline{f''} = -48/0.8 = -60`$ and $`E_a \approx \frac{0.8}{12} \times 60 \times h^2 = 4h^2`$. In Table 3 the actual error approaches $`4h^2`$ as $`n`$ grows, and at $`n = 64`$ the two agree to four significant figures.

> **Note:** The actual errors in Table 3 agree with $`4h^2 - \frac{8}{3}h^4`$ for every $`n`$ within roundoff. This is the case $`C = 4`$, $`D = -8/3`$ of the error expansion $`I = I_h + Ch^2 + Dh^4 + O(h^6)`$ given in Section 2.2 of the [Chapter 20 lecture note](../../Lecture-Notes/Ch20_Numerical-Integration-of-Functions/Concepts_Lecture.md), with no terms of order $`h^6`$ or higher. The ratio approaches 4 from below because of the $`h^4`$ term of opposite sign, which weighs more when $`n`$ is small.

### 3.4 Results Used by Later Cells

The function is reused throughout the notebook. The function `romberg` of Section 4 fills the first column of the Romberg table with values of `composite_trap` for $`n = 2^k`$ sub-intervals and then applies Richardson extrapolation. The output of the completed notebook is as follows.

```text
Romberg Table:
Level k   n=2^k      R[k,0]          R[k,1]          R[k,2]          R[k,3]    
----------------------------------------------------------------------------------
      0       1    0.1728000000                                                
      1       2    1.0688000000    1.3674666667                                
      2       4    1.4848000000    1.6234666667    1.6405333333                
      3       8    1.6008000000    1.6394666667    1.6405333333    1.6405333333

Best estimate R[3,3] = 1.640533333333
Exact                  = 1.640533333333
Error                  = 2.89e-15
```

- **First column:** `R[k,0]` equals the values for $`n = 1, 2, 4, 8`$ in Section 3.3. A wrong blank would have made the whole Romberg table wrong.
- **Second column:** `R[1,1]` $`= (4 \times 1.0688 - 0.1728)/3 = 1.3674667`$, which equals the composite Simpson's 1/3 value for two sub-intervals, 1.36746667, in Section 3 of the notebook. `R[2,1]` and `R[3,1]` also equal the Simpson's 1/3 values for four and eight sub-intervals, showing that extrapolating two trapezoidal values yields Simpson's 1/3 rule.
- **Third column:** From `R[2,2]` on, the exact value 1.6405333333 appears. As the note in Section 3.3 shows, the trapezoidal error of this function has only an $`h^2`$ term and an $`h^4`$ term, so two extrapolations remove both.

The function also serves as the reference method in the convergence plot of Section 2, the comparison with Simpson's rules in Section 3, the comparison with Gauss-Legendre quadrature in Section 5, and the heat-flux, arc-length, and thrust-integral applications of Section 8.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Composite trapezoidal rule | Applying the trapezoidal rule to $`n`$ sub-intervals and adding gives $`\frac{h}{2}(f_0 + 2\sum_{i=1}^{n-1} f_i + f_n)`$, in which interior nodes are added twice and the end nodes once. |
| Array implementation | The ends are `y[0]` and `y[-1]`, the interior sum is `np.sum(y[1:-1])`, and for $`n = 1`$ the interior sum is 0, so the single-segment formula results. |
| Order of convergence | The error ratio approaches 4 when $`n`$ doubles, so the error is proportional to $`h^2`$, and the log-log slope for the benchmark function was $`-1.98`$. |
| Error estimate | Substituting $`\overline{f''} = -60`$ into $`E_a \approx -\frac{b-a}{12}\overline{f''}h^2`$ gives $`4h^2`$, which approaches the actual error as $`n`$ grows. |
| Romberg integration | The first column of the Romberg table holds composite trapezoidal values, the once-extrapolated column holds composite Simpson's 1/3 values, and for the benchmark function the twice-extrapolated values are exact. |

---

<br>

## Self-Check Questions

1. **Coefficient of interior nodes:** Why are the function values at interior nodes multiplied by 2 in the composite trapezoidal rule?

   <details>
   <summary>Answer</summary>

   An interior node $`x_i`$ $`(1 \le i \le n-1)`$ is the right end of the left sub-interval $`[x_{i-1}, x_i]`$ and the left end of the right sub-interval $`[x_i, x_{i+1}]`$, so it is used once by each of the two trapezoids. Each end node belongs to only one sub-interval and is added once.

   </details>

2. **One sub-interval:** In the code written, what happens for $`n = 1`$, and what does the result equal?

   <details>
   <summary>Answer</summary>

   `x` holds the two points of $`[a, b]`$ and `y[1:-1]` is an empty array, so `np.sum(y[1:-1])` is 0. Hence $`I = \frac{h}{2}(f(a) + f(b))`$ with $`h = b - a`$, which is the single-segment trapezoidal rule. For the benchmark function it is 0.1728.

   </details>

3. **Error ratio:** If the error is proportional to $`h^2`$, by what factor does it shrink when $`n`$ doubles? Why was the ratio only 2.57 when $`n`$ went from 1 to 2 in Table 3?

   <details>
   <summary>Answer</summary>

   Halving $`h`$ shrinks the error to $`(1/2)^2 = 1/4`$, so the ratio is 4. When $`n`$ is small and $`h`$ is large, the $`h^4`$ term of the error expansion carries more weight, so the ratio deviates from the prediction of the leading $`h^2`$ term alone. For the benchmark function the $`h^4`$ term has the opposite sign, so the ratio is below 4.

   </details>

4. **Error estimate:** The error estimate $`E_a \approx -\frac{b-a}{12}\overline{f''}h^2`$ gives 2.56 at $`n = 1`$, while the actual error is 1.468. Why does the estimate become more accurate as $`n`$ grows?

   <details>
   <summary>Answer</summary>

   The estimate keeps only the leading $`h^2`$ term of the error expansion. When $`h`$ is large the neglected $`h^4`$ term ($`-\frac{8}{3}h^4 \approx -1.09`$) is large, but as $`h`$ shrinks this term decreases much faster than the $`h^2`$ term, so the estimate approaches the actual error.

   </details>

5. **First column of the Romberg table:** Write the first extrapolation `R[1,1]` of Romberg integration in terms of the trapezoidal values `R[0,0]` and `R[1,0]`, and compute its value for the benchmark function.

   <details>
   <summary>Answer</summary>

   `R[1,1]` $`= (4R_{1,0} - R_{0,0})/3 = (4 \times 1.0688 - 0.1728)/3 = 1.3674667`$. This value equals the composite Simpson's 1/3 value for two sub-intervals.

   </details>
