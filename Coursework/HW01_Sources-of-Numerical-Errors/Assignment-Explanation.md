# HW1 Assignment Explanation — Sources of Numerical Errors

[Repository Guide](../../README.md) | [English](Assignment-Explanation.md) | [한국어](Assignment-Explanation.ko.md)

> **Last Updated:** 2026-11-18
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 4

> **Learning Objectives:**
> 1. Compute Maclaurin partial sums of the exponential term by term and explain how the true percent relative error changes with the number of terms
> 2. Complete a series-evaluation loop that stops by the approximate percent relative error and the Scarborough criterion
> 3. Compare the accuracy of two expressions for the small root of a quadratic equation as the coefficient grows and explain the difference
> 4. Implement the forward, backward, and central differences and confirm numerically how their errors scale with the step size

---

## Table of Contents

- [1. Assignment Overview](#1-assignment-overview)
- [2. Requirements](#2-requirements)
- [3. Problem 1: Maclaurin Partial Sums of the Exponential](#3-problem-1-maclaurin-partial-sums-of-the-exponential)
  - [3.1 Approach](#31-approach)
  - [3.2 Code Written](#32-code-written)
  - [3.3 Results and Verification](#33-results-and-verification)
- [4. Problem 2: Series Evaluation with a Stopping Criterion](#4-problem-2-series-evaluation-with-a-stopping-criterion)
  - [4.1 Approach](#41-approach)
  - [4.2 Code Written](#42-code-written)
  - [4.3 Results and Verification](#43-results-and-verification)
- [5. Problem 3: Quadratic Roots and Catastrophic Cancellation](#5-problem-3-quadratic-roots-and-catastrophic-cancellation)
  - [5.1 Approach](#51-approach)
  - [5.2 Code Written](#52-code-written)
  - [5.3 Results and Verification](#53-results-and-verification)
- [6. Problem 4: Three Finite-Difference Formulas](#6-problem-4-three-finite-difference-formulas)
  - [6.1 Approach](#61-approach)
  - [6.2 Code Written](#62-code-written)
  - [6.3 Results and Verification](#63-results-and-verification)
- [7. Problem 5: Forward and Central Difference Errors versus Step Size](#7-problem-5-forward-and-central-difference-errors-versus-step-size)
  - [7.1 Approach](#71-approach)
  - [7.2 Code Written](#72-code-written)
  - [7.3 Results and Verification](#73-results-and-verification)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Assignment Overview

> **Try It First:** This document contains the answers and solutions. Solve the problems in the [assignment handout](Original/Assignment.md) yourself before reading on.

The assignment is to fill in the five blanks of the lab notebook `Lab02-SourceOfErrors.ipynb` distributed by the instructor and run it. The notebook has four parts: error definitions and metrics, truncation error and Taylor series, roundoff error and floating-point arithmetic, and finite differences with combined error analysis, and the blanks lie in the last three parts. The problem statements are in the [assignment handout](Original/Assignment.md), the distributed blank notebook is [Original/Lab02-SourceOfErrors.ipynb](Original/Lab02-SourceOfErrors.ipynb), and the notebook with the blanks filled in and run is [Solution/Lab02-SourceOfErrors_Completed.ipynb](Solution/Lab02-SourceOfErrors_Completed.ipynb).

| Problem | Notebook section | Blanks to fill |
|:-----|:-----|:-----|
| Problem 1 | 2.2 Maclaurin series for the exponential | `f_approx[3]`, `f_approx[4]`, `f_approx[5]` |
| Problem 2 | 2.3 Iterative method with a stopping criterion | The `sol` update in the function `exp_series` |
| Problem 3 | 3.5 Catastrophic cancellation | `x1_s`, `x2_stab` |
| Problem 4 | 4.3 Finite-difference example | `dfdx_forward`, `dfdx_backward`, `dfdx_central` |
| Problem 5 | 4.5 Comparing the difference methods | `E_fwd[i]`, `E_cen[i]` |

*Table 1. The five blanks that make up the assignment*

Error definitions, the Scarborough criterion, floating-point representation, and subtractive cancellation are covered in Sections 1-7 of the [Chapter 4 lecture note](../../Lecture-Notes/Ch04_Roundoff-and-Truncation-Errors/Concepts_Lecture.md) and in the [Chapter 4 lab note](../../Lecture-Notes/Ch04_Roundoff-and-Truncation-Errors/Concepts_Lab.md), and Taylor series, finite differences, and the total numerical error are covered in Sections 8-12 of the same lecture note. This document repeats only what the solutions need.

---

<br>

## 2. Requirements

The completed notebook must satisfy the following requirements. The rest of each cell (output format, loop structure, plots) stays as distributed, and only the blanks are filled in.

1. **Problem 1:** `f_approx[k]` must be the partial sum of the terms from order 0 through order $`k`$. It continues in the same way as the given `f_approx[0]`, `f_approx[1]`, and `f_approx[2]`.
2. **Problem 2:** In iteration $`k`$, the term of order $`k`$ must be added to the partial sum `sol`. The rest of the loop uses this value to compute the approximate percent relative error and to decide whether to stop, so the output must match the expected output printed in the notebook.
3. **Problem 3:** `x1_s` must be the large root for each value of $`b`$, and `x2_stab` must be the small root computed with the stable form `c / (a * x1)` given in the comment. The standard-formula value `x2_s` is already filled in.
4. **Problem 4:** The three variables must be the forward, backward, and central difference formulas of Table 1, using `fm`, `fi`, `fp`, and `h` computed in the preceding cell.
5. **Problem 5:** Each iteration must store the absolute errors of the forward and central differences at that iteration's step size `h` in `E_fwd[i]` and `E_cen[i]`. The true value is `true_val = dfdx(0.5)`.

---

<br>

## 3. Problem 1: Maclaurin Partial Sums of the Exponential

### 3.1 Approach

Let $`S_k`$ denote the partial sum stored in `f_approx[k]`. The partial sums satisfy the following recurrence.

$$
S_k = \sum_{i=0}^{k} \frac{x^i}{i!} = S_{k-1} + \frac{x^k}{k!}
$$

Each entry is therefore the previous entry plus one new term. The given line `f_approx[2] = f_approx[1] + x**2 / math.factorial(2)` follows this form, so the blanks add the terms for $`k = 3, 4, 5`$ in the same way. Taylor series and the order of approximation are covered in Sections 8 and 9 of the [Chapter 4 lecture note](../../Lecture-Notes/Ch04_Roundoff-and-Truncation-Errors/Concepts_Lecture.md).

### 3.2 Code Written

```python
f_approx[3] = f_approx[2] + x**3 / math.factorial(3)   # k=3: + x^3/3!
f_approx[4] = f_approx[3] + x**4 / math.factorial(4)   # k=4: + x^4/4!
f_approx[5] = f_approx[4] + x**5 / math.factorial(5)   # k=5: + x^5/5!
```

### 3.3 Results and Verification

```text
 Terms         Approximation    True Rel Err (%)
--------------------------------------------------
     1     1.000000000000000       39.3469340287
     2     1.500000000000000        9.0204010431
     3     1.625000000000000        1.4387677967
     4     1.645833333333333        0.1751622556
     5     1.648437500000000        0.0172115630
     6     1.648697916666667        0.0014164937
```

The true value is $`e^{0.5} = 1.648721270700128`$. Two checks confirm the result.

- **Values of the partial sums:** $`S_3 = 1.625 + 0.5^3/6 = 1.6458\overline{3}`$, $`S_4 = S_3 + 0.5^4/24 = 1.6484375`$, and $`S_5 = S_4 + 0.5^5/120 = 1.648697916\ldots`$, which match the output.
- **Bounds from the remainder:** Keeping six terms leaves the Taylor remainder $`R_5 = \frac{e^{\xi}}{6!} (0.5)^6`$ $`(0 < \xi < 0.5)`$ as the truncation error. Since $`e^{\xi}`$ lies between 1 and $`e^{0.5}`$, $`2.170 \times 10^{-5} < R_5 < 3.578 \times 10^{-5}`$, and the actual error $`e^{0.5} - S_5 = 2.335 \times 10^{-5}`$ lies within these bounds.

The factor by which each added term reduces the error is 0.229, 0.160, 0.122, 0.098, and 0.082, decreasing as more terms are added. This happens because the factorial in the denominator grows faster than the power in the numerator.

---

<br>

## 4. Problem 2: Series Evaluation with a Stopping Criterion

### 4.1 Approach

The function starts from `sol = 1.0` (the term of order 0) and iterates from $`k = 1`$, so the term to add in iteration $`k`$ is $`x^k / k!`$. The code after the blank computes the approximate percent relative error from the previous partial sum `sol_old` and the new partial sum `sol`.

$$
\varepsilon_a = \frac{|S_k - S_{k-1}|}{|S_k|} \times 100 = \frac{x^k / k!}{S_k} \times 100
$$

In other words, $`\varepsilon_a`$ is the share of the newly added term in the partial sum, and the loop stops once this value falls below `es`. The approximate percent relative error and the Scarborough criterion are covered in Sections 3 and 4 of the [Chapter 4 lecture note](../../Lecture-Notes/Ch04_Roundoff-and-Truncation-Errors/Concepts_Lecture.md).

### 4.2 Code Written

```python
        sol = sol + x**k / math.factorial(k)
```

### 4.3 Results and Verification

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

The output matches the expected output in the notebook character for character. The stopping criterion is `scarborough(4)` $`= 0.5 \times 10^{2-4} = 0.005\%`$. At $`k = 5`$, $`\varepsilon_a = 0.0158\% > 0.005\%`$, so the loop continues, and at $`k = 6`$, $`\varepsilon_a = 0.00132\% < 0.005\%`$, so it stops. The true percent relative error is then $`|1.6487212707 - 1.6487196181| / 1.6487212707 \times 100 \approx 1.0 \times 10^{-4}\%`$, which is smaller than the stopping criterion. Rounded to six significant figures, both the result and the true value become 1.64872, so more digits are correct than the four the Scarborough criterion guarantees.

> **Note:** `Terms  : 6` in the output is the order $`k`$ of the last term added. The partial sum contains seven terms, from order 0 through order 6, and the result $`1.648719618\ldots`$ is $`S_6`$, which is $`S_5`$ from Problem 1 plus $`0.5^6 / 6! = 2.170 \times 10^{-5}`$.

---

<br>

## 5. Problem 3: Quadratic Roots and Catastrophic Cancellation

### 5.1 Approach

The cell sets $`a = 1`$, $`b = -b_{\text{val}}`$, and $`c = 1`$, so $`-b`$ is positive. The two roots of the quadratic formula are as follows.

$$
x_1 = \frac{-b + \sqrt{b^2 - 4ac}}{2a}, \qquad x_2 = \frac{-b - \sqrt{b^2 - 4ac}}{2a}
$$

The numerator of $`x_1`$ adds two positive numbers, so no cancellation occurs. When $`b^2 \gg 4ac`$, however, $`\sqrt{b^2 - 4ac} \approx |b|`$, so the numerator of $`x_2`$ subtracts two nearly equal numbers and catastrophic cancellation occurs. The stable form comes from the relation between roots and coefficients, $`x_1 x_2 = c/a`$.

$$
x_2 = \frac{c}{a x_1}
$$

It applies only one division to the accurately computed $`x_1`$, so the relative error does not grow. The blank `x1_s` therefore takes the formula with the plus sign, and `x2_stab` takes the form in the comment. Subtractive cancellation is covered in Section 6.1 of the [Chapter 4 lecture note](../../Lecture-Notes/Ch04_Roundoff-and-Truncation-Errors/Concepts_Lecture.md) and Section 3.5 of the [Chapter 4 lab note](../../Lecture-Notes/Ch04_Roundoff-and-Truncation-Errors/Concepts_Lab.md).

### 5.2 Code Written

```python
    x1_s = (-b_ + disc) / (2*a_)
    x2_s = (-b_ - disc) / (2*a_)  # standard formula
    x2_stab = c_ / (a_ * x1_s)   # stable formula: c / (a * x1)
```

### 5.3 Results and Verification

```text
           b       x2 standard         x2 stable    rel err std (%)    rel err stable (%)
------------------------------------------------------------------------------------------
       1e+02      0.0100010002      0.0100010002       1.221123e-11          0.000000e+00
       1e+04      0.0001000000      0.0001000000       1.117663e-07          0.000000e+00
       1e+06      0.0000010000      0.0000010000       7.614492e-04          0.000000e+00
       1e+08      0.0000000075      0.0000000100       2.549419e+01          0.000000e+00
       1e+10      0.0000000000      0.0000000001       1.000000e+02          0.000000e+00
```

The output matches the output in the distributed notebook. The results are interpreted as follows.

- **Growing error of the standard formula:** Each time $`b`$ grows by a factor of 100, the relative error grows by roughly $`10^4`$, from $`1.2 \times 10^{-11}\%`$ to $`1.1 \times 10^{-7}\%`$ to $`7.6 \times 10^{-4}\%`$. The absolute rounding error of $`\sqrt{b^2 - 4}`$ grows with $`|b|`$, while the target $`x_2 \approx 1/|b|`$ shrinks.
- **The 100% error at $`b = 10^{10}`$:** Near $`b^2 = 10^{20}`$ the spacing of double-precision numbers is 16384, so $`b^2 - 4`$ rounds to $`b^2`$. Then `disc` equals $`|b|`$ exactly, `x2_s` becomes 0, and the relative error is 100%.
- **The zero error of the stable form:** The reference `x2_ref = 1.0 / x1_s` is the same expression as `x2_stab = c_ / (a_ * x1_s)` when $`a = c = 1`$, so the last column is exactly 0 by definition. The accuracy of the stable form rests not on this column but on the reasoning of Section 5.1: `x1_s` involves no cancellation and the division does not amplify the relative error.

> **Note:** In this cell $`b < 0`$, so the formula with the plus sign gives the large root. If $`b > 0`$, the expression $`-b - \sqrt{b^2 - 4ac}`$ becomes the same-sign addition, so in general one first computes the root that adds the square root with the same sign as $`-b`$ and then obtains the other root as $`c / (a x_1)`$.

> **Note:** The first example cell of Problem 3 in the assignment handout sets the true roots of $`x^2 - 100000.0001x + 1 = 0`$ to $`x_1 = 100000.00000999990`$ and $`x_2 = 0.000010000000100000`$, but the actual roots of this equation are $`x_1 = 100000.00009000\ldots`$ and $`x_2 = 9.999999991 \times 10^{-6}`$. The error $`8.00 \times 10^{-5}`$ that the cell prints for the standard-formula $`x_1`$ therefore comes from the misstated true value, and against the actual root the standard-formula $`x_1`$ is accurate to roundoff level. Against the actual root, the error in $`x_2`$ is $`3.9 \times 10^{-12}`$ for the standard formula and below $`10^{-20}`$ for the stable form, so the conclusion that the stable form is more accurate still holds.

---

<br>

## 6. Problem 4: Three Finite-Difference Formulas

### 6.1 Approach

The preceding cell computed `fm` $`= f(x-h)`$, `fi` $`= f(x)`$, and `fp` $`= f(x+h)`$, so the three formulas of Table 1 become differences of these values divided by the step size. The forward difference is $`(f(x+h) - f(x))/h`$, the backward difference is $`(f(x) - f(x-h))/h`$, and the central difference is $`(f(x+h) - f(x-h))/(2h)`$. Their derivations and error orders are covered in Section 10 of the [Chapter 4 lecture note](../../Lecture-Notes/Ch04_Roundoff-and-Truncation-Errors/Concepts_Lecture.md).

### 6.2 Code Written

```python
dfdx_forward  = (fp - fi) / h
dfdx_backward = (fi - fm) / h
dfdx_central  = (fp - fm) / (2*h)
```

### 6.3 Results and Verification

```text
      Method    Estimate    True Rel Err (%)
---------------------------------------------
     Forward    -1.45000             58.9041
    Backward    -0.55000             39.7260
     Central    -1.00000              9.5890
        True    -0.91250
```

By hand, $`(0.2 - 0.925)/0.5 = -1.45`$, $`(0.925 - 1.2)/0.5 = -0.55`$, and $`(0.2 - 1.2)/1 = -1.0`$, matching the output. The Taylor expansion explains the errors exactly. At $`x = 0.5`$, $`f''(0.5) = -1.75`$, $`f'''(0.5) = -2.1`$, and $`f^{(4)} = -2.4`$, so the error of the forward difference is as follows.

$$
\frac{f''}{2}h + \frac{f'''}{6}h^2 + \frac{f^{(4)}}{24}h^3 = -0.4375 - 0.0875 - 0.0125 = -0.5375
$$

This equals $`-1.45 - (-0.9125)`$. In the central difference the even-order terms cancel and only $`\frac{f'''}{6}h^2 = -0.0875`$ remains, and since $`f`$ is a quartic polynomial whose fifth derivative is 0, this is the exact error. The central value $`-1.0`$ also equals the average of the forward and backward values, $`(-1.45 - 0.55)/2`$. The first-order error terms of the two differences have equal size and opposite signs, so they cancel in the average.

---

<br>

## 7. Problem 5: Forward and Central Difference Errors versus Step Size

### 7.1 Approach

The loop records the current step size with `H[i] = h`, passes through the blanks, and finally shrinks the step size with `h /= 10`. The blanks therefore compute both differences with the current `h` and store the absolute difference from the true value. This is the same approach that the notebook's `step_size_analysis` uses for the central difference, and the forward difference divides the difference between $`f(0.5+h)`$ and $`f(0.5)`$ by $`h`$. The balance between truncation and roundoff errors and the optimal step size are covered in Section 12 of the [Chapter 4 lecture note](../../Lecture-Notes/Ch04_Roundoff-and-Truncation-Errors/Concepts_Lecture.md).

### 7.2 Code Written

```python
    E_fwd[i] = abs(true_val - (f(0.5 + h) - f(0.5)) / h)
    E_cen[i] = abs(true_val - (f(0.5 + h) - f(0.5 - h)) / (2*h))
```

### 7.3 Results and Verification

The cell plots the two errors on a log-log graph. The plotted values are listed in the table below.

| $`h`$ | Forward difference error | Central difference error |
|:-----|:-----|:-----|
| 1e+00 | 1.3250e+00 | 3.5000e-01 |
| 1e-01 | 9.1100e-02 | 3.5000e-03 |
| 1e-02 | 8.7851e-03 | 3.5000e-05 |
| 1e-03 | 8.7535e-04 | 3.5000e-07 |
| 1e-04 | 8.7503e-05 | 3.4998e-09 |
| 1e-05 | 8.7500e-06 | 3.3178e-11 |
| 1e-06 | 8.7497e-07 | 5.4228e-12 |
| 1e-07 | 8.6603e-08 | 5.4969e-10 |
| 1e-08 | 8.8872e-09 | 3.3361e-09 |
| 1e-09 | 3.5522e-08 | 1.9989e-08 |
| 1e-10 | 7.5501e-08 | 7.5501e-08 |
| 1e-11 | 3.4062e-06 | 3.4062e-06 |
| 1e-12 | 7.6961e-06 | 7.6961e-06 |
| 1e-13 | 1.0333e-04 | 1.0333e-04 |
| 1e-14 | 2.1171e-03 | 2.1171e-03 |

*Table 2. Absolute errors of the forward and central differences versus step size*

The central-difference column equals the errors printed by `step_size_analysis` in Section 4.4 of the notebook. The results are interpreted as follows.

- **Region dominated by truncation error:** The forward-difference error is close to $`|f''(0.5)| h / 2 = 0.875h`$ and shrinks by 10 when $`h`$ shrinks by 10 (slope 1), while the central-difference error equals $`|f'''(0.5)| h^2 / 6 = 0.35h^2`$ and shrinks by 100 (slope 2). The reference lines $`0.5h`$ and $`0.05h^2`$ in the plot are guides that show these slopes.
- **Step size with the smallest error:** The central difference reaches its minimum at $`h = 10^{-6}`$, close to the value $`\varepsilon_{\text{mach}}^{1/3} \approx 6.06 \times 10^{-6}`$ given in the notebook. The forward difference reaches its minimum at $`h = 10^{-8}`$; because its truncation error is first order in $`h`$, the balance point moves down to around $`\varepsilon_{\text{mach}}^{1/2} \approx 1.49 \times 10^{-8}`$.
- **Region dominated by roundoff error:** For smaller $`h`$ the roundoff error in the values of $`f`$ is divided by $`h`$, so the error grows again. For $`h \le 10^{-10}`$ the two methods show the same error values, and the truncation advantage of the central difference disappears.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Maclaurin partial sums | The partial sums build up as $`S_k = S_{k-1} + x^k / k!`$, and the true percent relative error of $`e^{0.5}`$ after six terms is about 0.0014%, within the bounds of the Taylor remainder. |
| Stopping criterion | Each iteration adds a new term, and the loop stops at $`k = 6`$ when $`\varepsilon_a`$ falls below the Scarborough criterion $`\varepsilon_s = 0.005\%`$, at which point the partial sum contains seven terms. |
| Catastrophic cancellation | When $`b^2 \gg 4ac`$, the small root from the standard formula subtracts two nearly equal numbers and its relative error grows, reaching 100% at $`b = 10^{10}`$ where it becomes 0. |
| Stable root computation | Computing the cancellation-free large root first and then the small root from the relation $`x_2 = c/(a x_1)`$ preserves accuracy. |
| Order of finite differences | The truncation error is $`O(h)`$ for the forward and backward differences and $`O(h^2)`$ for the central difference, which is the average of the forward and backward differences. |
| Optimal step size | Shrinking the step size reduces the truncation error but increases the roundoff error, and the smallest error occurred at $`h = 10^{-6}`$ for the central difference and $`10^{-8}`$ for the forward difference. |

---

<br>

## Self-Check Questions

1. **Terms in the output:** In the output of Problem 2, how many terms does `Terms  : 6` mean were added?

   <details>
   <summary>Answer</summary>

   Seven terms, from order 0 through order 6. The function starts from the order-0 term `sol = 1.0`, and the loop variable `k` is the order of the last term added, so stopping at `k = 6` means seven terms.

   </details>

2. **Positive b:** If $`b`$ had been positive in the cell of Problem 3, how should `x1_s` change to avoid cancellation?

   <details>
   <summary>Answer</summary>

   If $`b > 0`$, then $`-b`$ is negative, so the cancellation-free root is $`(-b - \sqrt{b^2 - 4ac})/(2a)`$. Compute this root first and obtain the other root as $`c/(a x_1)`$. In general, the square root is added with the same sign as $`-b`$.

   </details>

3. **Central difference as an average:** Show with a formula why the central difference in Problem 4 equals the average of the forward and backward differences.

   <details>
   <summary>Answer</summary>

   $`\frac{1}{2}\left(\frac{f(x+h) - f(x)}{h} + \frac{f(x) - f(x-h)}{h}\right) = \frac{f(x+h) - f(x-h)}{2h}`$. The $`f(x)`$ terms cancel and the central difference formula remains.

   </details>

4. **Different optimal step sizes:** In Problem 5, explain why the step size with the smallest forward-difference error is smaller than that of the central difference.

   <details>
   <summary>Answer</summary>

   Treating the total error as truncation plus roundoff, the forward difference behaves like $`Ch + \varepsilon/h`$ and the central difference like $`Ch^2 + \varepsilon/h`$. The first is minimized at $`h \sim \varepsilon^{1/2}`$ and the second at $`h \sim \varepsilon^{1/3}`$, and for $`\varepsilon \approx 2.2 \times 10^{-16}`$ the value $`\varepsilon^{1/2}`$ is smaller. The minima actually appeared at $`10^{-8}`$ and $`10^{-6}`$.

   </details>

5. **Bounds from the remainder:** In Problem 1, why must the truncation error of $`S_5`$ be larger than $`2.17 \times 10^{-5}`$ and smaller than $`3.58 \times 10^{-5}`$?

   <details>
   <summary>Answer</summary>

   The Taylor remainder is $`R_5 = e^{\xi} (0.5)^6 / 6!`$ $`(0 < \xi < 0.5)`$, and since $`1 < e^{\xi} < e^{0.5}`$, it lies between $`(0.5)^6 / 720 = 2.17 \times 10^{-5}`$ and $`e^{0.5} (0.5)^6 / 720 = 3.58 \times 10^{-5}`$. The actual error is $`2.34 \times 10^{-5}`$.

   </details>
