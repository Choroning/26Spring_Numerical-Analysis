# HW6 Assignment Explanation — Root-Finding Methods

[Repository Guide](../../README.md) | [English](Assignment-Explanation.md) | [한국어](Assignment-Explanation.ko.md)

> **Last Updated:** 2026-11-16
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 5, 6

> **Learning Objectives:**
> 1. Locate brackets that enclose roots by sign changes and explain why some roots are missed depending on the grid spacing
> 2. Implement the update formulas of bisection and false position and confirm the iteration counts and error bounds
> 3. Explain why false position can be slower than bisection in terms of the curvature of the function and a fixed endpoint
> 4. Implement the update formulas of the Newton-Raphson and secant methods and confirm quadratic and superlinear convergence from the behavior of the errors

---

## Table of Contents

- [1. Assignment Overview](#1-assignment-overview)
- [2. Requirements](#2-requirements)
- [3. Problem 1: Incremental Search](#3-problem-1-incremental-search)
  - [3.1 Approach](#31-approach)
  - [3.2 Code Written](#32-code-written)
  - [3.3 Results and Verification](#33-results-and-verification)
- [4. Problem 2: Bisection Method](#4-problem-2-bisection-method)
  - [4.1 Approach](#41-approach)
  - [4.2 Code Written](#42-code-written)
  - [4.3 Results and Verification](#43-results-and-verification)
- [5. Problem 3: False Position](#5-problem-3-false-position)
  - [5.1 Approach](#51-approach)
  - [5.2 Code Written](#52-code-written)
  - [5.3 Results and Verification](#53-results-and-verification)
- [6. Problem 4: Newton-Raphson Method](#6-problem-4-newton-raphson-method)
  - [6.1 Approach](#61-approach)
  - [6.2 Code Written](#62-code-written)
  - [6.3 Results and Verification](#63-results-and-verification)
- [7. Problem 5: Secant Method](#7-problem-5-secant-method)
  - [7.1 Approach](#71-approach)
  - [7.2 Code Written](#72-code-written)
  - [7.3 Results and Verification](#73-results-and-verification)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Assignment Overview

> **Try It First:** This document contains the answers and solutions. Solve the problems in the [assignment handout](Original/Assignment.md) yourself before reading on.

The assignment is to fill in the five blanks of the lab notebook `Lab07-RootSolving.ipynb` distributed by the instructor and run it. Using a savings problem and a bungee jumper problem as examples, the notebook implements incremental search, bisection, false position, Newton-Raphson, and the secant method in turn and compares them with `scipy.optimize`, and each of the five method functions contains one blank. The problem statements are in the [assignment handout](Original/Assignment.md), the distributed blank notebook is [Original/Lab07-RootSolving.ipynb](Original/Lab07-RootSolving.ipynb), and the notebook with the blanks filled in and run is [Solution/Lab07-RootSolving_Completed.ipynb](Solution/Lab07-RootSolving_Completed.ipynb).

| Problem | Notebook section | Blanks to fill |
|:-----|:-----|:-----|
| Problem 1 | 2 Incremental search | The `if` condition in the function `incsearch` |
| Problem 2 | 3.1 Bisection | `xr` and `xl` in the function `bisect` |
| Problem 3 | 4.1 False position | `xr` in the function `false_position` |
| Problem 4 | 5.1 Newton-Raphson | `x1` in the function `newtraph` |
| Problem 5 | 6.1 Secant method | `x2` in the function `secant` |

*Table 1. The five blanks that make up the assignment*

Incremental search, bisection, and false position are covered in the [Chapter 5 lecture note](../../Lecture-Notes/Ch05_Roots-Bracketing-Methods/Concepts_Lecture.md), and the Newton-Raphson and secant methods in the [Chapter 6 lecture note](../../Lecture-Notes/Ch06_Roots-Open-Methods/Concepts_Lecture.md). This document repeats only what the solutions need.

---

<br>

## 2. Requirements

The completed notebook must satisfy the following requirements. The rest of each cell (output format, loop structure, plots) stays as distributed, and only the blanks are filled in.

1. **Problem 1:** The `if` condition must be true when the function values at the two ends of the subinterval $`[x_k, x_{k+1}]`$ have different signs.
2. **Problem 2:** `xr` must be the midpoint of the current bracket, and when the sign test is false (the root is in the right half), the left end `xl` must move to the midpoint.
3. **Problem 3:** `xr` must be the x-intercept of the line connecting $`(x_l, f(x_l))`$ and $`(x_u, f(x_u))`$, using `fl` and `fu` computed at the top of the iteration.
4. **Problem 4:** `x1` must be the x-intercept of the tangent line at the current estimate `xi`, using the computed `fxi` and `fpxi`.
5. **Problem 5:** `x2` must be the x-intercept of the secant line through the two previous estimates `x0` and `x1`, using the computed `f0` and `f1`.

---

<br>

## 3. Problem 1: Incremental Search

### 3.1 Approach

The product of two real numbers is negative only when their signs differ, so $`f(x_k) \cdot f(x_{k+1}) < 0`$ means the sign changes across the subinterval. If $`f`$ is continuous, the Intermediate Value Theorem guarantees at least one root in that subinterval. Because the inequality is strict, a case where $`f`$ is exactly 0 at a grid point makes the product 0 and is not caught. Incremental search is covered in Section 5 of the [Chapter 5 lecture note](../../Lecture-Notes/Ch05_Roots-Bracketing-Methods/Concepts_Lecture.md).

### 3.2 Code Written

```python
        if func(x[k]) * func(x[k+1]) < 0:
```

### 3.3 Results and Verification

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

The output matches the output of the distributed notebook. Counting sign changes of $`f(x) = \sin(10x) + \cos(3x)`$ over 200,000 subintervals of $`[3, 6]`$ gives 9 roots, and $`n_s = 100`$ (spacing 0.0303) found all 9. Yet $`n_s = 50`$ (spacing 0.0612) finds only 5, fewer than the 7 of $`n_s = 20`$ (spacing 0.158). Where two roots lie close together, as at 4.21-4.27 and 5.15-5.21, a subinterval containing both has the same sign at its two ends and both roots are missed, and which pairs fall into one subinterval depends on where the grid points lie. The number of brackets found is therefore not monotone in $`n_s`$, and when roots may cluster, the spacing must be made small enough.

---

<br>

## 4. Problem 2: Bisection Method

### 4.1 Approach

The midpoint is $`x_r = (x_l + x_u)/2`$. In the sign test that follows, if $`f(x_l) \cdot f(x_r) < 0`$, the root is in the left half and `xu = xr`; otherwise the root is in the right half $`[x_r, x_u]`$, so the left end moves with `xl = xr`. Bisection and its convergence are covered in Sections 6 and 7 of the [Chapter 5 lecture note](../../Lecture-Notes/Ch05_Roots-Bracketing-Methods/Concepts_Lecture.md).

### 4.2 Code Written

```python
        xr = (xl + xu) / 2  # — midpoint
        ea = abs((xr - xr_old) / xr) if xr != 0 else abs(xr - xr_old)  # approx. relative error

        if func(xl) * func(xr) < 0:     #  — sign test: root in left half?
            xu = xr
        else:
            xl = xr  #      — root in right half
```

### 4.3 Results and Verification

```text
Root found : m = 142.73763308 kg
f(root)    : -5.203e-10
Rel. error : 6.525e-10
Iterations : 31
```

The 31-row iteration table matches the output of the distributed notebook, so see Problem 2 of the [assignment handout](Original/Assignment.md). The results are checked as follows.

- **First iteration:** `xr_old` starts at `xl = 50`, so the first midpoint 150 gives $`\varepsilon_a = |150 - 50|/150 = 0.6667`$, matching the first row of the table.
- **Iteration count:** In iteration $`n`$ the midpoint moves by $`200/2^n`$, so $`\varepsilon_a \approx 200/(2^n \times 142.74)`$. For this to fall below $`10^{-9}`$, $`n > \log_2(200/(142.74 \times 10^{-9})) = 30.38`$, so the loop stops at 31. The value at iteration 31, $`200/2^{31}/142.738 = 6.525 \times 10^{-10}`$, also matches the output.
- **Error bound:** After 31 iterations the true error must be at most $`200/2^{31} = 9.3 \times 10^{-8}`$, and the actual difference from the root 142.737633108 is $`2.5 \times 10^{-8}`$.

> **Note:** The algorithm in Section 3.1 of the notebook multiplies $`\varepsilon_a`$ by 100%, but `ea` in the code is a fraction without the factor 100. Hence `es=1e-9` means a relative error of $`10^{-9}`$, or $`10^{-7}\%`$ as a percentage. The code for false position, Newton-Raphson, and the secant method works the same way.

---

<br>

## 5. Problem 3: False Position

### 5.1 Approach

Setting $`y = 0`$ in the line $`y = f_u + \frac{f_l - f_u}{x_l - x_u}(x - x_u)`$ through $`(x_l, f_l)`$ and $`(x_u, f_u)`$ and solving for $`x`$ gives the formula of the notebook.

$$
x_r = x_u - f_u\,\frac{x_l - x_u}{f_l - f_u}
$$

The top of the iteration already computes `fl, fu = func(xl), func(xu)`, so the blank uses those values directly. False position and the one-sidedness problem are covered in Section 8 of the [Chapter 5 lecture note](../../Lecture-Notes/Ch05_Roots-Bracketing-Methods/Concepts_Lecture.md).

### 5.2 Code Written

```python
        xr = xu - fu * (xl - xu) / (fl - fu)
```

### 5.3 Results and Verification

```text
Bisection:      34 iterations,  final |error| = 3.312e-05
False Position: 42 iterations,  final |error| = 3.313e-05
```

The output matches the output of the distributed notebook. The results are interpreted as follows.

- **Why false position is slow:** The x-intercept of the first line is 205.56, to the right of the true root 142.74. On $`[50, 250]`$, $`f`$ is increasing and concave (its second derivative is negative), so the line connecting the endpoints lies below the curve, and the x-intercept always falls to the right of the root. The first eight iterations therefore all replace `xu`, while the left end `xl = 50` stays fixed. With one end fixed, the bracket closes in on the root slowly and takes more iterations than bisection. This is the one-sidedness mentioned in Section 4.2 of the notebook.
- **Why both final errors are equal:** The reference `m_true = 142.7376` is the actual root 142.737633108 rounded to four decimal places. The two differ by $`3.31 \times 10^{-5}`$, so the final errors of both methods show the rounding error of the reference, not the accuracy of convergence.

---

<br>

## 6. Problem 4: Newton-Raphson Method

### 6.1 Approach

Setting $`y = 0`$ in the tangent line $`y = f(x_i) + f'(x_i)(x - x_i)`$ at $`x_i`$ gives $`x_{i+1} = x_i - f(x_i)/f'(x_i)`$. The top of the iteration computes `fxi = f(xi)` and `fpxi = fp(xi)` and checks whether `fpxi` is near zero, so the blank only writes the update formula with these two values. The Newton-Raphson method and quadratic convergence are covered in Section 5 of the [Chapter 6 lecture note](../../Lecture-Notes/Ch06_Roots-Open-Methods/Concepts_Lecture.md).

### 6.2 Code Written

```python
        x1 = xi - fxi / fpxi  # — Newton-Raphson update
```

### 6.3 Results and Verification

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

The output matches the output of the distributed notebook. The first iteration can be checked by hand. At $`x_0 = 0`$, $`f(0) = -1`$ and $`f'(0) = 2`$, so $`x_1 = 0 - (-1)/2 = 0.5`$, and $`f(0.5) = 0.5 - e^{-0.5} = -0.10653`$.

The true error shrinks as $`6.7 \times 10^{-2}`$, $`8.3 \times 10^{-4}`$, $`1.3 \times 10^{-7}`$, $`2.8 \times 10^{-15}`$, so the number of correct digits roughly doubles each iteration. The coefficient of the quadratic convergence formula $`e_{i+1} \approx \frac{|f''(x^*)|}{2|f'(x^*)|} e_i^2`$ is $`0.181`$, from $`f''(x^*) = -e^{-x^*} = -0.5671`$ and $`f'(x^*) = 1 + e^{-x^*} = 1.5671`$. The actual ratios $`e_{i+1}/e_i^2`$ are $`0.185`$ and $`0.181`$, in agreement with this coefficient.

---

<br>

## 7. Problem 5: Secant Method

### 7.1 Approach

Replacing $`f'(x_1)`$ in the Newton-Raphson formula with the slope $`(f_0 - f_1)/(x_0 - x_1)`$ of the secant through $`(x_0, f_0)`$ and $`(x_1, f_1)`$ gives the following formula.

$$
x_2 = x_1 - f_1\,\frac{x_0 - x_1}{f_0 - f_1}
$$

Its form is the same as that of false position, but the secant method uses the two most recent points without a sign test, so the points do not need to bracket the root. The line `x0, x1 = x1, x2` at the end of the iteration shifts the points by one. The secant method is covered in Section 6 of the [Chapter 6 lecture note](../../Lecture-Notes/Ch06_Roots-Open-Methods/Concepts_Lecture.md).

### 7.2 Code Written

```python
        x2 = x1 - f1 * (x0 - x1) / (f0 - f1)
```

### 7.3 Results and Verification

```text
=== f(x) = x - e^{-x},  true root = 0.5671432904097838 ===
Newton-Raphson:  5 iters,  |error| = 1.110e-16
Secant:          6 iters,  |error| = 0.000e+00
```

The output matches the output of the distributed notebook. In the first iteration, $`x_0 = 0`$ and $`x_1 = 1`$ give $`f_0 = -1`$ and $`f_1 = 1 - e^{-1} = 0.63212`$, so $`x_2 = 1 - 0.63212 \times (0 - 1)/(-1 - 0.63212) = 0.61270`$. The true error then shrinks as $`4.6 \times 10^{-2}`$, $`3.3 \times 10^{-3}`$, $`2.7 \times 10^{-5}`$, $`1.6 \times 10^{-8}`$, $`7.9 \times 10^{-14}`$, 0, showing superlinear convergence in which the logarithm of the error grows by a factor of about 1.6-1.7 per iteration. It takes one more iteration than Newton-Raphson but needs no derivative.

Section 7 of the notebook compares the four completed functions on the bungee jumper problem with a tolerance of $`10^{-8}`$.

```text
Method                        Root       f(root)         Error   Iters  Type        
------------------------------------------------------------------------------------
Bisection             142.73763373     1.282e-08     3.373e-05      28  Bracketing  
False Position        142.73763447     2.776e-08     3.447e-05      34  Bracketing  
Newton-Raphson        142.73763311     7.105e-15     3.311e-05       4  Open        
Secant                142.73763311     0.000e+00     3.311e-05       7  Open        
scipy brentq          142.73763311     0.000e+00     3.311e-05       9  Hybrid      
```

The iteration counts, 28 for bisection, 34 for false position, 4 for Newton-Raphson, and 7 for the secant method, follow the order of the convergence orders. The Error column of about $`3.3 \times 10^{-5}`$ is the rounding error of the reference `m_true` seen in Section 5.3, and the open methods agree on the root itself to eight decimal places.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Incremental search | A subinterval whose end values have a negative product is recorded as a bracket, and two close roots are missed when they fall into one subinterval, so the number of brackets found is not monotone in the spacing. |
| Bisection | The bracket is halved at the midpoint, the error after $`n`$ iterations is at most $`(x_u - x_l)/2^n`$, and the bungee jumper problem stopped after 31 iterations at a relative error of $`10^{-9}`$. |
| False position | It uses the x-intercept of the line through the two endpoints, and when the curve is concave one endpoint stays fixed, which can make it slower than bisection. |
| Newton-Raphson | It updates with the tangent intercept $`x_i - f(x_i)/f'(x_i)`$ and shows quadratic convergence, with the error shrinking as $`e_{i+1} \approx C e_i^2`$. |
| Secant method | It replaces the derivative with the slope of the secant through the two previous points, needs no bracket, and converges superlinearly with order about 1.618. |
| Reference value for the error | An error measured against a rounded true value cannot fall below the rounding error of that reference, so the accuracy of the reference must be checked when judging convergence. |

---

<br>

## Self-Check Questions

1. **Missed roots:** In Problem 1, why did $`n_s = 50`$ find fewer brackets than $`n_s = 20`$?

   <details>
   <summary>Answer</summary>

   When two close roots fall into one subinterval, the function values at its two ends have the same sign and both roots are missed. Which pairs of roots fall into one subinterval depends on where the grid points lie, so more pairs can be grouped together at $`n_s = 50`$ even though its spacing is smaller.

   </details>

2. **Bisection iterations:** For the bracket $`[50, 250]`$ of Problem 2, how many iterations guarantee an absolute error of $`10^{-6}`$?

   <details>
   <summary>Answer</summary>

   From $`200/2^n \le 10^{-6}`$, $`n \ge \log_2(2 \times 10^8) = 27.6`$, so 28 iterations. The following cell of the notebook also prints 28.

   </details>

3. **One-sidedness:** In Problem 3, if $`f`$ were an increasing convex function, which endpoint would stay fixed in false position?

   <details>
   <summary>Answer</summary>

   If $`f`$ is convex, the line connecting the endpoints lies above the curve, so the x-intercept always falls to the left of the root. Then only `xl` changes, and the right end `xu` stays fixed.

   </details>

4. **Checking quadratic convergence:** In Problem 4, given $`e_2 = 8.3 \times 10^{-4}`$, predict $`e_3`$ with the quadratic convergence coefficient 0.181 and compare it with the actual value.

   <details>
   <summary>Answer</summary>

   $`e_3 \approx 0.181 \times (8.3 \times 10^{-4})^2 = 1.25 \times 10^{-7}`$. This agrees with the actual error $`|0.567143165035 - 0.567143290410| = 1.25 \times 10^{-7}`$.

   </details>

5. **Secant method and false position:** The update formulas of the two methods have the same form. Why does the secant method converge faster?

   <details>
   <summary>Answer</summary>

   False position chooses endpoints by the sign test so that the root stays bracketed, which can leave an old endpoint fixed. The secant method always uses the two most recent points, so the secant quickly approaches the tangent near the root and the method converges superlinearly. In exchange, there is no guarantee that the root stays bracketed, and the method can diverge.

   </details>
