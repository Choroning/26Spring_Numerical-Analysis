# HW7 Assignment Explanation — Optimization Methods

[Repository Guide](../../README.md) | [English](Assignment-Explanation.md) | [한국어](Assignment-Explanation.ko.md)

> **Last Updated:** 2026-11-17
>
> Applied Numerical Methods with Python for Engineers and Scientists, Chapra and Clough, Ch 7

> **Learning Objectives:**
> 1. Update the bracket and interior points of golden-section search and explain by a property of the golden ratio why an old interior point can be reused
> 2. Estimate a minimum by the vertex of the parabola through three points and compare its convergence speed with golden-section search
> 3. Complete the update formula of steepest descent, which moves along the negative gradient, and explain by floating-point limits why the iteration did not stop
> 4. Explain the convergence condition of steepest descent with a fixed step size by the eigenvalues of the Hessian

---

## Table of Contents

- [1. Assignment Overview](#1-assignment-overview)
- [2. Requirements](#2-requirements)
- [3. Problem 1: Golden-Section Search](#3-problem-1-golden-section-search)
  - [3.1 Approach](#31-approach)
  - [3.2 Code Written](#32-code-written)
  - [3.3 Results and Verification](#33-results-and-verification)
- [4. Problem 2: Parabolic Interpolation](#4-problem-2-parabolic-interpolation)
  - [4.1 Approach](#41-approach)
  - [4.2 Code Written](#42-code-written)
  - [4.3 Results and Verification](#43-results-and-verification)
- [5. Problem 3: Steepest Descent](#5-problem-3-steepest-descent)
  - [5.1 Approach](#51-approach)
  - [5.2 Code Written](#52-code-written)
  - [5.3 Results and Verification](#53-results-and-verification)
- [6. Additional Exercises in the Notebook](#6-additional-exercises-in-the-notebook)
  - [6.1 Exercise A: Local Minimum](#61-exercise-a-local-minimum)
  - [6.2 Exercise B: Steepest Descent with a Fixed Step Size](#62-exercise-b-steepest-descent-with-a-fixed-step-size)
- [Summary](#summary)
- [Self-Check Questions](#self-check-questions)

---

<br>

## 1. Assignment Overview

> **Try It First:** This document contains the answers and solutions. Solve the problems in the [assignment handout](Original/Assignment.md) yourself before reading on.

The assignment is to fill in the three blanks of the lab notebook `Lab08-Optimization.ipynb` distributed by the instructor and run it. Using the example function $`f(x) = x^2/10 - 2\sin(x)`$ and a two-dimensional quadratic, the notebook implements golden-section search, parabolic interpolation, and steepest descent before moving on to `scipy.optimize` and engineering applications, and each of the three method functions contains one blank. The notebook also has two Exercise cells with empty code cells, and since these cells are empty in the completed notebook, Section 6 solves them separately with the completed functions. The problem statements are in the [assignment handout](Original/Assignment.md), the distributed blank notebook is [Original/Lab08-Optimization.ipynb](Original/Lab08-Optimization.ipynb), and the notebook with the blanks filled in and run is [Solution/Lab08-Optimization_Completed.ipynb](Solution/Lab08-Optimization_Completed.ipynb).

| Problem | Notebook section | Blanks to fill |
|:-----|:-----|:-----|
| Problem 1 | 2.1-2.3 Golden-section search | `x1` and `x2` in the two branches of the function `goldmin` |
| Problem 2 | 3.1-3.3 Parabolic interpolation | `x4` in the function `parabolic` |
| Problem 3 | 4.1-4.4 Steepest descent | `x` in the function `steepest_descent` |

*Table 1. The three blanks that make up the assignment*

Golden-section search, parabolic interpolation, and steepest descent are covered in Sections 6, 7, and 10 of the [Chapter 7 lecture note](../../Lecture-Notes/Ch07_Optimization/Concepts_Lecture.md). This document repeats only what the solutions need.

---

<br>

## 2. Requirements

The completed notebook must satisfy the following requirements. The rest of each cell (output format, loop structure, plots) stays as distributed, and only the blanks are filled in.

1. **Problem 1:** After the bracket shrinks, only the one newly needed interior point must be computed. In the branch with $`f(x_1) < f(x_2)`$, the new right interior point `x1` is computed from the updated bracket, and in the other branch, the new left interior point `x2`.
2. **Problem 2:** `x4` must be the vertex of the parabola through the three points `x1`, `x2`, and `x3`, using the computed `num` and `den`.
3. **Problem 3:** `x` must be the point reached by moving from the current point opposite to the gradient `g` by the step size `alpha` chosen by the backtracking line search.

---

<br>

## 3. Problem 1: Golden-Section Search

### 3.1 Approach

If $`f(x_1) < f(x_2)`$, the new bracket is $`[x_2, x_u]`$, and the code has already moved the left end and the left interior point with `xl = x2` and `x2 = x1`. What remains is to compute the new right interior point $`x_1 = x_l + (\varphi - 1)(x_u - x_l)`$ in the updated bracket. In the other branch, the new left interior point $`x_2 = x_u - (\varphi - 1)(x_u - x_l)`$ is computed in the new bracket $`[x_l, x_1]`$. Both formulas must use the already updated `xl` or `xu`.

An old interior point can be reused as it is because $`r = \varphi - 1 \approx 0.618`$ satisfies $`r^2 = 1 - r`$. For example, shrinking $`[x_l, x_u]`$ to $`[x_l, x_1]`$ gives the new width $`rL`$ ($`L = x_u - x_l`$), and the new right interior point is $`x_l + r \cdot rL = x_l + (1 - r)L = x_u - rL`$, exactly the old $`x_2`$. The function is therefore evaluated only once per iteration. Golden-section search is covered in Section 6 of the [Chapter 7 lecture note](../../Lecture-Notes/Ch07_Optimization/Concepts_Lecture.md).

### 3.2 Code Written

```python
            x1 = xl + (phi - 1) * (xu - xl)     # — new right interior point
```

```python
            x2 = xu - (phi - 1) * (xu - xl)     # — new left interior point
```

### 3.3 Results and Verification

```text
Minimum found : x* = 1.4275517674
f(x*)         : -1.7757256531
True minimum  : 1.4275517797
|error|       : 1.229e-08
Iterations    : 45
```

The iteration table matches the output of the distributed notebook, so see Problem 1 of the [assignment handout](Original/Assignment.md). The results are checked as follows.

- **First step:** On $`[0, 4]`$, $`d = 0.618 \times 4 = 2.4721`$, so $`x_1 = 2.4721`$ and $`x_2 = 1.5279`$. Since $`f(x_1) = -0.630 \geq f(x_2) = -1.765`$, the new bracket is $`[0, 2.4721]`$, matching the interval width 2.472136 in the first row of the table. The new right interior point is $`0 + 0.618 \times 2.4721 = 1.5279`$, the same as the old $`x_2`$.
- **Iteration count:** In iteration $`i`$, `xint` is $`4 \times 0.618^{i-1}`$, so $`\varepsilon_a = 0.382 \times 4 \times 0.618^{i-1} / 1.4276`$. This formula reproduces $`7.568 \times 10^{-9}`$ at iteration 40 and $`1.787 \times 10^{-9}`$ at iteration 43, and at iteration 45 it becomes $`6.8 \times 10^{-10} < 10^{-9}`$, so the loop stops.
- **Why the error stops near $`10^{-8}`$:** Near the minimum, $`f(x) - f(x^*) \approx \frac{1}{2}f''(x^*)(x - x^*)^2`$, so when this difference falls below the roundoff of the function value, $`\varepsilon_{\text{mach}}|f(x^*)|`$, the function values at two points can no longer be told apart. With $`f''(x^*) = 0.2 + 2\sin(x^*) = 2.18`$ and $`|f(x^*)| = 1.776`$, the resolvable limit is $`|x - x^*| \approx \sqrt{2\varepsilon_{\text{mach}}|f^*|/f''} \approx 1.9 \times 10^{-8}`$, and the actual error $`1.229 \times 10^{-8}`$ is of the same size. For the same reason, the following cell reduces the tolerance to $`10^{-10}`$ and runs 49 iterations, yet the error stays the same.

---

<br>

## 4. Problem 2: Parabolic Interpolation

### 4.1 Approach

The code has already computed the numerator and the denominator of the vertex formula in Section 3.1 of the notebook as `num` and `den`, so $`x_4 = x_2 - \frac{1}{2} \cdot \text{num}/\text{den}`$. Whether the denominator is near zero is checked just before. Parabolic interpolation is covered in Section 7 of the [Chapter 7 lecture note](../../Lecture-Notes/Ch07_Optimization/Concepts_Lecture.md).

### 4.2 Code Written

```python
        x4 = x2 - 0.5 * num / den
```

### 4.3 Results and Verification

```text
Golden-Section:          49 iterations,  final |error| = 1.229e-08
Parabolic Interpolation:  9 iterations,  final |error| = 2.222e-09

Note: Parabolic converges faster when f is smooth near the minimum.

Use: curvature information of a function
```

The output matches the output of the distributed notebook. The results are checked as follows.

- **First step:** The vertex for the three points 0, 1, and 4 is 1.505535, the same value printed by the visualization cell in Section 3.1 of the notebook.
- **Convergence speed:** The error of parabolic interpolation shrinks as $`7.8 \times 10^{-2}`$, $`5.3 \times 10^{-2}`$, $`1.8 \times 10^{-3}`$, $`6.6 \times 10^{-5}`$, $`1.9 \times 10^{-6}`$, $`3.5 \times 10^{-9}`$, far faster than golden-section search, which shrinks by 0.618 per step. This is because it uses the curvature of the function.
- **Ninth iteration:** Eight estimates are recorded, yet the iteration count is 9. In the ninth iteration the function values at the three points became equal in floating point, so `den` was exactly 0, and the function stopped at the denominator check and returned the eighth estimate. The error not shrinking below about $`10^{-9}`$ is the same roundoff limit as in Section 3.3.

---

<br>

## 5. Problem 3: Steepest Descent

### 5.1 Approach

The blank transcribes step 3 of the algorithm, $`\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} - \alpha^{(k)} \mathbf{g}^{(k)}`$. The line just above calls `backtrack_linesearch`, which halves $`\alpha`$ until $`f(\mathbf{x} - \alpha\mathbf{g}) < f(\mathbf{x})`$, so the blank only moves opposite to the gradient by that $`\alpha`$. `x` and `g` are NumPy arrays, so one line updates both coordinates. Steepest descent is covered in Section 10 of the [Chapter 7 lecture note](../../Lecture-Notes/Ch07_Optimization/Concepts_Lecture.md).

### 5.2 Code Written

```python
        x = x - alpha * g                       # x - alpha * g  # — steepest descent step
```

### 5.3 Results and Verification

```text
Starting point    : [-2.  3.]
Minimizer found   : [-1.00000001  1.50000001]
True minimum      : [-1.   1.5]
f(minimizer)      : 0.7500000000
f(true min)       : 0.7500000000
|error|           : 2.107e-08
Iterations        : 200
```

The output matches the output of the distributed notebook. The distance to the minimum $`(-1, 1.5)`$ is $`2.1 \times 10^{-8}`$, yet the loop ran to the maximum of 200 iterations. Inspecting the iteration history of the completed function reveals the cause.

- The gradient norm decreases to $`3.1 \times 10^{-2}`$ at iteration 10, $`3.1 \times 10^{-5}`$ at iteration 30, and $`6.0 \times 10^{-8}`$ at iteration 48, and from iteration 50 on it stays at $`2.98 \times 10^{-8}`$. This value is larger than the tolerance $`10^{-8}`$, so the stopping condition is never met.
- After iteration 50, $`f(\mathbf{x}) - 0.75`$ is 0 in floating point, so no $`\alpha`$ satisfies $`f(\mathbf{x} - \alpha\mathbf{g}) < f(\mathbf{x})`$. The line search returns $`\alpha = 0.5^{50} = 8.9 \times 10^{-16}`$ after 50 halvings, and the point barely moves.

In other words, backtracking judges descent by differences of function values and cannot see an improvement smaller than the roundoff of $`f`$, so a gradient-norm tolerance set below this limit makes the iteration never end. The horizontal part of the notebook's gradient-norm plot from about iteration 50 is this phenomenon.

---

<br>

## 6. Additional Exercises in the Notebook

The code cells of the two Exercises are empty in the completed notebook. The code below was written for this explanation; it calls the completed `goldmin` and `steepest_descent` as they are.

### 6.1 Exercise A: Local Minimum

$`f(x) = x^2/10 - 2\sin(x)`$ is unimodal on $`[5, 8]`$, so `goldmin` can be applied to this interval directly.

```python
# Written for this explanation: the completed notebook leaves this cell empty.
xopt_loc, fopt_loc, ea_loc, n_loc, _ = goldmin(f1d, xl=5.0, xu=8.0, es=1e-9, maxit=100)
print(f'Local minimum : x* = {xopt_loc:.8f},  f(x*) = {fopt_loc:.8f}  ({n_loc} iterations)')
print(f'Global minimum: x* = {xopt:.8f},  f(x*) = {fopt:.8f}')
```

```text
Local minimum : x* = 7.06889124,  f(x*) = 3.58227359  (41 iterations)
Global minimum: x* = 1.42755177,  f(x*) = -1.77572565
```

The local minimum $`f(7.0689) = 3.582`$ is 5.36 larger than the global minimum $`-1.776`$. Golden-section search finds only the minimum inside the given interval, so if the starting interval does not contain the global minimum, it stays at a local minimum.

### 6.2 Exercise B: Steepest Descent with a Fixed Step Size

Following the hint, the line search is replaced with a constant, but to avoid editing the code inside `steepest_descent`, the name `backtrack_linesearch` is temporarily rebound to a function that returns the constant before the call.

```python
# Written for this explanation: the completed notebook leaves this cell empty.
backtrack_original = backtrack_linesearch
x0 = np.array([-2.0, 3.0])
for alpha_fixed in [0.01, 0.1, 0.5]:
    backtrack_linesearch = lambda f, x, grad, a=alpha_fixed: a   # constant step instead of backtracking
    xk, fk, nk, hk = steepest_descent(f2d, grad_f2d, x0, tol=1e-8, maxit=200)
    print(f'alpha = {alpha_fixed:<4}  iterations = {nk:>3}  '
          f'||grad|| = {np.linalg.norm(grad_f2d(xk)):.3e}  |x - x*| = {np.linalg.norm(xk - x2_true):.3e}')
backtrack_linesearch = backtrack_original
```

```text
alpha = 0.01  iterations = 200  ||grad|| = 2.969e-01  |x - x*| = 3.887e-01
alpha = 0.1   iterations = 200  ||grad|| = 1.723e-07  |x - x*| = 2.255e-07
alpha = 0.5   iterations = 200  ||grad|| = 2.038e+41  |x - x*| = 3.893e+40
```

$`f`$ is quadratic with the Hessian $`H = \begin{bmatrix} 4 & 2 \\ 2 & 2 \end{bmatrix}`$, so with a fixed step size the error is multiplied by $`I - \alpha H`$ at each iteration. The eigenvalues of $`H`$ are $`3 \pm \sqrt{5}`$, that is, 0.764 and 5.236, so the error reduction ratio is $`\max_i |1 - \alpha\lambda_i|`$.

| $`\alpha`$ | $`\max_i \vert 1 - \alpha\lambda_i \vert`$ | Result |
|:-----|:-----|:-----|
| 0.01 | 0.992 | Converges, but slowly: $`\lVert\mathbf{x} - \mathbf{x}^*\rVert = 0.39`$ even after 200 iterations |
| 0.1 | 0.924 | Converges to $`2.3 \times 10^{-7}`$ after 200 iterations |
| 0.5 | 1.618 | Diverges because the ratio exceeds 1 |

*Table 2. Error reduction ratio for each fixed step size*

A fixed step size converges only when $`\alpha < 2/\lambda_{\max} = 0.382`$. Too small a step is slow and too large a step diverges, which is why a line search that chooses the step size at each iteration, such as backtracking, is needed.

---

<br>

## Summary

| Concept | Key Summary |
|:-----|:-----|
| Golden-section search | Only one new interior point is computed in the updated bracket, and because $`r = \varphi - 1`$ satisfies $`r^2 = 1 - r`$, an old interior point becomes a new one as it is, so the function is evaluated once per iteration. |
| Accuracy limit of a minimum | Near a minimum, differences of function values are proportional to the square of the distance, so methods that compare function values can locate the minimum only to about $`\sqrt{\varepsilon_{\text{mach}}}`$. |
| Parabolic interpolation | It jumps to the vertex $`x_2 - \frac{1}{2}\text{num}/\text{den}`$ of the parabola through three points, and because it uses curvature, it stopped after 9 iterations in the example instead of the 49 of golden-section search. |
| Steepest descent | It updates with $`\mathbf{x} - \alpha\mathbf{g}`$ using the $`\alpha`$ chosen by backtracking, and when the improvement in the function value falls below the roundoff, the line search fails and the loop runs to the maximum iteration count. |
| Local minimum | A bracketing search finds only the minimum inside the given interval, so on $`[5, 8]`$ it finds the local minimum $`f(7.0689) = 3.582`$, larger than the global minimum. |
| Fixed step size | For a quadratic, the error reduction ratio is $`\max_i \vert 1 - \alpha\lambda_i \vert`$, and the method converges only when $`\alpha < 2/\lambda_{\max}`$. |

---

<br>

## Self-Check Questions

1. **Order of updates:** In the branch of Problem 1 with $`f(x_1) < f(x_2)`$, what happens if `x1 = xl + (phi - 1) * (xu - xl)` is computed before `xl = x2`?

   <details>
   <summary>Answer</summary>

   It is computed in the old bracket $`[x_l, x_u]`$, so the new `x1` becomes the old $`x_1`$, which is also the new `x2`, and the two interior points coincide. The comparison in the next iteration then loses its meaning and the ratio between the bracket and the interior points breaks, so neither the reuse property nor the convergence rate 0.618 holds.

   </details>

2. **Property of the golden ratio:** Verify with $`\varphi = (1 + \sqrt{5})/2`$ that $`r = \varphi - 1`$ satisfies $`r^2 = 1 - r`$.

   <details>
   <summary>Answer</summary>

   $`r = (\sqrt{5} - 1)/2`$, so $`r^2 = (6 - 2\sqrt{5})/4 = (3 - \sqrt{5})/2`$, and $`1 - r = (2 - \sqrt{5} + 1)/2 = (3 - \sqrt{5})/2`$. The two values are equal.

   </details>

3. **Risk of parabolic interpolation:** What problem arises in parabolic interpolation when the three points are nearly collinear?

   <details>
   <summary>Answer</summary>

   The denominator `den` approaches 0, so the vertex jumps very far away or cannot be computed. A parabola that is nearly a straight line has almost no curvature, so its vertex is not well defined. The notebook's function stops iterating when `abs(den) < 1e-30`.

   </details>

4. **An iteration that does not stop:** If the tolerance in Problem 3 had been `tol=1e-6`, around which iteration would the loop end?

   <details>
   <summary>Answer</summary>

   The gradient norm was $`3.1 \times 10^{-5}`$ at iteration 30 and shrank by about 0.7 per iteration, so it falls below $`10^{-6}`$ around iteration 40. This is before iteration 50, where the roundoff limit is reached, so the loop stops normally instead of running to the maximum. Running it indeed stops at iteration 40.

   </details>

5. **Upper limit of a fixed step size:** For the function of Exercise B, does $`\alpha = 0.3`$ converge?

   <details>
   <summary>Answer</summary>

   $`|1 - 0.3 \times 0.764| = 0.771`$ and $`|1 - 0.3 \times 5.236| = 0.571`$, so the reduction ratio is 0.771, smaller than 1, and the method converges. This is because 0.3 is below the limit $`2/\lambda_{\max} = 0.382`$.

   </details>
