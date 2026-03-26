# Chapter 18 Lecture --- Splines and Piecewise Interpolation

> **Last Updated:** 2026-03-26

---

<br>

## Table of Contents

1. [Introduction and Motivation](#1-introduction-and-motivation)
2. [Continuity Conditions](#2-continuity-conditions)
3. [Linear Spline](#3-linear-spline)
4. [Quadratic Spline](#4-quadratic-spline)
   - 4.1 [General Form](#41-general-form)
   - 4.2 [Deriving the Conditions](#42-deriving-the-conditions)
   - 4.3 [Example 18.2 --- Quadratic Spline](#43-example-182----quadratic-spline)
5. [Cubic Spline](#5-cubic-spline)
   - 5.1 [General Form](#51-general-form)
   - 5.2 [Deriving the Conditions (i -- iv)](#52-deriving-the-conditions-i----iv)
   - 5.3 [Eliminating Unknowns --- The Tridiagonal System for $c_i$](#53-eliminating-unknowns----the-tridiagonal-system-for-c_i)
   - 5.4 [End Conditions (Section 18.4.2)](#54-end-conditions-section-1842)
   - 5.5 [Example 18.3 --- Cubic Spline](#55-example-183----cubic-spline)
6. [Bilinear Interpolation](#6-bilinear-interpolation)
7. [Summary Table](#7-summary-table)

---

<br>

## 1. Introduction and Motivation

A **spline** is originally a physical drafting tool --- a long flexible strip fixed at a few control points whose tension creates a smooth curve passing through those points.

In numerical analysis, splines minimize oscillations by fitting **low-order polynomials** to data in a **piecewise** manner, rather than using a single high-order polynomial that can oscillate wildly between data points (Runge's phenomenon).

> **Key Insight:** A simple linear spline (connecting data points with straight line segments) can provide a **better approximation** than a single high-order polynomial, especially near the boundaries of the data range (Fig. 18.1).

The fundamental trade-off:

| Approach | Pros | Cons |
|---|---|---|
| Single high-order polynomial | Globally smooth | Oscillation (Runge's phenomenon) |
| Piecewise low-order splines | No oscillation, local control | Smoothness at knots must be enforced |

---

<br>

## 2. Continuity Conditions

When two piecewise polynomial segments meet at a point (called a **knot**), we must specify how smooth the connection should be. This is described by **$C^k$ continuity**.

- **$C^0$ continuity:** The positions (function values) of the curves are continuous at the knot.

$$f_1 = f_2$$

- **$C^1$ continuity:** The positions **and** 1st derivatives are continuous at the knot.

$$f_1 = f_2, \quad f_1' = f_2'$$

- **$C^2$ continuity:** The positions, 1st derivatives, **and** 2nd derivatives are continuous at the knot.

$$f_1 = f_2, \quad f_1' = f_2', \quad f_1'' = f_2''$$

- **$C^n$ continuity:** The 0th, 1st, 2nd, ..., $n$-th derivatives are all continuous at the knot.

> **Visual Interpretation:**
> - $C^0$: The curve is connected but may have sharp corners (kinks).
> - $C^1$: The curve is connected and tangent directions match --- no visible corners, but curvature may change abruptly.
> - $C^2$: The curve is connected, tangent directions match, and curvature is also continuous --- visually very smooth.

---

<br>

## 3. Linear Spline

**Requires:** $C^0$ Continuity

### Setup

Consider $n$ data points: $(x_1, f_1), (x_2, f_2), \dots, (x_n, f_n)$.

- There are $n - 1$ intervals.
- Each interval $i$ has its own spline function $S_i(x)$.

### Formula

For the $i$-th interval $[x_i,\; x_{i+1}]$:

$$S_i(x) = f_i + \frac{f_{i+1} - f_i}{x_{i+1} - x_i}(x - x_i)$$

This is simply **linear interpolation** between consecutive data points.

### Properties

- Satisfies $C^0$ continuity automatically (the segments share endpoints).
- **Not smooth at knots** --- the first derivative is generally discontinuous where two segments meet.

> **Limitation:** The lack of smoothness at knots motivates the use of higher-order polynomial splines (quadratic and cubic) to ensure smoothness.

---

<br>

## 4. Quadratic Spline

**Requires:** $C^1$ Continuity

### 4.1 General Form

Each interval $i$ uses a second-degree polynomial:

$$S_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2$$

Each segment has **3 unknowns** ($a_i$, $b_i$, $c_i$).

With $n$ data points and $n-1$ intervals, the total number of unknowns is $3(n-1)$.

### 4.2 Deriving the Conditions

We need $3(n-1)$ equations to determine all unknowns.

---

**Condition (i) --- Function values at left endpoints:**

$$S_i(x_i) = a_i = f_i, \quad \text{for } i = 1, 2, \dots, n-1$$

This gives $n - 1$ equations and immediately determines $a_i = f_i$.

---

**Condition (ii) --- $C^0$ continuity at interior knots $x_{i+1}$:**

$$S_i(x_{i+1}) = S_{i+1}(x_{i+1}) = f_{i+1}$$

Substituting the polynomial form and letting $h_i := x_{i+1} - x_i$:

$$f_i + b_i h_i + c_i h_i^2 = f_{i+1}, \quad \text{for } i = 1, 2, \dots, n-1$$

This gives $n - 1$ equations.

---

**Condition (iii) --- $C^1$ continuity at $(n-2)$ interior knots:**

The derivative of each segment is:

$$S_i'(x) = b_i + 2c_i(x - x_i)$$

At each interior knot $x_{i+1}$, equating the derivatives of adjacent segments:

$$S_i'(x_{i+1}) = S_{i+1}'(x_{i+1})$$

$$b_i + 2c_i h_i = b_{i+1}, \quad \text{for } i = 1, 2, \dots, n-2$$

This gives $n - 2$ equations.

---

**Equation count:**

| Source | Count |
|---|---|
| (i) Left endpoint values | $n - 1$ |
| (ii) $C^0$ continuity | $n - 1$ |
| (iii) $C^1$ continuity | $n - 2$ |
| **Total** | $3(n-1) - 1$ |

We are **one equation short**. Therefore, an additional condition is needed.

---

**Condition (iv) --- Extra condition:**

Choose $c_1 = 0$, which means the first segment becomes a **linear approximation**:

$$S_1(x) = f_1 + b_1(x - x_1)$$

This closes the system and allows sequential computation of all coefficients.

### 4.3 Example 18.2 --- Quadratic Spline

**Given:** 4 data points $(x_1, f_1),\, (x_2, f_2),\, (x_3, f_3),\, (x_4, f_4)$

- 3 intervals
- $3 \times 3 = 9$ unknowns (but $a_i = f_i$ reduces it to 6)

**Three spline functions (after applying $a_i = f_i$):**

$$S_1(x) = f_1 + b_1(x - x_1) + c_1(x - x_1)^2$$

$$S_2(x) = f_2 + b_2(x - x_2) + c_2(x - x_2)^2$$

$$S_3(x) = f_3 + b_3(x - x_3) + c_3(x - x_3)^2$$

**From $C^0$ continuity (condition ii):**

$$f_2 - f_1 = b_1 h_1 + c_1 h_1^2$$

$$f_3 - f_2 = b_2 h_2 + c_2 h_2^2$$

$$f_4 - f_3 = b_3 h_3 + c_3 h_3^2$$

**From $C^1$ continuity (condition iii):**

$$b_1 + 2c_1 h_1 - b_2 = 0$$

$$b_2 + 2c_2 h_2 - b_3 = 0$$

**From extra condition (iv):**

$$c_1 = 0$$

This gives **6 equations in 6 unknowns** ($b_1, c_1, b_2, c_2, b_3, c_3$).

**Matrix / vector form:**

$$\begin{pmatrix} f_2 - f_1 \\ f_3 - f_2 \\ f_4 - f_3 \\ 0 \\ 0 \\ 0 \end{pmatrix} = \begin{pmatrix} h_1 & h_1^2 & 0 & 0 & 0 & 0 \\ 0 & 0 & h_2 & h_2^2 & 0 & 0 \\ 0 & 0 & 0 & 0 & h_3 & h_3^2 \\ 1 & 2h_1 & -1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 2h_2 & -1 & 0 \\ 0 & 1 & 0 & 0 & 0 & 0 \end{pmatrix} \begin{pmatrix} b_1 \\ c_1 \\ b_2 \\ c_2 \\ b_3 \\ c_3 \end{pmatrix}$$

> **Note:** This is a linear system $\mathbf{A}\mathbf{x} = \mathbf{b}$ that can be solved using any standard technique (Gaussian elimination, LU decomposition, etc.).

---

<br>

## 5. Cubic Spline

**Requires:** $C^2$ Continuity

### 5.1 General Form

Each interval $i$ uses a third-degree polynomial:

$$S_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3$$

Each segment has **4 unknowns** ($a_i$, $b_i$, $c_i$, $d_i$).

With $n$ data points and $n-1$ intervals, the total number of unknowns is $4(n-1)$.

Key properties of cubic splines:

- Cubic curves are **popular in computer graphics**.
- Piecewise cubic polynomials allow for **$C^2$ continuity**.
- Nice **tradeoff** between computational cost and smoothness.

### 5.2 Deriving the Conditions (i -- iv)

Let $h_i := x_{i+1} - x_i$.

---

**Condition (i) --- Function values at left endpoints:**

$$S_i(x_i) = a_i = f_i, \quad \text{for } i = 1, 2, \dots, n-1$$

Gives $n - 1$ equations. Immediately: $a_i = f_i$.

---

**Condition (ii) --- $C^0$ continuity at knots:**

$$S_i(x_{i+1}) = f_{i+1}$$

$$f_i + b_i h_i + c_i h_i^2 + d_i h_i^3 = f_{i+1}, \quad \text{for } i = 1, 2, \dots, n-1 \quad \cdots (**)$$

Gives $n - 1$ equations.

---

**Condition (iii) --- $C^1$ continuity at interior knots:**

$$S_i'(x) = b_i + 2c_i(x - x_i) + 3d_i(x - x_i)^2$$

At $x_{i+1}$:

$$S_i'(x_{i+1}) = b_i + 2c_i h_i + 3d_i h_i^2 = S_{i+1}'(x_{i+1}) = b_{i+1}$$

$$\therefore \quad b_i + 2c_i h_i + 3d_i h_i^2 = b_{i+1}, \quad \text{for } i = 1, 2, \dots, n-2 \quad \cdots (*)$$

Gives $n - 2$ equations.

---

**Condition (iv) --- $C^2$ continuity at interior knots:**

$$S_i''(x) = 2c_i + 6d_i(x - x_i)$$

At $x_{i+1}$:

$$S_i''(x_{i+1}) = 2c_i + 6d_i h_i = S_{i+1}''(x_{i+1}) = 2c_{i+1}$$

$$\therefore \quad c_{i+1} = c_i + 3d_i h_i \quad \cdots (***)$$

$$\Longrightarrow \quad d_i = \frac{c_{i+1} - c_i}{3 h_i}$$

Gives $n - 2$ equations.

---

**Equation count so far:**

| Source | Count |
|---|---|
| (i) $a_i = f_i$ | $n - 1$ |
| (ii) $C^0$ continuity | $n - 1$ |
| (iii) $C^1$ continuity | $n - 2$ |
| (iv) $C^2$ continuity | $n - 2$ |
| **Total** | $4n - 6$ |

Unknowns: $4(n-1) = 4n - 4$.

We need **2 more equations** (the end conditions).

### 5.3 Eliminating Unknowns --- The Tridiagonal System for $c_i$

The key strategy is to express all unknowns in terms of $c_i$ alone, then solve a tridiagonal system.

**Step 1: Express $d_i$ in terms of $c_i$.**

From condition (iv):

$$d_i = \frac{c_{i+1} - c_i}{3h_i}$$

**Step 2: Express $b_i$ in terms of $c_i$.**

Substitute $d_i$ into the $C^0$ continuity equation $(**)$:

$$f_{i+1} = f_i + b_i h_i + c_i h_i^2 + \frac{c_{i+1} - c_i}{3} h_i^2$$

$$= f_i + b_i h_i + \frac{h_i^2}{3}(c_{i+1} + 2c_i)$$

Solving for $b_i$:

$$b_i = \frac{f_{i+1} - f_i}{h_i} - \frac{h_i}{3}(2c_i + c_{i+1})$$

**Step 3: Derive the tridiagonal equation.**

From the $C^1$ continuity condition $(*)$, substitute $d_i$ to get:

$$b_{i+1} = b_i + h_i(c_i + c_{i+1})$$

Now, reduce the index of the $b_i$ formula by 1:

$$b_{i-1} = \frac{f_i - f_{i-1}}{h_{i-1}} - \frac{h_{i-1}}{3}(2c_{i-1} + c_i)$$

and similarly:

$$b_i = b_{i-1} + h_{i-1}(c_{i-1} + c_i)$$

After algebraic manipulation (substituting the expressions for $b_i$ and $b_{i-1}$ and simplifying), the result is:

$$3\frac{f_{i+1} - f_i}{h_i} - 3\frac{f_i - f_{i-1}}{h_{i-1}} = h_{i-1}\,c_{i-1} + 2(h_{i-1} + h_i)\,c_i + h_i\,c_{i+1}$$

for $i = 2, 3, \dots, n-2$.

> **Note:** The left-hand side involves **divided differences**:
> $$3\,f[x_{i+1}, x_i] - 3\,f[x_i, x_{i-1}]$$
> where $f[x_{i+1}, x_i] = \dfrac{f_{i+1} - f_i}{h_i}$.

This gives $n - 3$ equations for the unknowns $c_1, c_2, \dots, c_{n-1}$ (which are $n - 1$ unknowns). Together with the 2 end conditions, we get $n - 1$ equations in $n - 1$ unknowns.

**Matrix form (with natural end conditions $c_1 = 0$, $c_n = 0$):**

$$\begin{pmatrix} 0 \\ 3(f[x_3,x_2] - f[x_2,x_1]) \\ 3(f[x_4,x_3] - f[x_3,x_2]) \\ \vdots \\ 3(f[x_n,x_{n-1}] - f[x_{n-1},x_{n-2}]) \\ 0 \end{pmatrix} = \begin{pmatrix} 1 & & & & \\ h_1 & 2(h_1+h_2) & h_2 & & \\ & h_2 & 2(h_2+h_3) & h_3 & \\ & & \ddots & \ddots & \ddots \\ & & & h_{n-2} & 2(h_{n-2}+h_{n-1}) & h_{n-1} \\ & & & & & 1 \end{pmatrix} \begin{pmatrix} c_1 \\ c_2 \\ c_3 \\ \vdots \\ c_{n-1} \\ c_n \end{pmatrix}$$

> **Important:** The coefficient matrix is **tridiagonal** and **diagonally dominant**, so it can be solved very efficiently using the **Thomas algorithm** in $O(n)$ operations.

**After solving for $c_i$, recover the other coefficients:**

$$d_i = \frac{c_{i+1} - c_i}{3h_i}, \quad i = 1, 2, \dots, n-1$$

$$b_i = \frac{f_{i+1} - f_i}{h_i} - \frac{h_i}{3}(2c_i + c_{i+1}), \quad i = 1, 2, \dots, n-1$$

### 5.4 End Conditions (Section 18.4.2)

Since the tridiagonal system has $n - 1$ unknowns ($c_1, \dots, c_{n-1}$) but only $n - 3$ interior equations, **2 additional boundary (end) conditions** are required:

| End Condition | Description | Equations |
|---|---|---|
| **Natural** | Second derivatives are zero at both endpoints | $c_1 = 0, \quad c_n = 0$ |
| **Clamped** | First derivatives $f_1'$ and $f_n'$ are specified at the endpoints | Modifies the first and last rows of the tridiagonal system |
| **Not-a-knot** | Continuity of the **3rd derivative** at the second and second-to-last knots | $S_1 = S_2$ (same cubic) and $S_{n-2} = S_{n-1}$ (same cubic) |

> **Natural spline** is the most common default choice. Setting $S_1''(x_1) = 0$ yields $c_1 = 0$, and setting $S_{n-1}''(x_n) = 0$ yields (via $(***)$) $c_n = 0$.
>
> It is **convenient to introduce** $c_n$ as an additional variable (even though the last segment is $S_{n-1}$) to make the matrix system uniform in size $n \times n$.

**Derivation for natural end conditions:**

At the left end:

$$S_1''(x_1) = 2c_1 + 6d_1(x_1 - x_1) = 2c_1 = 0 \implies c_1 = 0$$

At the right end:

$$S_{n-1}''(x_n) = 2c_{n-1} + 6d_{n-1}h_{n-1} = 0$$

Using $c_{i+1} = c_i + 3d_i h_i$ (condition $(***)$), this implies $c_n = 0$.

### 5.5 Example 18.3 --- Cubic Spline

**Given:** 4 data points $(x_1, f_1),\, (x_2, f_2),\, (x_3, f_3),\, (x_4, f_4)$

- 3 intervals
- $3 \times 4 = 12$ unknowns total

**Three spline functions (after $a_i = f_i$):**

$$S_1(x) = f_1 + b_1(x - x_1) + c_1(x - x_1)^2 + d_1(x - x_1)^3$$

$$S_2(x) = f_2 + b_2(x - x_2) + c_2(x - x_2)^2 + d_2(x - x_2)^3$$

$$S_3(x) = f_3 + b_3(x - x_3) + c_3(x - x_3)^2 + d_3(x - x_3)^3$$

**Derivatives:**

$$S_i'(x) = b_i + 2c_i(x - x_i) + 3d_i(x - x_i)^2$$

$$S_i''(x) = 2c_i + 6d_i(x - x_i)$$

**Tridiagonal system (with natural end conditions, $c_1 = 0$, $c_4 = 0$):**

$$\begin{pmatrix} 0 \\ 3(f[x_3, x_2] - f[x_2, x_1]) \\ 3(f[x_4, x_3] - f[x_3, x_2]) \\ 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ h_1 & 2(h_1 + h_2) & h_2 & 0 \\ 0 & h_2 & 2(h_2 + h_3) & h_3 \\ 0 & 0 & 0 & 1 \end{pmatrix} \begin{pmatrix} c_1 \\ c_2 \\ c_3 \\ c_4 \end{pmatrix}$$

**After solving for $c_1, c_2, c_3, c_4$, recover:**

$$b_i = \frac{f_{i+1} - f_i}{h_i} - \frac{h_i}{3}(2c_i + c_{i+1}), \quad i = 1, 2, 3$$

$$d_i = \frac{c_{i+1} - c_i}{3h_i}, \quad i = 1, 2, 3$$

```python
import numpy as np

def cubic_spline_natural(x, f):
    """
    Compute natural cubic spline coefficients.

    Parameters
    ----------
    x : array-like, shape (n,)
        Knot positions (sorted).
    f : array-like, shape (n,)
        Function values at knots.

    Returns
    -------
    a, b, c, d : arrays, each of shape (n-1,)
        Coefficients for S_i(t) = a_i + b_i*(t-x_i)
                                    + c_i*(t-x_i)^2
                                    + d_i*(t-x_i)^3
    """
    n = len(x)
    h = np.diff(x)                     # h[i] = x[i+1] - x[i]

    # --- Build the tridiagonal system for c ---
    # Size: n x n  (with c[0] = 0, c[n-1] = 0 as natural BCs)
    A = np.zeros((n, n))
    rhs = np.zeros(n)

    # Natural end conditions
    A[0, 0] = 1.0
    A[n-1, n-1] = 1.0

    # Interior equations  (i = 1, 2, ..., n-2)
    for i in range(1, n - 1):
        A[i, i-1] = h[i-1]
        A[i, i]   = 2.0 * (h[i-1] + h[i])
        A[i, i+1] = h[i]
        rhs[i] = 3.0 * ((f[i+1] - f[i]) / h[i]
                       - (f[i] - f[i-1]) / h[i-1])

    c = np.linalg.solve(A, rhs)

    # --- Recover b and d ---
    a = f[:-1].copy()
    b = np.zeros(n - 1)
    d = np.zeros(n - 1)

    for i in range(n - 1):
        d[i] = (c[i+1] - c[i]) / (3.0 * h[i])
        b[i] = (f[i+1] - f[i]) / h[i] \
             - h[i] * (2.0 * c[i] + c[i+1]) / 3.0

    c = c[:-1]   # only need c[0] ... c[n-2]

    return a, b, c, d


def evaluate_spline(x_knots, a, b, c, d, x_eval):
    """Evaluate the cubic spline at points x_eval."""
    results = np.zeros_like(x_eval, dtype=float)
    for j, xv in enumerate(x_eval):
        # Find the correct interval
        i = np.searchsorted(x_knots, xv, side='right') - 1
        i = np.clip(i, 0, len(a) - 1)
        dx = xv - x_knots[i]
        results[j] = a[i] + b[i]*dx + c[i]*dx**2 + d[i]*dx**3
    return results


# --- Example usage ---
if __name__ == "__main__":
    # Example 18.3: 4 data points
    x = np.array([3.0, 4.5, 7.0, 9.0])
    f = np.array([2.5, 1.0, 2.5, 0.5])

    a, b, c, d = cubic_spline_natural(x, f)

    print("Spline coefficients:")
    for i in range(len(a)):
        print(f"  S_{i+1}(x): a={a[i]:.6f}, b={b[i]:.6f}, "
              f"c={c[i]:.6f}, d={d[i]:.6f}")

    # Evaluate
    x_test = np.array([5.0])
    y_test = evaluate_spline(x, a, b, c, d, x_test)
    print(f"\nS(5.0) = {y_test[0]:.6f}")
```

---

<br>

## 6. Bilinear Interpolation

Bilinear interpolation extends piecewise interpolation to **two dimensions**. Given function values at the four corners of a rectangle, we can interpolate the value at any interior point $(x_i, y_i)$.

### Setup

Consider a rectangular cell with corners at $(x_1, y_1)$, $(x_2, y_1)$, $(x_1, y_2)$, $(x_2, y_2)$ with known function values $f(x_1, y_1)$, $f(x_2, y_1)$, $f(x_1, y_2)$, $f(x_2, y_2)$.

### Formula

$$f(x_i, y_i) = N_1(x_i, y_i)\,f(x_1, y_1) + N_2(x_i, y_i)\,f(x_2, y_2) + N_3(x_i, y_i)\,f(x_1, y_2) + N_4(x_i, y_i)\,f(x_2, y_1)$$

where the **bilinear basis functions** are:

$$N_1(x, y) = \frac{x - x_2}{x_1 - x_2} \cdot \frac{y - y_2}{y_1 - y_2}$$

$$N_2(x, y) = \frac{x - x_1}{x_2 - x_1} \cdot \frac{y - y_2}{y_1 - y_2}$$

$$N_3(x, y) = \frac{x - x_2}{x_1 - x_2} \cdot \frac{y - y_1}{y_2 - y_1}$$

$$N_4(x, y) = \frac{x - x_1}{x_2 - x_1} \cdot \frac{y - y_1}{y_2 - y_1}$$

> **Interpretation:** Each $N_k$ is a product of two 1D linear interpolation weights. $N_k$ equals 1 at its associated corner and 0 at all other corners.

```python
def bilinear_interpolation(x1, x2, y1, y2, f11, f21, f12, f22, xi, yi):
    """
    Bilinear interpolation on a rectangular cell.

    f11 = f(x1, y1), f21 = f(x2, y1)
    f12 = f(x1, y2), f22 = f(x2, y2)
    """
    N1 = ((x2 - xi) / (x2 - x1)) * ((y2 - yi) / (y2 - y1))   # corner (x1,y1)
    N2 = ((xi - x1) / (x2 - x1)) * ((y2 - yi) / (y2 - y1))   # corner (x2,y1)
    N3 = ((x2 - xi) / (x2 - x1)) * ((yi - y1) / (y2 - y1))   # corner (x1,y2)
    N4 = ((xi - x1) / (x2 - x1)) * ((yi - y1) / (y2 - y1))   # corner (x2,y2)

    return N1 * f11 + N2 * f21 + N3 * f12 + N4 * f22
```

---

<br>

## 7. Summary Table

| Spline Type | Polynomial Order | Continuity | Unknowns per Segment | Extra Conditions Needed | Key Equation |
|---|---|---|---|---|---|
| **Linear** | 1 | $C^0$ | 2 ($a_i, b_i$) | None | $S_i(x) = f_i + \frac{f_{i+1}-f_i}{h_i}(x-x_i)$ |
| **Quadratic** | 2 | $C^1$ | 3 ($a_i, b_i, c_i$) | 1 (set $c_1=0$) | $S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2$ |
| **Cubic** | 3 | $C^2$ | 4 ($a_i, b_i, c_i, d_i$) | 2 (end conditions) | $S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3$ |

| Property | Linear | Quadratic | Cubic |
|---|---|---|---|
| Smoothness at knots | Not smooth | 1st deriv. continuous | 1st & 2nd deriv. continuous |
| System to solve | None (direct formula) | Linear system | Tridiagonal system |
| Computational cost | $O(n)$ | $O(n)$ | $O(n)$ via Thomas algorithm |
| Typical application | Simple interpolation | Moderate smoothness | Computer graphics, engineering |
