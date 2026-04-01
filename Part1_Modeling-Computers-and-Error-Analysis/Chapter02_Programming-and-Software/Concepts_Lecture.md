# Chapter 2 Lecture — Python Fundamentals

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 2

> **Prerequisites**: [Programming Language] Basic programming. [Calculus] Derivatives.
>
> **Learning Objectives**:
> 1. Write structured programs for numerical computation
> 2. Apply modular programming principles
> 3. Debug and test numerical algorithms

---

<br>

## Table of Contents

- [1. Programming Environment](#1-programming-environment)
- [2. Python Basics](#2-python-basics)
  - [2.1 Variables and Assignment](#21-variables-and-assignment)
  - [2.2 Data Types](#22-data-types)
  - [2.3 NumPy Arrays](#23-numpy-arrays)
- [3. Basic Operations](#3-basic-operations)
  - [3.1 Element-wise Operations](#31-element-wise-operations)
  - [3.2 Matrix Operations](#32-matrix-operations)
- [4. Useful NumPy/SciPy Functions](#4-useful-numpyscipy-functions)
- [5. Plotting with Matplotlib](#5-plotting-with-matplotlib)
- [6. Control Structures](#6-control-structures)
  - [6.1 For Loop](#61-for-loop)
  - [6.2 While Loop](#62-while-loop)
  - [6.3 Conditional Statements](#63-conditional-statements)
- [7. Functions](#7-functions)
- [8. Example: Bungee Jumper in Python](#8-example-bungee-jumper-in-python)
  - [8.1 Implementing Euler's Method](#81-implementing-eulers-method)
  - [8.2 Plotting Numerical vs. Analytical Solution](#82-plotting-numerical-vs-analytical-solution)
  - [8.3 Using scipy.integrate.solve_ivp](#83-using-scipyintegratesolve_ivp)
- [Summary](#summary)

---

<br>

## 1. Programming Environment

This course uses **Python** as the primary programming language, together with the scientific computing ecosystem:

- **NumPy** — numerical arrays, linear algebra, and mathematical functions
- **SciPy** — advanced scientific computing (ODE solvers, optimization, interpolation, etc.)
- **Matplotlib** — 2D/3D plotting and visualization

| Feature | MATLAB | Python |
|:--------|:-------|:-------|
| License | Commercial (expensive) | Free and open-source |
| Array handling | Built-in matrix type | NumPy `ndarray` |
| Ecosystem | Toolboxes (paid add-ons) | Packages (free via pip/conda) |
| Syntax style | 1-indexed, end keyword | 0-indexed, indentation-based |
| Community | Engineering/academia | Broad (data science, web, ML, etc.) |

> **Note:** Most textbooks on numerical methods were historically written for MATLAB. The concepts and algorithms are identical regardless of the programming language — only the syntax differs. This course follows the Python edition.

---

<br>

## 2. Python Basics

### 2.1 Variables and Assignment

In Python, variables are created by simple assignment. There is no need to declare types explicitly — Python infers the type from the assigned value:

```python
x = 5        # int
y = 3.14     # float
name = "NA"  # string
```

Variable names are case-sensitive (`x` and `X` are different variables) and must start with a letter or underscore.

### 2.2 Data Types

Python provides several built-in data types relevant to scientific computing:

| Type | Example | Description |
|:-----|:--------|:------------|
| `int` | `5`, `-3` | Integer (arbitrary precision) |
| `float` | `3.14`, `1e-6` | Floating-point number (64-bit double) |
| `str` | `'hello'` | Text string |
| `list` | `[1, 2, 3]` | Mutable ordered sequence |
| `tuple` | `(1, 2, 3)` | Immutable ordered sequence |
| `dict` | `{'a': 1}` | Key-value mapping |

For numerical work, Python's built-in `list` is insufficient because it does not support vectorized arithmetic. This is why we use **NumPy arrays**.

### 2.3 NumPy Arrays

NumPy's `ndarray` is the fundamental data structure for numerical computing in Python. It provides fast, memory-efficient, multi-dimensional arrays with vectorized operations:

```python
import numpy as np

# 1D array (vector)
a = np.array([1, 2, 3])
print(a)       # [1 2 3]
print(a.shape) # (3,)

# 2D array (matrix)
A = np.array([[1, 2],
              [3, 4]])
print(A)
# [[1 2]
#  [3 4]]
print(A.shape) # (2, 2)
```

Key differences between Python lists and NumPy arrays:

| Feature | Python `list` | NumPy `ndarray` |
|:--------|:-------------|:----------------|
| Element types | Mixed allowed | Homogeneous (all same type) |
| Arithmetic | Concatenation (`+`) | Element-wise addition |
| Performance | Slow (interpreted loop) | Fast (compiled C backend) |
| Memory | Higher overhead | Compact, contiguous |

---

<br>

## 3. Basic Operations

### 3.1 Element-wise Operations

Standard arithmetic operators on NumPy arrays act **element-wise**:

```python
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

print(a + b)   # [5 7 9]     — element-wise addition
print(a - b)   # [-3 -3 -3]  — element-wise subtraction
print(a * b)   # [4 10 18]   — element-wise multiplication
print(a / b)   # [0.25 0.4 0.5] — element-wise division
print(a ** 2)  # [1 4 9]     — element-wise power
```

### 3.2 Matrix Operations

For matrix (linear algebra) operations, use the `@` operator or dedicated NumPy functions:

```python
A = np.array([[1, 2],
              [3, 4]])
B = np.array([[5, 6],
              [7, 8]])

# Matrix multiplication
C = A @ B            # preferred syntax (Python 3.5+)
C = np.dot(A, B)     # equivalent

# Transpose
print(A.T)
# [[1 3]
#  [2 4]]
```

> **Note:** Be careful to distinguish between `*` (element-wise multiplication) and `@` (matrix multiplication). This is one of the most common sources of bugs in numerical Python code.

---

<br>

## 4. Useful NumPy/SciPy Functions

The following functions are used frequently throughout this course:

**Array creation:**

| Function | Description | Example |
|:---------|:-----------|:--------|
| `np.linspace(start, stop, num)` | Evenly spaced points | `np.linspace(0, 1, 5)` $\to$ `[0, 0.25, 0.5, 0.75, 1]` |
| `np.arange(start, stop, step)` | Evenly spaced by step size | `np.arange(0, 1, 0.2)` $\to$ `[0, 0.2, 0.4, 0.6, 0.8]` |
| `np.zeros(n)` | Array of zeros | `np.zeros(3)` $\to$ `[0, 0, 0]` |
| `np.ones(n)` | Array of ones | `np.ones(3)` $\to$ `[1, 1, 1]` |
| `np.eye(n)` | $n \times n$ identity matrix | `np.eye(3)` $\to$ $I_3$ |

**Linear algebra (`np.linalg`):**

| Function | Description |
|:---------|:-----------|
| `np.linalg.solve(A, b)` | Solve $Ax = b$ for $x$ |
| `np.linalg.inv(A)` | Compute $A^{-1}$ (matrix inverse) |
| `np.linalg.det(A)` | Compute $\det(A)$ (determinant) |
| `np.linalg.norm(x)` | Compute vector/matrix norm |
| `np.linalg.eig(A)` | Compute eigenvalues and eigenvectors |

```python
A = np.array([[2, 1],
              [5, 3]])
b = np.array([4, 7])

x = np.linalg.solve(A, b)  # Solve Ax = b
print(x)                    # [5. -6.]

print(np.linalg.det(A))    # 1.0
print(np.linalg.inv(A))    # [[ 3. -1.] [-5.  2.]]
```

---

<br>

## 5. Plotting with Matplotlib

Matplotlib is the standard plotting library for Python. The `pyplot` interface provides MATLAB-like plotting commands:

```python
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 2 * np.pi, 100)
y = np.sin(x)

plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('sin(x)')
plt.title('Sine Function')
plt.grid(True)
plt.show()
```

Multiple curves can be plotted on the same figure:

```python
plt.plot(x, np.sin(x), label='sin(x)')
plt.plot(x, np.cos(x), label='cos(x)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Trigonometric Functions')
plt.legend()
plt.grid(True)
plt.show()
```

Common customizations include line styles (`'-'`, `'--'`, `':'`), colors (`'r'`, `'b'`, `'g'`), and markers (`'o'`, `'s'`, `'^'`).

---

<br>

## 6. Control Structures

### 6.1 For Loop

The `for` loop iterates over a sequence of values. The `range()` function generates integer sequences:

```python
# Sum of 1 to 10
total = 0
for i in range(1, 11):
    total += i
print(total)  # 55
```

`range(start, stop, step)` generates integers from `start` up to (but not including) `stop`.

### 6.2 While Loop

The `while` loop repeats as long as a condition is `True`:

```python
count = 10
while count > 0:
    print(count, end=' ')
    count -= 1
# Output: 10 9 8 7 6 5 4 3 2 1
```

While loops are especially useful when the number of iterations is not known in advance (e.g., iterating until convergence in a numerical algorithm).

### 6.3 Conditional Statements

The `if`/`elif`/`else` structure controls branching:

```python
x = 42

if x > 0:
    print('Positive')
elif x == 0:
    print('Zero')
else:
    print('Negative')
```

---

<br>

## 7. Functions

Functions encapsulate reusable blocks of code. They are defined with the `def` keyword:

```python
def myfunction(x):
    return x ** 2

print(myfunction(5))  # 25
```

Functions can return multiple values as a tuple:

```python
def circle_properties(r):
    area = 3.14159 * r ** 2
    circumference = 2 * 3.14159 * r
    return area, circumference

a, c = circle_properties(5)
print(f'Area = {a:.2f}, Circumference = {c:.2f}')
# Area = 78.54, Circumference = 31.42
```

Default arguments allow flexible function calls:

```python
def greet(name, greeting='Hello'):
    return f'{greeting}, {name}!'

print(greet('World'))            # Hello, World!
print(greet('World', 'Hi'))      # Hi, World!
```

---

<br>

## 8. Example: Bungee Jumper in Python

This section implements the bungee jumper problem from Chapter 1 in Python, demonstrating how to translate mathematical formulations into working code.

### 8.1 Implementing Euler's Method

Recall the ODE for a bungee jumper with linear drag:

$$\frac{dv}{dt} = g - \frac{c}{m}v$$

Euler's method approximation:

$$v(t_{i+1}) = v(t_i) + \left(g - \frac{c}{m}v(t_i)\right) \cdot h$$

```python
import numpy as np

# Parameters
g = 9.81      # m/s^2
m = 68.1      # kg
c = 12.5      # kg/s
h = 2.0       # step size (s)
t_end = 12.0  # final time (s)

# Initialize
t = np.arange(0, t_end + h, h)
v = np.zeros(len(t))
v[0] = 0  # initial condition

# Euler's method
for i in range(len(t) - 1):
    dvdt = g - (c / m) * v[i]
    v[i + 1] = v[i] + dvdt * h

print('t (s)  |  v (m/s)')
for ti, vi in zip(t, v):
    print(f'{ti:5.1f}  |  {vi:.4f}')
```

### 8.2 Plotting Numerical vs. Analytical Solution

The analytical solution (from Chapter 1) is:

$$v(t) = \frac{gm}{c}\left(1 - e^{-(c/m)t}\right)$$

```python
import matplotlib.pyplot as plt

# Analytical solution on a fine grid
t_fine = np.linspace(0, t_end, 200)
v_analytical = (g * m / c) * (1 - np.exp(-(c / m) * t_fine))

# Plot comparison
plt.plot(t_fine, v_analytical, 'b-', label='Analytical')
plt.plot(t, v, 'ro--', label=f'Euler (h={h}s)')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Bungee Jumper: Euler vs. Analytical')
plt.legend()
plt.grid(True)
plt.show()
```

### 8.3 Using scipy.integrate.solve_ivp

For more accurate ODE solving, SciPy provides `solve_ivp` which uses adaptive step-size methods:

```python
from scipy.integrate import solve_ivp

# Define the ODE as a function: dv/dt = f(t, v)
def bungee_ode(t, v):
    return g - (c / m) * v

# Solve
sol = solve_ivp(bungee_ode, [0, t_end], [0], t_eval=t_fine)

plt.plot(t_fine, v_analytical, 'b-', label='Analytical')
plt.plot(sol.t, sol.y[0], 'g--', label='solve_ivp (RK45)')
plt.plot(t, v, 'ro--', label=f'Euler (h={h}s)')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Bungee Jumper: Comparison of Methods')
plt.legend()
plt.grid(True)
plt.show()
```

> **[Python]** `scipy.integrate.solve_ivp` uses adaptive Runge-Kutta methods by default, which are much more accurate than Euler's method. We'll study these methods in detail in Chapter 22.

---

<br>

## Summary

| Topic | Key Point |
|:------|:----------|
| Programming Environment | Python + NumPy + SciPy + Matplotlib (free, open-source alternative to MATLAB) |
| Variables & Types | Dynamic typing; `int`, `float`, `str`, `list`, NumPy `ndarray` |
| NumPy Arrays | `np.array(...)` for vectors and matrices; homogeneous, fast, memory-efficient |
| Element-wise Operations | `+`, `-`, `*`, `/` act element-wise on arrays |
| Matrix Operations | `A @ B` for matrix multiply; `A.T` for transpose |
| Key Functions | `np.linspace`, `np.zeros`, `np.ones`, `np.eye`, `np.linalg.solve` |
| Plotting | `plt.plot(x, y)`, `plt.xlabel(...)`, `plt.title(...)`, `plt.show()` |
| Control Structures | `for`, `while`, `if`/`elif`/`else` |
| Functions | `def func(args): return result`; multiple return values via tuples |
| ODE Solving | Euler's method (manual) vs. `scipy.integrate.solve_ivp` (adaptive, accurate) |

---
