# Chapter 1 Lab — Python Fundamentals for Scientific Computing

> **Last Updated:** 2026-04-01
>
> Chapra & Canale, Numerical Methods for Engineers 8th Ed. Ch 1

> **Prerequisites**: [Programming Language] MATLAB/Python. No prior numerical analysis knowledge required. [Calculus] Basic derivatives and integrals.
>
> **Learning Objectives**:
> 1. Explain the role of mathematical modeling in engineering
> 2. Identify sources of error in numerical computations
> 3. Describe the trade-off between model complexity and accuracy

---

<br>

## Table of Contents

- [1. Variables and Assignment](#1-variables-and-assignment)
  - [1.1 Dynamic Typing](#11-dynamic-typing)
  - [1.2 Mathematical Constants and Formatted Output](#12-mathematical-constants-and-formatted-output)
- [2. NumPy Arrays, Vectors, and Matrices](#2-numpy-arrays-vectors-and-matrices)
  - [2.1 Why NumPy?](#21-why-numpy)
  - [2.2 Creating Arrays](#22-creating-arrays)
  - [2.3 Indexing and Slicing](#23-indexing-and-slicing)
  - [2.4 Sequence Generation: arange, linspace, logspace](#24-sequence-generation-arange-linspace-logspace)
- [3. Mathematical Operations](#3-mathematical-operations)
  - [3.1 Arithmetic Operators and Precedence](#31-arithmetic-operators-and-precedence)
  - [3.2 Element-wise vs. Matrix Multiplication](#32-element-wise-vs-matrix-multiplication)
  - [3.3 Dot Product, Cross Product, and Transpose](#33-dot-product-cross-product-and-transpose)
- [4. Built-in Functions](#4-built-in-functions)
  - [4.1 Common Built-in Functions](#41-common-built-in-functions)
  - [4.2 Banker's Rounding](#42-bankers-rounding)
  - [4.3 Application: Bungee Jumper Velocity](#43-application-bungee-jumper-velocity)
- [5. Visualization with Matplotlib](#5-visualization-with-matplotlib)
  - [5.1 Basic Line Plot](#51-basic-line-plot)
  - [5.2 Multiple Curves with Legend](#52-multiple-curves-with-legend)
- [6. Views vs. Copies in NumPy](#6-views-vs-copies-in-numpy)
  - [6.1 Views (Slicing)](#61-views-slicing)
  - [6.2 Copies (Fancy Indexing)](#62-copies-fancy-indexing)
- [Summary](#summary)

---

<br>

## 1. Variables and Assignment

### 1.1 Dynamic Typing

Python is a **dynamically typed** language, meaning you do not need to declare the type of a variable before using it. The type is determined automatically at runtime based on the value you assign. Variable names are **case-sensitive** (`myVar` and `myvar` are different variables).

You can check the type of any variable using the built-in `type()` function. Python also supports **multiple assignment** in a single line:

```python
a, b, c = 1, 2, 3  # assigns 1 to a, 2 to b, 3 to c
```

```python
a = 4
print(type(a))  # <class 'int'>

c = 2 + 4j
print(type(c))  # <class 'complex'>
```

> **[Python]** Python's dynamic typing means variables don't need explicit type declarations. The type is determined at runtime based on the assigned value. This is different from statically typed languages like C++ or Java.

### 1.2 Mathematical Constants and Formatted Output

Python's `math` module provides fundamental mathematical constants and functions. For physical constants used in scientific computing, the `scipy.constants` module is invaluable.

```python
import math
print(math.pi)   # 3.141592653589793
print(math.e)    # 2.718281828459045

# Formatted output using str.format()
# {0:7.4f} means: argument index 0, total width 7, 4 decimal places, fixed-point
print('{0:7.4f}'.format(math.pi))  # ' 3.1416'

# Physical constants from scipy
import scipy.constants as pc
print('{0:12.4e}'.format(pc.h))  # Planck's constant: 6.6261e-34
```

> **Note:** The format specifier `{0:7.4f}` breaks down as: `0` = first argument, `7` = minimum total width (including decimal point), `.4` = four decimal places, `f` = fixed-point notation. For scientific notation, use `e` instead of `f`.

---

<br>

## 2. NumPy Arrays, Vectors, and Matrices

### 2.1 Why NumPy?

**NumPy** (Numerical Python) is the foundational library for scientific computing in Python. Its key advantages are:

- **Speed**: NumPy arrays are backed by C code and stored in contiguous memory blocks, making operations orders of magnitude faster than Python lists
- **Vectorization**: Operations are applied element-wise to entire arrays without explicit Python loops
- **Broadcasting**: Allows arithmetic between arrays of different shapes following well-defined rules

> **[Python]** NumPy arrays are stored in contiguous memory blocks (unlike Python lists), enabling vectorized operations that run at C speed. Broadcasting allows operations between arrays of different shapes without explicit loops.

### 2.2 Creating Arrays

NumPy provides several functions for creating arrays:

| Function | Description | Example |
|:---------|:-----------|:--------|
| `np.array(list)` | Create array from a Python list | `np.array([1, 2, 3])` |
| `np.zeros(shape)` | Array of zeros | `np.zeros((3, 4))` |
| `np.ones(shape)` | Array of ones | `np.ones((2, 3))` |
| `np.eye(n)` | Identity matrix ($n \times n$) | `np.eye(3)` |
| `np.arange(start, stop, step)` | Evenly spaced values (stop excluded) | `np.arange(0, 1, 0.1)` |
| `np.linspace(start, stop, n)` | $n$ evenly spaced values (stop included) | `np.linspace(0, 1, 11)` |
| `np.logspace(start, stop, n)` | Log-spaced values from $10^{\text{start}}$ to $10^{\text{stop}}$ | `np.logspace(-1, 2, 4)` |

> **Note:** `np.matrix` is deprecated — always use `np.array` instead. 2D arrays support all matrix operations via the `@` operator.

```python
import numpy as np

# 1D array (vector)
data = np.array([12.2, 10.9, 13.6, 8.4, 11.1])
print(data.shape)  # (5,)

# 2D array (matrix)
A = np.array([[2, 4], [1, 3]])
print(A.shape)  # (2, 2)

# Special arrays
Z = np.zeros((5, 3))   # 5x3 matrix of zeros
O = np.ones((2, 3))    # 2x3 matrix of ones
I = np.eye(3)          # 3x3 identity matrix
```

### 2.3 Indexing and Slicing

NumPy uses **zero-based indexing**. For 2D arrays, indices are specified as `[row, col]`. The colon `:` selects entire rows or columns:

```python
A = np.array([[2, 4], [1, 3]])

print(A[0, 1])   # 4 (row 0, col 1)
print(A[0, :])   # [2, 4] (entire row 0)
print(A[:, 1])   # [4, 3] (entire column 1)
```

### 2.4 Sequence Generation: arange, linspace, logspace

These three functions generate sequences of numbers, but they differ in how the spacing is specified:

**`np.arange`** — specified by step size, endpoint **not** included:

```python
x = np.arange(0, 1, 0.1)  # [0.0, 0.1, 0.2, ..., 0.9] — 10 values, 1.0 NOT included
```

**`np.linspace`** — specified by number of points, endpoint **included** by default:

```python
x = np.linspace(0, 1, 11)  # [0.0, 0.1, 0.2, ..., 1.0] — 11 values, 1.0 included
x = np.linspace(0, 1, 10, endpoint=False)  # [0.0, 0.1, ..., 0.9] — excludes 1.0
```

**`np.logspace`** — logarithmically spaced values:

```python
# Default base=10: from 10^(-1) to 10^2
x = np.logspace(-1, 2, 4)  # [0.1, 1, 10, 100]

# Custom base: from 2^1 to 2^6
x = np.logspace(1, 6, 6, base=2.0)  # [2, 4, 8, 16, 32, 64]
```

---

<br>

## 3. Mathematical Operations

### 3.1 Arithmetic Operators and Precedence

Python provides the following arithmetic operators:

| Operator | Description | Example | Result |
|:---------|:-----------|:--------|:-------|
| `+` | Addition | `3 + 2` | `5` |
| `-` | Subtraction | `3 - 2` | `1` |
| `*` | Multiplication (element-wise for arrays) | `3 * 2` | `6` |
| `/` | Division | `7 / 2` | `3.5` |
| `//` | Floor division | `7 // 2` | `3` |
| `%` | Modulus (remainder) | `7 % 2` | `1` |
| `**` | Exponentiation | `2 ** 3` | `8` |
| `@` | Matrix multiplication | `A @ B` | Matrix product |

**Operator precedence** follows standard mathematical convention (PEMDAS). An important subtlety: **exponentiation is right-to-left associative**:

```python
x, y, z = 2, 3, 2
print(x ** y ** z)  # 2 ** (3 ** 2) = 2 ** 9 = 512, NOT (2 ** 3) ** 2 = 64
```

### 3.2 Element-wise vs. Matrix Multiplication

This is one of the most important distinctions in NumPy. The `*` operator performs **element-wise** multiplication, while `@` performs **matrix multiplication**:

```python
import numpy as np

# 1D arrays
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

print(a * b)  # [4, 10, 18] — element-wise: [1*4, 2*5, 3*6]
print(a @ b)  # 32 — dot product: 1*4 + 2*5 + 3*6

# 2D arrays (matrices)
A = np.array([[2, 4], [1, 3]])

print(A @ A)              # matrix multiplication
# [[8, 20],
#  [5, 13]]

print(np.multiply(A, A))  # element-wise multiplication
# [[4, 16],
#  [1,  9]]
```

### 3.3 Dot Product, Cross Product, and Transpose

```python
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])

# Dot product: sum of element-wise products
print(np.dot(a, b))    # 32 (same as a @ b for 1D arrays)

# Cross product: vector perpendicular to both inputs
print(np.cross(a, b))  # [-3, 6, -3]

# Transpose: swap rows and columns
A = np.array([[2, 4], [1, 3]])
print(A.T)
# [[2, 1],
#  [4, 3]]
```

---

<br>

## 4. Built-in Functions

### 4.1 Common Built-in Functions

Python provides several built-in functions that are essential for everyday computation:

| Function | Description | Example |
|:---------|:-----------|:--------|
| `abs(x)` | Absolute value (or magnitude for complex) | `abs(-6 + 4j)` → `7.211...` |
| `round(x, n)` | Round to $n$ decimal places | `round(3.14159, 2)` → `3.14` |
| `min(...)` | Minimum value | `min(3, 1, 4)` → `1` |
| `max(...)` | Maximum value | `max(3, 1, 4)` → `4` |
| `len(x)` | Length of a sequence | `len([1, 2, 3])` → `3` |
| `type(x)` | Type of an object | `type(3.14)` → `<class 'float'>` |
| `float(x)` | Convert to float | `float('3.14')` → `3.14` |
| `int(x)` | Convert to integer (truncates) | `int(3.7)` → `3` |
| `help(x)` | Display documentation | `help(abs)` |

```python
print(abs(-6 + 4j))  # 7.211102550927978 (magnitude: sqrt(36 + 16))
```

### 4.2 Banker's Rounding

Python's `round()` function uses **banker's rounding** (also called "round half to even"). When a value is exactly halfway between two integers, it rounds to the **nearest even number**:

```python
print(round(4.5))     # 4 (rounds to nearest even)
print(round(4.51))    # 5 (not exactly halfway, rounds normally)
print(round(0.5))     # 0 (rounds to nearest even)
print(round(1.5))     # 2 (rounds to nearest even)
print(round(2.5))     # 2 (rounds to nearest even)
print(round(3.5))     # 4 (rounds to nearest even)
```

> **Note:** Python uses "banker's rounding" (round half to even): `round(0.5) = 0`, `round(1.5) = 2`. This minimizes cumulative rounding bias in statistical computations.

### 4.3 Application: Bungee Jumper Velocity

Using the analytical solution for the bungee jumper with quadratic drag, we can compute the velocity vectorized across all time points:

$$v(t) = \sqrt{\frac{mg}{c_d}} \tanh\!\left(\sqrt{\frac{c_d\, g}{m}}\, t\right)$$

```python
import math
import numpy as np

m, g, cd = 68.1, 9.81, 0.25
tm = np.linspace(0, 20, 11)  # time from 0 to 20 s, 11 points

v = math.sqrt(m * g / cd) * np.tanh(math.sqrt(cd * g / m) * tm)
print(v)
# [ 0.          9.63...  17.85...  23.53...  27.29...  29.67...
#  31.09...  31.87...  32.26...  32.44...  32.52...]
```

> **[Python]** The `np.tanh()` function operates element-wise on the entire `tm` array, producing a result array of the same shape without any explicit loop. This is vectorization in action.

---

<br>

## 5. Visualization with Matplotlib

### 5.1 Basic Line Plot

**Matplotlib** is the standard plotting library for Python. The `pyplot` module provides a MATLAB-like interface for creating visualizations:

```python
import matplotlib.pyplot as plt
import numpy as np
import math

m, g, cd = 68.1, 9.81, 0.25
tm = np.linspace(0, 20, 11)
v = math.sqrt(m * g / cd) * np.tanh(math.sqrt(cd * g / m) * tm)

plt.plot(tm, v)
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Bungee Jumper Velocity')
plt.grid(True)
plt.show()
```

### 5.2 Multiple Curves with Legend

To compare different scenarios on the same plot, call `plt.plot()` multiple times and use `label` with `plt.legend()`:

```python
import matplotlib.pyplot as plt
import numpy as np
import math

m, g = 68.1, 9.81
tm = np.linspace(0, 20, 11)

cd1, cd2 = 0.25, 0.30
v1 = math.sqrt(m * g / cd1) * np.tanh(math.sqrt(cd1 * g / m) * tm)
v2 = math.sqrt(m * g / cd2) * np.tanh(math.sqrt(cd2 * g / m) * tm)

plt.plot(tm, v1, ls='-', marker='o', label=f'cd = {cd1}')
plt.plot(tm, v2, ls='--', marker='D', label=f'cd = {cd2}')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Bungee Jumper: Effect of Drag Coefficient')
plt.legend()
plt.grid(True)
plt.savefig('sample.png', bbox_inches='tight')
plt.show()
```

> **[Python]** The `ls` parameter controls line style (`'-'` solid, `'--'` dashed, `':'` dotted), `marker` controls data point markers (`'o'` circle, `'D'` diamond, `'s'` square), and `label` provides the text for the legend entry. `bbox_inches='tight'` removes excess whitespace when saving.

---

<br>

## 6. Views vs. Copies in NumPy

Understanding the difference between **views** and **copies** is critical in scientific computing. Modifying a view accidentally changes the original data, which can lead to subtle bugs.

### 6.1 Views (Slicing)

When you create a subset of an array using **basic slicing** (e.g., `x[1:3]`), NumPy returns a **view** — a window into the same underlying memory. Changes to the view affect the original array:

```python
import numpy as np

x = np.array([1, 2, 3, 4, 5])
y = x[1:3]       # y is a VIEW of x — shares the same memory
print(y)          # [2, 3]

x[2] = 99         # modify the original
print(y)          # [2, 99] — y changed too, because it's a view!

# Verify that y is a view
print(y.base)     # [1, 2, 99, 4, 5] — shows the original (base) array
```

### 6.2 Copies (Fancy Indexing)

When you use **fancy indexing** (indexing with a list or array of indices), NumPy returns a **copy** — a completely independent array:

```python
import numpy as np

x = np.array([1, 2, 3, 4, 5])
y = x[[1, 2]]    # y is a COPY — independent memory
print(y)          # [2, 3]

x[2] = 99         # modify the original
print(y)          # [2, 3] — y is unchanged, because it's a copy

# Verify that y is a copy
print(y.base is None)  # True — no base array, it's independent
```

> **Note:** This distinction is critical in scientific computing. Modifying a view accidentally changes the original data. Use `.copy()` explicitly when you need an independent copy: `y = x[1:3].copy()`.

---

<br>

## Summary

| Topic | Key Functions / Concepts |
|:------|:------------------------|
| Variables | `type()`, dynamic typing, `math.pi`, `scipy.constants` |
| NumPy Arrays | `np.array`, `np.zeros`, `np.ones`, `np.eye`, `np.linspace` |
| Math Ops | `*` (element-wise), `@` (matrix), `np.dot`, `np.cross` |
| Built-in | `abs`, `round` (banker's rounding), `min`, `max`, `help()` |
| Plotting | `plt.plot`, `plt.xlabel`, `plt.legend`, `plt.savefig` |
| View vs. Copy | Slicing → view (shared memory), Fancy indexing → copy (independent) |

---
